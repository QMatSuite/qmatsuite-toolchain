!
! Copyright (C) 2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------
MODULE  read_upf_new_module
  !-----------------------------------------------------
  !! this module contains the simplified code for reading
  !! pseudopotential files in either UPF v.2 or xml
  !
  USE xmltools
  USE upf_io,    ONLY: stdout
  USE upf_kinds, ONLY: dp
  USE pseudo_types, ONLY: pseudo_upf, pseudo_config
  !
  LOGICAL :: v2
  !! true if UPF v.2 version, false if new UPF with xml schema
  INTEGER :: iun
  !! unit for reading data
  !
  PUBLIC
  !
CONTAINS
  !
  !------------------------------------------------+
  SUBROUTINE read_upf_new (filename, upf, ierr)         !
    !---------------------------------------------+
    !! Reads pseudopotential in UPF format (either v.2 or upf_schema).
    !! Derived-type variable *upf* store in output the data read from file. 
    !! File *filename* is opened and closed inside the routine
    !
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: filename  
    !! i/o filename
    TYPE(pseudo_upf),INTENT(INOUT) :: upf
    !! the derived type storing the pseudo data
    !! INOUT because many variables are reset to default values in input
    INTEGER, INTENT(OUT) :: ierr
    !! ierr= -2 : UPF v.2
    !! ierr=  0 : xml schema
    !! ierr=1-4 : error reading PP file
    !! ierr= 81 : error opening PP file
    !
    CHARACTER(LEN=64) :: tag_buf  ! Fixed-length buffer for tag names (MinGW-safe)
    !
    iun = xml_open_file ( filename )
    if ( iun == -1 ) THEN
       ierr = 81
       go to 10
    end if
    call xmlr_opentag ( 'qe_pp:pseudo', IERR = ierr )
    if ( ierr == 0 ) then
       v2 =.false.
    else if ( ierr == 1 ) then
       rewind (iun) 
       call xmlr_opentag ( 'UPF', IERR = ierr )
       if ( ierr == 0 ) then
          v2 =.true.
          CALL get_attr ( 'version', upf%nv )
       end if
    else
       go to 10
    end if
    if ( ierr > 0 ) go to 10
    !
    ! The header sections differ a lot between UPF v.2 and UPF with schema
    !
    IF ( v2 ) THEN
       CALL read_pp_header_v2 ( upf )
    ELSE
       CALL read_pp_header_schema ( upf )
    END IF
    ! compatibility
    upf%is_gth = .false.
    upf%is_multiproj = .true.
    !
    ! From here on the format of v2 and schema do not differ much:
    ! the most frequent difference is capitalization of tags
    ! (see function capitalize_if_v2)
    !
    CALL read_pp_mesh ( upf )
    !
    allocate ( upf%rho_atc(upf%mesh) )
    !! FIXME: this is needed only if the nonlinear core correction is used,
    !! FIXME: but with PAW the pseudo-core charge is used also if no nlcc
    IF(upf%nlcc) then
       CALL capitalize_if_v2_into('pp_nlcc', v2, tag_buf)
       CALL xmlr_readtag( TRIM(tag_buf), &
            upf%rho_atc(:) )
    else
       upf%rho_atc(:) = 0.0_dp
    end if
    IF( .NOT. upf%tcoulombp) then
       allocate ( upf%vloc(upf%mesh) )
       CALL capitalize_if_v2_into('pp_local', v2, tag_buf)
       CALL xmlr_readtag( TRIM(tag_buf), &
            upf%vloc(:), ierr )
       !
       ! existing PP files may have pp_nlcc first, pp_local later,
       ! but also the other way round - check that everything was right
       !
       if ( ierr ==-10 ) ierr = 0
       if ( ierr /= 0 ) go to 10
    end if
    !
    CALL read_pp_semilocal ( upf, ierr )
    if ( ierr > 0 ) go to 10
    !
    CALL read_pp_nonlocal ( upf, ierr )
    if ( ierr > 0 ) go to 10
    !
    CALL read_pp_pswfc ( upf, ierr )
    if ( ierr > 0 ) go to 10
    !
    CALL read_pp_full_wfc ( upf, ierr )
    if ( ierr > 0 ) go to 10
    !
    ALLOCATE( upf%rho_at(1:upf%mesh) )
    CALL capitalize_if_v2_into('pp_rhoatom', v2, tag_buf)
    CALL xmlr_readtag( TRIM(tag_buf), &
         upf%rho_at(1:upf%mesh) )
    !
    CALL read_pp_metagga ( upf, ierr)
    if ( ierr > 0 ) go to 10
    !
    CALL read_pp_spinorb ( upf, ierr )
    if ( ierr > 0 ) go to 10
    !
    CALL read_pp_paw ( upf, ierr )
    if ( ierr > 0 ) go to 10
    !
    CALL read_pp_gipaw ( upf, ierr )
    if ( ierr > 0 ) go to 10
    !
    ! close initial tag, qe_pp:pseudo or UPF
    !
    CALL xmlr_closetag ( )
    !
    CALL xml_closefile ( )
    !
    ! normal return
    if ( v2 ) ierr = -2
    return
    !
    ! error return
10  call xml_closefile( )
    return

  END SUBROUTINE read_upf_new
  !
  SUBROUTINE capitalize_if_v2_into ( strin, v2_flag, strout )
    !
    ! Robust MinGW-safe version: writes capitalized string for UPF v.2, 
    ! same string otherwise, into fixed-length output buffer.
    ! (UPF v.2 uses capitalized tags, UPF with schema use lowercase)
    !
    ! This subroutine avoids allocatable CHARACTER returns that cause
    ! empty-tag bugs on MinGW gfortran.
    !
    USE upf_utils, ONLY: capital
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN)  :: strin
    LOGICAL,          INTENT(IN)  :: v2_flag
    CHARACTER(LEN=*), INTENT(OUT) :: strout
    !
    INTEGER :: n, lt, lout
    !
    ! Initialize output to empty
    strout = ''
    !
    ! Get trimmed input length
    lt = LEN_TRIM(strin)
    lout = LEN(strout)
    !
    ! Guard: detect empty input
    IF (lt == 0) THEN
       WRITE(stdout,'("FATAL: capitalize_if_v2_into called with empty string")')
       RETURN
    END IF
    !
    ! Guard: ensure output buffer is large enough
    IF (lout < lt) THEN
       WRITE(stdout,'("FATAL: capitalize_if_v2_into: output buffer too small")')
       WRITE(stdout,'("  Input length: ",I0,", buffer size: ",I0)') lt, lout
       ERROR STOP 'capitalize_if_v2_into: buffer overflow'
    END IF
    !
    IF ( v2_flag ) THEN
       ! UPF v.2: capitalize each character
       ! Use direct character assignment (NO concatenation)
       DO n = 1, lt
          strout(n:n) = capital(strin(n:n))
       END DO
       ! Pad remainder with blanks (already done by initialization)
    ELSE
       ! UPF schema: copy trimmed input
       strout(1:lt) = strin(1:lt)
       ! Remainder already blank from initialization
    END IF
    !
    ! Final safety check: ensure result is non-empty
    IF (LEN_TRIM(strout) == 0 .AND. lt > 0) THEN
       WRITE(stdout,'("FATAL: capitalize_if_v2_into produced empty string from non-empty input")')
       WRITE(stdout,'("  Input: [",A,"], v2_flag=",L1,", lt=",I0)') strin, v2_flag, lt
       ERROR STOP 'capitalize_if_v2_into: empty result'
    END IF
    !
  END SUBROUTINE capitalize_if_v2_into
  !--------------------------------------------------------
  SUBROUTINE read_pp_header_schema ( upf )
    !--------------------------------------------------------
    !
    IMPLICIT NONE
    TYPE(pseudo_upf), INTENT(INOUT) :: upf ! the pseudo data
    !
    CHARACTER(LEN=64) :: tag_buf
    !
    CALL capitalize_if_v2_into('pp_header', v2, tag_buf)
    CALL xmlr_opentag( TRIM(tag_buf) )
    !
    CALL xmlr_readtag( 'element', upf%psd )
    CALL xmlr_readtag( 'z_valence', upf%zp )
    CALL xmlr_readtag( 'type', upf%typ )
    CALL xmlr_readtag( 'functional', upf%dft )
    CALL xmlr_readtag( 'relativistic', upf%rel )
    CALL xmlr_readtag( 'is_ultrasoft', upf%tvanp )
    CALL xmlr_readtag( 'is_paw', upf%tpawp )
    CALL xmlr_readtag( 'is_coulomb', upf%tcoulombp )
    CALL xmlr_readtag( 'has_so', upf%has_so )
    CALL xmlr_readtag( 'has_wfc', upf%has_wfc )
    CALL xmlr_readtag( 'has_gipaw', upf%has_gipaw )
    CALL xmlr_readtag( 'paw_as_gipaw', upf%paw_as_gipaw)
    CALL xmlr_readtag( 'core_correction', upf%nlcc)
    CALL xmlr_readtag( 'with_metagga_info', upf%with_metagga_info )
    CALL xmlr_readtag( 'total_psenergy', upf%etotps )
    CALL xmlr_readtag( 'wfc_cutoff', upf%ecutwfc )
    CALL xmlr_readtag( 'rho_cutoff', upf%ecutrho )
    CALL xmlr_readtag( 'l_max', upf%lmax )
    CALL xmlr_readtag( 'l_max_rho', upf%lmax_rho )
    CALL xmlr_readtag( 'l_local', upf%lloc )
    CALL xmlr_readtag( 'mesh_size', upf%mesh )
    CALL xmlr_readtag( 'number_of_wfc', upf%nwfc )
    CALL xmlr_readtag( 'number_of_proj', upf%nbeta )
    !
    CALL xmlr_closetag( )
    !
  END SUBROUTINE read_pp_header_schema
  !
  !--------------------------------------------------------
  SUBROUTINE read_pp_header_v2 ( upf )
    !--------------------------------------------------------
    !
    IMPLICIT NONE
    TYPE(pseudo_upf), INTENT(INOUT) :: upf ! the pseudo data
    !
    CHARACTER(LEN=1) :: dummy
    CHARACTER(LEN=64) :: tag_buf
    !
    ! Robust MinGW-safe tag handling: use fixed-length buffer
    CALL capitalize_if_v2_into('pp_header', v2, tag_buf)
    CALL xmlr_readtag ( TRIM(tag_buf), dummy )
    CALL get_attr ('generated', upf%generated)
    CALL get_attr ('author', upf%author)
    CALL get_attr ('date', upf%date)
    CALL get_attr ('comment', upf%comment)
    CALL get_attr ('element', upf%psd)
    CALL get_attr ('pseudo_type', upf%typ)
    CALL get_attr ('relativistic', upf%rel)
    CALL get_attr ('is_ultrasoft', upf%tvanp)
    CALL get_attr ('is_paw', upf%tpawp)
    CALL get_attr ('is_coulomb', upf%tcoulombp)
    CALL get_attr ('has_so', upf%has_so)
    CALL get_attr ('has_wfc', upf%has_wfc)
    CALL get_attr ('has_gipaw', upf%has_gipaw)
    CALL get_attr ('paw_as_gipaw', upf%paw_as_gipaw)
    CALL get_attr ('core_correction', upf%nlcc)
    CALL get_attr( 'with_metagga_info', upf%with_metagga_info )
    CALL get_attr ('functional', upf%dft)
    CALL get_attr ('z_valence', upf%zp)
    CALL get_attr ('total_psenergy', upf%etotps)
    CALL get_attr ('wfc_cutoff', upf%ecutwfc)
    CALL get_attr ('rho_cutoff', upf%ecutrho)
    CALL get_attr ('l_max', upf%lmax)
    CALL get_attr ('l_max_rho', upf%lmax_rho)
    CALL get_attr ('l_local', upf%lloc)
    CALL get_attr ('mesh_size', upf%mesh)
    CALL get_attr ('number_of_wfc', upf%nwfc)
    CALL get_attr ('number_of_proj', upf%nbeta )
    !
  END SUBROUTINE read_pp_header_v2
  !
  !--------------------------------------------------------
  SUBROUTINE read_pp_mesh ( upf )
    !--------------------------------------------------------
    !
    IMPLICIT NONE
    TYPE(pseudo_upf),INTENT(INOUT) :: upf ! the pseudo data
    integer :: mesh
    CHARACTER(LEN=64) :: tag_buf
    !
    CALL capitalize_if_v2_into('pp_mesh', v2, tag_buf)
    CALL xmlr_opentag( TRIM(tag_buf) )
    CALL get_attr ( 'mesh', mesh )
    if ( mesh == 0 ) THEN
#if defined (__debug)
       WRITE(stdout,'("read_pp_mesh: mesh size missing, using the one in header")'
#else
       continue
#endif
    else if ( mesh /= upf%mesh ) THEN
       WRITE(stdout,'("read_pp_mesh; mismatch in mesh size, discarding the one in header")')
       upf%mesh = mesh
    end if
    CALL get_attr ( 'dx'  , upf%dx   )
    CALL get_attr ( 'xmin', upf%xmin )
    CALL get_attr ( 'rmax', upf%rmax )
    CALL get_attr ( 'zmesh', upf%zmesh )
    allocate ( upf%r(1:upf%mesh) )
    CALL capitalize_if_v2_into('pp_r', v2, tag_buf)
    CALL xmlr_readtag( TRIM(tag_buf), upf%r(1:upf%mesh) )
    ! Debug: verify r array
    WRITE(stdout,'("DEBUG read_pp_mesh: r array size=",I0,", mesh=",I0)') &
         SIZE(upf%r), upf%mesh
    FLUSH(stdout)
    IF (SIZE(upf%r) == 0) THEN
       WRITE(stdout,'("FATAL read_pp_mesh: r array is empty")')
       ERROR STOP 'read_pp_mesh: empty r array'
    END IF
    allocate ( upf%rab(1:upf%mesh) )
    CALL capitalize_if_v2_into('pp_rab', v2, tag_buf)
    CALL xmlr_readtag( TRIM(tag_buf), upf%rab(1:upf%mesh) )
    ! Debug: verify rab array
    WRITE(stdout,'("DEBUG read_pp_mesh: rab array size=",I0,", mesh=",I0)') &
         SIZE(upf%rab), upf%mesh
    FLUSH(stdout)
    IF (SIZE(upf%rab) == 0) THEN
       WRITE(stdout,'("FATAL read_pp_mesh: rab array is empty")')
       ERROR STOP 'read_pp_mesh: empty rab array'
    END IF
    !
    CALL xmlr_closetag( ) ! end pp_mesh
    !
  END SUBROUTINE read_pp_mesh
  !
  !--------------------------------------------------------
  SUBROUTINE read_pp_semilocal ( upf, ierr )
    !--------------------------------------------------------
    !
    IMPLICIT NONE
    TYPE(pseudo_upf),INTENT(INOUT) :: upf ! the pseudo data
    INTEGER, INTENT(INOUT) :: ierr
    !
    INTEGER :: nb, ind, l, j
    CHARACTER(LEN=9) :: tag
    CHARACTER(LEN=64) :: tag_buf
    CHARACTER(LEN=32) :: nbuf  ! For MinGW-safe tag construction
    real(dp), allocatable :: vnl(:)
    !
    IF ( upf%typ == "SL" ) THEN
       !
       IF ( upf%has_so ) then
          ALLOCATE(upf%vnl(upf%mesh,0:upf%lmax,2))
       else
          ALLOCATE(upf%vnl(upf%mesh,0:upf%lmax,1))
       end if
       allocate ( vnl(1:upf%mesh) )
       CALL capitalize_if_v2_into('pp_semilocal', v2, tag_buf)
       CALL xmlr_opentag( TRIM(tag_buf) )       
       !
       tag = 'vnl'
       DO nb = 1,upf%nbeta
          IF ( v2 ) THEN
             ! NOTA BENE: v2 format follows available PP files, written 
             ! using original write_upf_v2; not FoX-based write_upf_v2
             IF ( nb - 1 == upf%lloc ) CYCLE
             ! MinGW-safe tag construction
             WRITE(nbuf,'(I0)') nb-1
             tag = 'PP_VNL.' // TRIM(nbuf)
          END IF
          CALL xmlr_readtag( tag, vnl, ierr )
          if ( ierr /= 0 ) then
             WRITE(stdout,'("read_pp_semiloca: error reading SL PPs")')
             return
          end if
          CALL get_attr ( 'l', l)
          ind = 1
          IF ( upf%has_so ) then
             CALL get_attr ( 'j', j)
             IF ( l > 0 .AND. ABS(j-l-0.5_dp) < 0.001_dp ) ind = 2
             ! FIXME: what about spin-orbit case for v.2 upf?
             if ( v2 ) then
                WRITE(stdout,'("read_pp_semilocal: check spinorbit case")')
                ierr = 1
                return
             end if
          END IF
          upf%vnl(:,l,ind) = vnl(:)
       END DO
       deallocate ( vnl )
       !
       CALL xmlr_closetag( ) ! end pp_semilocal
       !
    END IF
    !
  END SUBROUTINE read_pp_semilocal
  !
  !--------------------------------------------------------
  SUBROUTINE read_pp_nonlocal ( upf, ierr )
    !--------------------------------------------------------
    !
    IMPLICIT NONE
    TYPE(pseudo_upf),INTENT(INOUT) :: upf ! the pseudo data
    INTEGER, INTENT(INOUT) :: ierr
    !
    LOGICAL :: isnull
    INTEGER :: nb, ind, l, l_, ln, lm, mb, nmb
    CHARACTER(LEN=15) :: tag
    CHARACTER(LEN=64) :: tag_buf
    CHARACTER(LEN=32) :: nbuf, nbuf1, nbuf2, nbuf3  ! For MinGW-safe tag construction
    REAL(dp), ALLOCATABLE :: aux(:)
    !
    nb = upf%nbeta
    IF ( nb == 0 ) nb = 1
    ALLOCATE (upf%beta(upf%mesh,nb) )
    ALLOCATE (upf%els_beta(nb), &
              upf%lll(nb),      &
              upf%kbeta(nb),    &
              upf%rcut(nb),     &
              upf%rcutus(nb),   &
              upf%dion(nb,nb),  &
              upf%qqq(nb,nb)    )
    !
    IF (upf%has_so) ALLOCATE( upf%jjj(upf%nbeta)) 
    !
    IF ( upf%nbeta == 0 ) THEN
       upf%nqf = 0
       upf%nqlc= 0
       upf%kkbeta = 0
       upf%qqq_eps=-1.0_dp
       RETURN
    END IF
    !
    CALL capitalize_if_v2_into('pp_nonlocal', v2, tag_buf)
    CALL xmlr_opentag( TRIM(tag_buf) )
    !
    DO nb = 1,upf%nbeta
       !
       IF ( v2 ) THEN
          ! MinGW-safe tag construction: use direct WRITE instead of i2c()
          ! This avoids allocatable character issues on MinGW gfortran
          WRITE(nbuf,'(I0)') nb
          tag = 'PP_BETA.' // TRIM(nbuf)
          ! Safety check: ensure tag is non-empty
          IF (LEN_TRIM(tag) == 0) THEN
             WRITE(stdout,'("FATAL read_pp_nonlocal: nb=",I0,", constructed tag is empty")') nb
             ERROR STOP 'read_pp_nonlocal: empty tag after construction'
          END IF
          ! Optional debug output (controlled by QE_UPF_DEBUG env var)
          ! Note: GET_ENVIRONMENT_VARIABLE is Fortran 2003, use simpler check
          ! For now, always print debug info if needed (can be disabled later)
          ! IF (LEN_TRIM(GET_ENVIRONMENT_VARIABLE('QE_UPF_DEBUG')) > 0) THEN
          !    WRITE(stdout,'("DEBUG read_pp_nonlocal: nb=",I0,", constructed tag=[",A,"]")') &
          !         nb, trim(tag)
          !    FLUSH(stdout)
          ! END IF
       ELSE
          tag = 'pp_beta'
       END IF
       ! Debug: print tag being read
       WRITE(stdout,'("DEBUG read_pp_nonlocal: nb=",I0,", reading tag=[",A,"]")') nb, trim(tag)
       FLUSH(stdout)
       CALL xmlr_readtag( tag, upf%beta(1:upf%mesh,nb) )
       ! Debug: print array size after reading
       WRITE(stdout,'("DEBUG read_pp_nonlocal: nb=",I0,", beta array size=",I0)') &
            nb, SIZE(upf%beta(1:upf%mesh,nb))
       FLUSH(stdout)
       ! Note: attributes should be available in attrlist after xmlr_readtag
       ! get_attr will read from attrlist which was populated by xmlr_opentag (called inside xmlr_readtag)
       ! Robust check: ensure array is non-empty
       IF (SIZE(upf%beta(1:upf%mesh,nb)) == 0) THEN
          WRITE(stdout,'("FATAL read_pp_nonlocal: nb=",I0,", tag=[",A,"] returned empty array")') &
               nb, trim(tag)
          WRITE(stdout,'("  mesh=",I0,", upf%mesh=",I0)') upf%mesh, upf%mesh
          ERROR STOP 'read_pp_nonlocal: empty beta array'
       END IF
       CALL get_attr('index', mb)
       ! not-so-strict test: index is absent or incorrect in some UPF v.2 files
       IF ( .NOT. v2 .AND. nb /= mb ) then
          write(stdout,'("read_pp_nonlocal: mismatch")')
          ierr = nb
          return
       end if
       CALL get_attr('label', upf%els_beta(nb))
       CALL get_attr('angular_momentum', upf%lll(nb))
       ! Debug: print parsed angular_momentum
       WRITE(stdout,'("DEBUG read_pp_nonlocal: nb=",I0,", angular_momentum=",I0)') &
            nb, upf%lll(nb)
       FLUSH(stdout)
       ! Robust check: ensure angular_momentum is valid (0 <= l <= lmax)
       IF (upf%lll(nb) < 0 .OR. upf%lll(nb) > upf%lmax) THEN
          WRITE(stdout,'("FATAL read_pp_nonlocal: nb=",I0,", invalid angular_momentum=",I0,", lmax=",I0)') &
               nb, upf%lll(nb), upf%lmax
          ERROR STOP 'read_pp_nonlocal: invalid angular_momentum'
       END IF
       IF ( .NOT. v2 .AND. upf%has_so ) &
            CALL get_attr('tot_ang_mom', upf%jjj(nb))
       CALL get_attr('cutoff_radius_index', upf%kbeta(nb))
       CALL get_attr('cutoff_radius', upf%rcut(nb))
       CALL get_attr('ultrasoft_cutoff_radius', upf%rcutus(nb))
       !
       ! Old version of UPF PPs v.2 contained an error in the tag.
       ! To be able to read the old PPs we need the following
       ! Copied from read_upf_v2.f90 :: read_upf_nonlocal
       IF ( upf%rcutus(nb) .EQ. 0._DP ) &
          CALL get_attr('norm_conserving_radius', upf%rcutus(nb))
       !
    END DO
    !
    ! Debug: verify all beta functions were read correctly
    WRITE(stdout,'("DEBUG read_pp_nonlocal: Finished reading ",I0," beta functions")') upf%nbeta
    DO nb = 1, upf%nbeta
       WRITE(stdout,'("  beta(",I0,"): l=",I0,", size=",I0,", kbeta=",I0)') &
            nb, upf%lll(nb), SIZE(upf%beta(1:upf%mesh,nb)), upf%kbeta(nb)
       IF (SIZE(upf%beta(1:upf%mesh,nb)) == 0) THEN
          WRITE(stdout,'("FATAL: beta(",I0,") array is empty!")') nb
          ERROR STOP 'read_pp_nonlocal: empty beta array after reading'
       END IF
    END DO
    FLUSH(stdout)
    !
    ! pp_dij (D_lm matrix)
    !
    CALL capitalize_if_v2_into('pp_dij', v2, tag_buf)
    CALL xmlr_readtag ( TRIM(tag_buf), upf%dion )
    !
    ! pp_augmentation
    !
    IF (upf%tvanp .or. upf%tpawp) THEN
       CALL capitalize_if_v2_into('pp_augmentation', v2, tag_buf)
       CALL xmlr_opentag( TRIM(tag_buf) )
       !
       IF ( v2 ) THEN
          CALL get_attr ( 'q_with_l', upf%q_with_l )
          CALL get_attr ( 'nqf', upf%nqf )
          CALL get_attr ( 'nqlc', upf%nqlc )
          IF (upf%tpawp) THEN
             CALL get_attr ( 'shape', upf%paw%augshape )
             CALL get_attr ( 'cutoff_r', upf%paw%raug )
             CALL get_attr ( 'cutoff_r_index', upf%paw%iraug )
             CALL get_attr ( 'augmentation_epsilon', upf%qqq_eps )
             CALL get_attr ( 'l_max_aug', upf%paw%lmax_aug )
          ENDIF
       ELSE
          CALL xmlr_readtag( 'q_with_l', upf%q_with_l )
          CALL xmlr_readtag( 'nqf', upf%nqf )
          CALL xmlr_readtag( 'nqlc', upf%nqlc )
          IF (upf%tpawp) THEN
             CALL xmlr_readtag( 'shape', upf%paw%augshape )
             CALL xmlr_readtag( 'cutoff_r', upf%paw%raug )
             CALL xmlr_readtag( 'cutoff_r_index', upf%paw%iraug )
             CALL xmlr_readtag( 'augmentation_epsilon', upf%qqq_eps )
             CALL xmlr_readtag( 'l_max_aug', upf%paw%lmax_aug )
          ENDIF
       ENDIF
       !
       CALL capitalize_if_v2_into('pp_q', v2, tag_buf)
       CALL xmlr_readtag( TRIM(tag_buf), upf%qqq )
       !
       IF ( upf%tpawp ) THEN
          ALLOCATE ( upf%paw%augmom(1:upf%nbeta,1:upf%nbeta,0:2*upf%lmax) )
          CALL capitalize_if_v2_into('pp_multipoles', v2, tag_buf)
          CALL xmlr_readtag( TRIM(tag_buf), upf%paw%augmom )
       ENDIF
       !
       ! read polinomial coefficients for Q_ij expansion at small radius
       !
       IF ( upf%nqlc == 0 ) upf%nqlc = 2*upf%lmax+1
       ALLOCATE( upf%rinner( upf%nqlc ) )
       IF ( v2 .AND. upf%nqf > 0) THEN
          ALLOCATE ( upf%qfcoef(upf%nqf, upf%nqlc, upf%nbeta, upf%nbeta) )
          CALL xmlr_opentag('PP_QFCOEF')
          READ(iun,*) upf%qfcoef
          CALL xmlr_closetag ()
          CALL xmlr_readtag('PP_RINNER',upf%rinner)
       ELSE IF ( upf%nqf == 0 ) THEN
          ALLOCATE( upf%qfcoef(1,1,1,1) )
          upf%qfcoef =0.0_dp
       ENDIF
       !
       ! Read augmentation charge Q_ij
       !
       IF( upf%q_with_l ) THEN
          ALLOCATE( upf%qfuncl(upf%mesh,upf%nbeta*(upf%nbeta+1)/2,0:2*upf%lmax) )
          upf%qfuncl(:,:,:) = 0.0_dp
          ! NOTE: it would be wiser to dimension qfuncl as (:,:,0:upf%lmax)
          ! and store the q_l(r) with index l=L/2 (see loop_on_l below)
          ! This would save some storage and avoid "holes" in the array
          ! that may be a source of trouble if not initialized to zero 
       ELSE
          ALLOCATE ( upf%qfunc(upf%mesh,upf%nbeta*(upf%nbeta+1)/2) )
          upf%qfunc (:,:) = 0.0_dp
       END IF
       ALLOCATE ( aux(upf%mesh) )
       loop_on_nb: DO nb = 1,upf%nbeta
          ln = upf%lll(nb)
          loop_on_mb: DO mb = nb,upf%nbeta
             lm = upf%lll(mb)
             IF( upf%q_with_l ) THEN
                loop_on_l: DO l = abs(ln-lm),ln+lm,2 ! only even terms
                   isnull = .FALSE. 
                   IF( upf%tpawp ) isnull = (abs(upf%paw%augmom(nb,mb,l)) < upf%qqq_eps)
                   IF(isnull) CYCLE loop_on_l
                   IF ( v2 ) THEN
                      ! MinGW-safe tag construction
                      WRITE(nbuf1,'(I0)') nb
                      WRITE(nbuf2,'(I0)') mb
                      WRITE(nbuf3,'(I0)') l
                      tag = 'PP_QIJL.' // TRIM(nbuf1) // '.' // TRIM(nbuf2) // '.' // TRIM(nbuf3)
                   ELSE
                      tag = 'pp_qijl'
                   END IF
                   CALL xmlr_readtag( tag, aux )
                   CALL get_attr ('composite_index', nmb)
                   IF ( nmb /= mb*(mb-1)/2 + nb ) then
                      write(stdout,'("read_pp_nonlocal: mismatch")')
                      ierr = 1
                      return
                   end if
                   CALL get_attr ('angular_momentum', l_)
                   IF ( l /= l_ ) then
                      write(stdout,'("read_pp_nonlocal: mismatch")')
                      ierr = 2
                      return
                   end if
                   upf%qfuncl(:,nmb,l) = aux(:)
                   IF (upf%tpawp) upf%qfuncl(upf%paw%iraug+1:,nmb,l) = 0._DP
                ENDDO loop_on_l
             ELSE
                isnull = .FALSE. 
                IF  ( upf%tpawp ) isnull = ( abs(upf%qqq(nb,mb)) < upf%qqq_eps )
                IF (isnull) CYCLE loop_on_mb
                IF ( v2 ) THEN
                   ! MinGW-safe tag construction
                   WRITE(nbuf1,'(I0)') nb
                   WRITE(nbuf2,'(I0)') mb
                   tag = 'PP_QIJ.' // TRIM(nbuf1) // '.' // TRIM(nbuf2)
                ELSE
                   tag = 'pp_qij'
                END IF
                CALL xmlr_readtag( tag, aux )
                CALL get_attr ('composite_index', nmb)
                IF ( nmb /= mb*(mb-1)/2 + nb ) then
                   write(stdout,'("read_pp_nonlocal: mismatch")')
                   ierr = 3
                   return
                end if
                upf%qfunc(:,nmb) = aux(:)
                !
             ENDIF
          ENDDO loop_on_mb
       ENDDO  loop_on_nb
       !
       DEALLOCATE (aux)
       CALL xmlr_closetag( ) ! end pp_augmentation
       !
    END IF
    CALL xmlr_closetag( ) ! end pp_nonlocal
    !
    ! Maximum radius of beta projector: outer radius to integrate
    upf%kkbeta = MAXVAL(upf%kbeta(1:upf%nbeta))
    ! For PAW, augmentation charge may extend a bit further:
    IF(upf%tpawp) upf%kkbeta = MAX(upf%kkbeta, upf%paw%iraug)
    !
  END SUBROUTINE read_pp_nonlocal
  !
  !--------------------------------------------------------
  SUBROUTINE read_pp_pswfc ( upf, ierr )
    !--------------------------------------------------------
    !
    IMPLICIT NONE
    TYPE(pseudo_upf),INTENT(INOUT) :: upf ! the pseudo data
    INTEGER, INTENT(INOUT) :: ierr
    !
    INTEGER :: nw, ind, l
    CHARACTER(LEN=9) :: tag
    CHARACTER(LEN=64) :: tag_buf
    CHARACTER(LEN=32) :: nbuf  ! For MinGW-safe tag construction
    !
    allocate ( upf%chi(1:upf%mesh,upf%nwfc) )
    allocate ( upf%els(upf%nwfc), &
                upf%oc(upf%nwfc), &
                upf%lchi(upf%nwfc), &
                upf%nchi(upf%nwfc), &
                upf%rcut_chi(upf%nwfc), &
                upf%rcutus_chi(upf%nwfc), &
                upf%epseu(upf%nwfc) )
    IF ( upf%has_so ) allocate ( upf%jchi(upf%nwfc) )
    !
    CALL capitalize_if_v2_into('pp_pswfc', v2, tag_buf)
    CALL xmlr_opentag( TRIM(tag_buf) )
    DO nw=1,upf%nwfc
       IF ( v2 ) THEN
          ! MinGW-safe tag construction: use direct WRITE instead of i2c()
          WRITE(nbuf,'(I0)') nw
          tag = 'PP_CHI.' // TRIM(nbuf)
          ! Safety check: ensure tag is non-empty
          IF (LEN_TRIM(tag) == 0) THEN
             WRITE(stdout,'("FATAL read_pp_pswfc: nw=",I0,", constructed tag is empty")') nw
             ERROR STOP 'read_pp_pswfc: empty tag after construction'
          END IF
          ! Optional debug output (controlled by QE_UPF_DEBUG env var)
          ! Note: GET_ENVIRONMENT_VARIABLE is Fortran 2003, use simpler check
          ! For now, always print debug info if needed (can be disabled later)
          ! IF (LEN_TRIM(GET_ENVIRONMENT_VARIABLE('QE_UPF_DEBUG')) > 0) THEN
          !    WRITE(stdout,'("DEBUG read_pp_pswfc: nw=",I0,", constructed tag=[",A,"]")') &
          !         nw, trim(tag)
          !    FLUSH(stdout)
          ! END IF
       ELSE
          tag = 'pp_chi'
       END IF
       ! Debug: print tag being read
       WRITE(stdout,'("DEBUG read_pp_pswfc: nw=",I0,", reading tag=[",A,"]")') nw, trim(tag)
       FLUSH(stdout)
       CALL xmlr_readtag( tag, upf%chi(1:upf%mesh,nw) )
       ! Debug: print array size after reading
       WRITE(stdout,'("DEBUG read_pp_pswfc: nw=",I0,", chi array size=",I0)') &
            nw, SIZE(upf%chi(1:upf%mesh,nw))
       FLUSH(stdout)
       ! Robust check: ensure array is non-empty
       IF (SIZE(upf%chi(1:upf%mesh,nw)) == 0) THEN
          WRITE(stdout,'("FATAL read_pp_pswfc: nw=",I0,", tag=[",A,"] returned empty array")') &
               nw, trim(tag)
          WRITE(stdout,'("  mesh=",I0,", upf%mesh=",I0)') upf%mesh, upf%mesh
          ERROR STOP 'read_pp_pswfc: empty chi array'
       END IF
       call get_attr('index', ind)
       ! not-so-strict test: index is absent or incorrect in some UPF v.2 files
       if ( .NOT. v2 .AND. ind /= nw ) then
          write(stdout,'("read_pp_pswfc: mismatch reading PSWFC")')
          ierr = nw
          return
       end if
       call get_attr( 'label', upf%els(nw) )
       call get_attr( 'l', upf%lchi(nw) )
       ! Debug: print parsed l value
       WRITE(stdout,'("DEBUG read_pp_pswfc: nw=",I0,", l=",I0)') nw, upf%lchi(nw)
       FLUSH(stdout)
       ! Robust check: ensure l is valid (0 <= l <= lmax)
       IF (upf%lchi(nw) < 0 .OR. upf%lchi(nw) > upf%lmax) THEN
          WRITE(stdout,'("FATAL read_pp_pswfc: nw=",I0,", invalid l=",I0,", lmax=",I0)') &
               nw, upf%lchi(nw), upf%lmax
          ERROR STOP 'read_pp_pswfc: invalid l'
       END IF
       IF ( .not. v2 .and. upf%has_so ) call get_attr( 'jchi', upf%jchi(nw) )
       call get_attr( 'occupation', upf%oc(nw) )
       call get_attr( 'n', upf%nchi(nw) )
       call get_attr( 'pseudo_energy', upf%epseu(nw) )
       call get_attr( 'cutoff_radius', upf%rcut_chi(nw) )
       call get_attr( 'ultrasoft_cutoff_radius', upf%rcutus_chi(nw) )
    END DO
    !
    ! Debug: verify all chi functions were read correctly
    WRITE(stdout,'("DEBUG read_pp_pswfc: Finished reading ",I0," chi functions")') upf%nwfc
    DO nw = 1, upf%nwfc
       WRITE(stdout,'("  chi(",I0,"): l=",I0,", size=",I0)') &
            nw, upf%lchi(nw), SIZE(upf%chi(1:upf%mesh,nw))
       IF (SIZE(upf%chi(1:upf%mesh,nw)) == 0) THEN
          WRITE(stdout,'("FATAL: chi(",I0,") array is empty!")') nw
          ERROR STOP 'read_pp_pswfc: empty chi array after reading'
       END IF
    END DO
    FLUSH(stdout)
    !
    CALL xmlr_closetag( ) ! end pp_pswfc
    !
  END SUBROUTINE read_pp_pswfc
  !
  !--------------------------------------------------------
  SUBROUTINE read_pp_full_wfc ( upf, ierr )
    !--------------------------------------------------------
    !
    IMPLICIT NONE
    TYPE(pseudo_upf),INTENT(INOUT) :: upf ! the pseudo data
    INTEGER, INTENT(INOUT) :: ierr
    !
    INTEGER :: nb, mb
    CHARACTER(LEN=15) :: tag
    CHARACTER(LEN=64) :: tag_buf
    CHARACTER(LEN=32) :: nbuf  ! For MinGW-safe tag construction
    !
    IF ( upf%has_wfc ) THEN
       !
       ALLOCATE (upf%aewfc(1:upf%mesh,upf%nbeta) )
       CALL capitalize_if_v2_into('pp_full_wfc', v2, tag_buf)
       CALL xmlr_opentag( TRIM(tag_buf) )
       !
       DO nb = 1, upf%nbeta
          IF ( v2 ) THEN
             ! MinGW-safe tag construction
             WRITE(nbuf,'(I0)') nb
             tag = 'PP_AEWFC.' // TRIM(nbuf)
          ELSE
             tag = 'pp_aewfc'
          END IF
          CALL xmlr_readtag( tag, upf%aewfc(1:upf%mesh,nb) )
          CALL get_attr ('index',mb)
          ! not-so-strict test (and two more below):
          ! index may be absent or incorrect in some UPF v.2 files
          IF ( .NOT. v2 .AND. nb /= mb ) THEN
             WRITE(stdout,'("read_pp_full_wfc: mismatch")')
             ierr = 1
             return
          END IF
       END DO
       !
       IF ( upf%has_so .AND. upf%tpawp ) THEN
          ALLOCATE (upf%paw%aewfc_rel(1:upf%mesh,upf%nbeta) )
          DO nb = 1, upf%nbeta
             IF ( v2 ) THEN
                ! MinGW-safe tag construction
                WRITE(nbuf,'(I0)') nb
                tag = 'PP_AEWFC_REL.' // TRIM(nbuf)
             ELSE
                tag = 'pp_aewfc_rel'
             END IF
             CALL xmlr_readtag(tag, upf%paw%aewfc_rel(1:upf%mesh,nb) )
             CALL get_attr ('index',mb)
             IF ( .NOT. v2 .AND. nb /= mb ) THEN
                WRITE(stdout,'("read_pp_full_wfc: mismatch")')
                ierr = 2
                return
             END IF
          END DO
       END IF
       !
       ALLOCATE (upf%pswfc(1:upf%mesh,upf%nbeta) )
       DO nb = 1, upf%nbeta
          IF ( v2 ) THEN
             ! MinGW-safe tag construction
             WRITE(nbuf,'(I0)') nb
             tag = 'PP_PSWFC.' // TRIM(nbuf)
          ELSE
             tag = 'pp_pswfc'
          END IF
          CALL xmlr_readtag(tag, upf%pswfc(1:upf%mesh,nb) )
          CALL get_attr ('index',mb)
          IF ( .NOT. v2 .AND. nb /= mb )  THEN
             WRITE(stdout,'("read_pp_full_wfc: mismatch")')
             ierr = 3
             return
          END IF
       END DO
       !
       CALL xmlr_closetag( )
       !
    END IF
    !
  END SUBROUTINE read_pp_full_wfc
  !
  !--------------------------------------------------------
  SUBROUTINE read_pp_metagga ( upf, ierr )
    !--------------------------------------------------------
    !
    IMPLICIT NONE
    TYPE(pseudo_upf),INTENT(INOUT) :: upf ! the pseudo data
    INTEGER, INTENT(INOUT) :: ierr
    !
    CHARACTER(LEN=64) :: tag_buf
    !
    ierr = 0
    if ( .NOT. upf%with_metagga_info ) RETURN
    !
    allocate ( upf%tau_core(upf%mesh) )
    allocate ( upf%tau_atom(upf%mesh) )
    CALL capitalize_if_v2_into('pp_taumod', v2, tag_buf)
    CALL xmlr_readtag( TRIM(tag_buf), upf%tau_core(:) )
    CALL capitalize_if_v2_into('pp_tauatom', v2, tag_buf)
    CALL xmlr_readtag( TRIM(tag_buf), upf%tau_atom(:) )
    !
  END SUBROUTINE read_pp_metagga
  !
  !--------------------------------------------------------
  SUBROUTINE read_pp_spinorb ( upf, ierr )
    !--------------------------------------------------------
    !
    IMPLICIT NONE
    TYPE(pseudo_upf),INTENT(INOUT) :: upf ! the pseudo data
    INTEGER, INTENT(INOUT) :: ierr
    INTEGER :: nw, nb, nn
    CHARACTER(LEN=1) :: dummy
    CHARACTER(LEN=32) :: nbuf  ! For MinGW-safe tag construction
    !
    IF ( .NOT. v2 .OR. .NOT. upf%has_so ) RETURN
    !
    CALL xmlr_opentag( 'PP_SPIN_ORB' )
    DO nw = 1,upf%nwfc
       ! MinGW-safe tag construction
       WRITE(nbuf,'(I0)') nw
       CALL xmlr_readtag( 'PP_RELWFC.' // TRIM(nbuf), dummy )
       CALL get_attr( 'index' , nb )
       ! not-so-strict test: index absent or incorrect in some UPF v.2 files
       IF ( .NOT. v2 .AND. nb /= nw ) THEN
          WRITE(stdout,'("read_pp_spinor: mismatch")')
          ierr = 1
          return
       end if
       CALL get_attr( 'nn',    nn ) ! obsolete
       CALL get_attr( 'jchi',  upf%jchi(nw) )
       !
       ! the following data is already known and was not read in old versions
       ! of UPF-reading code. upf%oc is actually missing in some UPF files:
       ! reading it here may spoil the value read earlier and break DFT+U
       !
       ! CALL get_attr( 'lchi',  upf%lchi(nw) )
       ! CALL get_attr( 'els',   upf%els(nw) )
       ! CALL get_attr( 'oc',    upf%oc(nw) )
    ENDDO
    !
    DO nb = 1,upf%nbeta
       ! MinGW-safe tag construction
       WRITE(nbuf,'(I0)') nb
       CALL xmlr_readtag( 'PP_RELBETA.' // TRIM(nbuf), dummy, ierr )
       !
       ! existing PP files may have pp_relbeta first, pp_relwfc later,
       ! but also the other way round
       !
       if ( ierr > 0 ) return
       CALL get_attr( 'index' , nw )
       IF ( .NOT.v2 .AND. nb /= nw ) THEN
          WRITE(stdout,'("read_pp_spinorb: mismatch")')
          ierr = 2
       END IF
       CALL get_attr( 'lll',  upf%lll(nb) )
       CALL get_attr( 'jjj',  upf%jjj(nb) )
    ENDDO
    CALL xmlr_closetag () ! end pp_spin_orb
    !
  END SUBROUTINE read_pp_spinorb
  !
  !--------------------------------------------------------
  SUBROUTINE read_pp_paw ( upf, ierr )
    !--------------------------------------------------------
    !
    IMPLICIT NONE
    TYPE(pseudo_upf),INTENT(INOUT) :: upf ! the pseudo data
    INTEGER, INTENT(INOUT) :: ierr
    INTEGER :: nb, mb
    CHARACTER(LEN=64) :: tag_buf
    !
    IF ( .NOT. upf%tpawp ) RETURN
    !
    CALL capitalize_if_v2_into('pp_paw', v2, tag_buf)
    CALL xmlr_opentag( TRIM(tag_buf) )
    CALL get_attr ('paw_data_format', upf%paw_data_format)
    CALL get_attr ('core_energy', upf%paw%core_energy) 
    ! Full occupation (not only > 0 ones)
    ALLOCATE (upf%paw%oc(upf%nbeta) )
    ALLOCATE (upf%paw%ae_rho_atc(upf%mesh) )
    ALLOCATE (upf%paw%ae_vloc(upf%mesh) )
    CALL capitalize_if_v2_into('pp_occupations', v2, tag_buf)
    CALL xmlr_readtag( TRIM(tag_buf), &
         upf%paw%oc(1:upf%nbeta) )
    ! All-electron core charge
    CALL capitalize_if_v2_into('pp_ae_nlcc', v2, tag_buf)
    CALL xmlr_readtag( TRIM(tag_buf), &
         upf%paw%ae_rho_atc(1:upf%mesh) )
    ! All-electron local potential
    CALL capitalize_if_v2_into('pp_ae_vloc', v2, tag_buf)
    CALL xmlr_readtag( TRIM(tag_buf), &
         upf%paw%ae_vloc(1:upf%mesh) )
    CALL xmlr_closetag () ! end pp_paw
    !
    ALLOCATE(upf%paw%pfunc(upf%mesh, upf%nbeta,upf%nbeta) )
    upf%paw%pfunc(:,:,:) = 0._dp
    IF (upf%has_so) THEN
       ALLOCATE(upf%paw%pfunc_rel(upf%mesh, upf%nbeta,upf%nbeta) )
       upf%paw%pfunc_rel(:,:,:) = 0._dp
    ENDIF
    DO nb=1,upf%nbeta
       DO mb=1,nb
          upf%paw%pfunc (1:upf%mesh, nb, mb) = &
               upf%aewfc(1:upf%mesh, nb) * upf%aewfc(1:upf%mesh, mb)
          IF (upf%has_so) THEN
             upf%paw%pfunc_rel (1:upf%paw%iraug, nb, mb) =  &
                  upf%paw%aewfc_rel(1:upf%paw%iraug, nb) *   &
                  upf%paw%aewfc_rel(1:upf%paw%iraug, mb)
!
!    The small component is added to pfunc. pfunc_rel is useful only
!    to add a small magnetic contribution
!
             upf%paw%pfunc (1:upf%paw%iraug, nb, mb) = &
                        upf%paw%pfunc (1:upf%paw%iraug, nb, mb) + &
                        upf%paw%pfunc_rel (1:upf%paw%iraug, nb, mb)
          ENDIF
          upf%paw%pfunc(upf%paw%iraug+1:,nb,mb) = 0._dp
          !
          upf%paw%pfunc (1:upf%mesh, mb, nb) = upf%paw%pfunc (1:upf%mesh, nb, mb)
          IF (upf%has_so) upf%paw%pfunc_rel (1:upf%mesh, mb, nb) =  &
               upf%paw%pfunc_rel (1:upf%mesh, nb, mb)
       ENDDO
    ENDDO
    !
    ! Pseudo wavefunctions (not only the ones for oc > 0)
    ! All-electron wavefunctions
    ALLOCATE(upf%paw%ptfunc(upf%mesh, upf%nbeta,upf%nbeta) )
    upf%paw%ptfunc(:,:,:) = 0._dp
    DO nb=1,upf%nbeta
       DO mb=1,upf%nbeta
          upf%paw%ptfunc (1:upf%mesh, nb, mb) = &
               upf%pswfc(1:upf%mesh, nb) * upf%pswfc(1:upf%mesh, mb)
          upf%paw%ptfunc(upf%paw%iraug+1:,nb,mb) = 0._dp
          !
          upf%paw%ptfunc (1:upf%mesh, mb, nb) = upf%paw%ptfunc (1:upf%mesh, nb, mb)
       ENDDO
    ENDDO
    !
  END SUBROUTINE read_pp_paw
  !--------------------------------------------------------
  SUBROUTINE read_pp_gipaw ( upf, ierr )
    !--------------------------------------------------------
    !
    IMPLICIT NONE
    TYPE(pseudo_upf),INTENT(INOUT) :: upf ! the pseudo data
    INTEGER, INTENT(INOUT) :: ierr
    !
    INTEGER :: nb, mb
    CHARACTER(LEN=24) :: tag
    CHARACTER(LEN=64) :: tag_buf
    CHARACTER(LEN=32) :: nbuf  ! For MinGW-safe tag construction
    !
    IF (.NOT. upf%has_gipaw) RETURN
    !
    CALL capitalize_if_v2_into('pp_gipaw', v2, tag_buf)
    CALL xmlr_opentag( TRIM(tag_buf) )
    CALL get_attr ('gipaw_data_format', upf%gipaw_data_format ) 
    IF ( v2 ) THEN
       CALL xmlr_opentag( 'PP_GIPAW_CORE_ORBITALS', IERR=ierr )
       ! ierr= 0 for case <PP_GIPAW_CORE_ORBITALS>...</PP_GIPAW_CORE_ORBITALS>
       ! ierr=-1 for case <PP_GIPAW_CORE_ORBITALS ... />
       CALL get_attr ('number_of_core_orbitals', upf%gipaw_ncore_orbitals)
    ELSE
       CALL xmlr_readtag ('number_of_core_orbitals', upf%gipaw_ncore_orbitals) 
       IF ( .NOT. upf%paw_as_gipaw) & 
          CALL xmlr_readtag( 'number_of_valence_orbitals', upf%gipaw_wfs_nchannels)  
    END IF
    ALLOCATE ( upf%gipaw_core_orbital(upf%mesh,upf%gipaw_ncore_orbitals) )
    ALLOCATE ( upf%gipaw_core_orbital_n(upf%gipaw_ncore_orbitals) )
    ALLOCATE ( upf%gipaw_core_orbital_el(upf%gipaw_ncore_orbitals) )
    ALLOCATE ( upf%gipaw_core_orbital_l(upf%gipaw_ncore_orbitals) )
    DO nb = 1,upf%gipaw_ncore_orbitals
       IF ( v2 ) THEN
          ! MinGW-safe tag construction
          WRITE(nbuf,'(I0)') nb
          tag = "PP_GIPAW_CORE_ORBITAL." // TRIM(nbuf)
       ELSE
          tag = 'pp_gipaw_core_orbital'
       END IF
       CALL xmlr_readtag( tag, upf%gipaw_core_orbital(1:upf%mesh,nb) )
       CALL get_attr ('index', mb)
       IF ( nb /= mb ) THEN
          WRITE(stdout,'("read_pp_gipaw: mismatch")')
          ierr = 1
          return
       END IF
       CALL get_attr ('label', upf%gipaw_core_orbital_el(nb) )
       CALL get_attr ('n', upf%gipaw_core_orbital_n(nb) )
       CALL get_attr ('l', upf%gipaw_core_orbital_l(nb) )
    END DO
    ! close only for case <PP_GIPAW_CORE_ORBITALS> ... </PP_GIPAW_CORE_ORBITALS>
    IF ( v2 .AND. ierr == 0  ) CALL xmlr_closetag ( )
    !
    IF ( upf%paw_as_gipaw) THEN
       !
       !    PAW as GIPAW case: all-electron and pseudo-orbitals not read here
       !
       upf%gipaw_wfs_nchannels = upf%nbeta
       ALLOCATE ( upf%gipaw_wfs_el(upf%gipaw_wfs_nchannels) )
       ALLOCATE ( upf%gipaw_wfs_ll(upf%gipaw_wfs_nchannels) )
       ALLOCATE ( upf%gipaw_wfs_rcut(upf%gipaw_wfs_nchannels) )
       ALLOCATE ( upf%gipaw_wfs_rcutus(upf%gipaw_wfs_nchannels) )
       ALLOCATE ( upf%gipaw_wfs_ae(upf%mesh,upf%gipaw_wfs_nchannels) )
       ALLOCATE ( upf%gipaw_wfs_ps(upf%mesh,upf%gipaw_wfs_nchannels) )
       DO nb = 1,upf%gipaw_wfs_nchannels
          upf%gipaw_wfs_el(nb) = upf%els_beta(nb)
          upf%gipaw_wfs_ll(nb) = upf%lll(nb)
          upf%gipaw_wfs_ae(:,nb) = upf%aewfc(:,nb)
       ENDDO
       DO nb = 1,upf%gipaw_wfs_nchannels
          upf%gipaw_wfs_ps(:,nb) = upf%pswfc(:,nb) 
       ENDDO
       ALLOCATE ( upf%gipaw_vlocal_ae(upf%mesh) )
       ALLOCATE ( upf%gipaw_vlocal_ps(upf%mesh) )
       upf%gipaw_vlocal_ae(:)= upf%paw%ae_vloc(:)  
       upf%gipaw_vlocal_ps(:)= upf%vloc(:)
       DO nb = 1,upf%gipaw_wfs_nchannels
          upf%gipaw_wfs_rcut(nb)=upf%rcut(nb)
          upf%gipaw_wfs_rcutus(nb)=upf%rcutus(nb)
       ENDDO
       !
    ELSE
       !
       ! Read valence all-electron and pseudo orbitals
       !
       IF ( v2 ) THEN
          CALL xmlr_opentag( 'PP_GIPAW_ORBITALS' )
          CALL get_attr( 'number_of_valence_orbitals', &
               upf%gipaw_wfs_nchannels )
       END IF
       ALLOCATE ( upf%gipaw_wfs_el(upf%gipaw_wfs_nchannels) )
       ALLOCATE ( upf%gipaw_wfs_ll(upf%gipaw_wfs_nchannels) )
       ALLOCATE ( upf%gipaw_wfs_rcut(upf%gipaw_wfs_nchannels) )
       ALLOCATE ( upf%gipaw_wfs_rcutus(upf%gipaw_wfs_nchannels) )
       ALLOCATE ( upf%gipaw_wfs_ae(upf%mesh,upf%gipaw_wfs_nchannels) )
       ALLOCATE ( upf%gipaw_wfs_ps(upf%mesh,upf%gipaw_wfs_nchannels) )
       DO nb = 1,upf%gipaw_wfs_nchannels
          IF ( v2 ) THEN
             ! MinGW-safe tag construction
             WRITE(nbuf,'(I0)') nb
             tag = "PP_GIPAW_ORBITAL." // TRIM(nbuf)
          ELSE
             tag = 'pp_gipaw_orbital'
          END IF
          CALL xmlr_opentag( tag )
          CALL get_attr ('index', mb)
          IF ( nb /= mb ) THEN
             WRITE(stdout,'("read_pp_gipaw: mismatch")')
             ierr = 2
             return
          end if
          CALL get_attr ('label', upf%gipaw_wfs_el(nb) )
          CALL get_attr ('l',     upf%gipaw_wfs_ll(nb) )
          CALL get_attr ('cutoff_radius', upf%gipaw_wfs_rcut(nb) )
          CALL get_attr ('ultrasoft_cutoff_radius', upf%gipaw_wfs_rcutus(nb) )
          CALL capitalize_if_v2_into('pp_gipaw_wfs_ae', v2, tag_buf)
          CALL xmlr_readtag( TRIM(tag_buf), &
               upf%gipaw_wfs_ae(1:upf%mesh,nb) )
          CALL capitalize_if_v2_into('pp_gipaw_wfs_ps', v2, tag_buf)
          CALL xmlr_readtag( TRIM(tag_buf),&
               upf%gipaw_wfs_ps(1:upf%mesh,nb) )
          CALL xmlr_closetag ()
       END DO
       IF ( v2 ) CALL xmlr_closetag( )
       !
       ! Read all-electron and pseudo local potentials
       !
       ALLOCATE ( upf%gipaw_vlocal_ae(upf%mesh) )
       ALLOCATE ( upf%gipaw_vlocal_ps(upf%mesh) )
       CALL capitalize_if_v2_into('pp_gipaw_vlocal', v2, tag_buf)
       CALL xmlr_opentag( TRIM(tag_buf) )
       CALL capitalize_if_v2_into('pp_gipaw_vlocal_ae', v2, tag_buf)
       CALL xmlr_readtag( TRIM(tag_buf), &
            upf%gipaw_vlocal_ae(1:upf%mesh) )
       CALL capitalize_if_v2_into('pp_gipaw_vlocal_ps', v2, tag_buf)
       CALL xmlr_readtag( TRIM(tag_buf), &
            upf%gipaw_vlocal_ps(1:upf%mesh) )
       CALL xmlr_closetag ()
    END IF
    CALL xmlr_closetag () ! end pp_gipaw
    !
  END SUBROUTINE read_pp_gipaw
  !
END MODULE read_upf_new_module
