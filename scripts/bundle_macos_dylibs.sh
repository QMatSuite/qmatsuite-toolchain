#!/usr/bin/env bash
set -euo pipefail

STAGE_DIR="${1:?Usage: bundle_macos_dylibs.sh <stage-dir>}"
BIN_DIR="$STAGE_DIR/bin"
LIB_DIR="$STAGE_DIR/lib"

mkdir -p "$LIB_DIR"

if [ ! -d "$BIN_DIR" ]; then
  echo "ERROR: bin directory not found: $BIN_DIR"
  exit 1
fi

shopt -s nullglob
BINARIES=("$BIN_DIR"/*.x)
if [ "${#BINARIES[@]}" -eq 0 ]; then
  echo "ERROR: no executables found in $BIN_DIR"
  exit 1
fi

get_rpaths() {
  local target="$1"
  otool -l "$target" | awk '
    $1=="cmd" && $2=="LC_RPATH" {in_rpath_cmd=1; next}
    in_rpath_cmd && $1=="path" {print $2; in_rpath_cmd=0}
  '
}

expand_special_path() {
  local target="$1"
  local path="$2"
  local target_dir

  target_dir="$(cd "$(dirname "$target")" && pwd)"

  case "$path" in
    @loader_path|@loader_path/*)
      # @loader_path or @loader_path/subdir → resolve relative to the target directory.
      local suffix="${path#@loader_path}"
      suffix="${suffix#/}"
      if [ -n "$suffix" ]; then
        echo "$target_dir/$suffix"
      else
        echo "$target_dir"
      fi
      return 0
      ;;
    @executable_path|@executable_path/*)
      # Good enough for build-time resolution: treat executable dir as target dir.
      local suffix="${path#@executable_path}"
      suffix="${suffix#/}"
      if [ -n "$suffix" ]; then
        echo "$target_dir/$suffix"
      else
        echo "$target_dir"
      fi
      return 0
      ;;
  esac

  case "$path" in
    /*)
      # Absolute path, leave as-is.
      echo "$path"
      ;;
    *)
      # Relative path: interpret as relative to target directory.
      echo "$target_dir/$path"
      ;;
  esac
}

resolve_dep_path() {
  local target="$1"
  local dep="$2"

  # Existing absolute path
  if [[ "$dep" = /* ]] && [ -f "$dep" ]; then
    echo "$dep"
    return 0
  fi

  # @loader_path and @executable_path
  if [[ "$dep" == @loader_path/* || "$dep" == @executable_path/* ]]; then
    local expanded
    expanded="$(expand_special_path "$target" "$dep")"
    if [ -f "$expanded" ]; then
      echo "$expanded"
      return 0
    fi
    return 1
  fi

  # @rpath resolution
  if [[ "$dep" == @rpath/* ]]; then
    local name="${dep#@rpath/}"
    local rpath expanded_rpath candidate

    while IFS= read -r rpath; do
      [ -n "$rpath" ] || continue
      expanded_rpath="$(expand_special_path "$target" "$rpath")"
      candidate="$expanded_rpath/$name"
      if [ -f "$candidate" ]; then
        echo "$candidate"
        return 0
      fi
    done < <(get_rpaths "$target")

    return 1
  fi

  # Relative dependency path
  if [[ "$dep" != /* ]]; then
    local expanded
    expanded="$(expand_special_path "$target" "$dep")"
    if [ -f "$expanded" ]; then
      echo "$expanded"
      return 0
    fi
  fi

  return 1
}

collect_nonsystem_deps() {
  local binary="$1"
  otool -L "$binary" | tail -n +2 | awk '{print $1}' | while read -r dep; do
    case "$dep" in
      /System/*|/usr/lib/*) continue ;;
    esac

    local resolved
    if resolved="$(resolve_dep_path "$binary" "$dep")"; then
      echo "$resolved"
    else
      if [[ "$dep" = /* ]]; then
        echo "WARNING: unresolved absolute dep $dep for $binary" >&2
      elif [[ "$dep" == @rpath/* ]]; then
        echo "WARNING: unresolved @rpath dep $dep for $binary" >&2
      else
        echo "WARNING: unresolved dep $dep for $binary" >&2
      fi
    fi
  done
}

DYLIBS_TO_BUNDLE=$(mktemp)

for exe in "${BINARIES[@]}"; do
  collect_nonsystem_deps "$exe" >> "$DYLIBS_TO_BUNDLE"
done

PREV_COUNT=0
while true; do
  sort -u "$DYLIBS_TO_BUNDLE" -o "$DYLIBS_TO_BUNDLE"
  CURR_COUNT=$(wc -l < "$DYLIBS_TO_BUNDLE" | tr -d ' ')
  if [ "$CURR_COUNT" -eq "$PREV_COUNT" ]; then
    break
  fi
  PREV_COUNT=$CURR_COUNT
  SNAPSHOT=$(mktemp)
  cp "$DYLIBS_TO_BUNDLE" "$SNAPSHOT"
  while IFS= read -r dylib; do
    if [ -f "$dylib" ]; then
      collect_nonsystem_deps "$dylib" >> "$DYLIBS_TO_BUNDLE"
    fi
  done < "$SNAPSHOT"
  rm -f "$SNAPSHOT"
done
sort -u "$DYLIBS_TO_BUNDLE" -o "$DYLIBS_TO_BUNDLE"

echo "=== Non-system dylibs to bundle ==="
cat "$DYLIBS_TO_BUNDLE"
echo "=== Total: $(wc -l < "$DYLIBS_TO_BUNDLE" | tr -d ' ') ==="

while IFS= read -r dylib; do
  if [ -f "$dylib" ]; then
    cp -fL "$dylib" "$LIB_DIR/"
    echo "Bundled: $(basename "$dylib")"
  else
    echo "WARNING: $dylib not found, skipping"
  fi
done < "$DYLIBS_TO_BUNDLE"

rewrite_deps() {
  local target="$1"
  local rel_prefix="$2"
  otool -L "$target" | tail -n +2 | awk '{print $1}' | while read -r dep; do
    case "$dep" in
      /System/*|/usr/lib/*) continue ;;
    esac
    local base
    base=$(basename "$dep")
    if [ -f "$LIB_DIR/$base" ]; then
      local new_path="$rel_prefix/$base"
      if [ "$dep" != "$new_path" ]; then
        install_name_tool -change "$dep" "$new_path" "$target"
      fi
    fi
  done
}

for exe in "${BINARIES[@]}"; do
  rewrite_deps "$exe" "@executable_path/../lib"
  # Remove Homebrew-style rpaths that should no longer be needed after bundling.
  while IFS= read -r rpath; do
    case "$rpath" in
      /opt/homebrew/*|/usr/local/*)
        install_name_tool -delete_rpath "$rpath" "$exe" || true
        ;;
    esac
  done < <(get_rpaths "$exe")
done

LIBS=("$LIB_DIR"/*.dylib)
for dylib in "${LIBS[@]}"; do
  base=$(basename "$dylib")
  install_name_tool -id "@loader_path/$base" "$dylib"
  rewrite_deps "$dylib" "@loader_path"
  # Remove Homebrew-style rpaths from bundled dylibs as well.
  while IFS= read -r rpath; do
    case "$rpath" in
      /opt/homebrew/*|/usr/local/*)
        install_name_tool -delete_rpath "$rpath" "$dylib" || true
        ;;
    esac
  done < <(get_rpaths "$dylib")
done

for dylib in "${LIBS[@]}"; do
  codesign --force --sign - "$dylib"
done
for exe in "${BINARIES[@]}"; do
  codesign --force --sign - "$exe"
done

echo ""
echo "=== Portability verification ==="
FAIL=0
for item in "${BINARIES[@]}" "${LIBS[@]}"; do
  [ -f "$item" ] || continue
  BREW_DEPS=$(otool -L "$item" | grep -c "/opt/homebrew\|/usr/local" || true)
  RPATH_DEPS=$(otool -L "$item" | grep -c "@rpath/" || true)
  BREW_RPATHS=$(otool -l "$item" | grep -c "/opt/homebrew\|/usr/local" || true)

  if [ "$BREW_DEPS" -gt 0 ] || [ "$RPATH_DEPS" -gt 0 ] || [ "$BREW_RPATHS" -gt 0 ]; then
    echo "FAIL: $(basename "$item") still has non-portable linkage:"
    if [ "$BREW_DEPS" -gt 0 ] || [ "$RPATH_DEPS" -gt 0 ]; then
      echo "  otool -L:"
      otool -L "$item" | sed 's/^/    /'
    fi
    if [ "$BREW_RPATHS" -gt 0 ]; then
      echo "  otool -l (LC_RPATH entries with Homebrew-style paths):"
      otool -l "$item" | awk '
        $1=="cmd" && $2=="LC_RPATH" {in_rpath_cmd=1; next}
        in_rpath_cmd && $1=="path" {print "    " $0; in_rpath_cmd=0}
      ' | grep "/opt/homebrew\|/usr/local" || true
    fi
    FAIL=1
  fi
done

if [ "$FAIL" -eq 0 ]; then
  echo "All binaries and dylibs are self-contained ✓"
else
  echo "ERROR: Portability check failed!"
  exit 1
fi

rm -f "$DYLIBS_TO_BUNDLE"

echo ""
echo "=== Bundle summary ==="
echo "Executables: $(ls "$BIN_DIR"/*.x | wc -l | tr -d ' ')"
echo "Bundled dylibs: $(ls "$LIB_DIR"/*.dylib 2>/dev/null | wc -l | tr -d ' ')"
echo "Total bin/ size: $(du -sh "$BIN_DIR" | cut -f1)"
echo "Total lib/ size: $(du -sh "$LIB_DIR" | cut -f1)"
