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

collect_nonsystem_deps() {
  local binary="$1"
  otool -L "$binary" | tail -n +2 | awk '{print $1}' | while read -r dep; do
    case "$dep" in
      /System/*|/usr/lib/*|@*) continue ;;
    esac
    echo "$dep"
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
      /System/*|/usr/lib/*|@*) continue ;;
    esac
    local base
    base=$(basename "$dep")
    if [ -f "$LIB_DIR/$base" ]; then
      install_name_tool -change "$dep" "$rel_prefix/$base" "$target"
    fi
  done
}

for exe in "${BINARIES[@]}"; do
  rewrite_deps "$exe" "@executable_path/../lib"
done

LIBS=("$LIB_DIR"/*.dylib)
for dylib in "${LIBS[@]}"; do
  base=$(basename "$dylib")
  install_name_tool -id "@loader_path/$base" "$dylib"
  rewrite_deps "$dylib" "@loader_path"
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
  if [ "$BREW_DEPS" -gt 0 ]; then
    echo "FAIL: $(basename "$item") still links to non-bundled libs:"
    otool -L "$item" | grep "/opt/homebrew\|/usr/local"
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
