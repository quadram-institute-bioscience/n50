#!/bin/bash

VER="$1"
VER_REGEX="${VER//./\\.}"
NEW="$2"

if [[ -z $VER ]]; then
  echo "Usage: $0 <current_version> [<new_version>]"
  echo "  If <new_version> is omitted, the script lists files containing <current_version>."
  exit 1
else
  echo "Current version: $VER"
fi

if [[ -z $NEW ]]; then
  echo "Find mode: looking for files containing $VER"
  find . \( -name "*.md" -o -name "*.c" -o -name "*.sh" -o -name "*.txt" \) -exec grep -l "$VER_REGEX" {} +
else
  echo "Replace mode: replacing $VER with $NEW"
  find . \( -name "*.md" -o -name "*.c" -o -name "*.sh" -o -name "*.txt" \) -exec sed -i '.bakver' "s/$VER_REGEX/$NEW/g" {} +
fi

