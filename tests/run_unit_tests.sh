#!/bin/bash
# Unit tests of the closures in src/sqCfdem/common/spheroidGeometry.H, against independent
# references. Needs only a C++14 compiler (no OpenFOAM):   bash tests/run_unit_tests.sh
set -e
here=$(cd "$(dirname "$0")" && pwd); top=$(dirname "$here"); b=$(mktemp -d)
for t in test_geometry test_closures test_place; do
    g++ -O2 -std=c++14 -Wall -Wno-unused-variable -I"$here/unit/shim" -I"$top/src/sqCfdem/common" \
        "$here/unit/$t.cpp" -o "$b/$t"
    echo "== $t"; "$b/$t"
done
python3 "$here/unit/test_case_checks.py"
rm -rf "$b"; echo "ALL UNIT TESTS PASSED"
