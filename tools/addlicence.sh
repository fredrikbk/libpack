#!/bin/sh

cd `git rev-parse --show-toplevel`

c_source_files=`git ls-files | grep -v ddtbench | grep -v hrtimer | grep '\(\.c\|\.cpp\|\.h\|\.hpp\)$'`
fortran_source_files=`git ls-files | grep -v ddtbench | grep -v hrtimer | grep '\(\.f\|\.f77\|\.f90\|\.F\|\.F77\|\.F90\)$'`

# TODO this only works if the first line is a good identifier
copyright_line=`head -1 tools/c_license_header`

for f in `grep -L "$copyright_line" $c_source_files`
do
    echo "Adding license header to $f"
    dest="$f.with_header"
    cat tools/c_license_header $f > $dest
    mv $dest $f
done

# TODO this only works if the first line is a good identifier
copyright_line=`head -1 tools/f77_license_header`

for f in `grep -L "$copyright_line" $fortran_source_files`
do
    echo "Adding license header to $f"
    dest="$f.with_header"
    cat tools/f77_license_header $f > $dest
    mv $dest $f
done


