#! /usr/bin/env bash

# Get the directory where the script is located. Should be the
# repository root directoy.
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
SUITE=""

usage() {
    echo "Usage: ./runtests.sh [unit|fastfunc|kauffman|nowak]"
    echo ""
    echo "  Run the test suite specified by the provided argument."
    echo "  Omitting the argument is the same as specifying 'unit'."
}

runtests() {
    echo "Running test suite: $1"
    # Run the tests based on the code in the working directory, not just what
    # has been committed.
    cd "$DIR/test" && julia --color=yes --check-bounds=yes -e "include(\"../src/NKLandscapes.jl\"); include(\"${1}.jl\");"
}

if [ "$#" -eq 0 ]; then
    SUITE="unit"
else
    SUITE="$1"
fi

if [ "$#" -eq 1 ]; then
    if [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
        usage
        exit 0
    fi
fi

if [ "$#" -gt 1 ]; then
    usage
    exit 1
fi

runtests "$SUITE"

exit 0

