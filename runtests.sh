#! /usr/bin/env bash

# Get the directory where the script is located. Should be the
# repository root directoy.
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

usage() {
    echo "Usage: ./runtests.sh {unit,fastfunc,fastnets,kauffman,nowak} [{all,nightlies,releases}]"
    echo ""
    echo "  Run the test suite specified by the first argument"
    echo "  in the environment (via Vagrant) specified by the"
    echo "  second argument. Omitting the second argument runs"
    echo "  the tests locally."
    echo ""
    echo "  To run tests on Vagrant the VM must be started before"
    echo "  the tests are run."
}

runtests() {
    echo "Environment: $1"
    vagrant ssh "$2" -c "/vagrant/runtests.sh \"$1\""
}

if [ "$#" -eq 0 ]; then
    usage
    exit 1
fi

if [ "$#" -eq 1 ]; then
    if [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
        usage
    else
        # Run the tests based on the code in the working directory, not just
        # what has been committed.
        cd "$DIR/test" && julia -e "include(\"../src/NKLandscapes.jl\"); include(\"${1}.jl\");"
    fi
fi

if [ "$#" -eq 2 ]; then
    if [ "$2" == "nightlies" ] || [ "$2" == "releases" ]; then
        runtests "$1" "$2"
    elif [ "$2" == "all" ]; then
        runtests "$1" "nightlies"
        runtests "$1" "releases"
    else
        usage
        exit 1
    fi
fi

if [ "$#" -gt 2 ]; then
    usage
    exit 1
fi

exit 0

