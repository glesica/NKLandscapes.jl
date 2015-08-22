#! /usr/bin/env bash

# Get the directory where the script is located. Should be the
# repository root directoy.
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

usage() {
    echo "./runtests.sh {kauffman,nowak} [{all,nightlies,releases}]"
}

runtests() {
    echo "Environment: $1"
    vagrant ssh "$2" -c "/vagrant/runtests.sh $1"
}

if [ "$#" -eq 2 ]; then
    if [ "$2" == "nightlies" ] || [ "$2" == "releases" ]; then
        runtests $2
        exit 0
    elif [ "$2" == "all" ]; then
        runtests "nightlies"
        runtests "releases"
        exit 0
    else
        usage
        exit 1
    fi
fi

# Run the tests based on the code in the working directory, not just
# what has been committed.
cd "$DIR/test" && julia -e "include(\"../src/NK.jl\"); include(\"${1}.jl\");"
