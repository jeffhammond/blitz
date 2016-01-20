#! /bin/sh

if test ! -d .git && test ! -f src/globals.cpp ; then
    echo You really need to run this script in the top-level blitz directory
    exit 1
fi

set -x

autoreconf -ivf
