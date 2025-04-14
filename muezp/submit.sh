#!/bin/bash

if [ ! "$#" = 1 ]; then
    >&2 echo "usage: $(basename "$0") <incoming_energy>"
    exit 1
fi

hep_sub run.sh -argu "$1" %{ProcId} -n 10000 \
    -o run_1.log -e run_2.log
