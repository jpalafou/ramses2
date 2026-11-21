#!/bin/sh
set -eu

namelist=$(mktemp)
trap 'rm -f "$namelist"' EXIT

for level in 6 7; do
    sed -e "s/^levelmin=.*/levelmin=$level/" \
        -e "s/^levelmax=.*/levelmax=$level/" \
        namelist/sinus1d.nml > "$namelist"

    for executable in ramses1d ramses_1dof1d ramses_2dof1d; do
        "bin/$executable" "$namelist" > "out/${executable}_level${level}.out"
    done

    sed -i '' "s/scheme='muscl'/scheme='p3'\\
time_integrator='rk4'/" "$namelist"
    bin/ramses_2dof1d "$namelist" > "out/ramses_2dof1d_level${level}_p3.out"
done
