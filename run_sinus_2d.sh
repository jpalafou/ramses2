#!/bin/sh
set -eu

root=$(pwd)
namelist=$(mktemp)
trap 'rm -f "$namelist"' EXIT

run() (
    executable=$1
    output=$2
    mkdir -p "$output"
    cd "$output"
    "$root/bin/$executable" "$namelist" > run.out
    echo "Finished $output"
)

for level in 4 5 6; do
    sed -e "s/^levelmin=.*/levelmin=$level/" \
        -e "s/^levelmax=.*/levelmax=$level/" \
        -e 's/^ngridmax=.*/ngridmax=10000/' \
        namelist/sinus2d.nml > "$namelist"

    for executable in ramses2d ramses_1dof2d ramses_2dof2d; do
        run "$executable" "out/${executable}_level${level}"
    done

    sed -i '' "s/scheme='muscl'/scheme='p3'\\
time_integrator='rk4'/" "$namelist"
    run ramses_2dof2d "out/ramses_2dof2d_p3_level${level}"
done
