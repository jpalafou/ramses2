bin/ramses_old1d namelist/sinus1d.nml > sinus1d_old.out
mv sinus1d_old.out jp-workspace/runs/sinus1d_old

bin/ramses_old2d namelist/sinus2d.nml > sinus2d_old.out
mv sinus2d_old.out jp-workspace/runs/sinus2d_old
rm -r jp-workspace/runs/sinus2d_old/output_00001
rm -r jp-workspace/runs/sinus2d_old/output_00002
mv output_00001 jp-workspace/runs/sinus2d_old/
mv output_00002 jp-workspace/runs/sinus2d_old/
