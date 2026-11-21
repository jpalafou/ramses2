bin/ramses1d namelist/sinus1d.nml > ramses1d_level6.out
bin/ramses_1dof1d namelist/sinus1d.nml > ramses_1dof1d_level6.out
bin/ramses_2dof1d namelist/sinus1d.nml > ramses_2dof1d_level6.out

sed -i '' 's/levelmin=6/levelmin=7/g' namelist/sinus1d.nml
sed -i '' 's/levelmax=6/levelmax=7/g' namelist/sinus1d.nml
bin/ramses1d namelist/sinus1d.nml > ramses1d_level7.out
bin/ramses_1dof1d namelist/sinus1d.nml > ramses_1dof1d_level7.out
bin/ramses_2dof1d namelist/sinus1d.nml > ramses_2dof1d_level7.out

sed -i '' 's/levelmin=7/levelmin=8/g' namelist/sinus1d.nml
sed -i '' 's/levelmax=7/levelmax=8/g' namelist/sinus1d.nml
bin/ramses1d namelist/sinus1d.nml > ramses1d_level8.out
bin/ramses_1dof1d namelist/sinus1d.nml > ramses_1dof1d_level8.out
bin/ramses_2dof1d namelist/sinus1d.nml > ramses_2dof1d_level8.out

sed -i '' 's/levelmin=8/levelmin=9/g' namelist/sinus1d.nml
sed -i '' 's/levelmax=8/levelmax=9/g' namelist/sinus1d.nml
bin/ramses1d namelist/sinus1d.nml > ramses1d_level9.out
bin/ramses_1dof1d namelist/sinus1d.nml > ramses_1dof1d_level9.out
bin/ramses_2dof1d namelist/sinus1d.nml > ramses_2dof1d_level9.out

sed -i '' 's/levelmin=9/levelmin=6/g' namelist/sinus1d.nml
sed -i '' 's/levelmax=9/levelmax=6/g' namelist/sinus1d.nml