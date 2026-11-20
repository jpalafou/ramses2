git checkout develop
cd bin

make clean
make NDIM=1 HYDRO=1 INIT=SINUS

make clean
make NDIM=2 HYDRO=1 INIT=SINUS

cd ..
git checkout high-order
cd bin

make clean
make NDIM=1 HYDRO=1 INIT=SINUS DOFS_PER_CELL=1 EXEC=ramses_1dof

make clean
make NDIM=1 HYDRO=1 INIT=SINUS DOFS_PER_CELL=2 EXEC=ramses_2dof

make clean
make NDIM=2 HYDRO=1 INIT=SINUS DOFS_PER_CELL=1 EXEC=ramses_1dof

make clean
make NDIM=2 HYDRO=1 INIT=SINUS DOFS_PER_CELL=2 EXEC=ramses_2dof
make clean

cd ..