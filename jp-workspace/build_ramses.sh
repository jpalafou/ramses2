cd bin

make clean
make NDIM=1 HYDRO=1 INIT=SINUS EXEC=ramses_old

make clean
make NDIM=2 HYDRO=1 INIT=SINUS EXEC=ramses_old

make clean
make NDIM=3 HYDRO=1 INIT=SINUS EXEC=ramses_old
make clean

cd ..