gcc -m64 -O3 -c ./boysGaussJacobi.c
gcc -m64 -O3 -c ./boysGaussJacobiTst.c
gcc -m64 -O3 -o ./boysGaussJacobiTst.exe boysGaussJacobiTst.o boysGaussJacobi.o -lm
gcc -m64 -O3 -shared ./boysGaussJacobi.c -o ./../.build/boysGaussJacobi.lib
.\boysGaussJacobiTst.exe > ./test.out