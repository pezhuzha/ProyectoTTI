g++ tests/main.cpp src/matrix.cpp include/matrix.h src/JPL_Eph_DE430.cpp include/JPL_Eph_DE430.h src/GLOBAL.cpp include/SAT_Const.h include/GLOBAL.h src/Cheb3D.cpp include/Cheb3D.h -lm -std=c++23 -o bin/main.out
cd bin
./main.out