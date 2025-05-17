g++ tests/EKF_GEOS3.cpp src/*.cpp include/*.h -lm -std=c++23 -o bin/main.out
cd bin
./main.out