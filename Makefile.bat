g++ tests/EKF_GEOS3.cpp src/*.cpp include/*.h  -lm -std=c++23 -o bin/main.exe
cd bin
main.exe
pause
