# Proyecto TTI

# Programa principal
g++ tests/EKF_GEOS3.cpp src/*.cpp include/*.h  -lm -std=c++23 -o bin/main.exe
cd bin
main.exe
pause

# Tests
g++ tests/tests.cpp src/*.cpp include/*.h -lm -std=c++23 -o bin/tests.exe
cd bin
tests.exe
pause