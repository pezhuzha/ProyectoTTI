g++ tests/tests.cpp src/*.cpp include/*.h -pipe -lm -std=c++23 -o bin/tests.out
cd bin
./tests.out