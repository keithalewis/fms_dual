CC = g++
CXXFLAGS += -D_DEBUG -g -Wall -std=c++2b 
LDFLAGE += -lm

fms_dual: fms_dual.cpp 

.PHONY: test clean

test: fms_dual
	./fms_dual

clean:
	rem fms_dual
