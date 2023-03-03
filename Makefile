CXXFLAGS += -D_DEBUG -g -Wall -std=c++2b 

fms_dual: fms_dual.cpp 

.PHONY: test clean

test: fms_dual
	./fms_dual

clean:
	rm fms_dual
