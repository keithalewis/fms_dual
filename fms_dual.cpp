// fms_dual.cpp - test dual numbers
#include <cassert>
#include <stdexcept>
#include "fms_dual.h"

using namespace fms;

int main()
{
	try {
		int test_dual = dual<double>::test();
		
		{
			dual d(0);
			//assert(d == inv(d));
		}
	}
	catch (const std::exception& ex) {
		puts(ex.what());
	}

	return 0;
}
