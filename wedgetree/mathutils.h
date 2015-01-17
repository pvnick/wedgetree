#ifndef HEADER_E4D805E540984D638A9F803DF28B880B
#define HEADER_E4D805E540984D638A9F803DF28B880B

#include <cmath>

namespace MathUtils {
	inline double fastdist(double x, double y) {
		return (x - y) * (x - y);
	}
}

#endif