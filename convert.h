#ifndef __CONVERT_H__
#define __CONVERT_H__

#include "array.h"

struct coord {
	float x, y, z;
};

struct coordAndSide {
	int side, i, j, k;
};

coord getCenter (int i, int j, int k, int side, const fort::array <1, float> &r, 
		const fort::array<3, float> &x, const fort::array<3, float>&y, const fort::array<3, float> &z ) {
	float rrUp = r (k);
	float rrDown = r(k-1);
	float xUp, yUp, zUp, xDown, yDown, zDown;
	xUp = rrUp * x(i, j, side);
	yUp = rrUp * y(i, j, side);
	zUp = rrUp * z(i, j, side);
	xDown = rrDown * x(i-1, j-1, side);
	yDown = rrDown * y(i-1, j-1, side);
	zDown = rrDown * z(i-1, j-1, side);

	return coord{0.5f*(xUp + xDown), 0.5f*(yUp + yDown), 0.5f*(zUp + zDown)};
}

coordAndSide getCoordAndSide(float x, float y, float z, const fort::array <1, float> &r, int Nr, int n) {
	int side = 0, i, j, k;
	float maxXY, maxYZ, maxZX, radVect, fi, psi;

	maxXY = std::max(std::abs(x), std::abs(y));
	maxYZ = std::max(std::abs(y), std::abs(z));
	maxZX = std::max(std::abs(z), std::abs(x));

	double pi = std::atan(1)*4;

	radVect = std::sqrt(x*x + y*y + z*z);
	k = std::lower_bound (&r(0), &r(Nr) + 1, radVect) - &r(0);
	if (k < 1) 
		k = 1;
	if (k > Nr)
		k = Nr; 

	if (-x >= maxYZ ) {
		side = 1;
		fi = std::atan2(-y,-x);
		psi = std::atan2(z,-x);
	}

	if ( z >= maxXY ) { 
		side = 2;
		fi = std::atan2(-y,z);
		psi = std::atan2(x,z);
	}

	if ( y >= maxZX ) {
		side = 3;
		fi = std::atan2(-x,y);
		psi = std::atan2(z,y);
	}

	if ( -y >= maxZX ) {
		side = 4;
		fi = std::atan2(x,-y);
		psi = std::atan2(z,-y);
	}

	if ( x >= maxYZ ) {
		side = 5;
		fi = std::atan2(y,x);
		psi = std::atan2(z,x);
	}

	if ( -z >= maxXY ) {
		side = 6;
		fi = std::atan2(-y,-z);
		psi = std::atan2(-x,-z);
	}

	i = std::ceil(2*n*fi/pi + n*0.5);
	j = std::ceil(2*n*psi/pi + n*0.5);

	return coordAndSide {side, i, j, k};

}



#endif