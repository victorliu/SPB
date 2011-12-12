#include <cmath>
#include <cstdlib>

double max3(const double a, const double b, const double c){
	double ret = a;
	if(b > ret){ ret = b; }
	if(c > ret){ ret = c; }
	return ret;
}

double hypot3(const double x, const double y, const double z){
	double xabs = fabs(x);
	double yabs = fabs(y);
	double zabs = fabs(z);
	double w = xabs;
	if(yabs > w){ w = yabs; }
	if(zabs > w){ w = zabs; }

	if(w == 0.0){
		return (0.0);
	}else{
		xabs /= w;
		yabs /= w;
		zabs /= w;
		return w * sqrt(xabs * xabs + yabs * yabs + zabs * zabs);
	}
}

double frand(){
	return (double)rand() / (double)RAND_MAX - 0.5;
}
