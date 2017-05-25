
#ifndef _rsbeziercurve_h
#define _rsbeziercurve_h

#include "COVector3D.h"

class RSBezierCurve {
public:
	RSBezierCurve();


private:
	int pointCount_;
	int degree_;
	
	COVector3D* points_;
	double* coefficients_;
};


#endif