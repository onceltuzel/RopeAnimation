//
//  Created by Oncel Tuzel.
//  Copyright © 2004 Oncel Tuzel. All rights reserved.
//


#ifndef _rscubicbspline_h
#define _rscubicbspline_h

#include "COVector3D.h"
#include "COMatrix.h"

class RSCubicBSpline {
public:
	RSCubicBSpline(int pointCount);
	RSCubicBSpline();
	~RSCubicBSpline();

	int GetMultipliers(double s, COVector& multipliers);

private:
	void Init();

	int pointCount_;
	int segmentCount_;
	double knotDistance_;

	COMatrix baseFunctions_;

	COVector svector_;
	COVector multipliers_;

};

#endif