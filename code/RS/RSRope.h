//
//  Created by Oncel Tuzel.
//  Copyright © 2004 Oncel Tuzel. All rights reserved.
//


#ifndef _rsrope_h
#define _rsrope_h

#include "RSCubicBSpline.h"

class RSRope {
public:
	RSRope(int pointCount);
	~RSRope();

	RSCubicBSpline* GetSpline();
	void GetPosition(double s, COVector3D& pos);
	void GetVelocity(double s, COVector3D& vel);
	COVector3D* GetControlPoint(int index);
	void operator=(RSRope& rope);
//private:

	int pointCount_;

	COVector3D* pos_;
	COVector3D* vel_;

	RSCubicBSpline* spline_;

	COVector multipliers_;

	void Init();
};


#endif