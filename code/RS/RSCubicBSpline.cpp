//
//  Created by Oncel Tuzel.
//  Copyright © 2004 Oncel Tuzel. All rights reserved.
//



#include "RSCubicBSpline.h"
#include "RSCubicBSplineBF.h"

void RSCubicBSpline::Init()
{
	baseFunctions_ = RSCubicBSplineBF::Instance()->GetBaseFunctions();

	svector_.Resize(4);
	multipliers_.Resize(4);

}

RSCubicBSpline::RSCubicBSpline()
{
	Init();
}

RSCubicBSpline::RSCubicBSpline(int pointCount) 
{
	Init();
	
	pointCount_ = pointCount;
	segmentCount_ = pointCount - 3;
	knotDistance_ = 1.0 / segmentCount_;
}


RSCubicBSpline::~RSCubicBSpline()
{
	
}


int RSCubicBSpline::GetMultipliers(double s, COVector& multipliers)
{
	int segment;
	double rs;

	segment = (int) (s / knotDistance_);

	if (segment>=segmentCount_)
		segment = segmentCount_-1;

	rs = s/knotDistance_ - segment;

	svector_(0) = 1;
	svector_(1) = rs;
	svector_(2) = rs*rs;
	svector_(3) = rs*rs*rs;

	multipliers =	baseFunctions_*svector_;

	return segment;
}
