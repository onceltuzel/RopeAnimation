//
//  Created by Oncel Tuzel.
//  Copyright © 2004 Oncel Tuzel. All rights reserved.
//

#include "RSRope.h"



RSRope::RSRope(int pointCount)
{
	pointCount_ = pointCount;
	Init();
}

RSRope::~RSRope()
{
	delete[] pos_;
	delete[] vel_;
	delete spline_;
}

void RSRope::Init()
{
	pos_ = new COVector3D[pointCount_];
	vel_ = new COVector3D[pointCount_];
	spline_ = new RSCubicBSpline(pointCount_);

	for (int i=0 ; i< pointCount_ ; i++)
	{
		//pos_[i].Set(0, ((pointCount_-1-i)*1.0)/(pointCount_-1), 0);
		pos_[i].Set(((i-1)*1.0)/(pointCount_-1), 0.4, 0);
	}
}


RSCubicBSpline* RSRope::GetSpline()
{
	return spline_;
}

void RSRope::GetPosition(double s, COVector3D& pos)
{
	int i;
	int segment;

	segment = spline_->GetMultipliers(s, multipliers_);

	pos.Set(0,0,0);
	for (i=0 ; i<4 ; i++)
	{
		pos += pos_[segment+i]*multipliers_(i);
	}
}

void RSRope::GetVelocity(double s, COVector3D& vel)
{
	int i;
	int segment;

	segment = spline_->GetMultipliers(s, multipliers_);

	vel.Set(0,0,0);
	for (i=0 ; i<4 ; i++)
	{
		vel += vel_[segment+i]*multipliers_(i);
	}
}

void RSRope::operator=(RSRope& rope)
{
	memcpy(pos_, rope.pos_, sizeof(COVector3D)*pointCount_);
	memcpy(vel_, rope.vel_, sizeof(COVector3D)*pointCount_);
}


COVector3D* RSRope::GetControlPoint(int index)
{
	return &pos_[index];
}
