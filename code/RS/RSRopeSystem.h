//
//  Created by Oncel Tuzel.
//  Copyright © 2004 Oncel Tuzel. All rights reserved.
//


#ifndef _rsropesystem_h
#define _rsropesystem_h

#include "RSRope.h"
#include <vector>

using namespace std;

class RSRopeSystem {
public:
	RSRopeSystem(int pointCount);
	~RSRopeSystem();
	void Init();

	void AddRope(RSRope* rope);
	RSRope* GetRope(int index);

	void Iterate(double seconds);

	COVector3D* GetFixedPoint(); 

private:
	void ConvertToRaw(RSRope* rope);
	void ConvertFromRaw(RSRope* rope);

	void CalculateForces();

	void ApplyGravity();
	void ApplyViscosity();
	void ApplyRopePotential();
	void ApplyCollision();
	void ApplyFixedPointSpring();
	void ApplyInternalRopeSpring();
	void ApplyForce(double s, COVector3D& force);

	void Solve(double seconds);
	void CalculateAcceleration();
	void Integrate(double seconds);

	int pointCount_;
	int segmentCount_;

	COMatrix baseDIntegral_;
	COVector baseIntegral_;
	
	COMatrix massMatrix_;
	COMatrix inverseMassMatrix_;

	COVector forceVector_;
	COVector accVector_;
	COVector velVector_;
	COVector posVector_;

	vector<RSRope*> ropes_;

	COVector3D fixedPoint_;

	RSRope* rope_;

	COVector multipliers_;

	COVector tposVector_[4];
	COVector tvelVector_[4];
	RSRope* trope_;

};


#endif