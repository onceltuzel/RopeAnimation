//
//  Created by Oncel Tuzel.
//  Copyright © 2004 Oncel Tuzel. All rights reserved.
//

#include "RSRopeSystem.h"
#include "RSCubicBSplineBF.h"

#define ITERATION_COUNT 100
#define VISCOSITY_CONSTANT -0.1

#define FIXED_POINT_S_CONSTANT 10000.0
#define FIXED_POINT_D_CONSTANT 10.0

#define INTERNAL_ROPE_S_CONSTANT 100000.0
#define INTERNAL_ROPE_D_CONSTANT 100.0

#define COLLISION_S_CONSTANT 1000.0
#define COLLISION_D_CONSTANT 200.0
#define FRICTION_CONSTANT -0.1


#define THRESHOLD 0.001

RSRopeSystem::RSRopeSystem(int pointCount)
{
	pointCount_ = pointCount;
	segmentCount_ = pointCount_ - 3;

	Init();
	trope_ = new RSRope(pointCount_);
}

RSRopeSystem::~RSRopeSystem()
{
	delete trope_;
}


void RSRopeSystem::Init()
{
	int i;
	baseDIntegral_ = RSCubicBSplineBF::Instance()->GetBaseDIntegral();
	baseIntegral_ = RSCubicBSplineBF::Instance()->GetBaseIntegral();

	COMatrix massCoorMatrix_;
	massCoorMatrix_.Resize(pointCount_, pointCount_, 0.0);

	COMatrix tmpMassCoordMatrix;
	for (i=0 ; i<segmentCount_ ; i++)
	{
		tmpMassCoordMatrix.Resize(pointCount_, pointCount_, 0.0);
		tmpMassCoordMatrix.SetSubMatrix(i,i,baseDIntegral_);
		massCoorMatrix_ += tmpMassCoordMatrix;
	}

	massMatrix_.Resize(pointCount_*3, pointCount_*3, 0.0);
	massMatrix_.SetSubMatrix(0,0,massCoorMatrix_);
	massMatrix_.SetSubMatrix(pointCount_,pointCount_,massCoorMatrix_);
	massMatrix_.SetSubMatrix(2*pointCount_,2*pointCount_,massCoorMatrix_);

	inverseMassMatrix_ = massMatrix_.Inverse();

	forceVector_.Resize(pointCount_*3, 0);
	accVector_.Resize(pointCount_*3, 0);
	velVector_.Resize(pointCount_*3, 0);
	posVector_.Resize(pointCount_*3, 0);

	fixedPoint_.Set(0.0,0.4,0);

	for (i=0 ; i<11 ; i++)
	{
		for (int j=0 ; j<11 ; j++)
		{
			printf("%.4f ", massMatrix_(i,j));
		}
		printf("\n");
	}
}

void RSRopeSystem::AddRope(RSRope* rope)
{
	ropes_.push_back(rope);
}

RSRope* RSRopeSystem::GetRope(int index)
{
	if (ropes_.size() > index)
		return ropes_[index];
	else
		return NULL;
}

COVector3D* RSRopeSystem::GetFixedPoint()
{ 
	return &fixedPoint_;
}
void RSRopeSystem::Iterate(double seconds)
{

	int i, j;
	
	if (seconds>0.05)
		return;

	for (i=0 ; i<ITERATION_COUNT ; i++)
	{
		for (j=0 ; j<ropes_.size() ; j++)
		{
			rope_ = ropes_[j];
			trope_ = rope_;
			Solve(seconds/ITERATION_COUNT);
		}
	}

	printf("%f\n", seconds);
}

void RSRopeSystem::ConvertToRaw(RSRope* rope)
{
	int i;

	for (i=0 ;i<pointCount_ ; i++)
	{
		posVector_(i) = rope->pos_[i].x;
		posVector_(pointCount_ + i) = rope->pos_[i].y;
		posVector_(2*pointCount_ + i) = rope->pos_[i].z;

		velVector_(i) = rope->vel_[i].x;
		velVector_(pointCount_ + i) = rope->vel_[i].y;
		velVector_(2*pointCount_ + i) = rope->vel_[i].z;
	}	
}

void RSRopeSystem::ConvertFromRaw(RSRope* rope)
{
	int i;

	for (i=0 ;i<pointCount_ ; i++)
	{
		rope->pos_[i].x = posVector_(i);
		rope->pos_[i].y = posVector_(pointCount_ + i);
		rope->pos_[i].z = posVector_(2*pointCount_ + i);

		rope->vel_[i].x = velVector_(i);
		rope->vel_[i].y = velVector_(pointCount_ + i);
		rope->vel_[i].z = velVector_(2*pointCount_ + i);
	}	
}


void RSRopeSystem::CalculateForces()
{
	accVector_.Clear();
	forceVector_.Clear();

	ApplyGravity();
	ApplyFixedPointSpring();
	ApplyInternalRopeSpring();
	ApplyViscosity();
	ApplyCollision();
}


void RSRopeSystem::Solve(double seconds)
{


	ConvertToRaw(rope_);
	CalculateForces();
	CalculateAcceleration();
	Integrate(seconds);		
	ConvertFromRaw(rope_);

/*
	int i;

	ConvertToRaw(rope_);

	for(i=0 ; i< 4 ; i++)
	{
		CalculateForces();
		CalculateAcceleration();
		if (i==3)
			Integrate(seconds);
		else
			Integrate(seconds/2);
		
		ConvertFromRaw(rope_);

		tposVector_[i] = velVector_*seconds;
		tvelVector_[i] = accVector_*seconds;
	}

	ConvertToRaw(trope_);

	posVector_ += (tposVector_[0]*(1.0/6.0)) + 
					  (tposVector_[1]*(1.0/3.0)) + 
					  (tposVector_[2]*(1.0/3.0)) + 
					  (tposVector_[3]*(1.0/6.0));

	velVector_ += (tvelVector_[0]*(1.0/6.0)) + 
					  (tvelVector_[1]*(1.0/3.0)) + 
					  (tvelVector_[2]*(1.0/3.0)) + 
					  (tvelVector_[3]*(1.0/6.0));

	ConvertFromRaw(rope_);
*/
}

void RSRopeSystem::CalculateAcceleration()
{
	accVector_ = inverseMassMatrix_*forceVector_;
}

//Euler integration
void RSRopeSystem::Integrate(double seconds)
{
	velVector_ += accVector_*seconds;
	posVector_ += velVector_*seconds;
}


void RSRopeSystem::ApplyGravity()
{
	int i, j;

	for (i=0 ;i<segmentCount_ ; i++)
	{
		for (j=0 ; j<4 ; j++)
		{
			forceVector_(pointCount_ + i + j) += -baseIntegral_(j);
		}
	}	
}

void RSRopeSystem::ApplyViscosity()
{
	forceVector_ += ((massMatrix_ * velVector_) * VISCOSITY_CONSTANT);
}


void RSRopeSystem::ApplyFixedPointSpring()
{
	COVector3D pos;
	COVector3D force;
	COVector3D dir;
	COVector3D vel;
	double distance;

	rope_->GetPosition(0.0, pos);
	rope_->GetVelocity(0.0, vel);

	dir = fixedPoint_ - pos;
	distance = dir.length();

	if (distance > THRESHOLD)
	{
		COVector3D sf = dir * FIXED_POINT_S_CONSTANT;
		COVector3D sff = vel * FIXED_POINT_D_CONSTANT;
		
		force = sf-sff;
		
		ApplyForce(0.0, force);
	}
}

void RSRopeSystem::ApplyCollision()
{
	int i;
	COVector3D pos, vel, force;
	double s;
	double dist;
	
	for (i=0 ; i<=20 ; i++)
	{
		s = i/20.0f;

		rope_->GetPosition(s, pos);
		rope_->GetVelocity(s, vel);

		if (pos.y < 0.012 && vel.y < 0)
		{
			//Apply collision response
			dist = fabs(pos.y);
			
			COVector3D sf(0,1,0);
			COVector3D sff(0,1,0);
			sf *= (dist * 	COLLISION_S_CONSTANT);
			sff *= (vel.y * COLLISION_D_CONSTANT);

			force = sf-sff;
			ApplyForce(s, force);
		}

		if (pos.y < 0.013)
		{
			//Apply friction
			vel.y = 0;
			vel.unitize();
			force = vel * FRICTION_CONSTANT;
			ApplyForce(s, force);
		}

	}
}

void RSRopeSystem::ApplyInternalRopeSpring()
{
	int i;
	double s0, s1;
	COVector3D pos0, pos1, dir, force0, force1, vel0, vel1, veldif;
	double l0 = 0.8 / 20;
	double distance;

	for (i=0 ; i<20 ; i++)
	{
		s0 = i/20.0;
		s1 = (i+1)/20.0;

		rope_->GetPosition(s0, pos0);
		rope_->GetPosition(s1, pos1);
		dir = pos1 - pos0;
		distance = dir.length();

		if (distance>l0+THRESHOLD)
		{
			dir.unitize();
			rope_->GetVelocity(s0, vel0);
			rope_->GetVelocity(s1, vel1);
			veldif = vel0 - vel1;

			COVector3D sf = dir * ((distance-l0) * 	INTERNAL_ROPE_S_CONSTANT);
			COVector3D sff = veldif * INTERNAL_ROPE_D_CONSTANT;
			//COVector3D sff = dir * (veldif.dot(dir) * INTERNAL_ROPE_D_CONSTANT);

			force0 = sf-sff;
			force1 = force0*-1;


			ApplyForce(s0, force0);
			ApplyForce(s1, force1);
			
		}
		else if (distance<l0-THRESHOLD)
		{
		}
	}

}


void RSRopeSystem::ApplyForce(double s, COVector3D& force)
{
	int segment;

	segment = rope_->GetSpline()->GetMultipliers(s, multipliers_);

	for (int i=0 ; i<4 ;i++)
	{
		forceVector_(segment + i) += force.x*multipliers_(i);
		forceVector_(pointCount_ + segment + i) += force.y*multipliers_(i);
		forceVector_(2*pointCount_ +segment +i) += force.z*multipliers_(i);
	}

	
}