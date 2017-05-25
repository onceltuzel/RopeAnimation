//
//  Created by Oncel Tuzel.
//  Copyright © 2004 Oncel Tuzel. All rights reserved.
//

#include "RSCubicBSplineBF.h"
#include "COPolynomialFunctions.h"

RSCubicBSplineBF* RSCubicBSplineBF::instance_ = 0;

RSCubicBSplineBF* RSCubicBSplineBF::Instance() {
    if (instance_ == 0) {
        instance_ = new RSCubicBSplineBF;
    }
    return instance_;
}

RSCubicBSplineBF::RSCubicBSplineBF()
{
	Init();
}

void RSCubicBSplineBF::Init()
{
	tvector_.Resize(4, 0);

	baseFunctions_.Resize(4,4,0.0);
	
	baseFunctions_(0, 0) = 1;
	baseFunctions_(1, 0) = 4;
	baseFunctions_(2, 0) = 1;
	baseFunctions_(3, 0) = 0;

	baseFunctions_(0, 1) = -3;
	baseFunctions_(1, 1) = 0;
	baseFunctions_(2, 1) = 3;
	baseFunctions_(3, 1) = 0;

	baseFunctions_(0, 2) = 3;
	baseFunctions_(1, 2) = -6;
	baseFunctions_(2, 2) = 3;
	baseFunctions_(3, 2) = 0;	

	baseFunctions_(0, 3) = -1;
	baseFunctions_(1, 3) = 3;
	baseFunctions_(2, 3) = -3;
	baseFunctions_(3, 3) = 1;

	baseFunctions_/=6;
	
	ConstructIntegrals();
}


void RSCubicBSplineBF::ConstructIntegrals()
{
	int i, j;
	COVector p;
	COVector p1;
	COVector p2;
	COVector pintegral;

	baseIntegral_.Resize(4, 0);
	for (i=0 ; i<4 ; i++)
	{
		p = baseFunctions_.GetRow(i);
		pintegral = COPolynomialFunctions::Integral(p);
		baseIntegral_(i) = COPolynomialFunctions::Evaluate(pintegral,1) - 
								 COPolynomialFunctions::Evaluate(pintegral,0);
	}

	baseDIntegral_.Resize(4, 4, 0);
	for (i=0 ; i<4 ; i++)
	{
		for (j=0 ; j<4 ; j++)
		{
			p1 = baseFunctions_.GetRow(i);
			p2 = baseFunctions_.GetRow(j);

			p = COPolynomialFunctions::Multiply(p1,p2);
			pintegral = COPolynomialFunctions::Integral(p);
			baseDIntegral_(i,j) = COPolynomialFunctions::Evaluate(pintegral,1) - 
										 COPolynomialFunctions::Evaluate(pintegral,0);
		}
	}
}

COVector RSCubicBSplineBF::GetBaseIntegral()
{
	return baseIntegral_;
}

COMatrix RSCubicBSplineBF::GetBaseDIntegral()
{
	return baseDIntegral_;
}

COMatrix RSCubicBSplineBF::GetBaseFunctions()
{
	return baseFunctions_;
}

COVector RSCubicBSplineBF::GetCoefficients(double s)
{
	tvector_(0) = 1;
	tvector_(1) = s;
	tvector_(2) = s*s;
	tvector_(3) = s*s*s;

	return baseFunctions_*tvector_;
}

