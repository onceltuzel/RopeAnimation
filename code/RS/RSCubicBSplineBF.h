//
//  Created by Oncel Tuzel.
//  Copyright © 2004 Oncel Tuzel. All rights reserved.
//


#ifndef _rscubicbsplinebf_h
#define _rscubicbsplinebf_h

#include "COVector3D.h"
#include "COMatrix.h"

class RSCubicBSplineBF {
public:
	static RSCubicBSplineBF* Instance();
	~RSCubicBSplineBF();

	COVector GetBaseIntegral();
	COMatrix GetBaseDIntegral();
	COMatrix GetBaseFunctions();
	COVector GetCoefficients(double s);

protected:
	RSCubicBSplineBF();
 
	static RSCubicBSplineBF* instance_;
	
private:
	void Init();
	void ConstructIntegrals();

	COMatrix baseFunctions_;

	COVector baseIntegral_;
	COMatrix baseDIntegral_;

	COVector tvector_;
};

#endif