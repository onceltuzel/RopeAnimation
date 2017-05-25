//
//  Created by Oncel Tuzel.
//  Copyright © 2004 Oncel Tuzel. All rights reserved.
//

#ifndef _copolynomialfunctions_h
#define _copolynomialfunctions_h

#include "COVector.h"

class COPolynomialFunctions {
public:
	COPolynomialFunctions();

	static COVector Multiply(COVector& p1, COVector& p2);
	static COVector Integral(COVector& p);
	static double   Evaluate(COVector& p, double t);

private:	

};


#endif