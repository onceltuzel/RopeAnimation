//
//  Created by Oncel Tuzel.
//  Copyright © 2004 Oncel Tuzel. All rights reserved.
//

#include "COPolynomialFunctions.h"


COVector COPolynomialFunctions::Multiply(COVector& p1, COVector& p2)
{
	int i, j;
	COVector poly;
	int d1, d2;

	d1 = p1.GetLength();
	d2 = p2.GetLength();

	poly.Resize(d1+d2-1, 0);

	for (i=0 ; i<d1 ; i++)
	{
		for (j=0 ; j<d2 ; j++)
		{
			poly(i+j) += p1(i)*p2(j);
		}
	}

	return poly;
}


double COPolynomialFunctions::Evaluate(COVector& p, double t)
{
	int i;
	int d = p.GetLength();
	double res = 0.0;
	double mt = 1.0;
	for (i=0 ; i<d ; i++)
	{
		res += mt * p(i);
		mt = mt*t;
	}

	return res;
}


COVector COPolynomialFunctions::Integral(COVector& p)
{
	int i;
	COVector poly;
	int d = p.GetLength();

	poly.Resize(d+1, 0);

	for (i=0 ; i<d ; i++)
	{
		poly(i+1) = p(i)/(i+1);
	}

	return poly;
}
