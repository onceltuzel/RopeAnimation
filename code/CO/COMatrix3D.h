//
//  Created by Oncel Tuzel.
//  Copyright © 2004 Oncel Tuzel. All rights reserved.
//


#ifndef _comatrix3d_h
#define _comatrix3d_h

#include <math.h>

class COMatrix3D {
public:
   
   COMatrix3D()
   {
	  int i, j;

	  for (i=0 ; i<3 ; i++)
		 for (j=0 ; j<3 ; j++)
			data_[i][j] = 0;

	  for (i=0 ; i<3 ; i++)
		 data_[i][j] = 1;
   }
   
   
   const COMatrix3D& operator=(const COMatrix3D& in)
   {
	  int i, j;

	  for (i=0 ; i<3 ; i++)
		 for (j=0 ; j<3 ; j++)
			data_[i][j] = in.data_[i][j];
   }

   const COMatrix3D operator*(const COMatrix3D& in)
   {
   }

   inline double operator()(int i, int j) const
   {
	  return data_[i-1][j-1];
   }
   
   inline double& operator()(int i, int j)
   {
	  return data_[i-1][j-1];
   }


   double data_[3][3];

};
#endif