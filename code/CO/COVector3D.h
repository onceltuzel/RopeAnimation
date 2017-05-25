//
//  Created by Oncel Tuzel.
//  Copyright © 2004 Oncel Tuzel. All rights reserved.
//


#ifndef _covector3d_h
#define _covector3d_h

#define COVVector3DThreshold 0.00001

#include <math.h>



class COVector3D {
public:
   double x;
   double y;
   double z;

   double*	 data_[3];


   COVector3D()
   {
	  x = 0;
	  y = 0;
	  z = 0;

	  data_[0] = &x;
	  data_[1] = &y;
	  data_[2] = &z;
   }

   COVector3D(double px, double py, double pz)
   {
	  x = px;
	  y = py;
	  z = pz;

	  data_[0] = &x;
	  data_[1] = &y;
	  data_[2] = &z;
   }

	inline void Set(double px, double py, double pz)
   {
	  x = px;
	  y = py;
	  z = pz;
	}
	
   inline const int operator==(const COVector3D& in)
   {
	  if (in.x == x && in.y == y && in.z == z)
		 return 1;
	  else
		 return 0;
   }


   inline void operator=(const COVector3D& in)
   {
	  x = in.x;
	  y = in.y;
	  z = in.z;   
   }

   inline void operator-=(const COVector3D& in)
   {
	  x = x - in.x;
	  y = y - in.y;
	  z = z - in.z;   
   }

   inline void operator+=(const COVector3D& in)
   {
	  x = x + in.x;
	  y = y + in.y;
	  z = z + in.z;   
   }

   inline void operator/=(double in)
   {
	  x = x/in;
	  y = y/in;
	  z = z/in;
   }

   inline void operator*=(double in)
   {
	  x = x*in;
	  y = y*in;
	  z = z*in;
   }


   inline const COVector3D operator-(const COVector3D& in)
   {
	  COVector3D out;

	  out.x = x - in.x;
	  out.y = y - in.y;
	  out.z = z - in.z;
   
	  return out;
   }

   inline const COVector3D operator+(const COVector3D& in)
   {
	  COVector3D out;

	  out.x = x + in.x;
	  out.y = y + in.y;
	  out.z = z + in.z;
   
	  return out;
   }

   inline const COVector3D operator/(double in)
   {
	  COVector3D out;

	  out.x = x/in;
	  out.y = y/in;
	  out.z = z/in;

	  return out;
   }

   inline const COVector3D operator*(double in)
   {
	  COVector3D out;

	  out.x = x*in;
	  out.y = y*in;
	  out.z = z*in;

	  return out;
   }

   inline const COVector3D operator*(const COVector3D& in)
   {
	  COVector3D out;

	  out.x = y*in.z - z*in.y;
	  out.y = z*in.x - x*in.z;
	  out.z = x*in.y - y*in.x;

	  return out;
   }

   inline const double dot(const COVector3D& in)
   {
	  return x*in.x + y*in.y + z*in.z;
   }

   inline const double distance(const COVector3D& in)
   {
	  COVector3D dist = this->operator-(in);

	  return dist.length();
   }

   inline const double length()
   {
	  return (double) sqrt(x*x + y*y + z*z);
   }


   inline const COVector3D unit()
   {
	  COVector3D out;
	  double total;

	  total = length() + COVVector3DThreshold;

	  out.x = x / total;
	  out.y = y / total;
	  out.z = z / total;

	  return out;
   }

   inline const double unitize()
   {
	  double total;

	  total = length() + COVVector3DThreshold;

	  x = x / total;
	  y = y / total;
	  z = z / total;

	  return total;
   }
   

};
#endif