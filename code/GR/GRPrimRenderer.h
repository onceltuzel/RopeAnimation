//
//  Created by Oncel Tuzel.
//  Copyright © 2004 Oncel Tuzel. All rights reserved.
//


#ifndef _grprimrenderer_h
#define _grprimrenderer_h
#include <stdlib.h>
#include "COVector3D.h"
#include "GL/glut.h"

class GRPrimRenderer {
public:
   void drawUnitPlane();
   void drawUnitCube();

   void drawPlane(COVector3D& p1,
			  	  COVector3D& p2,							   			 
				  COVector3D& p3,						 
				  COVector3D& p4);

   void drawPlane(COVector3D& p1,
			  	  COVector3D& p2,							   			 
				  COVector3D& p3,	 
				  COVector3D& p4,
				  COVector3D& n);

private:


};

#endif