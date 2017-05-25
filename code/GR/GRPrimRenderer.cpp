//
//  Created by Oncel Tuzel.
//  Copyright © 2004 Oncel Tuzel. All rights reserved.
//


#include "GRPrimRenderer.h"


void GRPrimRenderer::drawUnitPlane()
{

   drawPlane(COVector3D(0.5, 0.5, 0.0),
			 COVector3D(0.5, -0.5, 0.0),
			 COVector3D(-0.5, -0.5, 0.0),
	         COVector3D(-0.5, 0.5, 0.0)
			 );

}


void GRPrimRenderer::drawUnitCube()
{
   drawPlane(
			 COVector3D(0.5, 0.5, -0.5),
			 COVector3D(0.5, -0.5, -0.5),
			 COVector3D(-0.5, -0.5, -0.5),	         
			 COVector3D(-0.5, 0.5, -0.5),
			 COVector3D(0.0, 0.0, -1.0)
			 );


   drawPlane(
	         COVector3D(-0.5, 0.5, 0.5),
			 COVector3D(-0.5, -0.5, 0.5),
			 COVector3D(0.5, -0.5, 0.5),
			 COVector3D(0.5, 0.5, 0.5),
			 COVector3D(0.0, 0.0, 1.0)
			 );
//

   drawPlane(
	         COVector3D(-0.5, -0.5, 0.5),
			 COVector3D(-0.5, -0.5, -0.5),
			 COVector3D(0.5, -0.5, -0.5),
			 COVector3D(0.5, -0.5, 0.5),
			 COVector3D(0.0, -1.0, 0.0)
			 );

   drawPlane(
			 COVector3D(0.5, 0.5, 0.5),
			 COVector3D(0.5, 0.5, -0.5),
			 COVector3D(-0.5, 0.5, -0.5),
	         COVector3D(-0.5, 0.5, 0.5),
			 COVector3D(0.0, 1.0, 0.0)
			 );

//

   drawPlane(
			 COVector3D(-0.5, 0.5, 0.5),
			 COVector3D(-0.5, 0.5, -0.5),
			 COVector3D(-0.5, -0.5, -0.5),
	         COVector3D(-0.5, -0.5, 0.5),
			 COVector3D(-1.0, 0.0, 0.0)
			 );

   drawPlane(
	         COVector3D( 0.5, -0.5, 0.5),
			 COVector3D( 0.5, -0.5, -0.5),
			 COVector3D( 0.5, 0.5, -0.5),
			 COVector3D( 0.5, 0.5, 0.5),
			 COVector3D( 1.0, 0.0, 0.0)
			 );


}

void GRPrimRenderer::drawPlane(COVector3D& p1,
							   COVector3D& p2,
							   COVector3D& p3,
							   COVector3D& p4)
{
   COVector3D v1, v2;
   COVector3D n;

   v1 = p2 - p1;
   v2 = p3 - p1;

   n = v1 * v2;
   n = n.unit();

   glPushMatrix();

   glBegin(GL_TRIANGLES);

   glNormal3f( n.x , n.y, n.z );

   glVertex3f(p1.x , p1.y, p1.z);
   glVertex3f(p2.x , p2.y, p2.z);
   glVertex3f(p3.x , p3.y, p3.z);
   
   glVertex3f(p3.x , p3.y, p3.z);
   glVertex3f(p4.x , p4.y, p4.z);
   glVertex3f(p1.x , p1.y, p1.z);
   
   glEnd();

   glPopMatrix();
}


void GRPrimRenderer::drawPlane(COVector3D& p1,
							   COVector3D& p2,
							   COVector3D& p3,
							   COVector3D& p4,
							   COVector3D& n)
{
   COVector3D un = n.unit();

   glPushMatrix();

   glBegin(GL_TRIANGLES);

   glNormal3f( un.x , un.y, un.z );

   glVertex3f(p1.x , p1.y, p1.z);
   glVertex3f(p2.x , p2.y, p2.z);
   glVertex3f(p3.x , p3.y, p3.z);
   
   glVertex3f(p3.x , p3.y, p3.z);
   glVertex3f(p4.x , p4.y, p4.z);
   glVertex3f(p1.x , p1.y, p1.z);
   
   glEnd();



   glPopMatrix();
}