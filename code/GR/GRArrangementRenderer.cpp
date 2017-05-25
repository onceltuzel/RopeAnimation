
//
//  Created by Oncel Tuzel.
//  Copyright © 2004 Oncel Tuzel. All rights reserved.
//

#include "GRArrangementRenderer.h"
#include "CMFunctions.h"


void GRArrangementRenderer::setRenderStatus(int status)
{
   status_ = status;
}
int  GRArrangementRenderer::getRenderStatus()
{
   return status_;
}


void GRArrangementRenderer::sortPoints(TVectorVector& points, CMVector3D& pNormal)
{
   int i;
   int n = points.size();
   CMVector3D mean(0,0,0);
   TVectorVector sortedPoints;

   float* fPoints  = new float[n];
   int*	  iPoints  = new int[n];

   for (i=0 ; i < n ; i++)
   {
	  mean = mean + points[i];
	  sortedPoints.push_back(mean);
   }

   mean = mean/i;

   CMVector3D v1 = points[0] - mean;
   float l1 = v1.length();
   
   for (i=0 ; i<n ; i++)
   {
	  CMVector3D v2 = points[i] - mean;
	  float l2 = v2.length();
	  
	  fPoints[i] = (v1.dot(v2) / (l1*l2));
	  iPoints[i] = i;

	  CMVector3D cross = v1*v2;
	  cross.unitize();

	  if (pNormal.dot(cross) >= 0)
		 fPoints[i]=5.0f+fPoints[i];
	  else
		 fPoints[i]=-5.0-fPoints[i];
   }

   CMFunctions::sort(fPoints, n, iPoints);

   for (i=0 ; i<n ; i++)
   {
	  sortedPoints[i] = points[iPoints[i]];
   }
   
   for (i=0 ; i<n ; i++)
   {
	  points[i] = sortedPoints[i];
   }


   delete [] fPoints;
   delete [] iPoints;

}


void GRArrangementRenderer::setArrangement(ARMainSystem* arrangement)
{
   status_ = 1;
   neighborNo_ = 0;
   arrangement_ = arrangement;
   activeCell_ = (ARCell*) arrangement_->levelArrangements_[3][0];
   if (activeCell_ != activeCell_->subFaces_[neighborNo_]->superFaces_[0])
	  neighborCell_ = (ARCell*) activeCell_->subFaces_[neighborNo_]->superFaces_[0];
   else
	  neighborCell_ = (ARCell*) activeCell_->subFaces_[neighborNo_]->superFaces_[1];
}

void GRArrangementRenderer::incrementArrangement(int n)
{
   if (n+neighborNo_>=0 && n+neighborNo_<activeCell_->subFaces_.size())
   {
	  neighborNo_ += n;
   }
   if (activeCell_ != activeCell_->subFaces_[neighborNo_]->superFaces_[0])
	  neighborCell_ = (ARCell*) activeCell_->subFaces_[neighborNo_]->superFaces_[0];
   else
	  neighborCell_ = (ARCell*) activeCell_->subFaces_[neighborNo_]->superFaces_[1];
}



void GRArrangementRenderer::change()
{
   activeCell_ = neighborCell_;
   neighborNo_ = 0;

   if (activeCell_ != activeCell_->subFaces_[neighborNo_]->superFaces_[0])
	  neighborCell_ = (ARCell*) activeCell_->subFaces_[neighborNo_]->superFaces_[0];
   else
	  neighborCell_ = (ARCell*) activeCell_->subFaces_[neighborNo_]->superFaces_[1];
}

void GRArrangementRenderer::render()
{
   glPushMatrix();

   renderCell(activeCell_, 1);
   renderCell(neighborCell_, 0);

   glPopMatrix();

}

void GRArrangementRenderer::renderCell(ARCell* cell, int st)
{
   int i;
   ARPlane* plane;
   CMVector3D position;
   int sign;

   for (i=0 ; i<cell->subFaces_.size() ; i++)
   {
	  plane = (ARPlane*) cell->subFaces_[i];
	  position = cell->position_ - plane->position_;
	  if (position.dot(plane->plane_.normal_)>0)
		  sign = -1;
	  else
		  sign = 1;


	  renderPlane(plane, activeCell_->id_, st, sign);
   }
}


void GRArrangementRenderer::renderPlane(ARPlane* plane, int* id, int st, int sign)
{
   int i;
   AREdge* edge;
   TEdgeVectorPointer edges;
   TVectorVector points;
/*
   points.push_back(CMVector3D(1, 0, 0));
   points.push_back(CMVector3D(0.3, 0.7, 0));
   points.push_back(CMVector3D(0.0, 1, 0));
   points.push_back(CMVector3D(-0.3, 0.7, 0));
   points.push_back(CMVector3D(-0.3, -0.7, 0));
   points.push_back(CMVector3D(0.3, -0.7, 0));

   CMVector3D x(0, 0, -1);
   sortPoints(points, x);

   CMVector3D a[10];

   for (i=0 ; i<points.size() ; i++)
	  a[i] = points[i];
*/
/*
   int sign=0;

   for (i=0 ; i< MAX_PLANE ; i++)
   {
      if (plane->id_[i] != id[i])
	  {
		 sign = id[i];
		 break;
	  }
	  
   }
*/
   CMVector3D normal = plane->plane_.normal_ * -sign;

   edges = (TEdgeVectorPointer) &(plane->subFaces_);   

   CMVector3D v[2];
   int t = edges->size();
   for (i=0 ; i<edges->size(); i++)
   {
      edge = (*edges)[i];

	  edge->getVertices(v);

	  points.push_back(v[0]);
	  points.push_back(v[1]);
   }


   if (!status_ || !st)
   {
	  if (st)
		 glColor3f(0,0,0);
	  else
		 glColor3f(1,1,0);
	  //glDisable(GL_DEPTH_TEST);
	  glDisable(GL_LIGHTING);
      glLineWidth(2);
	 
	  glBegin(GL_LINES);
   
	  for (i=0 ; i<points.size() ; i++)
	  {
		 glVertex3f(points[i].x , points[i].y, points[i].z);
		 //glVertex3f(points[i+1].x , points[i+1].y, points[i+1].z);
	  }

	  glEnd();
	  glEnable(GL_LIGHTING);
	  glEnable(GL_DEPTH_TEST);

      glLineWidth(1);
	  glColor3f(1,1,1);

   }
   else
   {
      sortPoints(points, normal);

	  //glDisable(GL_LIGHTING);
	  glBegin(GL_POLYGON);

	  glNormal3f( -normal.x , -normal.y, -normal.z );

	  for (i=0 ; i<points.size() ; i++)
	  {
		 glVertex3f(points[i].x , points[i].y, points[i].z);

	  }

	  glEnd();
	  //glEnable(GL_LIGHTING);
   }

   points.clear();
}

