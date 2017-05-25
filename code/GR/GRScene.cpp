//
//  Created by Oncel Tuzel.
//  Copyright © 2004 Oncel Tuzel. All rights reserved.
//


#include "GRScene.h"
#include <math.h>


#define PAN_STEP  1
#define TILT_STEP 1
#define MOVE_STEP 1
#define GRID_STEP 5
#define POINT_STEP 0.5
#define GRID_RANGE 100
#define ROPE_SEGMENT  0.01
#define ROPE_POINT_COUNT  11



GRScene::GRScene()
{
}


GRScene::GRScene(GRWindow* window)
{
   window_= window;
   drawGrid_ = 0;
   drawCoordinates_ = 0;
	drawControlPoints_ = 0;
}


GRScene::~GRScene()
{
	delete ropeSystem_;
}



void GRScene::initCamera()
{
   glMatrixMode(GL_MODELVIEW); 	  

   gluLookAt(0.0, 0.3, -1.2,
			    0.0, 0.3, 0.0,			 
				 0.0, 1.0, 0.0);

}

void GRScene::initGLStates()
{
   glEnable(GL_DEPTH_TEST);

}


void GRScene::setLights() 
{
   //Lights
   GLfloat diffuseLight0[4] = {0.3, 0.3f, 0.3f, 1.0};
   GLfloat lightDirection0[4] = {1.0f, 1.0f, -1.0f, 0.0f};

   GLfloat diffuseLight1[4] = {0.5, 0.5f, 0.5f, 1.0};
   GLfloat lightDirection1[4] = {-1.0f, 1.0f, 1.0f, 0.0f};

   GLfloat diffuseLight2[4] = {0.2, 0.2f, 0.2f, 1.0};
   GLfloat lightDirection2[4] = {0.0f, 1.0f, 0.0f, 0.0f};


   GLfloat ambientLight[4] = {0.3, 0.3f, 0.3f, 1.0};
   
   
   glLightfv(GL_LIGHT0, GL_AMBIENT, ambientLight);

   glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseLight0);
   glLightfv(GL_LIGHT0, GL_POSITION, lightDirection0);

   glLightfv(GL_LIGHT1, GL_DIFFUSE, diffuseLight1);
   glLightfv(GL_LIGHT1, GL_POSITION, lightDirection1);

   glLightfv(GL_LIGHT2, GL_DIFFUSE, diffuseLight2);
   glLightfv(GL_LIGHT2, GL_POSITION, lightDirection2);
   
   //glEnable(GL_LIGHT0);
   //glEnable(GL_LIGHT1);
   //glEnable(GL_LIGHT2);
   glDisable(GL_LIGHTING);
   
}


void GRScene::init()
{
	QueryPerformanceFrequency(&frequency_);
	
   initGLStates();
   initCamera();
	initRopeSystem();
   
}



void GRScene::reshape(int width, int height)
{
   window_->reshape(width, height);

   glViewport(0, 0, window_->width_, window_->height_);


   glMatrixMode(GL_PROJECTION);
   glLoadIdentity(); 
   
   gluPerspective( 40.0,
			          ((double)window_->width_)/window_->height_,				   
						 0.5, 				   
						 10000.0);

   glMatrixMode(GL_MODELVIEW); 

}


void GRScene::keyboardInput(unsigned char key, int x, int y)
{
   switch (key)
   {
	  case 'o':		  
		  rotateCamera(ETILT, -TILT_STEP);
		  break;
	  case 'l':
		  rotateCamera(ETILT, TILT_STEP);
		  break;
	  case 'e':		 
		  moveCamera(MOVE_STEP);		 
		  break;
	  case 'd':		 
		  moveCamera(-MOVE_STEP);		 
		  break;
	  case 'f':		 
		  rotateCamera(EPAN, PAN_STEP);		 
		  break;
	  case 's': 		 
		  rotateCamera(EPAN, -PAN_STEP);		 
		  break;
	  case 'g':		 
		  drawGrid_ = !drawGrid_;		 
		  break;
	  case 'c':		 
		  drawCoordinates_ = !drawCoordinates_;		 
		  break;
	  case 'p':		 
		  drawControlPoints_ = !drawControlPoints_;		 
		  break;
	  case 27 :
		 exit(0);
		 break;
   }

}
void GRScene::specialInput(int key, int x, int y)
{
	COVector3D* fixedPoint = ropeSystem_->GetFixedPoint();
	RSRope* rope = ropeSystem_->GetRope(0);

	switch (key)
	{
	case GLUT_KEY_LEFT:
		fixedPoint->x+=0.02;
		break;
	case GLUT_KEY_RIGHT:
		fixedPoint->x-=0.02;
		break;
	case GLUT_KEY_UP:
		fixedPoint->z+=0.02;
		break;
	case GLUT_KEY_DOWN:	  
		fixedPoint->z-=0.02;
		break;
	case GLUT_KEY_PAGE_UP:
		fixedPoint->y+=0.02;
		break;
	case GLUT_KEY_PAGE_DOWN:
		fixedPoint->y-=0.02;
		break;
	}
}



void GRScene::getCameraParams()
{
	double mat[16];

	glGetDoublev( GL_MODELVIEW_MATRIX , mat);

	cameraLeft_.x = mat[0];
	cameraLeft_.y = mat[4];
	cameraLeft_.z = mat[8];

	cameraUp_.x = mat[1];
	cameraUp_.y = mat[5];
	cameraUp_.z = mat[9];

	cameraDir_.x = mat[2];
	cameraDir_.y = mat[6];
	cameraDir_.z = mat[10];
}


void GRScene::rotateCamera(ERotationType type, double degree)
{

	double mat[16];
	glGetDoublev( GL_MODELVIEW_MATRIX , mat);

	glLoadIdentity();

	switch (type)
	{
	case EPAN:
		glRotatef(degree, 0, 1, 0);
		break;
	case ETILT:
		glRotatef(degree, 1, 0, 0);
		break;
	}
	
	glMultMatrixd(mat);
	
}

void GRScene::moveCamera(double distance)
{
	double mat[16];
	glGetDoublev( GL_MODELVIEW_MATRIX , mat);

	glLoadIdentity();
	glTranslatef(0, 0, distance);
	glMultMatrixd(mat);

}





void GRScene::setShadingParams()
{
   
   glPolygonMode( GL_FRONT_AND_BACK ,
 			         GL_FILL );
				  
   glEnable(GL_CULL_FACE);
   glCullFace(GL_BACK);
   
}

void GRScene::drawGroundPlane()
{

	glPushMatrix();

	glBegin(GL_QUADS);
	glColor3d(1.0,0.1,0.1);
	glVertex3d(10.0,0.0,-10.0);
	glVertex3d(-10.0,0.0,-10.0);
	glColor3d(0.0,0.0,0.0);
	glVertex3d(-10.0,0.0,20.0);
	glVertex3d(10.0,0.0,20.0);
	glEnd();

	glPopMatrix();
	
}


void GRScene::drawCoordinates()
{
   int i;


   glLineWidth(1);
   glPointSize(2);

   glColor3f(1,0,0);
   glBegin(GL_POINTS);
   for (i=-GRID_RANGE ; i<0 ; i+=POINT_STEP)
   {
	  glVertex3f(0, i, 0);
   }
   glEnd();
   glBegin(GL_LINES);
   glVertex3f(0, 0 , 0);
   glVertex3f(0, GRID_RANGE , 0);
   glEnd();

   glColor3f(0,1,0);
   glBegin(GL_POINTS);
   for (i=-GRID_RANGE ; i<0 ; i+=POINT_STEP)
   {
	  glVertex3f(i, 0, 0);
   }
   glEnd();
   glBegin(GL_LINES);
   glVertex3f(0, 0 , 0);
   glVertex3f(GRID_RANGE , 0, 0);
   glEnd();

   glColor3f(0,0,1);
   glBegin(GL_POINTS);
   for (i=-GRID_RANGE ; i<0 ; i+=POINT_STEP)
   {
	  glVertex3f(0, 0, i);
   }
   glEnd();

   glBegin(GL_LINES);
   glVertex3f(0 , 0, 0);
   glVertex3f(0, 0, GRID_RANGE);
   glEnd();

   glColor3f(1,1,1);

   glLineWidth(1);
   glPointSize(1);


}

void GRScene::drawGrid()
{
   int i;

   glLineWidth(1);
   glPointSize(2);

   glColor3f(1,1,1);
   glBegin(GL_LINES);

   for (i=-GRID_RANGE ; i<= GRID_RANGE ; i+=GRID_STEP)
   {
	  if (i!=0)
	  {
		 glVertex3f(-GRID_RANGE, 0 , i);
		 glVertex3f(GRID_RANGE, 0 , i);
		 glVertex3f(i, 0, -GRID_RANGE);
		 glVertex3f(i, 0, GRID_RANGE);
	  }
   }

   glEnd();


}


void GRScene::display()
{
	
   setShadingParams();
   setLights();
   
   glClearColor(0.0f, 0.0f, 0.0f, 1.0f);

   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);   

   glPushMatrix();

   if (drawGrid_)
	  drawGrid();

   if (drawCoordinates_)
	  drawCoordinates();


	drawGroundPlane();
	drawRopes();
	drawShadows();
   
   glPopMatrix();
	
   glutSwapBuffers();
}


void GRScene::drawRopes()
{
	int i, j;
	double s;
	RSRope* rope;
	COVector3D point,*ctrlPoint;


   glLineWidth(3);
   glPointSize(5);

	j = 0;
	while (rope = ropeSystem_->GetRope(j++))
	{
		if (drawControlPoints_)
		{
			glColor3f(0,0,1);
			glBegin(GL_POINTS);	
			for (i=0 ; i<ROPE_POINT_COUNT ; i++)
			{
				ctrlPoint = rope->GetControlPoint(i);
				glVertex3d(ctrlPoint->x, ctrlPoint->y, ctrlPoint->z);
			}
			glEnd();
		}

		glColor3f(1.0,1.0,0.0);
		glBegin(GL_LINE_STRIP);	
		for (s=0 ; s<=1.0 ; s+=ROPE_SEGMENT)
		{
			rope->GetPosition(s, point);
			glVertex3d(point.x, point.y,point.z);
		}
		s = 1.0;
		rope->GetPosition(s, point);
		glVertex3d(point.x, point.y,point.z);
		glEnd();
	}

}

void GRScene::drawShadows()
{
	int j;
	double s;
	RSRope* rope;
	COVector3D point;


   glLineWidth(4);

	j = 0;
	while (rope = ropeSystem_->GetRope(j++))
	{
		glColor3f(0.01,0.01,0.01);
		glBegin(GL_LINE_STRIP);	
		for (s=0 ; s<=1.0 ; s+=ROPE_SEGMENT)
		{
			rope->GetPosition(s, point);
			glVertex3d(point.x, 0.01,point.z);
		}
		s = 1.0;
		rope->GetPosition(s, point);
		glVertex3d(point.x, 0.01,point.z);
		glEnd();
	}

}


void GRScene::initRopeSystem()
{
	RSRope* rope;
	
	ropeSystem_ = new RSRopeSystem(ROPE_POINT_COUNT);
	rope = new RSRope(ROPE_POINT_COUNT); 
	ropeSystem_->AddRope(rope);

	QueryPerformanceCounter(&lastUpdateTime_);

}

void GRScene::move()
{
	LARGE_INTEGER currentTime;

	QueryPerformanceCounter(&currentTime);	
	ropeSystem_->Iterate(((double)(currentTime.QuadPart - lastUpdateTime_.QuadPart))/(frequency_.QuadPart));
	lastUpdateTime_ = currentTime;
}