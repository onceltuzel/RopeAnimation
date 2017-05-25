//
//  Created by Oncel Tuzel.
//  Copyright © 2004 Oncel Tuzel. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include "GL/glut.h"
#include <windows.h>

#include "GRPrimRenderer.h"
#include "GRWindow.h"

#include "RSRopeSystem.h"

#ifndef _grscene_h
#define _grscene_h

class GRScene {
public:
   GRScene();
   GRScene(GRWindow* window);
   ~GRScene();


   void init();

   void display();
	void move();

   void reshape(int width, int height);

   void keyboardInput(unsigned char key, int x, int y);
   void specialInput(int key, int x, int y);

private:

   void initGLStates();
   void initCamera();
	void initRopeSystem();
   void setShadingParams();
   void setLights(); 


   enum ERotationType {
	   EPAN,
	   ETILT
   };

   void getCameraParams();
   void rotateCamera(ERotationType type, double degree);
   void moveCamera(double distance);

   void drawGrid();
   void drawCoordinates();
	void drawRopes();
	void drawShadows();
	void drawGroundPlane();
   
private:

   COVector3D cameraUp_;
   COVector3D cameraDir_;
   COVector3D cameraLeft_;

   //Geometric Primitive Renderer
   GRPrimRenderer primRenderer_;

   //Arrangement Renderer

   //Window Object
   GRWindow*   window_;

   int drawGrid_;
   int drawCoordinates_;
   int drawControlPoints_; 

	RSRopeSystem*	ropeSystem_;

	LARGE_INTEGER lastUpdateTime_;
	LARGE_INTEGER frequency_;
};



#endif