//
//  Created by Oncel Tuzel.
//  Copyright © 2004 Oncel Tuzel. All rights reserved.
//

#include <stdlib.h>
#include <stdio.h>
#include "GL/glut.h"

#include "GRScene.h"
#include "GRWindow.h"

GRScene* scene;
GRWindow* window;


void display()
{
   scene->display();	
}
 

void idle()
{
	scene->move();
   glutPostRedisplay();
}


void reshape(int width, int height)
{
   scene->reshape(width, height);

}


void keyboard(unsigned char key, int x, int y)
{
   scene->keyboardInput(key, x, y);
}

void special(int key, int x, int y)
{
   scene->specialInput(key, x, y);
}



void init(int argc, char** argv)
{
   //Arrangement construction


   //Window Initialization
   window = new GRWindow;
   window->init(argc, argv);

   //Scene initialization
   scene = new GRScene(window);
   scene->init();

   //Callback Initialization
   glutDisplayFunc(display);   
   glutReshapeFunc(reshape);
   glutIdleFunc(idle);

   //Input Controllers
   glutKeyboardFunc(keyboard);
   glutSpecialFunc(special);

}


void destroy()
{
   delete scene;
   delete window;
}


int main(int argc, char** argv)
{
   init(argc, argv);
   glutMainLoop();  
   destroy();


   return 0;
}