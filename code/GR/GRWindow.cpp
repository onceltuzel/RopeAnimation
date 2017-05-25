//
//  Created by Oncel Tuzel.
//  Copyright © 2004 Oncel Tuzel. All rights reserved.
//


#include "GRWindow.h"


GRWindow::GRWindow()
{
   width_ = 640;
   height_ = 480;

}
GRWindow::~GRWindow()
{

}

void GRWindow::init(int argc, char** argv)
{
   glutInit(&argc, argv);
   glutInitWindowPosition(0,0);
   glutInitWindowSize(width_, height_);

   glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);

   glutCreateWindow("Rope Animation");
}


void GRWindow::reshape(int width, int height)
{
   width_ = width;
   height_ = height;   
}
