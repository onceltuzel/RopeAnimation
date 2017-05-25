//
//  Created by Oncel Tuzel.
//  Copyright © 2004 Oncel Tuzel. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include "GL/glut.h"


#ifndef _grwindow_h
#define _grwindow_h

class GRWindow {
public:
   GRWindow();
   ~GRWindow();


   void init(int argc, char** argv);
   void reshape(int width, int height);

   int width_;
   int height_;

private:
  

   
};



#endif