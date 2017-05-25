
//
//  Created by Oncel Tuzel.
//  Copyright © 2004 Oncel Tuzel. All rights reserved.
//


#ifndef _grarrangementrenderer_h
#define _grarrangementrenderer_h

#include "CMDataTypes.h"
#include "GRPrimRenderer.h"
#include "ARMainSystem.h"

class GRArrangementRenderer {
public:
   void incrementArrangement(int n);

   void setArrangement(ARMainSystem* arrangement);
   void change();

   void render();
   void renderPlane(ARPlane* plane, int* id, int st, int sign);
   void renderCell(ARCell* cell, int st);

   void setRenderStatus(int status);
   int  getRenderStatus();

private:
   void sortPoints(TVectorVector& points, CMVector3D& n);
   
private:
   ARMainSystem* arrangement_;
   
   GRPrimRenderer primRenderer_;

   ARCell* activeCell_;
   ARCell* neighborCell_;


   ARPlane* activePlane_;

   int neighborNo_;
   int status_;
};

#endif