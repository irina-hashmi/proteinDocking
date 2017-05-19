/*
Author: Irina Hashmi (ihashmi@masonlive.gmu.edu)
Advisor: Prof. Amarda Shehu
Time and place: Oct 28, 2011, Department of Computer Science, George Mason University, VA, USA
Keywords: Baricenter
Goal: Triangle Class
Reference: None
Compile:
Run:
*/

#ifndef TRIANGLE_H_INCLUDED
#define TRIANGLE_H_INCLUDED

#include <cstring>
#include <iomanip>
#include <vector>
#include "../utils/Vector3.h"
#include "../math/mymath.h"
#include "../utils/definitions.h"
#include "../math/Random.h"


typedef float real;


namespace ihashmi{
  class Triangle{
    public:
      //constructor
      Triangle();
      Triangle(real x1, real y1, real z1,
                        real x2, real y2, real z2,
                        real x3, real y3, real z3,
                        real score, int shape);
      //Destructor
      ~Triangle();
      void CalculateCOM();
      real Getx1();
      real Gety1();
      real Getz1();

      real Getx2();
      real Gety2();
      real Getz2();

      real Getx3();
      real Gety3();
      real Getz3();

      Vector3 GetA();
      Vector3 GetB();
      Vector3 GetC();
      real* ReturnCOM();
      friend class Configuration;

      bool PairwiseTriangle(Triangle& tr);
      void TransformTriangle(real* rotationMatrix, Vector3& finalTranslation);
      bool CheckShapeComplementarity(Triangle& tr);
      int ReturnShape();
      int FindANeighborTriangle(int triangleId, std::vector<Triangle>& triangleVect);

    private:
      real x1, y1, z1;
      real x2, y2, z2;
      real x3, y3, z3;
      Vector3 A, B, C;
      real v1[3], v2[3], v3[3];
      real centerOfMass[3];
      int shape;    //shape of the triangle

  };
}

#endif // TRIANGLE_H_INCLUDED
