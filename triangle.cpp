/*

Author: Irina Hashmi (ihashmi@masonlive.gmu.edu)
Advisor: Prof. Amarda Shehu
Time and place: Oct, 2013, Department of Computer Science, George Mason University, VA, USA
Goal: Triangle class implementation
Compile:
Run:
Copyright (c) 2013 Irina Hashmi. All rights reserved.

*/

#include "triangle.h"
#include "../utils/matchtriangle.h"
#include "../utils/transformation.h"
#include "../utils/definitions.h"
#include "../utils/matrix/matrix.h"


using namespace ihashmi;
typedef float real;

namespace ihashmi{

  /* ********* Class Implementation ****** */
  Triangle::Triangle(){

  }

  Triangle::Triangle(real x1_, real y1_, real z1_,
                     real x2_, real y2_, real z2_,
                     real x3_, real y3_, real z3_,
                     real score, int shape_){
    x1 = x1_;
    y1 = y1_;
    z1 = z1_;
    x2 = x2_;
    y2 = y2_;
    z2 = z2_;
    x3 = x3_;
    y3 = y3_;
    z3 = z3_;
    //init vertices
    v1[0] = x1_;
    v1[1] = y1_;
    v1[2] = z1_;

    v2[0] = x2_;
    v2[1] = y2_;
    v2[2] = z2_;

    v3[0] = x3_;
    v3[1] = y3_;
    v3[2] = z3_;

    //update vectors
    A.updateX(x1_);
    A.updateY(y1_);
    A.updateZ(z1_);

    B.updateX(x2_);
    B.updateY(y2_);
    B.updateZ(z2_);

    C.updateX(x3_);
    C.updateY(y3_);
    C.updateZ(z3_);

    shape = shape_;
    CalculateCOM();
  }

  Triangle::~Triangle(){

  }

  void Triangle::CalculateCOM(){
    centerOfMass[0] = (x1 + y1 + z1)/3;
    centerOfMass[1] = (x2 + y2 + z2)/3;
    centerOfMass[2] = (x3 + y3 + z3)/3;

  }

  //vertex1
  real Triangle::Getx1(){
    return x1;
  }

   real Triangle::Gety1(){
    return y1;
  }

   real Triangle::Getz1(){
    return z1;
  }

  //vertex2
  real Triangle::Getx2(){
    return x2;
  }

   real Triangle::Gety2(){
    return y2;
  }

   real Triangle::Getz2(){
    return z2;
  }

  //vertex3
  real Triangle::Getx3(){
    return x3;
  }

   real Triangle::Gety3(){
    return y3;
  }

   real Triangle::Getz3(){
    return z3;
  }

  real* Triangle::ReturnCOM(){
    return centerOfMass;
  }

  Vector3 Triangle::GetA(){
    return A;
  }
  Vector3 Triangle::GetB(){
    return B;
  }
  Vector3 Triangle::GetC(){
    return C;
  }

  int Triangle::ReturnShape(){
      return shape;
  }

  bool Triangle::CheckShapeComplementarity(Triangle& tr){

      if (
          (ReturnShape() == CAP  && tr.ReturnShape() == PIT)
      ||  (ReturnShape() == PIT  && tr.ReturnShape() == CAP)
      ||  (ReturnShape() == CAP  && tr.ReturnShape() == BELT)
      ||  (ReturnShape() == BELT && tr.ReturnShape() == CAP)
      ||  (ReturnShape() == BELT && tr.ReturnShape() == PIT)
      ||  (ReturnShape() == PIT  && tr.ReturnShape() == BELT)){
        return true;

      }
      else{
        return false;
      }


  }

  bool Triangle::PairwiseTriangle(Triangle& tr){
      bool flag = ifFit(x1, y1, z1,
                     x2, y2, z2,
                     x3, y3, z3,
                     tr.x1, tr.y1, tr.z1,
                     tr.x2, tr.y2, tr.z2,
                     tr.x3, tr.y3, tr.z3);
      return flag;
  }

  //take three non collinear points from a triangle A, B, C
  //output: quat, rotation matrix and translation vector component wrt global
  //If you apply this rotation and translation A = 0,0,0, B = dist,0,0 and
  //C=point on xy plane
  void Triangle::TransformTriangle(real* rotationMatrix, Vector3& translationVector){
    //cout << "inside transformTriangle " << endl;
    int i;
    vector<real> q;
    real elt[4];
    Vector3 orig = GetA();
    Vector3 x = GetB();
    Vector3 xy = GetC();

    /* calculates the rotation matrix */
    Vector3 r1 = x-orig;
    Vector3 xyminusorigin = xy - orig;

    Vector3 r3 = xyminusorigin&r1;
    Vector3 r2 = r3 &r1;

    /* calculates the unit vector */
    r1 = r1.getUnitVector();
    r2 = r2.getUnitVector();
    r3 = r3.getUnitVector();

    Quaternion(q, r1, r2, r3);
    real mag = q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3];
    mag = sqrt(mag);
    /* normalize the quaternion value */
    q[0] /= mag;
    q[1] /= mag;
    q[2] /= mag;
    q[3] /= mag;
    elt[0] = q[1], elt[1] = q[2], elt[2] = q[3], elt[3] = q[0];
    //rotation component to be applied on orig
    Vector3 rotateOrig = TransformVector(rotationMatrix, elt, orig);
    // The translation component is the negation of the orig vector
    translationVector =  rotateOrig * (-1.0);
    //cout << translationVector << endl;
    //cout << "END TRANSFORM TRI ************************************\n";

  }

  //finds a triangle wihtin 7A from the base chain
  int Triangle::FindANeighborTriangle(int triangleId,vector<Triangle>& triangleVect){

    real* baseBaricenter = this->ReturnCOM();
    int totalTry =  triangleVect.size();
    int newTriangleId;
    bool f = false;
    while (totalTry){

      int r = IntUrandom(triangleVect.size());

      if (r!= triangleId){
        real* newBaricenter = triangleVect[r].ReturnCOM();
         // if baricenter is within 5A
        real b = Distance(baseBaricenter[0], baseBaricenter[1], baseBaricenter[2],
              newBaricenter[0], newBaricenter[1], newBaricenter[2]);
        if( b  < (PERTURBBARICENTERTHRESHOLD - RMP_EPSILON)){
          newTriangleId = r;
          totalTry = 0;
          f = true;
        }
        else{
          totalTry--;
        }
      }
    }//while

    if (f)
      return newTriangleId;
    else
      return -1;
  }//function



  /*
  Input:
  referenceRotationMatrix = unit1 reference frame rotation
  moveRotationMatrix = unit2 move frame rotation
  referenceTranslationVector3 = unit1 reference frame translation
  moveTranslationVector3 = unit2 move frame translation

  Output:
  finalRotationMatrix = final rotation between traingle 1 and triangle 2
  finalTranslationVector3 = final translation between traingle 1 and triangle 2
  */
//  void Triangle::TransformationBetweenTwoTriangles(real* referenceRotationMatrix, real *moveRotationMatrix, real* finalRotationMatrix,
//                  Vector3& referenceTranslationVector3, Vector3& moveTranslationVector3, Vector3& finalTranslationVector3){
//
//
//	  //if two traingles are P(reference) and Q(move)
//	  //final rot from Q to P -> R_PQ = R_WP.transpose() * R_WQ
//	  real referenceRotationMatrixTranspose[9];
//      //1. R_WP.transpose()
//	  TransposeMatrix(referenceRotationMatrixTranspose, referenceRotationMatrix);
//	  //2. R_WP.transpose() * R_WQ
//	  MultMatrix(finalRotationMatrix,
//				 referenceRotationMatrixTranspose,
//				 moveRotationMatrix);
//
//	  //final translation
//	  //T = R_WP * A - R_WQ*D
//	  //T_PQ = R_WP.transpose * T
//	  Vector3 distanceOrigReferenceMove =
//					moveTranslationVector3 - referenceTranslationVector3;
//	  finalTranslationVector3 =
//	  MultMatrixVector(referenceRotationMatrixTranspose,
//					   distanceOrigReferenceMove);
//
//
//
//  }


//  void Triangle::TransformTriangle1(real* rotationMatrix, Vector3& translationVector, real *elt){
//    cout << "inside transformTriangle " << endl;
//    int i;
//    vector<real> q;
//    real elt[4];
//    Vector3 orig = GetA();
//    Vector3 x = GetB();
//    Vector3 xy = GetC();
//
//    /* calculates the rotation matrix */
//    Vector3 r1 = x-orig;
//    Vector3 xyminusorigin = xy - orig;
//
//    Vector3 r3 = xyminusorigin&r1;
//    Vector3 r2 = r3 &r1;
//
//    /* calculates the unit vector */
//    r1 = r1.getUnitVector();
//    r2 = r2.getUnitVector();
//    r3 = r3.getUnitVector();
//
//    Quaternion(q, r1, r2, r3);
//    real mag = q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3];
//    mag = sqrt(mag);
//    /* normalize the quaternion value */
//    q[0] /= mag;
//    q[1] /= mag;
//    q[2] /= mag;
//    q[3] /= mag;
//    elt[0] = q[1], elt[1] = q[2], elt[2] = q[3], elt[3] = q[0];
//    //rotation component to be applied on orig
//    Vector3 rotateOrig = TransformVector(rotationMatrix, elt, orig);
//    // The translation component is the negation of the orig vector
//    translationVector =  rotateOrig * (-1.0);
//    //cout << translationVector << endl;
//    //cout << "END TRANSFORM TRI ************************************\n";
//
//  }

}
