/*
Author: Irina Hashmi (ihashmi@masonlive.gmu.edu)

Advisor: Prof. Amarda Shehu

Time and place: April, 2011, Department of Computer Science, George Mason University, VA, USA

Keywords: Lexicographic ordering

Goal: Take three points of a triangle order them according to the distance from origin

Reference: None

Note: Nne

Change log: FILL THIS UP AS IT CHANGES

Who calls: HashTriangle.cpp

*/

#ifndef LEXICOGRAPHIC_H_INCLUDED
#define LEXICOGRAPHIC_H_INCLUDED
#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>

using namespace std;
typedef float real;

int* sort_func(int *sorted_index, real d1, real d2, real d3)
{
   	std::vector< std::vector<real> > v( 3, std::vector<real>( 2 ) );

    v[0][0] = d1;
    v[0][1] = 1;
    v[1][0] = d2;
    v[1][1] = 2;
    v[2][0] = d3;
    v[2][1] = 3;


    sort(v.begin(), v.end());
    sorted_index[0] = v[0][1];
    sorted_index[1] = v[1][1];
    sorted_index[2] = v[2][1];
}

void copy_arr(real* result, real* a, real* b, real* c){
  result[0]= a[0];
  result[1]= a[1];
  result[2]= a[2];

  result[3]= b[0];
  result[4]= b[1];
  result[5]= b[2];

  result[6]= c[0];
  result[7]= c[1];
  result[8]= c[2];

}

real* lexicographic_order(real*result, real* p1, real* p2, real* p3){
  real dist1, dist2, dist3;
  int sorted_index[3];

  dist1 = sqrtf(p1[0]*p1[0] + p1[1]*p1[1] + p1[2]*p1[2]);
  dist2 = sqrtf(p2[0]*p2[0] + p2[1]*p2[1] + p2[2]*p2[2]);
  dist3 = sqrtf(p3[0]*p3[0] + p3[1]*p3[1] + p3[2]*p3[2]);

  sort_func(sorted_index, dist1, dist2, dist3);


  if (sorted_index[0]==1 and sorted_index[1]==2 and sorted_index[2]==3){
    copy_arr(result, p1, p2, p3);
  }
  else if (sorted_index[0]==1 and sorted_index[1]==3 and sorted_index[2]==2){
    copy_arr(result, p1, p3, p2);
  }
  else if (sorted_index[0]==2 and sorted_index[1]==1 and sorted_index[2]==3){
    copy_arr(result, p2, p1, p3);
  }
  else if (sorted_index[0]==2 and sorted_index[1]==3 and sorted_index[2]==1){
    copy_arr(result, p2, p3, p1);
  }
  else if (sorted_index[0]==3 and sorted_index[1]==2 and sorted_index[2]==1){
    copy_arr(result, p3, p2, p1);
  }
  else {
    copy_arr(result, p3, p1, p2);
  }

  return result;
}



#endif // LEXICOGRAPHIC_H_INCLUDED
