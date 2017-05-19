#ifndef MYMATH_H_INCLUDED
#define MYMATH_H_INCLUDED

/*

Author: Irina Hashmi (ihashmi@masonlive.gmu.edu)
Advisor: Prof. Amarda Shehu
Time and place: March, 2011, Department of Computer Science, George Mason University, VA, USA
Keywords:
Goal:
Reference:
Change log: FILL THIS UP AS CODE CHANGES
Note:

*/

#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>

//using namespace std;
typedef float real;

// dist = (ax-bx)^2 + (ay-by)^2 + (az-bz)^2
inline real Distance(real rx, real ry, real rz, real lx, real ly, real lz){
  return (
          sqrt((rx-lx)*(rx-lx) + (ry-ly)*(ry-ly)  + (rz-lz)*(rz-lz))
        );
}

inline real DistanceSqr(real rx, real ry, real rz, real lx, real ly, real lz){
  return ((rx-lx)*(rx-lx) + (ry-ly)*(ry-ly)  + (rz-lz)*(rz-lz));
}

inline unsigned long Choose(int n , int k){
    if (k > n) {
        return 0;
    }
    unsigned long r = 1;
    for (unsigned long d = 1; d <= k; ++d) {
        r *= n--;
        r /= d;
    }
    return r;
}


inline void Combination(int n, int r, std::vector<int>&combination) {

  std::vector<bool> v(n);
  std::fill(v.begin() + r, v.end(), true);
  do {
       for (int i = 0; i < n; ++i) {
           if (!v[i]) {
               combination.push_back(i+1);
           }
       }

  }while (std::next_permutation(v.begin(), v.end()));
  //return 0;
}

//count votes between 1 and 0
inline int CountOccurence(std::vector<int> v){
	int ones=0, zeros=0;
	for(int i = 0 ; i < v.size(); i++){
		if (v[i]==1){
			ones++;
		}
		else{
			zeros++;
		}
	}
	if (ones>zeros){
		return 1;
	}
	else if(zeros>ones) {
		return 0;
	}
	else{
		return 1;
	}
}




#endif // MYMATH_H_INCLUDED
