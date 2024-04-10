/*
Sparse implementation demo 
for the numerical methods lecture

Lukasz Kuszner 
2019
*/

#ifndef _my_sparse
#define _my_sparse

#include<iostream>
#include<map>
#include<utility>
#include"matrix.h"

class sparse : matrix_base {
 public:
  std::map<std::pair<int,int>, double> v;
  sparse(int n) {
    s = n;
    set_init_values();
  }
  void setvalue(int i, int j, double x);
  double getvalue(int i, int j) const;
  void set_init_values();
  friend std::ostream& operator<< (std::ostream& stream, const sparse& m);
};

#endif
