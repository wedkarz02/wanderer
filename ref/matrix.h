/*
Sparse implementation demo 
for the numerical methods lecture

Lukasz Kuszner 
2019
*/

#ifndef _my_matxix
#define _my_matxix

#include<iostream>

class matrix_base {
 protected:
  int s;
 public:
  virtual  void setvalue(int i, int j, double x) = 0; 
  virtual  double getvalue(int i, int j) const = 0;
  int size() const {return s;}

  virtual void set_init_values() {}
  void serialize (std::ostream& stream) const;
};

class matrix : matrix_base {
 private:
  double ** v = nullptr;
 public:
  virtual void set_init_values();
  matrix(int n) {
    s = n;
    v = new double*[n];
    for (int i=0; i<n; i++)
      v[i] = new double[n];
    set_init_values();
  }

  void setvalue(int i, int j, double x) {
    v[i][j] = x;
  }

  double getvalue(int i, int j) const {
    return v[i][j];
  }
  
  friend std::ostream& operator<< (std::ostream& stream, const matrix& m) {
    m.serialize(stream);
    return stream;
  }
};

#endif
