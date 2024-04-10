#include "sparse.h"

void sparse::setvalue(int i, int j, double x) {
  auto p = std::make_pair(i,j);
  v[p]=x;
}

double sparse::getvalue(int i, int j) const {
  auto p = std::make_pair(i,j);
  if (v.find(p)!=v.end())
      return v.at(p);
  else
    return 0;
}

void sparse::set_init_values() {
  setvalue(0,0,1);
  setvalue(size()-1,size()-1, 1);
  
  for (int w=1; w<size()-1; w++) {
    setvalue(w,w-1,-0.5);	
    setvalue(w,w,1);
    setvalue(w,w+1,-0.5);      
  }
  return;
}

std::ostream& operator<< (std::ostream& stream, const sparse& m) {
  m.serialize(stream);
  return stream;
}
