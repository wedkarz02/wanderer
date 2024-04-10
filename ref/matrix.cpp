#include "matrix.h"

void matrix::set_init_values() {
  for (int w=0; w<size(); w++) {
    for (int k=0; k<size(); k++) {
      setvalue(w,k,0);
    }
  }
  
  setvalue(0,0,1);
  setvalue(size()-1,size()-1, 1);
  
  for (int w=1; w<size()-1; w++) {
    setvalue(w,w-1,-0.5);	
    setvalue(w,w,1);
    setvalue(w,w+1,-0.5);      
  }
  return;
}

void matrix_base::serialize (std::ostream& stream) const {
  for (int i=0; i<size(); i++) {
    for (int j=0; j<size(); j++)
      stream << getvalue(i,j) << " ";
      stream << std::endl;
  }
}

