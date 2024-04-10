/*
Sparse implementation demo 
for the numerical methods lecture

Lukasz Kuszner 
2019
*/

#include "matrix.h"
#include "sparse.h"
#include <iostream>
#include <vector>
#include <cstring>


typedef std::vector<double> VecD;

void printV (VecD X) {
  for (auto elem: X){
    std::cout << elem <<  '\t';
  }
  std::cout << std::endl;
}

void Jacobi(matrix &M, VecD &B, VecD &X, int iter) {
  for (int i=0; i<iter; i++){   
    VecD Xtemp(X);
    for (int w=0; w<X.size(); w++){
      X[w]=B[w];
      for (int k=0; k<X.size(); k++){
	if (w!=k)
	  X[w]-=Xtemp[k]*M.getvalue(w,k);	  
      }
      X[w]/=M.getvalue(w,w);
    }
  }
}

void Jacobi(sparse &M, VecD &B, VecD &X, int iter) {
  for (int i=0; i<iter; i++){
    VecD Xtemp(X);
    for (int w=0; w<X.size(); w++){
      X[w]=B[w];
    }
    for (auto elem: M.v) {
      if (elem.first.first!=elem.first.second)
	X[elem.first.first]-=Xtemp[elem.first.second]*elem.second;
    }
    
    for (auto elem:M.v) {
      if (elem.first.first==elem.first.second)
	X[elem.first.first]/=elem.second;
    }
  }
}

void setB(VecD &B) {
  B[0]=1;
  for (int w=1; w<B.size(); w++) {
    B[w]=0;
  }
}



int main(int argc, char **argv) {

  if (argc==1) {
    int n=6;
    matrix M(n);
    VecD B(n);
    setB(B);
    std::cout << M;
    printV(B);
    VecD X1(n);
  
    std::cout <<std::endl;
  
    printV(X1);
    std::cout <<std::endl;
    Jacobi(M, B, X1, 20);
    std::cout <<std::endl;
    printV(X1);

    sparse S(n);
    VecD X2(n);
    std::cout << S;
    Jacobi(S, B, X2, 20);
    printV(X2);
  } 
  else if (argc==3) {
    int size_m = atoi(argv[2]);
    bool is_sparse = (0==strcmp(argv[1], "sparse"));
    VecD X(size_m);

    VecD B(size_m);
    setB(B);
    
    if (is_sparse) {
      sparse S(size_m);
      Jacobi(S, B, X, 100);
    } else {
      matrix M(size_m);
      Jacobi(M, B, X, 100);
    }
    double sum=0;
    for (auto x : X)
      sum+=x;
    std::cout << sum << std::endl;
  }
  else 
    {
      std::cout << "usage:" <<  std::endl;
      std::cout << "runme sparsë <size of linear system>" <<  std::endl;
      std::cout << "or:" <<  std::endl;
      std::cout << "runme normal <size of linear system>" <<  std::endl;
      std::cout << "like:" <<  std::endl;
      std::cout << "runme sparse 12000" <<  std::endl;
    }
  return 0;
}
