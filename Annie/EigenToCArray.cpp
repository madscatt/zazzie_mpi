#include <Eigen/Dense>
#include <iostream>
using namespace Eigen;
using namespace std;
int main()
{
  ArrayXXf  m(2,2);
  
  // assign some values coefficient by coefficient
  m(0,0) = 1.0; m(0,1) = 2.0;
  m(1,0) = 3.0; m(1,1) = m(0,1) + m(1,0);

  // print values to standard output
  cout << m << endl << endl;
 
ArrayXXf  t(2,2);
  // using the comma-initializer is also allowed
  t << 1.0,2.0,
       3.0,4.0;
     
  // print values to standard output
  cout << t << endl;

int foo[2][2];
int n; 
int k;
for(n=0;n<2;n++){
  for(k=0;k<2;k++){
     foo[n][k] = t(n,k);

     cout << foo[n][k] << endl;
}
}

}
