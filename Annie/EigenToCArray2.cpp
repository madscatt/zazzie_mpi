#include <Eigen/Dense>
#include <iostream>
#include <vector>
//#include <array>

using namespace Eigen;
using namespace std;
int main()
{
	std::vector<float> x ;
	std::vector<float> y ;
	std::vector<float> z ;

	int x_values[5] = {0, 1, 2, 3, 4};
	int y_values[4] = {0, 1, 2, 3};
	int z_values[3] = {0, 1, 2};


	int count_x = sizeof(x_values)/sizeof(*x_values);
	int count_y = sizeof(y_values)/sizeof(*y_values);
	int count_z = sizeof(z_values)/sizeof(*z_values);

	int count; 
	for(count = 0; count!=count_x; ++count){
	x.push_back(x_values[count]);
	cout << x[count] << endl;
	}

cout << "" << endl;

	for(count = 0; count!=count_y; ++count){
	y.push_back(y_values[count]);
	cout << y[count] << endl;
	}

cout << "" << endl;

	for(count = 0; count!=count_z; ++count){
	z.push_back(z_values[count]);
	cout << z[count] << endl;
	}
	

	Eigen::Array<int,Eigen::Dynamic,Eigen::Dynamic> x_Eigen ; 
        Eigen::Array<int,Eigen::Dynamic,Eigen::Dynamic> y_Eigen ; 
        Eigen::Array<int,Eigen::Dynamic,Eigen::Dynamic> z_Eigen ; 

	for(count = 0; count!=count_x; ++count){
	x_Eigen(count) = x[count];
	 cout << x_Eigen(count) << endl;
	}
cout << x_Eigen << endl;
cout << "" << endl;

	for(count = 0; count!=count_y; ++count){
	y_Eigen(count) = y[count];
	 
	}
cout << y_Eigen << endl;
cout << "" << endl;

	for(count = 0; count!=count_z; ++count){
	z_Eigen(count) = z[count];
	 
	}
cout << z_Eigen << endl;

	
}
