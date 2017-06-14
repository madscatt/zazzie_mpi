#include <Eigen/Dense>
#include <iostream>
#include <vector>

using namespace Eigen;
using namespace std;
int main()
{

//Defining vectors

	std::vector<float> x ;
	std::vector<float> y ;
	std::vector<float> z ;

	int x_values[5] = {0, 1, 2, 3, 4};
	int y_values[4] = {0, 1, 2, 3};
	int z_values[3] = {0, 1, 2};

//Finding size of the arrays

	int count_x = sizeof(x_values)/sizeof(*x_values);
	int count_y = sizeof(y_values)/sizeof(*y_values);
	int count_z = sizeof(z_values)/sizeof(*z_values);

//Populating vectors

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
cout << "" << endl;

//Defining Eigen Arrays

	Eigen::Array<int,Eigen::Dynamic,Eigen::Dynamic> x_Eigen; 
        Eigen::Array<int,Eigen::Dynamic,Eigen::Dynamic> y_Eigen ; 
        Eigen::Array<int,Eigen::Dynamic,Eigen::Dynamic> z_Eigen ; 

	x_Eigen.setZero(count_x, 1);
        y_Eigen.setZero(count_y, 1);
        z_Eigen.setZero(count_z, 1);

//Populating Eigen Arrays

	for(count = 0; count!=count_x; ++count){
	x_Eigen(count) = x[count];
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

//Defining C Arrays

	int xc[count_x];
	int yc[count_y];
	int zc[count_z];

//Populating C Arrays

	for(count = 0; count!=count_x; ++count){
	xc[count] = x_Eigen(count);
	cout << xc[count] << endl;
	}
cout << "" << endl;

	for(count = 0; count!=count_y; ++count){
	yc[count] = y_Eigen(count);
 	cout << yc[count] << endl;
	}
cout << "" << endl;

	for(count = 0; count!=count_z; ++count){
	zc[count] = z_Eigen(count);
	cout << zc[count] << endl;
	}

}
