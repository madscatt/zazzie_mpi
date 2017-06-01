#include <iostream>
#include <vector>
#include <string>
#include <fstream> 

using namespace std;

int main(){

cout << "Hello World" << endl;

const double PI = 3.1415926535;

char myGrade = 'A';

bool isHappy = true;

int myAge = 1;

float favNum = 3.14159;

double otherfavNum = 1.23255436565;

cout << "Favorite Number: " << favNum << endl;

cout << "Size of int " << sizeof(myAge) <<endl;

int five = 5;

cout << "5++ = " << five++ << endl;
cout << "++5 = " << ++five << endl;
cout << "5-- = " << five-- << endl;
cout << "--5 = " << --five << endl;

int age = 70;
int ageAtLastExam = 16;
bool isNotIntoxicated = true;

if((age >=1) && (age < 16)){
	
	cout << "You can't drive" << endl;

} else if(! isNotIntoxicated){

	cout << "You can't drive" << endl;

} else if(age >= 80 && ((age > 100) || ((age - ageAtLastExam) > 5))){

	cout << "You can't drive" << endl;

} else {
	
	cout << "You can drive" << endl;

}

int greetingOption = 2;

switch(greetingOption){

	case 1:
		cout << "bonjour" <<endl;
		break;
	case 2:
		cout << "hola" << endl;
		break;
	case 3: cout << "hello" << endl;
		break;

}

int mynums[5];

int nonums[5] = {1, 2, 3, 4, 5};

cout << "nonum1: " << nonums[0] << endl;

char myName[2][5] = {{'d', 'e', 'r', 'e', 'k'}, {'a','b','c','d','e'}};

cout<< "a" << myName[1][0] <<endl;

string numberGuessed;
int intNumberGuessed = 0;

do {

	cout << "Guess between 1 and 10: " << endl;
	
	getline(cin, numberGuessed);

	intNumberGuessed = stoi(numberGuessed);
	
	cout << intNumberGuessed << endl;

} while(intNumberGuessed != 4);

cout << "you win" << endl;


return 0;


}




































