// programmingAssignment2.cpp : This file contains the 'main' function. Program execution begins and ends there.

#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <Eigen/Core>

using namespace std;
using namespace Eigen; 

// constants
#define H  0.1 
#define U 2



// initalize the array
template <typename Derived>
void initalizeArray( DenseBase<Derived>& a, int size) {
	for (int i = 0; i < size; i++) {
		a(0, i) = sin(M_PI * (i + 0.5) * H);
	}
}
//index shifting macro as recommended by Carl
int indexShift(int N, int index) {
	int newIndex = (index + N) % N;
	return newIndex; 
}

// UW2 flux integral evaulation
template <typename Derived>
void UW2_Flux_Integral(DenseBase<Derived>& data, DenseBase<Derived>&  result, int size) {
	for (int i = 0; i < size; i++) {
		int indexI = indexShift(size,i);
		int indexIM1 = indexShift(size, i - 1);
		int indexIM2 = indexShift(size, i - 2);

		double tI = data[indexI];
		double tIM1 = data[indexIM1];
		double tIM2 = data[indexIM2];

		double CVAverage = ( 1.0 ) / (2 * H) * (3 * tI - 4 * tIM1 + tIM2);
		result(0,i) = CVAverage;
	}
}

//UWB3 flux integral evaulation
template <typename Derived>
void UWB3_Flux_Integral(DenseBase<Derived>& data, DenseBase<Derived>& result, int size) {
	for (int i = 0; i < size; i++) {
		int indexI = indexShift(size, i);
		int indexIP1 = indexShift(size, i + 1);
		int indexIM1 = indexShift(size, i - 1);
		int indexIM2 = indexShift(size, i - 2);

		double tI = data[indexI];
		double tIP1 = data[indexIP1];
		double tIM1 = data[indexIM1];
		double tIM2 = data[indexIM2];

		double CVAverage = (1.0 ) / (6 * H) * (2*tIP1+ 3*tI- 6*tIM1 + tIM2);
		result(0, i) = CVAverage;
	 }
}

int main()
{
	const int numCells = 1.0 / H;
	Eigen::Matrix<double, 1, numCells> cellArray;
	Eigen::Matrix<double, 1, numCells> CVArray;

	// initializing the array
	initalizeArray(cellArray, numCells);
	std::cout << cellArray << endl;

	// cv for UW2
	UW2_Flux_Integral(cellArray, CVArray, numCells);
	cout << CVArray << endl;




}
