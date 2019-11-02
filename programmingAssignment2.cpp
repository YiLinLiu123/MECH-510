// programmingAssignment2.cpp : This file contains the 'main' function. Program execution begins and ends there.

#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <Eigen/Core>

using namespace std;
using namespace Eigen; 

// constants
#define H  0.025
#define U 2.0



// initalize the array
void initalizeArray(RowVectorXd& a, int size) {
	for (int i = 0; i < size; i++) {
		double x = (i + 0.5) * H - 1;
		a(i) = sin(M_PI * x);
		//proper flux integral
		//double x_start = (i * H - 1);
		//a(i) = (-cos(M_PI * (x_start+H)) + cos(M_PI * x_start)) / (M_PI * H);
	}
}
//index shifting macro as recommended by Carl
int indexShift(int N, int index) {
	int newIndex = (index + N) % N;
	return newIndex; 
}

// UW2 flux integral evaulation

void UW2_Flux_Integral(RowVectorXd& data, RowVectorXd&  result, int size) {
	for (int i = 0; i < size; i++) {
		int indexI = indexShift(size,i);
		int indexIM1 = indexShift(size, i - 1);
		int indexIM2 = indexShift(size, i - 2);

		double tI = data(indexI);
		double tIM1 = data(indexIM1);
		double tIM2 = data(indexIM2);

		double CVAverage = ( U ) / (2.0 * H) * (3.0 * tI - 4.0 * tIM1 + tIM2);
		result(i) = CVAverage;
	}
}

// C2 flux integral evaulation

void C2_Flux_Integral(RowVectorXd& data, RowVectorXd& result, int size) {
	for (int i = 0; i < size; i++) {
		int indexI = indexShift(size, i);
		int indexIM1 = indexShift(size, i - 1);
		int indexIP1 = indexShift(size, i +1);

		double tI = data(indexI);
		double tIM1 = data(indexIM1);
		double tIP1 = data(indexIP1);

		double CVAverage = (U) / (2.0 * H) * (tIP1- tIM1);
		result(i) = CVAverage;
	}
}

//UWB3 flux integral evaulation
// note i am calculating U*(Ti+1/2- Ti-1/2)/(deltaX)
void UWB3_Flux_Integral(RowVectorXd& data, RowVectorXd& result,  int size) {
	for (int i = 0; i < size; i++) {
		int indexI = indexShift(size, i);
		int indexIP1 = indexShift(size, i + 1);
		int indexIM1 = indexShift(size, i - 1);
		int indexIM2 = indexShift(size, i - 2);

		double tI = data[indexI];
		double tIP1 = data[indexIP1];
		double tIM1 = data[indexIM1];
		double tIM2 = data[indexIM2];

		double CVAverage = (U) / (6.0 * H) * (2.0*tIP1+ 3.0*tI- 6.0*tIM1 + tIM2);
		result(i) = CVAverage;
	 }
}


// compute flux L2 norm
double flux_L2_Norm(RowVectorXd& CV, int size, RowVectorXd & errorVector, RowVectorXd& correctVector) {
	double cumError = 0; 
	for (int i = 0; i < size; i++) {
		double x = (i + 1.0 / 2.0) * H-1;
		double correct = U*M_PI * cos(M_PI*x);
		double error = abs(CV[i] - correct);
		errorVector(i) = error;
		correctVector(i) = correct; 
		//printf("error: %f \n", error);
		cumError += error*error; 
	}
	double L2Norm = pow(cumError / size, 0.5);
	return L2Norm;
}

template <typename Derived>
void RK2(RowVectorXd& data, RowVectorXd& updatedResult, void flux_Integral(RowVectorXd& data, RowVectorXd& result, int size), int size)
{
	// calculate the space flux
	Eigen::RowVectorXd fluxIntegral (size);
	flux_Integral(data, fluxIntegral, size);

	// constants
	double timeStep = 2; // need to change this later
	RowVectorXd w_1(size);
	w_1 = data + timeStep * data; 
	updatedResult = data + timeStep / 2 * (w_1 + data);

}


void RK4(RowVectorXd& data, RowVectorXd& updatedResult, void flux_Integral(RowVectorXd& data, RowVectorXd& result, int size), int size)
{
	// calculate the space flux
	Eigen::RowVectorXd fluxIntegral (size);
	flux_Integral(data, fluxIntegral, size);

	// constants
	double timeStep = 2; // need to change this later
	Eigen::RowVectorXd w_1(size);
	Eigen::RowVectorXd w_2(size);
	Eigen::RowVectorXd w_3(size);
	w_1 = data + timeStep/2 * data;
	w_2 = data + timeStep / 2 *w_1;
	w_3 = data + timeStep / 2 * w_3;
	updatedResult = data + timeStep / 6 * (data + 2 * w_1 + 2 * w_2 + w_3);

}


int main()
{
	const int numCells = 2.0 / H;
	RowVectorXd cellArray(numCells);
	RowVectorXd CV_UW2_Array(numCells);
	RowVectorXd CV_UWB3_Array(numCells);
	RowVectorXd errorVector(numCells); 
	RowVectorXd correctVector(numCells);

	// initializing the array
	initalizeArray(cellArray, numCells);
	printf("Initial Values: \n");
	std::cout << cellArray << endl;

	// cv for UW2
	UW2_Flux_Integral(cellArray, CV_UW2_Array, numCells);
	printf("CV_UW2_Values: \n");
	cout << CV_UW2_Array << endl;
	printf("L2 Norm for UW2 where N = %d: %f \n", numCells, flux_L2_Norm(CV_UW2_Array, numCells, errorVector, correctVector));
	cout << errorVector << endl; 

	// cv for UWB3
	UWB3_Flux_Integral(cellArray, CV_UWB3_Array, numCells);
	printf("CV_UWB3_Values: \n");
	cout << CV_UWB3_Array << endl;
	printf("L2 Norm for UWB3 where N = %d: %f \n", numCells, flux_L2_Norm(CV_UWB3_Array, numCells, errorVector, correctVector));
	cout << errorVector << endl;
	printf("CorrectVector \n");
	cout << correctVector << endl;
}
