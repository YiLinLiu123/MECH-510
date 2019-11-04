// programmingAssignment2.cpp : This file contains the 'main' function. Program execution begins and ends there.

#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <Eigen/Core>

using namespace std;
using namespace Eigen; 

// constants
const double H = 0.025; //step value in space, I used h cause thats what we use in math classes. 
const double U = 2.0; 



void displayVector(RowVectorXd& a, int size) {
	for (int i = 0; i < size; i++) {
		cout << a[i] << ", ";
	}
	cout << endl;
}


// initalize the array
void initalizeArray(RowVectorXd& a, int size) {
	for (int i = 0; i < size; i++) {
		//double x = (i + 0.5) * H - 1;
		//a(i) = sin(M_PI * x);
		//proper flux integral
		double x_start = (i * H - 1);
		double x_end = x_start + H;
		a(i) = (-cos(M_PI * (x_end)) + cos(M_PI * x_start)) / (M_PI * H);
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

		double term1 = 3.0 * tI;
		double term2 = tIM1 * 4.0;
		double constant = -1.0*U / (2.0 * H);
		double CVAverage =  constant*(term1-term2+tIM2);
		result(i) = CVAverage;
	}
}

// C2 flux integral evaulatio
void C2_Flux_Integral(RowVectorXd& data, RowVectorXd& result, int size) {
	for (int i = 0; i < size; i++) {
		int indexIM1 = indexShift(size, i-1);
		int indexIP1 = indexShift(size, i+1);

		double tIM1 = data(indexIM1);
		double tIP1 = data(indexIP1);

		double CVAverage = (-1.0*U) / (2.0 * H) * (tIP1- tIM1);
		result(i) = CVAverage;
	}
}

//UWB3 flux integral evaulation
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

		
		double term1 = 2 * tIP1;
		double term2 = 3* tI;
		double term3 = 6 * tIM1;
		double constant = -1.0*U / (6.0 * H);

		double CVAverage = constant * (2.0*tIP1+ 3.0*tI- 6.0*tIM1 + tIM2);
		double CVAverage2 = constant * (term1 + term2 - term3 + tIM2);
		result(i) = CVAverage2;
	 }
}


// compute flux L2 norm
double flux_L2_Norm(RowVectorXd& CV, int size, RowVectorXd & errorVector, RowVectorXd& correctVector, double (*correctValueFunction)(double x) ){
	double cumError = 0; 
	for (int i = 0; i < size; i++) {
		double x = (i + 1.0 / 2.0) * H-1;
		//double x = (i) * H - 1;
		double correct = correctValueFunction(x);
		double guess = CV(i);
		double error = guess - correct;
		errorVector(i) = error;
		correctVector(i) = correct; 
		//printf("error: %f \n", error);
		cumError += error*error; 
		cumError += 0;
	}
	double L2Norm = pow(cumError / size, 0.5);
	return L2Norm;
}

void RK2(RowVectorXd& data, RowVectorXd& updatedResult, void (*flux_Integral)(RowVectorXd& data, RowVectorXd& result, int size), int size, double timeStep)
{
	// vectors 
	Eigen::RowVectorXd fluxIntegral (size);
	Eigen::RowVectorXd fluxIntegral1(size);
	RowVectorXd w_1(size);

	//steps:
	flux_Integral(data, fluxIntegral, size);
	//cout << "W_n Flux: " << fluxIntegral << endl;
	w_1 = data + timeStep/2.0 * fluxIntegral;  
	//cout << "W_1: " << w_1 << endl;

	flux_Integral(w_1, fluxIntegral1, size);
	//cout << "W_1 Flux: " << fluxIntegral1 << endl;
	updatedResult = data + timeStep * (fluxIntegral1);
}




void RK4(RowVectorXd& data, RowVectorXd& updatedResult, void flux_Integral(RowVectorXd& data, RowVectorXd& result, int size), int size, double timeStep)
{
	// vectors
	Eigen::RowVectorXd fluxIntegral(size);
	Eigen::RowVectorXd w_1(size);
	Eigen::RowVectorXd fluxIntegral1(size);
	Eigen::RowVectorXd w_2(size);
	Eigen::RowVectorXd fluxIntegral2(size);
	Eigen::RowVectorXd w_3(size);
	Eigen::RowVectorXd fluxIntegral3(size);

	// steps
	flux_Integral(data, fluxIntegral, size);
	w_1 = data + timeStep / 2.0 * fluxIntegral;

	flux_Integral(w_1, fluxIntegral1, size);
	w_2 = data + timeStep / 2.0 * fluxIntegral1;


	flux_Integral(w_2, fluxIntegral2, size);
	w_3 = data + timeStep * fluxIntegral2;

	flux_Integral(w_3, fluxIntegral3, size);
	updatedResult = data + timeStep / 6.0 * (fluxIntegral + 2.0 * fluxIntegral1 + 2.0 * fluxIntegral2 + fluxIntegral3);

}

// made the same thing except has a lot of print statements for debugging. 
void RK4Debug(RowVectorXd& data, RowVectorXd& updatedResult, void flux_Integral(RowVectorXd& data, RowVectorXd& result, int size), int size, double timeStep)
{
	// vectors
	Eigen::RowVectorXd fluxIntegral(size);
	Eigen::RowVectorXd w_1(size);
	Eigen::RowVectorXd fluxIntegral1(size);
	Eigen::RowVectorXd w_2(size);
	Eigen::RowVectorXd fluxIntegral2(size);
	Eigen::RowVectorXd w_3(size);
	Eigen::RowVectorXd fluxIntegral3(size);

	// steps
	flux_Integral(data, fluxIntegral, size);
	cout << "flux integral: " <<fluxIntegral << endl;
	w_1 = data + timeStep / 2.0 * fluxIntegral;
	cout << "w1: " << endl << w_1 << endl;

	flux_Integral(w_1, fluxIntegral1, size);
	cout << "flux integral1: " << endl << fluxIntegral1 << endl;
	w_2 = data + timeStep / 2.0 * fluxIntegral1;
	cout << "w2: " << endl << w_2 << endl;


	flux_Integral(w_2, fluxIntegral2, size);
	cout << "flux integral2: " << endl << fluxIntegral2 << endl;
	w_3 = data + timeStep * fluxIntegral2;
	cout << "w3: " << endl << w_3 << endl;

	flux_Integral(w_3, fluxIntegral3, size);
	cout << "flux integral3: " << endl << fluxIntegral3 << endl;
	updatedResult = data + timeStep / 6.0 * (fluxIntegral + 2.0 * fluxIntegral1 + 2.0 * fluxIntegral2 + fluxIntegral3);

}

double partOneRightValue(double x) {
	return -1.0*U * M_PI * cos(M_PI * x);
}

void partOne()
{
	const int numCells = 2.0 / H;
	RowVectorXd cellArray(numCells);

	RowVectorXd CV_UW2_Array(numCells);
	RowVectorXd CV_C2_Array(numCells);
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
	printf("L2 Norm for UW2 where N = %d: %f \n", numCells, flux_L2_Norm(CV_UW2_Array, numCells, errorVector, correctVector, partOneRightValue));
	displayVector(errorVector, numCells);

	// cv for C2
	C2_Flux_Integral(cellArray, CV_C2_Array, numCells);
	printf("CV_C2_Values: \n");
	cout << CV_C2_Array << endl;
	printf("L2 Norm for C2 where N = %d: %f \n", numCells, flux_L2_Norm(CV_C2_Array, numCells, errorVector, correctVector, partOneRightValue));
	displayVector(errorVector, numCells);
}


void exactFluxIntegral(RowVectorXd& data, RowVectorXd& result, int size, double x) {
	// we know the exact solution
	// solve over 0 to .5
	result(0) = -U*M_PI*cos(M_PI*x);

}


double partThreeRightValue(double x) {
	return  sin(M_PI * x);
}
void partThree_RK2_And_UW2() {
	double targetTime = 2.0; 
	double timeStep = 0.05/(0.2/H)*0.8; 
	//double timeStep = 0.006*0.8; // for H = 0.025
	int steps = ceil(2.0 / timeStep);
	timeStep = 2.0 / steps;
	const int numCells = 2.0 / H;
	double time = 0; 
	
	RowVectorXd cellArray(numCells);
	RowVectorXd correct(numCells);
	RowVectorXd error(numCells);

	initalizeArray(cellArray,numCells);
	cout << "time step: " << timeStep << endl;
	cout << "dx: " << H << endl; 
	cout << "initial values: " <<cellArray << endl;

	for (int i = 0; i < steps; i++) {
		time = time +timeStep;
		RK2(cellArray, cellArray, UW2_Flux_Integral, numCells,timeStep);
	}
	cout << "final time: " << time << endl;
	cout << "final iteration: " << endl;
	displayVector(cellArray, numCells);
	cout <<"error norm : "<<flux_L2_Norm(cellArray, numCells, error, correct, partThreeRightValue) << endl;
	cout << "error vector : " << endl; 
	displayVector(error, numCells);
	//cout << "correct: "<< correct << endl;


}

void partThree_RK4_And_C2() {
    double timeStep = 0.28/(.2/H)*0.8; // for C2
	//double timeStep = 0.05 / (0.2 / H) * 0.8; // for UW2
	// double timeStep = 0.16 / (0.2 / H) * 0.8; // for UWB3
	//double constant = 0.7071067;
	//double timeStep = 0.01*constant*constant;
	double CFL = timeStep / H * U;
	int steps = ceil(2.0/timeStep);// happens that we need to ceil it
	timeStep = 2.0 / steps; 
	const int numCells = 2.0 / H;
	double time = 0;

	RowVectorXd cellArray(numCells);
	RowVectorXd correct(numCells);
	RowVectorXd error(numCells);

	initalizeArray(cellArray, numCells);
	cout << "time step: " << timeStep << endl;
	cout << "dx: " << H << endl;
	cout << "CFL: " << CFL << endl;
	cout << "initial values: " << cellArray << endl;

	for (int i = 0; i <steps; i++) {
		time = time + timeStep;
		RK4(cellArray, cellArray, C2_Flux_Integral, numCells, timeStep);
	}
	cout << "final time: " << time << endl;
	cout << "final iteration: " << endl;
	displayVector(cellArray, numCells);
	cout << "error norm : " << flux_L2_Norm(cellArray, numCells, error, correct, partThreeRightValue) << endl;
	cout << "error vector : " << endl;
	displayVector(error, numCells);
	cout << "correct: " << correct << endl;


}

int main()
{
	partOne();
	partThree_RK2_And_UW2();
	partThree_RK4_And_C2();

}
