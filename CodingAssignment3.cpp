// CodingAssignment3.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <Eigen>
using namespace std;

//helper functions
double itoY(double i, double dy, double size) {
	double Y = (size - i + 0.5) * dy;
	return Y;
}

double jtoX(double j, double dx, double size) {
	double X = (j - 0.5) * dx;
	return X;
}

double U_PartA(double x, double y) {
	double U0 = 1.0;
	double u = U0 *y* sin(M_PI * x);
	return u;
}

double V_PartA(double x, double y) {
	double V0 = 1.0;
	double v = V0 * x * cos(M_PI * y);
	return v;
}

int initializeCV(Eigen::MatrixXd& T,int size){

}

// calculate the flux integral of the matrix
// think in the coordinate system of the matrix, just take care to convert between the mesh and matrix definition. 
int EnergyFluxIntegral(Eigen::MatrixXd& T, Eigen::MatrixXd& FluxInts, double Re, double Pr, double dx, double dy, int size,
	double (*U)(double x,double y), double(*V)(double x, double y), double (*itoY)(double i, double dy, double size), double (*jtoX)(double j, double dx, double size)){
	// assumes we are using ghost cells 
	for (int i = 1; i < size - 1; i++)
	{
		// calculating the y values 
		double y = itoY(i, dy, size);
		double yAbove = y + dy;
		double yBelow = y - dy;

		for (int j = 1; j < size - 1; j++) {
			//calculating the x -coordinates
			double x = jtoX(j, dx, size);
			double xRight = x + dx; 
			double xLeft = x - dx;

			// getting all the cell values that I need
			double currentCell = T(i, j);
			double aboveCell = T(i - 1, j);
			double belowCell = T(i + 1, j);
			double rightCell = T(i, j + 1);
			double leftCell = T(i, j - 1);
			
			//calculating velocities
			double u_Right = U(xRight, y);
			double u_Left = U(xLeft, y);
			double v_Above = V(x, yAbove);
			double v_Below = V(x, yBelow);

			//seperatly calculating all the flux integral terms
			double dxTerm = -1.0/(2 * dx) * (u_Right * rightCell - u_Left * leftCell);
			double dx2Term = 1.0/(Re * Pr* dx*dx) * (rightCell - 2*currentCell + leftCell);
			double dyTerm = -1.0 / (2 * dy) * (v_Above * aboveCell - v_Below * belowCell);
			double dy2Term = 1.0 / (Re * Pr * dy * dy) * (aboveCell - 2 * currentCell + belowCell);

			FluxInts(i, j) = dxTerm + dx2Term + dyTerm + dy2Term;


		}
	}
	  
	return 1;
}


int main()
{
	int cellLength = 10+2;// plus to account for ghost cells.   
	Eigen::MatrixXd T(cellLength, cellLength);

    std::cout << "Hello World!\n";
}


