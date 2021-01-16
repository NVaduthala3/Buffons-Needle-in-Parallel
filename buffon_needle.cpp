#include <iostream>
#include <random>
#include <chrono>
#include <cmath>
#include<iomanip>
#include "mpi.h""

const double pi_real = 2 * acos(0.0);
using namespace std;
double random_num(double a, double b);
int buffon_needle_2D(double X, double Y, double len, int needle_amount);

int main(int argc, char** argv)
{
    double X = 1.0, Y = 1.0, len = 0.9; //len should be <=a and b
    int crossed_count, crossed_total;
    int master_branch = 0;             //master process will calculate pi approximation
    double prob, error, pi_approx;
    int pr_number, pr_rank;
    int needle_amount = 100000, needle_total;

    MPI::Init(argc, argv);

    pr_number = MPI::COMM_WORLD.Get_size();
    pr_rank = MPI::COMM_WORLD.Get_rank();
    
    if (pr_rank == master_branch)
    {
        cout << "The number of processes is: " << pr_number << endl;
        cout << "The cell width is:    " << X << endl;
        cout << "The cell length is:" << Y << endl;
        cout << "The needle length is: " << len << endl;
    }

    crossed_count = buffon_needle_2D(X, Y, len, needle_amount);

    MPI::COMM_WORLD.Reduce(&crossed_count, &crossed_total, 1, MPI::INT, MPI::SUM, master_branch);
    
    if (pr_rank == master_branch)
    {
        needle_total = pr_number * needle_amount;
        prob = static_cast<double>(crossed_total) / static_cast<double>(needle_total);
        if (crossed_total == 0)
        {
            pi_approx = 0; //avoids "divide-by-zero" error
        }
        else
        {
            pi_approx = (2.0*len*(X + Y)-len*len)/(X*Y*prob);
        }
        error = abs(pi_approx - pi_real);

        cout << endl;
        cout << "The number of trials is: " << needle_total << endl;
        cout << "The number of crossed needles is: " << crossed_total << endl;
        cout << "The calculated probability is: " << prob << endl;
        cout << "The pi approximation value is: " << pi_approx << endl;
        cout << "The error between the real and approximate pi values is: " << error << endl;
    }

    MPI::Finalize();

    return 0;
}

double random_num(double a, double b)
{
	random_device ran;  
	mt19937 generator(ran()); 
	uniform_real_distribution<> dist(a, b);
	return dist(generator);
}

int buffon_needle_2D(double X, double Y, double len, int needle_amount)
{
	double theta;
	int crossed_count = 0;
	double x1, x2, y1, y2;

	for (int i= 1; i<= needle_amount; i++)
	{
		//Set location of one end of needle
		x1 = X * random_num(0.0, 1.0); 
		y1 = Y * random_num(0.0, 1.0); 
		theta = 2.0*pi_real*random_num(0.0, 1.0); 

		//Set location of other end of needle
		x2 = x1+len*cos(theta);
		y2 = y1+len*sin(theta);

		// Increment if end point crosses a line
		if ((x2 <= 0.0) || (y2 <= 0.0) || (X <= x2) || (Y <= y2)) //cell is the region [0, X] x [0, Y]
		{
			crossed_count = crossed_count + 1;
		}
	}

	return crossed_count;
}
