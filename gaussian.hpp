#ifndef GAUSSIAN_HPP
#define GAUSSIAN_HPP

// C Includes
#include <cmath>

// CPP Includes
#include <vector>

using namespace std;

const double pi = 3.141592653589793;

// Multivariate Gaussain function
class Gaussian
{
private:
public:
	// Member variables
	int k;								// Number of dimensions that the multivariate gaussian takes as input
	double gaussain_denominator;		// stores the value of (2*pi)^(k/2)
	vector<double> mu;					// mean of the gaussian
	vector<double> sigma;				// diagonal values of covariance matrix of the gaussian
	vector<double> sigma_inv;			// Inverse of the diagonal sigma
	double det_sigma;					// Determinant of the sigma

	Gaussian(int n_dimension, vector<double> &mu_inp, vector<double> &sigma_inp);
	~Gaussian() {}
	// Member Functions
	void update_mu(vector<double> &mu_inp);
	void update_sigma(vector<double> &sigma_inp);
	void update_parameters(vector<double> &mu_inp, vector<double> &sigma_inp);
	double evaluate(double<vector> &x);	// Output of the gaussian for the input vector x
};

#endif
