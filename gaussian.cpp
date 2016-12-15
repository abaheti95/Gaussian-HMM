#include "gaussian.hpp"

Gaussian::Gaussian(int n_dimension, vector<double> &mu_inp, vector<double> &sigma_inp)
{
	k = n_dimension;
	update_parameters(mu_inp, sigma_inp);
}

void Gaussian::update_mu(vector<double> &mu_inp)
{
	mu = mu_inp;
}

void Gaussian::update_sigma(vector<double> &sigma_inp)
{
	sigma = sigma_inp;
	sigma_inv = sigma;
	det_sigma = 1.0;
	for(int i = 0; i < k; i++) {
		det_sigma *= sigma[i];
		sigma_inv = 1.0 / sigma_inv[i];
	}
	gaussain_denominator = pow(2*pi, k/2.0) * sqrt(det_sigma);
}

void Gaussian::update_parameters(vector<double> &mu_inp, vector<double> &sigma_inp)
{
	update_mu(mu_inp);
	update_sigma(sigma_inp);
}

/*
 * Input: Set of k-dimensional data points which will model the gaussian
 * Compute: maximum likelihood estimates of mu and sigma for the given data points
*/
void Gaussian::tune_parameters(vector<vector<double> &X)
{
	int N_points = (int)X.size();
	// mu is just the mean of all data points
	// zero the mu
	vector<double> dummy_mu(k, 0.0);
	for(int i = 0; i < N_points; i++) {
		assert(X[i].size() == k && "Some data point has dimensional != k");
		for(int j = 0; j < k; j++) {
			dummy_mu[j] += X[i][j];
		}
	}
	// Normalize mu
	for(int j = 0; j < k; j++) {
		dummy_mu[j] /= (double)(N_points);
	}
	update_mu(dummy_mu);

	// sigma is simply the diagonal covariances 
	// zero the sigma
	vector<double> dummy_sigma(k, 0.0);
	for(int i = 0; i < N_points; i++) {
		for(int j = 0; j < k; j++) {
			dummy_sigma[i][j] += pow((X[i][j] - dummy_mu[j]), 2.0);
		}
	}
	// Normalize sigma
	for(int j = 0; j < k; j++) {
		dummy_sigma[j] /= (double)(N_points);
	}
	update_sigma(dummy_sigma);
}

double Gaussian::evaluate(double<vector> &x)
{
	double summation = 0.0;
	for(int i = 0; i < k; i++) {
		summation += pow(x[i] - mu[i], 2.0) * sigma_inv[i];
	}
	return exp(-0.5 * summation) / gaussain_denominator;
}


