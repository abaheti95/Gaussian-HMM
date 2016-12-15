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

double Gaussian::evaluate(double<vector> &x)
{
	double summation = 0.0;
	for(int i = 0; i < k; i++) {
		summation += pow(x[i] - mu[i], 2.0) * sigma_inv[i];
	}
	return exp(-0.5 * summation) / gaussain_denominator;
}


