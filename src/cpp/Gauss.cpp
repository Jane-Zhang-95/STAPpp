#include "Gauss.h"

double GaussianQuadrature::weights_1p[1] = { 2.0 };
double GaussianQuadrature::points_1p[1] = { 0.0 };

double GaussianQuadrature::weights_2p[2] = { 1.0, 1.0 };
double GaussianQuadrature::points_2p[2] = { -0.57735026919, 0.57735026919 };

double GaussianQuadrature::weights_3p[3] = { 0.555555555556, 0.888888888889, 0.555555555556 };
double GaussianQuadrature::points_3p[3] = { -0.77459666924, 0.0, 0.77459666924 };