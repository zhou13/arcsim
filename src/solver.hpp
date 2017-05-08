#ifndef SOLVER_HPP
#define SOLVER_HPP

#include "sparse.hpp"
#include "vectors.hpp"

std::vector<double> eigen_linear_solve(const SpMat<double>& A, const std::vector<double>& b);

#endif
