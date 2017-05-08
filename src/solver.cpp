/*
  Copyright ©2017 The Regents of the University of California
  (Regents). All Rights Reserved. Permission to use, copy, modify, and
  distribute this software and its documentation for educational,
  research, and not-for-profit purposes, without fee and without a
  signed licensing agreement, is hereby granted, provided that the
  above copyright notice, this paragraph and the following two
  paragraphs appear in all copies, modifications, and
  distributions. Contact The Office of Technology Licensing, UC
  Berkeley, 2150 Shattuck Avenue, Suite 510, Berkeley, CA 94720-1620,
  (510) 643-7201, for commercial licensing opportunities.

  IN NO EVENT SHALL REGENTS BE LIABLE TO ANY PARTY FOR DIRECT,
  INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING
  LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS
  DOCUMENTATION, EVEN IF REGENTS HAS BEEN ADVISED OF THE POSSIBILITY
  OF SUCH DAMAGE.

  REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
  FOR A PARTICULAR PURPOSE. THE SOFTWARE AND ACCOMPANYING
  DOCUMENTATION, IF ANY, PROVIDED HEREUNDER IS PROVIDED "AS
  IS". REGENTS HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT,
  UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
*/

#include "solver.hpp"
#include "timer.hpp"

#include <Eigen/CholmodSupport>
// #include <Eigen/SparseLU>

std::vector<double> eigen_linear_solve(const SpMat<double>& A, const std::vector<double>& b)
{
    TimerGuard tg("eigen_linear_solve");
    typedef Eigen::Triplet<double> ET;

    // assert the dimension is correct
    size_t nn = A.n;
    assert(A.n == A.m);
    assert(nn == b.size());

    Eigen::SparseMatrix<double> eA(nn, nn);
    Eigen::Map<const Eigen::VectorXd> eb(b.data(), nn);
    std::vector<ET> coeff;

    // build the eigen matrix using triplet
    for (size_t i = 0; i < nn; ++i) {
        size_t m = A.rows[i].indices.size();
        for (size_t jj = 0; jj < m; ++jj) {
            size_t j = A.rows[i].indices[jj];
            double entry = A.rows[i].entries[jj];
            coeff.push_back(ET(i, j, entry));
        }
    }
    eA.setFromTriplets(coeff.begin(), coeff.end());

    Eigen::CholmodSimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(eA);

    // Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    // solver.analyzePattern(eA);
    // solver.factorize(eA);

    Eigen::VectorXd x = solver.solve(eb);
    return std::vector<double>(x.data(), x.data() + nn);
}
