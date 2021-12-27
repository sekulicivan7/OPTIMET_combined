// (C) University College London 2017
// This file is part of Optimet, licensed under the terms of the GNU Public License
//
// Optimet is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Optimet is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Optimet. If not, see <http://www.gnu.org/licenses/>.

#ifndef OPTIMET_PRECONDITIONNED_MATRIX_SOLVER_H
#define OPTIMET_PRECONDITIONNED_MATRIX_SOLVER_H

#include "PreconditionedMatrix.h"
#include "Solver.h"
#include "Types.h"
#include <Eigen/Dense>


namespace optimet {
namespace solver {

//! Use an actual matrix, and Eigen's Householder QR method
class PreconditionedMatrix : public AbstractSolver {
public:
  PreconditionedMatrix(std::shared_ptr<Geometry> geometry,
                       std::shared_ptr<Excitation const> incWave,
                       mpi::Communicator const &communicator = mpi::Communicator())
      : AbstractSolver(geometry, incWave, communicator) {
  #ifndef OPTIMET_MPI
    update();
   #endif
  }


  PreconditionedMatrix(Run const &run)
      : PreconditionedMatrix(run.geometry, run.excitation, run.communicator) {}

  void solve(Vector<t_complex> &X_sca_, Vector<t_complex> &X_int_, Vector<t_complex> &X_sca_SH,

              Vector<t_complex> &X_int_SH, std::vector<double *> CGcoeff) const override {
  
    // fundamental frequency
    
    Matrix<t_complex> TmatrixFF, RgQmatrixFF; 
    int nMax = geometry->nMax();
    int pMax = nMax * (nMax + 2);
    TmatrixFF = S.block(0 ,0 , 2*pMax, 2*pMax);
    RgQmatrixFF = S.block(0 , 2*pMax , 2*pMax, 2*pMax);
 
    X_sca_ = Q;
    
    unprecondition(X_sca_, X_int_, TmatrixFF, RgQmatrixFF);

    Vector<t_complex> K, K1;
    Matrix<t_complex> TmatrixSH, RgQmatrixSH;
    
    // SH frequency
    if(incWave->SH_cond){

    int nMaxS = geometry->nMaxS();
    int pMax = nMaxS * (nMaxS + 2);
    TmatrixSH = V.block(0 ,0 , 2*pMax, 2*pMax);
    RgQmatrixSH = V.block(0 , 2*pMax , 2*pMax, 2*pMax);

    K = source_vectorSH(*geometry, incWave, X_int_, X_sca_, TmatrixSH);

    K1 = source_vectorSHarb1(*geometry, incWave, X_int_, X_sca_);

    X_sca_SH = K;
    
    unprecondition_SH(X_sca_SH, X_int_SH, K1, RgQmatrixSH);
    }
  }
  
    

  void update() override {
  
    Q = source_vector(*geometry, incWave);

    S = getTRgQmatrix_FF_parr(*geometry, incWave);
    if(incWave->SH_cond)
    V = getTRgQmatrix_SH_parr(*geometry, incWave);

  }
  

protected:
  
  Matrix<t_complex> S;
  
  Vector<t_complex> Q;
  
  Matrix<t_complex> V;
  
  Vector<t_complex> K;

  //! Unpreconditions the result of preconditioned computation
  void unprecondition(Vector<t_complex> &X_sca_, Vector<t_complex> &X_int_, Matrix<t_complex> &Tmat, Matrix<t_complex> &RgQ) const {
    X_sca_ = AbstractSolver::convertIndirect(X_sca_, Tmat);
    X_int_ = AbstractSolver::solveInternal(X_sca_, RgQ);
    }
     
    
    void unprecondition_SH(Vector<t_complex> &X_sca_SH, Vector<t_complex> &X_int_SH, Vector<t_complex> &K1, Matrix<t_complex> &RgQ) const {
    X_sca_SH = AbstractSolver::convertIndirect_SH_outer(X_sca_SH);
    X_int_SH = AbstractSolver::solveInternal_SH(X_sca_SH, K1, RgQ);
    
    }

   
};
}
}

#endif
