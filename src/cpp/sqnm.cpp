/**
 * @file sqnm.cpp
 * @author Moritz Gubler (moritz.gubler@unibas.ch)
 * @brief Implementation of the SQNM method. More informations about the algorithm can be found here: https://aip.scitation.org/doi/10.1063/1.4905665
 * @date 2022-07-13
 * 
 */
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <memory>
#include "historylist.hpp"
#include "sqnm.hpp"

using namespace Eigen;
using namespace std;
using namespace hlist_space;
using namespace sqnm_space;

/**
 * @brief Construct a new SQNM::SQNM object using default parameters
 * 
 * @param ndim_ number of dimensions of target function
 * @param nhistx_ Maximal number of steps that will be stored in the history list. Use a value between 3 and 20. Must be <= than ndim_.
 * @param alpha_ initial step size. default is 1.0. For systems with hard bonds (e.g. C-C) use a value between and 1.0 and
 */
SQNM::SQNM(int ndim_, int nhistx_, double alpha_) {
  ndim = ndim_;
  nhistx = nhistx_;
  alpha = alpha_;
  xlist = make_unique<HistoryList>(ndim, nhistx);
  flist = make_unique<HistoryList>(ndim, nhistx);
}

/**
 * @brief Construct a new SQNM::SQNM object using custom parameters.
 * 
 * @param ndim_ number of dimensions of target function
 * @param nhistx_ Maximal number of steps that will be stored in the history list. Use a value between 3 and 20. Must be <= than ndim_.
 * @param alpha_ initial step size. default is 1.0. For systems with hard bonds (e.g. C-C) use a value between and 1.0 and
 * @param alpha0_  * @param alpha0 Lower limit on the step size. 1.e-2 is the default.
 * @param eps_subsp_ Lower limit on linear dependencies of basis vectors in history list. Default 1.e-4.
 */
SQNM::SQNM(int ndim_, int nhistx_, double alpha_, double alpha0_, double eps_subsp_) {
  ndim = ndim_;
  nhistx = nhistx_;
  alpha = alpha_;
  xlist = make_unique<HistoryList>(ndim, nhistx);
  flist = make_unique<HistoryList>(ndim, nhistx);
  alpha0 = alpha0_;
  eps_subsp = eps_subsp_;
}

/**
 * @brief Calculates new coordinates that are closer to local minimum that the current coordinates. This function should be used the following way:
 * 1. calculate f(x) and the derivative.
 * 2. call the step function.
 * 3. add return value of step function to x.
 * 4. repeat.
 * 
 * @param x Current position vector
 * @param f_of_x value of the target function evaluated at position x.
 * @param df_dx derivative of the target function evaluated at x.
 * @return VectorXd displacent that can be added to x in order to get new improved coordinates.
 */
VectorXd SQNM::step(VectorXd &x, double &f_of_x, VectorXd &df_dx) {
  nhist = xlist->add(x);
  flist->add(df_dx);
  if (nhist == 0) { // initial and first step
    //cout << "first iteration of sqnm" << endl;
    this->dir_of_descent = - alpha * df_dx;

  } else {
    double gainratio = calc_gainratio(f_of_x);
    adjust_stepsize(gainratio);
    //cout << "gainratio, stepsize " << gainratio << " " << alpha << endl;
    //cout << "gainratio " << gainratio << " "<< alpha << endl;
    MatrixXd S = calc_ovrlp();
    //cout << "S " << S << endl;
    esolve.compute(S);
    VectorXd S_eval = esolve.eigenvalues();
    MatrixXd S_evec = esolve.eigenvectors();
    //std::cout << xlist->normalized_difflist(0,nhist-1) << "\n";

    // compute eq 8
    //println("eq. 8");
    int dim_subsp = 0;
    for (int i = 0; i < S_eval.size(); i++){
      if (S_eval(i) / S_eval(0) > eps_subsp)
      {
        dim_subsp+=1;
      }
    }
    // remove dimensions from subspace
    for (int i = 0; i < dim_subsp; i++)
    {
      S_evec.col(i) = S_evec.col(nhist - dim_subsp + i);
      S_eval(i) = S_eval(nhist - dim_subsp + i);
    }
    
    MatrixXd dr_subsp(ndim, dim_subsp);
    dr_subsp.setZero();
    //cout << "dimsp " << dim_subsp << " nhist " << nhist << endl;
    //println("start loop");
    for (int i = 0; i < dim_subsp; i++) {
      for (int ihist = 0; ihist < nhist; ihist++){
        dr_subsp.col(i) += S_evec(ihist, i) * xlist->normalized_difflist.col(ihist);
      }
      dr_subsp.col(i) /= sqrt(S_eval(i));
    }

    // compute eq. 11
    //println("eq. 11");
    MatrixXd df_subsp(ndim, dim_subsp);
    df_subsp.setZero();
    for (int i = 0; i < dim_subsp; i++) {
      for (int ihist = 0; ihist < nhist; ihist++){
        df_subsp.col(i) += S_evec(ihist, i) * flist->difflist.col(ihist) / xlist->difflist.col(ihist).norm();
      }
      df_subsp.col(i) /= sqrt(S_eval(i));
    }
    // compute eq. 13
    //println("eq. 13");
    h_subsp = .5 * (df_subsp.transpose() * dr_subsp + dr_subsp.transpose() * df_subsp);
    esolve.compute(h_subsp);
    h_eval = esolve.eigenvalues();
    h_evec_subsp = esolve.eigenvectors();

    // compute eq. 15
    //println("eq. 15");
    h_evec.resize(ndim, dim_subsp);
    h_evec.setZero();
    for (int i = 0; i < dim_subsp; i++){
      for (int k = 0; k < dim_subsp; k++){
        h_evec.col(i) += h_evec_subsp(k, i) * dr_subsp.col(k);
      }
    }

    // compute residues (eq. 20)
    //println("eq. 20");
    res.resize(dim_subsp);
    for (int j = 0; j < dim_subsp; j++){
      res_temp = - h_eval(j) * h_evec.col(j);
      for (int k = 0; k < dim_subsp; k++){
        res_temp += h_evec_subsp(k, j) * df_subsp.col(k);
      }
      res(j) = res_temp.norm();
    }

    // modify eigenvalues (eq. 18)
    //println("eq. 18");
    for (int idim = 0; idim < dim_subsp; idim++){
      h_eval(idim) = sqrt(pow(h_eval(idim), 2) + pow(res(idim), 2));
    }
    
    // decompose gradient (eq. 16)
    //println("eq. 16");
    dir_of_descent = df_dx;
    for (int i = 0; i < dim_subsp; i++){
      //cout << h_evec.col(0).size() <<  " " << df_dx.size() << " " << dir_of_descent.size() << endl;
      //cout << h_evec << endl;
      //cout << h_evec.col(i).dot(df_dx) << endl;
      dir_of_descent -= h_evec.col(i).dot(df_dx) * h_evec.col(i);
    }
    dir_of_descent *= alpha;

    // appy preconditioning to subspace gradient (eq. 21)
    //println("eq. 21");
    //dir_of_descent.setZero();
    for (int idim = 0; idim < dim_subsp; idim++)
    {
      dir_of_descent += (df_dx.dot(h_evec.col(idim)) / h_eval(idim)) * h_evec.col(idim);
    }
    dir_of_descent *= -1.0;
    
  }
  prev_f = f_of_x;
  prev_df_dx = df_dx;
  return this->dir_of_descent;
}

double SQNM::calc_gainratio(double &f){
  //cout << "de " << f - prev_f << endl;
  double gr = (f - prev_f) / ( .5 * this->dir_of_descent.dot(prev_df_dx));
  return gr;
}

void SQNM::adjust_stepsize(double &gainratio){
  if ( gainratio < 0.5 ) alpha = max(alpha * 0.65, alpha0);
  else if(gainratio > 1.05) alpha = alpha * 1.05;
}


MatrixXd SQNM::calc_ovrlp(){
  MatrixXd S = xlist->normalized_difflist.block(0,0, ndim, nhist).transpose() * xlist->normalized_difflist.block(0,0, ndim, nhist);
  return S;
}
