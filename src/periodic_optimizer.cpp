/**
 * @file periodic_optimizer.cpp
 * @author Moritz Gubler (moritz.gubler@unibas.ch)
 * @brief Implementation of the vc-sqnm method. More informations about the algorithm can be found here: https://arxiv.org/abs/2206.07339
 * @date 2022-07-13
 * 
 */

#include <eigen3/Eigen/Dense>
#include <iostream>
#include "sqnm.hpp"
#include "periodic_optimizer.hpp"

using namespace Eigen;
using namespace std;

/**
 * @brief Construct a new periodic optimizer::periodic optimizer object for fixed cell optimization with default parameters.
 * 
 * @param nat number of atoms
 */
periodic_optimizer::periodic_optimizer(int &nat)
{
  this->nat = nat;
  this->ndim = 3*nat;
  this->opt_lattice = false;
  this->opt = make_unique<sqnm_space::SQNM>(ndim, n_hist_max, initial_step_size);
}

/**
 * @brief Construct a new periodic optimizer::periodic optimizer object for fixed cell optimization with custom parameters.
 * 
 * @param nat number of atoms
 * @param initial_step_size initial step size. default is 1.0. For systems with hard bonds (e.g. C-C) use a value between and 1.0 and
 * 2.5. If a system only contains weaker bonds a value up to 5.0 may speed up the convergence.
 * @param nhist_max Maximal number of steps that will be stored in the history list. Use a value between 3 and 20. Must be <= than 3*nat.
 * @param alpha0 Lower limit on the step size. 1.e-2 is the default.
 * @param eps_subsp Lower limit on linear dependencies of basis vectors in history list. Default 1.e-4.
 */
periodic_optimizer::periodic_optimizer(int &nat, double initial_step_size, int nhist_max, double alpha0, double eps_subsp)
{
  this->nat = nat;
  this->ndim = 3*nat;
  this->initial_step_size = initial_step_size;
  this->n_hist_max = nhist_max;
  this->opt_lattice = false;
  this->opt = make_unique<sqnm_space::SQNM>(ndim, n_hist_max, initial_step_size, alpha0, eps_subsp);
}

/**
 * @brief Construct a new periodic optimizer::periodic optimizer object for variable cell shape optimization with default parameters.
 * 
 * @param nat number of atoms
 * @param lat_a first lattice vector
 * @param lat_b second lattice vector
 * @param lat_c third lattice vector
 */
periodic_optimizer::periodic_optimizer(int &nat, Vector3d &lat_a, Vector3d &lat_b, Vector3d &lat_c)
{
  setupInitialLattice(nat, lat_a, lat_b, lat_c);
  this->opt = make_unique<sqnm_space::SQNM>(ndim, n_hist_max, initial_step_size);
}

/**
 * @brief Construct a new periodic optimizer::periodic optimizer object for variable cell shape optimization with custom parameters.
 * 
 * @param nat number of atoms
 * @param lat_a first lattice vector
 * @param lat_b second lattice vector
 * @param lat_c third lattice vector
* @param initial_step_size initial step size. default is 1.0. For systems with hard bonds (e.g. C-C) use a value between and 1.0 and
 * 2.5. If a system only contains weaker bonds a value up to 5.0 may speed up the convergence.
 * @param nhist_max Maximal number of steps that will be stored in the history list. Use a value between 3 and 20. Must be <= than 3*nat.
 * @param lattice_weight weight / size of the supercell that is used to transform lattice derivatives. Use a value between 1 and 2. Default is 2.
 * @param alpha0 Lower limit on the step size. 1.e-2 is the default.
 * @param eps_subsp Lower limit on linear dependencies of basis vectors in history list. Default 1.e-4.
 */
periodic_optimizer::periodic_optimizer(int &nat, Vector3d &lat_a, Vector3d &lat_b, Vector3d &lat_c, double initial_step_size, int nhist_max, double lattice_weight, double alpha0, double eps_subsp)
{
  this->w = lattice_weight;
  this->n_hist_max = nhist_max;
  this->initial_step_size = initial_step_size;
  setupInitialLattice(nat, lat_a, lat_b, lat_c);
  this->opt = make_unique<sqnm_space::SQNM>(ndim, n_hist_max, initial_step_size, alpha0, eps_subsp);
}

/**
 * @brief Calculates new atomic coordinates that are closer to the local minimum. Fixed cell optimization. This function should be used the following way:
 * 1. calculate energies and forces at positions r.
 * 2. call the step function to update positions r.
 * 3. repeat.
 * 
 * @param r Input: atomic coordinates, dimension(3, nat). Output: improved coordinates that are calculated based on forces from this and previous iterations.
 * @param energy Potential energy of the system in it's current state
 * @param f Forces of the system in it's current state. dimension(3, nat)
 */
void periodic_optimizer::step(MatrixXd &r, double &energy, MatrixXd &f){
  VectorXd pos_all = Map<VectorXd>(r.data(), 3*nat);
  VectorXd force_all = - Map<VectorXd>(f.data(), 3*nat);
  pos_all += opt->step(pos_all, energy, force_all);
  r = Map<MatrixXd>(pos_all.data(), 3, nat);
  
}

/**
 * @brief Calculates new atomic coordinates that are closer to the local minimum. Variable cell shape optimization. This function should be used the following way:
 * 1. calculate energies, forces and stress tensor at positions r and lattice vectors a, b, c.
 * 2. call the step function to update positions r and lattice vectors.
 * 3. repeat.
 * 
 * @param r Input: atomic coordinates, dimension(3, nat). Output: improved coordinates that are calculated based on forces from this and previous iterations.
 * @param energy Potential energy of the system in it's current state
 * @param f Forces of the system in it's current state. dimension(3, nat)
 * @param lat_a first lattice vector
 * @param lat_b second lattice vector
 * @param lat_c third lattice vector
 * @param stress stress tensor of the system in it' current state.
 */
void periodic_optimizer::step(MatrixXd &r, double &energy, MatrixXd &f, Vector3d &lat_a, Vector3d &lat_b, Vector3d &lat_c, Matrix3d &stress){
  Matrix3d alat;
  MatrixXd alat_tilde;
  alat.col(0) = lat_a;
  alat.col(1) = lat_b;
  alat.col(2) = lat_c;

  //cout << "transform atom coordinates" << endl;
  //  calculate transformed coordinates
  MatrixXd q(3, nat);
  MatrixXd dq(3, nat);
  q = initial_latttice * alat.inverse() * r;
  dq = - alat * this->initial_latttice_inv * f;

  //cout << "transform lattice vectors" << endl;
  // transform lattice vectors
  alat_tilde = alat * lattice_transformer;
  //cout << calc_lattice_derivatices(stress, alat).array() / da.array() << endl;
  MatrixXd dalat = calc_lattice_derivatices(stress, alat) * lattice_transformer_inv;
  VectorXd qall = combine_matrices(q, alat_tilde);
  VectorXd dqall = combine_matrices(dq, dalat);
  
  // 
  //cout << "update coordinates" << endl;
  qall += this->opt->step(qall, energy, dqall);
  //cout << "update done" << endl;
  //cout << qall(1) << endl;

  split_matrices(q, alat_tilde, qall);
  alat = alat_tilde * lattice_transformer_inv;
  r = alat * this->initial_latttice_inv * q;
  lat_a = alat.col(0);
  lat_b = alat.col(1);
  lat_c = alat.col(2);
}

VectorXd periodic_optimizer::combine_matrices(MatrixXd a, MatrixXd b){
  VectorXd result(this->ndim);
  for (int i = 0; i < 3*nat; i++)
  {
    result(i) = a(i);
  }
  int j = 0;
  for (int i = 3*nat; i < ndim; i++)
  {
    result(i) = b(j);
    ++j; 
  }
  return result;
}

void periodic_optimizer::split_matrices(MatrixXd &a, MatrixXd &b, VectorXd &c){
  for (int i = 0; i < 3*nat; i++)
  {
    a(i) = c(i);
  }
  int j = 0;
  for (int i = 3*nat; i < ndim; i++)
  {
    b(j) = c(i);
    j++;
  }
  
}

Matrix3d periodic_optimizer::calc_lattice_derivatices(Matrix3d &stress, Matrix3d &alat){
  Matrix3d da = - stress * alat.inverse().transpose() * alat.determinant();
  return da;
}

void periodic_optimizer::setupInitialLattice(int &nat, Eigen::Vector3d &lat_a, Eigen::Vector3d &lat_b, Eigen::Vector3d &lat_c)
{
  this->nat = nat;
  this->ndim = 3*nat + 9;
  this->initial_latttice.col(0) = lat_a;
  this->initial_latttice.col(1) = lat_b;
  this->initial_latttice.col(2) = lat_c;
  this->initial_latttice_inv = initial_latttice.inverse();
  lattice_transformer.setZero();
  for (int i = 0; i < 3; i++)
  {
    lattice_transformer(i, i) = 1.0 / (initial_latttice.col(i).norm());
  }
  lattice_transformer = lattice_transformer * (w * sqrt(nat));
  lattice_transformer_inv = lattice_transformer.inverse();
  this->opt_lattice = true;
}
