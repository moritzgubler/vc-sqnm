#include <eigen3/Eigen/Dense>
#include <iostream>
#include "sqnm.hpp"
#include "periodic_optimizer.hpp"

using namespace Eigen;
using namespace std;

periodic_optimizer::periodic_optimizer(int &nat)
{
  this->nat = nat;
  this->ndim = 3*nat;
  this->opt_lattice = false;
  is_initialized = false;
}

periodic_optimizer::periodic_optimizer(int &nat, Vector3d &lat_a, Vector3d &lat_b, Vector3d &lat_c)
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
  is_initialized = false;
}

void periodic_optimizer::initialize(){
  this->opt = make_unique<sqnm_space::SQNM>(ndim, n_hist_max, initial_step_size);
  is_initialized = true;
}

void periodic_optimizer::step(MatrixXd &r, double &energy, MatrixXd &f){
  if (!is_initialized){
    cerr << "periodic_optimizer needs to be initialized before step is called." << "\n";
  }
  VectorXd pos_all = Map<VectorXd>(r.data(), 3*nat);
  VectorXd force_all = - Map<VectorXd>(f.data(), 3*nat);
  pos_all += opt->step(pos_all, energy, force_all);
  r = Map<MatrixXd>(pos_all.data(), 3, nat);
  
}

void periodic_optimizer::step(MatrixXd &r, double &energy, MatrixXd &f, Vector3d &lat_a, Vector3d &lat_b, Vector3d &lat_c, Matrix3d &stress){
    if (!is_initialized){
    cerr << "periodic_optimizer needs to be initialized before step is called." << "\n";
  }
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
