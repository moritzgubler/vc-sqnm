#ifndef PERIODIC_OPTIMIZER_HPP
#define PERIODIC_OPTIMIZER_HPP
#include <eigen3/Eigen/Dense>
#include <iostream>
#include "sqnm.hpp"

class periodic_optimizer
{
  private:
  Eigen::Matrix3d initial_latttice;
  Eigen::Matrix3d initial_latttice_inv;
  Eigen::Matrix3d lattice_transformer;
  Eigen::Matrix3d lattice_transformer_inv;
  std::unique_ptr<sqnm_space::SQNM> opt;
  int nat;
  int ndim;
  bool opt_lattice;
  bool is_initialized;
  double initial_step_size = 1;
  int n_hist_max = 10;
  double w = 2.0;
  public:

  void set_w(double w){
    if (is_initialized) return;
    this->w = w;
  }

  double get_w(){
    return w;
  }

  void set_n_hist_max(int nhistx){
    if (is_initialized) return;
    n_hist_max = nhistx;
  }

  int get_n_hist_max(){
    return n_hist_max;
  }

  void set_initial_step_size(double alpha_in){
    if (is_initialized) return;
    initial_step_size = alpha_in;
  }

  double get_initial_step_size(){
    return initial_step_size;
  }

  periodic_optimizer(int &nat);

  periodic_optimizer(int &nat, Eigen::Vector3d &lat_a, Eigen::Vector3d &lat_b, Eigen::Vector3d &lat_c);

  void initialize();

  void step(Eigen::MatrixXd &r, double &energy, Eigen::MatrixXd &f);

  void step(Eigen::MatrixXd &r, double &energy, Eigen::MatrixXd &f, Eigen::Vector3d &lat_a, Eigen::Vector3d &lat_b, Eigen::Vector3d &lat_c, Eigen::Matrix3d &stress);

  private:
  Eigen::VectorXd combine_matrices(Eigen::MatrixXd a, Eigen::MatrixXd b);

  void split_matrices(Eigen::MatrixXd &a, Eigen::MatrixXd &b, Eigen::VectorXd &c);

  Eigen::Matrix3d calc_lattice_derivatices(Eigen::Matrix3d &stress, Eigen::Matrix3d &alat);
};
#endif
