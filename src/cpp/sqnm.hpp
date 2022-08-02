#ifndef SQNM_HPP
#define SQNM_HPP
#include <Eigen/Dense>
#include <iostream>
#include <memory>
#include "historylist.hpp"

namespace sqnm_space
{
  class SQNM {
    private:
    int ndim;
    int nhistx;
    double eps_subsp = 1.e-4;
    double alpha0 = 1.e-2;
    std::unique_ptr<hlist_space::HistoryList> xlist;
    std::unique_ptr<hlist_space::HistoryList> flist;
    double alpha;
    Eigen::VectorXd dir_of_descent;
    double prev_f;
    Eigen::VectorXd prev_df_dx;
    Eigen::MatrixXd h_subsp;
    Eigen::MatrixXd h_evec_subsp;
    Eigen::MatrixXd h_evec;
    Eigen::VectorXd h_eval;
    Eigen::VectorXd res;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esolve;
    Eigen::VectorXd res_temp;
    int nhist = 0;

    public:
    SQNM(int ndim_, int nhistx_, double alpha_);
    SQNM(int ndim_, int nhistx_, double alpha_, double alpha0_, double eps_subsp_);

    Eigen::VectorXd step(Eigen::VectorXd &x, double &f_of_x, Eigen::VectorXd &df_dx);

    private:
    double calc_gainratio(double &f);

    void adjust_stepsize(double &gainratio);

    Eigen::MatrixXd calc_ovrlp();

  };
}
#endif