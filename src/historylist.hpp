#ifndef HISTORYLIST_HPP
#define HISTORYLIST_HPP
#include <eigen3/Eigen/Dense>

namespace hlist_space{
class HistoryList {
  public: 
  HistoryList(int &ndim_, int &nhistx_);
  int add(Eigen::VectorXd &x_in);

  Eigen::MatrixXd hlist;
  Eigen::MatrixXd difflist;
  Eigen::MatrixXd normalized_difflist;

  private:
  int i_count = 0;
  int nhistx;
  int ndim;
  Eigen::VectorXd oldElem;

};
}
#endif