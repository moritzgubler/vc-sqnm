#include "historylist.hpp"
#include <eigen3/Eigen/Dense>

using namespace Eigen;
using namespace std;
using namespace hlist_space;

HistoryList::HistoryList(int &ndim_, int &nhistx_) {
  nhistx = nhistx_;
  ndim = ndim_;
  hlist.resize(ndim, nhistx);
  difflist.resize(ndim, nhistx);
  normalized_difflist.resize(ndim, nhistx);
  i_count = 0;
}

int HistoryList::add(VectorXd &x_in) {
  if (i_count < nhistx){
    hlist.col(i_count) = x_in;
    i_count += 1;
    for (int i = 1; i < i_count; i++){
      difflist.col(i-1) = hlist.col(i) - hlist.col(i - 1);
      normalized_difflist.col(i-1) = difflist.col(i-1) / difflist.col(i-1).norm();
    }
    return i_count - 1;
  } else {
    oldElem = hlist.col(0);
    for (int i = 0; i < nhistx - 1; i++){
      hlist.col(i) = hlist.col(i+1);
    }
    hlist.col(nhistx - 1) = x_in;
    // calculate difference list
    difflist.col(0) = hlist.col(0) - oldElem;
    normalized_difflist.col(0) = difflist.col(0) / difflist.col(0).norm();
    for (int i = 1; i < nhistx; i++){
      difflist.col(i) = hlist.col(i) - hlist.col(i - 1);
      normalized_difflist.col(i) = difflist.col(i) / difflist.col(i).norm();
    }
    return nhistx;
  }
}