#include <eigen3/Eigen/Dense>
#include <iostream>
#include <string>
#include <fstream>
#include "periodic_optimizer.hpp"

using namespace Eigen;
using namespace std;

#ifdef __cplusplus
extern"C" {
    #endif
    void forces_bazant(const int * , double [], double [], double *, double [], double []);
    #ifdef __cplusplus
}
#endif

void ff_wrapper(int &nat, MatrixXd &r, Vector3d &lat_a, Vector3d &lat_b, Vector3d &lat_c, double &etot, MatrixXd &f, Matrix3d &stress){
  double *pos = r.data();
  double *force = f.data();
  double lat[9];
  lat[0] = lat_a(0);
  lat[1] = lat_a(1);
  lat[2] = lat_a(2);
  lat[3] = lat_b(0);
  lat[4] = lat_b(1);
  lat[5] = lat_b(2);
  lat[6] = lat_c(0);
  lat[7] = lat_c(1);
  lat[8] = lat_c(2);

  double stressvec[9];
  forces_bazant(&nat, lat, pos, &etot, force, stressvec);
  stress = Map<Matrix3d>(stressvec);
  f = Map<MatrixXd>(force, 3, nat);
}

void read_ascii(string &fname, int &nat, MatrixXd &r, Vector3d &lat_a, Vector3d &lat_b, Vector3d &lat_c){
  std::ifstream ascii_stream(fname);
  string comment;
  getline(ascii_stream, comment);
  //cout << comment << endl;
  lat_a.setZero();
  lat_b.setZero();
  ascii_stream >> lat_a(0) >> lat_b(0) >> lat_b(1);
  ascii_stream >> lat_c(0) >>  lat_c(1) >> lat_c(2);
  lat_a /= 0.52917721067;
  lat_b /= 0.52917721067;
  lat_c /= 0.52917721067;
  for (int i = 0; i < nat; i++)
  {
    ascii_stream >> r(0, i) >> r(1, i) >> r(2, i) >> comment;
    r.col(i) /= 0.52917721067;
  }
  ascii_stream.close();
}

int main(int argc, char **argv) {

  int nat = 8;
  //string fname = argv[1];
  //cout << fname << endl;
  string fname = "../../test/test.ascii";
  MatrixXd r(3, nat);
  MatrixXd f(3, nat);
  Vector3d lat_a, lat_b, lat_c;
  Matrix3d stress;
  double epot;

  read_ascii(fname, nat, r, lat_a, lat_b, lat_c);
  

  // call this constructor for variable cell shape optimization
  periodic_optimizer test(nat, lat_a, lat_b, lat_c, 2.0, 10, 2.0, 1.e-2, 1.e-4);

  // call this constructor for fixed cell optimization:
  //periodic_optimizer test(nat, 2.0, 10, 1.e-2, 1.e-4);

  for (int i = 0; i < 30; i++)
  {
    ff_wrapper(nat, r, lat_a, lat_b, lat_c, epot, f, stress);
    cout << i << " energy " << epot << "\n";

    // call this step function for variable cell shape optimization
    test.step(r, epot, f, lat_a, lat_b, lat_c, stress);
    // call this step function for fixed cell optimization
    //test.step(r, epot, f);
  }
  
  return 0;
}
