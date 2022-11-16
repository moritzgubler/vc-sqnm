#include <Eigen/Dense>
#include <iostream>
#include <string>
#include <fstream>
#include "periodic_optimizer.hpp"

#ifdef __cplusplus
extern"C" {
    #endif
    void forces_bazant(const int * , double [], double [], double *, double [], double []);
    #ifdef __cplusplus
}
#endif

void ff_wrapper(int &nat, Eigen::MatrixXd &r, Eigen::Vector3d &lat_a, Eigen::Vector3d &lat_b, Eigen::Vector3d &lat_c, double &etot, Eigen::MatrixXd &f, Eigen::Matrix3d &stress){
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
  stress = Eigen::Map<Eigen::Matrix3d>(stressvec);
  f = Eigen::Map<Eigen::MatrixXd>(force, 3, nat);
}

void read_ascii(std::string &fname, int &nat, Eigen::MatrixXd &r, Eigen::Vector3d &lat_a, Eigen::Vector3d &lat_b, Eigen::Vector3d &lat_c){
  std::ifstream ascii_stream(fname);
  std::string comment;
  std::getline(ascii_stream, comment);
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

/*
 * This program creates an optimizer object and does 30 steps of a geometry optimization.
 * Detailed documentation is in the periodic_optimizer.hpp file. 
 * To integrate the optimizer in your code you simply need to include the historylist.hpp, sqnm.hpp and the
 * periodic_optimizer.hpp file.
*/
int main(int argc, char **argv) {

   // number of atoms
  int nat = 8;
  // filename of the ascii file that contains input coordinates. use test/test.ascii for example.
  std::string fname = argv[1];
  // positions of atoms
  Eigen::MatrixXd r(3, nat);
  // forces
  Eigen::MatrixXd f(3, nat);
  // lattice vectors
  Eigen::Vector3d lat_a, lat_b, lat_c;
  // stress tensor
  Eigen::Matrix3d stress;
  // potential energy
  double epot;

  // read the input coordinates.
  read_ascii(fname, nat, r, lat_a, lat_b, lat_c);
  

  // call this constructor for variable cell shape optimization
  PES_optimizer::periodic_optimizer test(nat, lat_a, lat_b, lat_c, -.1, 10, 2.0, 1.e-2, 1.e-3);

  // call this constructor for fixed cell optimization:
  //periodic_optimizer test(nat, 2.0, 10, 1.e-2, 1.e-4);

  // do 30 geometry optimization iterations.
  for (int i = 0; i < 30; i++)
  {
    // calculate energy, forces and stress tensor
    ff_wrapper(nat, r, lat_a, lat_b, lat_c, epot, f, stress);
    std::cout << i << " energy " << epot << " forcenorm " << f.norm() << "\n";

    /* call this step function for variable cell shape optimization
     * the step function calculates a new geometry that is closer to the local minima.
     * your job is to calculate enery and forces of the new guess of the step function and provide it 
     * later to the step function.
    */
    test.step(r, epot, f, lat_a, lat_b, lat_c, stress);
    // call this step function for fixed cell optimization
    //test.step(r, epot, f);
  }

  std::cout << "Current energy: " << epot<< "\n";
  // The SQNM method can estimate a lower bound of the ground state energy.
  // Here it is shown how you can access this estimate. The estimate is only accurate when 
  // the geometry optimization is reasonably well converged.
  std::cout << "Estimated lower bound of ground state energy: " << test.lower_bound()<< "\n";
  std::cout << "Energy uncertainty: " << epot - test.lower_bound()<< "\n";

  return 0;
}
