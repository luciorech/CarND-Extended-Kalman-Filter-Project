#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse << 0,0,0,0;
  if (estimations.size() == 0 || estimations.size() != ground_truth.size()) {
    std::cerr << "ERROR: mismatched sizes for estimations and ground truth vectors."
              << std::endl;
    return rmse;
  }

  for (int i = 0; i < estimations.size(); ++i) {
    VectorXd diff = estimations[i] - ground_truth[i];
    VectorXd diff2 = diff.array() * diff.array();
    rmse = diff2 + rmse;
  }
  
  rmse = rmse.array() / estimations.size();
  rmse = rmse.array().sqrt();
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  MatrixXd Hj(3,4);

  //recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);
	
  //compute the Jacobian matrix
  float c1 = px * px + py* py;
  float c2 = sqrt(c1);
  float c3 = c1 * c2;

  //check division by zero
  static double epsilon = 1e-10;
  if (fabs(c1) < epsilon) {
    std::cerr << "CalculateJacobian() - Error - Division by Zero" << std::endl;
    Hj << 1, 1, 1, 1,
          1, 1, 1, 1,
          1, 1, 1, 1;
  } else {
    Hj << (px/c2), (py/c2), 0, 0,
         -(py/c1), (px/c1), 0, 0,
          py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;
  }
  return Hj;
}
