#include "kalman_filter.h"
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                           MatrixXd &H_in, MatrixXd &R_kf_in, MatrixXd &R_ekf_in,
                           MatrixXd &Q_in) :
  x_(x_in), P_(P_in), F_(F_in), H_(H_in), R_kf_(R_kf_in),
  R_ekf_(R_ekf_in), Q_(Q_in)
{
}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Predict() {
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd PHt = P_ * Ht;
  MatrixXd S = H_ * PHt + R_kf_;
  MatrixXd Si = S.inverse();
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  // If cannot calculate jacobian, skip update
  MatrixXd Hj = Tools::CalculateJacobian(x_);
  if (Hj.rows() == 0) return;

  double px = x_[0];
  double py = x_[1];
  double vx = x_[2];
  double vy = x_[3];
  
  double rho_hx = sqrt((px * px) + (py * py));
  double phi_hx = atan2(py, px);

  std::cout << "px = " << px << "\n";
  std::cout << "py = " << py << "\n";
  std::cout << "rho_hx = " << rho_hx << "\n";
  std::cout << "phi_hx = " << phi_hx << "\n";
  std::cout << "rho = " << z[0] << "\n";
  std::cout << "phi = " << z[1] << "\n";

  VectorXd hx(3);
  hx << rho_hx, phi_hx, ((px * vx) + (py * vy)) / rho_hx;
  VectorXd y = z - hx;
  while (y[1] < -M_PI) y[1] += 2 * M_PI;
  while (y[1] > M_PI) y[1] -= 2 * M_PI;

  MatrixXd Ht = Hj.transpose();
  MatrixXd PHt = P_ * Ht;
  MatrixXd S = Hj * PHt + R_ekf_;
  MatrixXd Si = S.inverse();
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * Hj) * P_;
}
