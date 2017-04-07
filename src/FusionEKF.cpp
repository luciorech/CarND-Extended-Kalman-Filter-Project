#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <cmath>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 * Uses noise_ax = 9 and noise_ay = 9 for Q matrix.
 */
FusionEKF::FusionEKF() :
  ekf_(nullptr),
  previous_timestamp_(0),
  R_laser_(MatrixXd(2, 2)),
  R_radar_(MatrixXd(3, 3)),
  H_laser_(MatrixXd(2, 4)),
  Hj_(MatrixXd(3, 4)),
  P_(MatrixXd(4, 4)),
  noise_ax_(9),
  noise_ay_(9)
{
  R_laser_ << 0.0225, 0,
              0, 0.0225;
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;
  P_ << 1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1000, 0,
        0, 0, 0, 1000;
}

FusionEKF::~FusionEKF() {
  if (ekf_) delete ekf_;
}

bool FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

  static double epsilon = 1e-8;

  // A measurement cannot be at origin (no matter if lidar or radar)
  // If so, we simply discard it
  if (fabs(measurement_pack.raw_measurements_[0]) < epsilon &&
      (fabs(measurement_pack.raw_measurements_[1]) < epsilon ||
       measurement_pack.sensor_type_ == MeasurementPackage::RADAR))
  {
    std::cout << "Skipping invalid measurement" << std::endl;
    return false;
  }

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!ekf_) {
    cout << "EKF: " << endl;
    VectorXd x(4);
    MatrixXd F(4, 4);
    F << 1, 0, 0, 0,
         0, 1, 0, 0,
         0, 0, 1, 0,
         0, 0, 0, 1;
    MatrixXd Q(4, 4);
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      double rho = measurement_pack.raw_measurements_[0];
      double phi = measurement_pack.raw_measurements_[1];
      double px = rho * cos(phi);
      double py = rho * sin(phi);
      x << px, py, 0, 0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      x << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }
    previous_timestamp_ = measurement_pack.timestamp_;
    ekf_ = new KalmanFilter(x, P_, F, H_laser_, R_laser_, R_radar_, Q);
    return true;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
  // If no time elapsed since last measurement, there's no need (nor reason)
  // to redo the prediction
  if (measurement_pack.timestamp_ > previous_timestamp_) {
    double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
    previous_timestamp_ = measurement_pack.timestamp_;
    ekf_->F()(0, 2) = dt;
    ekf_->F()(1, 3) = dt;
     
    double dt_2 = dt * dt;
    double dt_3 = dt_2 * dt;
    double dt_4 = dt_3 * dt;
    double dt_3_2 = dt_3 / 2;
    double dt_4_4 = dt_4 / 4;
    double dt_3_2_noise_ax = dt_3_2 * noise_ax_;
    double dt_3_2_noise_ay = dt_3_2 * noise_ay_;

    ekf_->Q() << (dt_4_4 * noise_ax_), 0, dt_3_2_noise_ax, 0,
                 0, (dt_4_4 * noise_ay_), 0, dt_3_2_noise_ay,
                 dt_3_2_noise_ax, 0, (dt_2 * noise_ax_), 0,
                 0, dt_3_2_noise_ay, 0, (dt_2 * noise_ay_);
    ekf_->Predict();
  }

  /*****************************************************************************
   *  Update
   ****************************************************************************/
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    ekf_->UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    ekf_->Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_->x() << endl;
  cout << "P_ = " << ekf_->P() << endl;
  return true;
}
