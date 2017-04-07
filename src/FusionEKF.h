#ifndef FusionEKF_H_
#define FusionEKF_H_

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>
#include "kalman_filter.h"

class FusionEKF {
public:
  /**
  * Constructor.
  */
  FusionEKF();

  /**
  * Destructor.
  */
  virtual ~FusionEKF();

  /**
  * Run the whole flow of the Kalman Filter from here.
  */
  bool ProcessMeasurement(const MeasurementPackage &measurement_pack);

  inline KalmanFilter* ekf() const { return ekf_; }
  
private:
  /*
   * Kalman Filter update and prediction math lives in here.
   */
  KalmanFilter *ekf_;

  // previous timestamp
  unsigned long long previous_timestamp_;

  // Laser measurement noise covariance
  Eigen::MatrixXd R_laser_;

  // Radar measurement noise covariance
  Eigen::MatrixXd R_radar_;

  // Measurament function laser
  Eigen::MatrixXd H_laser_;

  // Measurement jacobian radar
  Eigen::MatrixXd Hj_;

  // State covariance
  Eigen::MatrixXd P_;

  // noise components
  double noise_ax_;
  double noise_ay_;
};

#endif /* FusionEKF_H_ */
