#ifndef KALMAN_FILTER_H_
#define KALMAN_FILTER_H_
#include "Eigen/Dense"

class KalmanFilter {
public:
  /**
   * Constructor
   * @param x_in Initial state
   * @param P_in Initial state covariance
   * @param F_in Transition matrix
   * @param H_in Measurement matrix
   * @param R_kf_in Measurement covariance matrix for kalman filter
   * @param R_ekf_in Measurement covariance matrix for extended kalman filter   
   * @param Q_in Process covariance matrix   
   */
  KalmanFilter(Eigen::VectorXd &x_in, Eigen::MatrixXd &P_in, Eigen::MatrixXd &F_in,
               Eigen::MatrixXd &H_in, Eigen::MatrixXd &R_kf_in, Eigen::MatrixXd &R_ekf_in,
               Eigen::MatrixXd &Q_in);

  /**
   * Destructor
   */
  virtual ~KalmanFilter();

  /**
   * Init Initializes Kalman filter
   */
  void Init();

  /**
   * Prediction Predicts the state and the state covariance
   * using the process model
   * @param delta_T Time between k and k+1 in s
   */
  void Predict();

  /**
   * Updates the state by using standard Kalman Filter equations
   * @param z The measurement at k+1
   */
  void Update(const Eigen::VectorXd &z);

  /**
   * Updates the state by using Extended Kalman Filter equations
   * @param z The measurement at k+1
   */
  void UpdateEKF(const Eigen::VectorXd &z);

  inline Eigen::MatrixXd& F() { return F_; }
  inline Eigen::MatrixXd& P() { return P_; }
  inline Eigen::MatrixXd& Q() { return Q_; }
  inline Eigen::VectorXd& x() { return x_; }
  
private:
  // state vector
  Eigen::VectorXd x_;

  // state covariance matrix
  Eigen::MatrixXd P_;

  // state transistion matrix
  Eigen::MatrixXd F_;

  // process covariance matrix
  Eigen::MatrixXd Q_;

  // measurement matrix
  Eigen::MatrixXd H_;

  // measurement covariance matrix
  Eigen::MatrixXd R_ekf_;
  Eigen::MatrixXd R_kf_;
};

#endif /* KALMAN_FILTER_H_ */
