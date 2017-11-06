/* Class for an Unscented Kalman Filter
 * performs the following steps:
 *
 * Prediction
 * ---------------------------------------
 * 1. Generate Sigma Points
 * 2. Predict Sigma Points
 * 3. Predict Mean and Covariance
 *
 * Update
 * ---------------------------------------
 * 1. Predict Measurement
 * 2. Update State
 *
 * */

// newest version
#include "ukf.h"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;
using std::cout;
using std::endl;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);
  x_.fill(0);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);
  P_.fill(0);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.9;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.
  Hint: one or more values initialized above might be wildly off...
  */
  is_initialized_ = false;
  n_x_ = 5;
  n_aug_ = 7;
  P_ = Eigen::MatrixXd::Identity(n_x_, n_x_);
  //P_(0, 0) = 0.15;
  //P_(1,1) = 0.15;

  lambda_ = 3 - n_aug_;

  Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_aug_.fill(0);

  x_aug_ = VectorXd(n_aug_);
  x_aug_.fill(0);

  P_aug_ = MatrixXd(n_aug_, n_aug_);

  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  Xsig_pred_.fill(0);

  weights_ = VectorXd(2 * n_aug_ +1);
  weights_.fill(0);
  // calculate weights
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (int i = 1; i < 2 * n_aug_ + 1; ++i) {
    weights_(i) = 1. / (2 * (lambda_ + n_aug_));
  }
  cout << "Weights \n" << weights_ << endl;

  R_Radar_ = MatrixXd(3,3);
  R_Radar_ << pow(std_radr_, 2), 0, 0,
              0, pow(std_radphi_, 2), 0,
              0, 0, pow(std_radrd_, 2);
  R_Lidar_ = MatrixXd(2,2);
  R_Lidar_ << pow(std_laspx_, 2), 0,
              0, pow(std_laspy_, 2);

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(const MeasurementPackage& meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  // initializing state vector and process covariance matrix with first measurement
  if (!is_initialized_) {
    InitState(meas_package);
    return;
  }
  cout << "Verarbeite " << meas_package.sensor_type_ << " (0=Laser, 1=Radar) Messung" << endl;

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
    //std::cout <<"Verarbeite Radar Messung\n" << x_ << std::endl;
    double delta_t = (meas_package.timestamp_ - time_us_) / 1e6;
    time_us_ = meas_package.timestamp_;
    Prediction(delta_t);
    UpdateRadar(meas_package);
    //std::cout <<"x nach radar update:\n" << x_ << std::endl;
  } else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
    //std::cout <<"Verarbeite Laser Messung\n" << x_ << std::endl;
    double delta_t = (meas_package.timestamp_ - time_us_) / 1e6;
    time_us_ = meas_package.timestamp_;
    Prediction(delta_t);
    UpdateLidar(meas_package);
    //std::cout <<"x nach laser update:\n" << x_ << std::endl;
  }


}

void UKF::InitState(const MeasurementPackage& meas_package) {
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    double rho = meas_package.raw_measurements_(0);
    double phi = meas_package.raw_measurements_(1);

    double px = rho * cos(phi);
    double py = rho * sin(phi);
    //initialize position from radar measurement
    x_.head(2) << px, py;
    // TODO: convert radar speed to state speed
  } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    //initialize position from lidar measurement
    double px = meas_package.raw_measurements_(0);
    double py = meas_package.raw_measurements_(1);

    x_.head(2) << px, py;
  }
  x_.tail(3) << 0, 0, 0;

  time_us_ = meas_package.timestamp_;
  is_initialized_ = true;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {

  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  AugmentState();

  GenerateSigmaPoints();

  PredictSigmaPoints(delta_t);

  PredictMeanCovariance();
}


void UKF::AugmentState() {
  // augmentation is needed because process noise has a non-linear effect on the state

  // augment state vector and covariance matrix to account for process noise
  x_aug_ << x_, 0, 0;
  //std::cout << "x_aug_\n" << x_aug_ << std::endl;
  P_aug_.fill(0);
  P_aug_.topLeftCorner(n_x_, n_x_) = P_;
  P_aug_(5, 5) = pow(std_a_, 2);
  P_aug_(6, 6) = pow(std_yawdd_, 2);
  std::cout << "P_aug_\n" << P_aug_ << std::endl;
}

void UKF::GenerateSigmaPoints() {
  Xsig_aug_.col(0) = x_aug_;

  // create square root matrix
  MatrixXd sqP_aug = P_aug_.llt().matrixL();
  //std::cout << "sqP_aug\n" << sqP_aug << std::endl;

  for (int i = 1; i < n_aug_ + 1; ++i) {
    Xsig_aug_.col(i)          = x_aug_ + sqrt(lambda_ + n_aug_) * sqP_aug.col(i - 1);
    Xsig_aug_.col(i + n_aug_) = x_aug_ - sqrt(lambda_ + n_aug_) * sqP_aug.col(i - 1);

  }
  std::cout << "Xsig_aug_\n" << Xsig_aug_ << std::endl;
}

// Generated sigma points are mapped from k to k+1 through the non-linear
// motion model
void UKF::PredictSigmaPoints(double delta_t) {
  Xsig_pred_.fill(0);
  for (int k = 0; k < 2 * n_aug_ + 1; ++k) {
    double vel     = Xsig_aug_.col(k)(2);
    double psi     = Xsig_aug_.col(k)(3);
    double psi_dot = Xsig_aug_.col(k)(4);
    double nu_a    = Xsig_aug_.col(k)(5);
    double nu_psi  = Xsig_aug_.col(k)(6);

    // TODO check equations!!!
    //std::cout << "Xsig_pred_:\n" << Xsig_pred_ << std::endl;

    Xsig_pred_.col(k) = Xsig_aug_.col(k).head(5);
    // case for driving curves
    if (fabs(psi_dot) >= 1e-3) {
      Xsig_pred_.col(k)(0) += vel / psi_dot * (sin(psi + psi_dot * delta_t) - sin(psi))
                              + 0.5 * pow(delta_t, 2) * cos(psi) * nu_a;
      Xsig_pred_.col(k)(1) += vel / psi_dot * (-cos(psi + psi_dot * delta_t) + cos(psi))
                              + 0.5 * pow(delta_t, 2) * sin(psi) * nu_a;
    }
    // case for driving a straight line, avoid division by zero
    else if (fabs(psi_dot) < 1e-3) {
      Xsig_pred_.col(k)(0) += vel * cos(psi) * delta_t
                              + 0.5 * pow(delta_t, 2) * cos(psi) * nu_a;
      Xsig_pred_.col(k)(1) += vel * sin(psi) * delta_t
                              + 0.5 * pow(delta_t, 2) * sin(psi) * nu_a;
    }

    Xsig_pred_.col(k)(2) += delta_t * nu_a;
    Xsig_pred_.col(k)(3) += psi_dot * delta_t
                            + 0.5 * pow(delta_t, 2) * nu_psi;
    Xsig_pred_.col(k)(4) += delta_t * nu_psi;
  }
  std::cout << "Xsig_pred_\n" << Xsig_pred_ << std::endl;
}

// The state vector mean and the process covariance at k+1 is estimated
// based on the sigma points which have been mapped to k+1 with the motion model
void UKF::PredictMeanCovariance() {
  // estimate state vector mean based on sigma points
  x_.fill(0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    x_ += weights_(i) * Xsig_pred_.col(i);
  }

  std::cout <<"Predicted Mean\n" << x_ << std::endl;

  // estimate covariance based on sigma points
  P_.fill(0);
  for (int j = 0; j < 2 * n_aug_ + 1; ++j) {
    VectorXd x_diff = Xsig_pred_.col(j) - x_;
    // normalize angles
    x_diff(3) = NormalizeAngle(x_diff(3));

    P_ += weights_(j) * x_diff * x_diff.transpose();
  }
  std::cout << "Predicted covariance\n" << P_ << std::endl;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(const MeasurementPackage& meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  // no sigma points needed since the transformation to measurement space
  // is linear in case of lidar
  MatrixXd H(n_z_laser_, n_x_);
  H << 1, 0, 0, 0, 0,
       0, 1, 0, 0, 0;
  VectorXd z_pred = H * x_;
  VectorXd z = meas_package.raw_measurements_;
  VectorXd y = z - z_pred;
  MatrixXd S(n_z_laser_, n_z_laser_);
  S = H * P_ * H.transpose();
  S += R_Lidar_;
  MatrixXd K = P_ * H.transpose() * S.inverse();

  x_ +=  K * y;
  MatrixXd I = MatrixXd::Identity(n_x_, n_x_);
  P_ = (I - K * H) * P_;
  NIS_lidar_ = (z - z_pred).transpose() * S.inverse() * (z - z_pred);
  cout << "x_ nach Laser Update\n" << x_ << endl;
  cout << "P_ nach Laser Update\n" << P_ << endl;
  std::cout << "Lidar NIS is\n" << NIS_lidar_ << std::endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(const MeasurementPackage& meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  // Predict Measurement
  //###############################################
  // map sigma points to measurement space
  MatrixXd Zsig(n_z_radar_, 2 * n_aug_ + 1);
  Zsig.fill(0);
  SimgaPointsToRadarSpace(Zsig);

  // map predictions to measurement space
  VectorXd z_pred(n_z_radar_);
  z_pred.fill(0);
  StateToMeasurementSpace(z_pred, Zsig);

  MatrixXd S(n_z_radar_, n_z_radar_);
  S.fill(0);
  RadarSpaceCovariance(S, Zsig, z_pred);

  // Update State
  //###############################################
  // calculate Cross-correlation Matrix
  MatrixXd Tc(n_x_, n_z_radar_);
  //std::cout << "Tc UpdateRadar\n" << Tc << std::endl << std::endl;
  Tc.fill(0);
  //std::cout << "Tc UpdateRadar\n" << Tc << std::endl << std::endl;
  CalculateCcMatrix(Tc, Zsig, z_pred, true);

  VectorXd z(n_z_radar_);
  z = meas_package.raw_measurements_;
  UpdateState(z, Tc, S, z_pred, n_z_radar_, true);

  // calculate NIS
  // TODO normalize angles???
  NIS_radar_ = (z - z_pred).transpose() * S.inverse() * (z - z_pred);
  std::cout << "Radar NIS is\n" << NIS_radar_ << std::endl;
}


void UKF::SimgaPointsToRadarSpace(MatrixXd& Zsig) {
  // map matrix for sigma points to measurement space
  double px, py, vel, psi;
  for (int j = 0; j < 2 * n_aug_ + 1; ++j) {
    px = Xsig_pred_.col(j)(0);
    py = Xsig_pred_.col(j)(1);
    vel = Xsig_pred_.col(j)(2);
    psi = Xsig_pred_.col(j)(3);

    Zsig.col(j)(0) = sqrt(pow(px, 2) + pow(py, 2));
    Zsig.col(j)(1) = atan2(py, px);
    Zsig.col(j)(2) = (px * cos(psi) * vel + py * sin(psi) * vel) / sqrt(pow(px, 2) + pow(py, 2));
  }
  //std::cout << "Xsig_pred\n" << Xsig_pred_ << std::endl;
}

void UKF::StateToMeasurementSpace(VectorXd &z_pred, MatrixXd &Zsig) {
  //calculate mean predicted measurement
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    z_pred += weights_(i) * Zsig.col(i);
  }
  //std::cout << "z_pred\n" << z_pred << std::endl << std::endl;
}

void UKF::RadarSpaceCovariance(MatrixXd& S, MatrixXd& Zsig, VectorXd& z_pred) {
  //calculate measurement covariance matrix S
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // normalizing angle
    z_diff(1) = NormalizeAngle(z_diff(1));

    S += weights_(i) * z_diff * z_diff.transpose();
  }
  S += R_Radar_;
  //std::cout << "S\n" << S << std::endl << std::endl;
}

void UKF::CalculateCcMatrix(MatrixXd& Tc, MatrixXd& Zsig, VectorXd& z_pred, bool isRadar) {

  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    VectorXd z_diff = Zsig.col(i) - z_pred;

    x_diff(3) = NormalizeAngle(x_diff(3));

    if (isRadar) {
      z_diff(1) = NormalizeAngle(z_diff(1));
    }
//    std::cout << "Tc\n" << Tc << std::endl << std::endl;
//    std::cout << "x_diff\n" << x_diff << std::endl << std::endl;
//    std::cout << "z_diff\n" << z_diff << std::endl << std::endl;
    Tc += weights_(i) * x_diff * z_diff.transpose();
  }
}

void UKF::UpdateState(VectorXd& z, MatrixXd& Tc,
                      MatrixXd& S, VectorXd& z_pred,
                      int n_z, bool isRadar)
{
  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  VectorXd z_diff = z - z_pred;

  //std::cout << "Xsig_pred_\n" << Xsig_pred_ << std::endl << std::endl;
  //std::cout << "Tc\n" << Tc << std::endl << std::endl;

  //std::cout << "K\n" << K << std::endl;
  //update state mean and covariance matrix

  if (isRadar) {
    z_diff(1) = NormalizeAngle(z_diff(1));
  }
  //std::cout << "x_\n" << x_ << std::endl;
  x_ += K * z_diff;
  P_ -= K * S * K.transpose();
  //std::cout << "x_\n" << x_ << std::endl;
}

double UKF::NormalizeAngle(double angle)
{
  //std::cout << angle << std::endl;
  while (angle > M_PI) {
    std::cout << "Angle too big:\n" << angle << std::endl;
    angle -= 2 * M_PI;
    std::cout << "Angle corrected:\n" << angle << std::endl;
  }

  while (angle < -M_PI) {
    std::cout << "Angle too small:\n" << angle << std::endl;
    angle += 2 * M_PI;
    std::cout << "Angle corrected:\n" << angle << std::endl;
  }
  //std::cout << "angle is normalized\n\n\n\n";
  return angle;
}

