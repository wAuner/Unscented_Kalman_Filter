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

#include "ukf.h"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

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

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 5;

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
  lambda_ = 3 - n_x_;
  Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  x_aug_ = VectorXd(n_aug_);
  P_aug_ = MatrixXd(n_aug_, n_aug_);
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  weights_ = VectorXd(2 * n_aug_ +1);

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
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
  double delta_t = (meas_package.timestamp_ - time_us_) / 1e6;
  time_us_ = meas_package.timestamp_;

  Prediction(delta_t);
  if (use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  } else if (use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER) {
    UpdateLidar(meas_package);
  }

}

void UKF::InitState(MeasurementPackage& meas_package) {
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    double rho = meas_package.raw_measurements_(0);
    double phi = meas_package.raw_measurements_(1);
    double rho_dot = meas_package.raw_measurements_(2);

    //initialize position from radar measurement
    x_.head(2) << rho * cos(phi), rho * sin(phi);
    // TODO: convert radar speed to state speed
    x_.tail(3) << 0, 0, 0;

  } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    //initialize position from lidar measurement
    x_.head(2) << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1);
    // TODO: initialize remaining state values
    x_.tail(3) << 0, 0, 0;
  }

  time_us_ = meas_package.timestamp_;
  is_initialized_ = true;

  // cap px if it is too small
  if (fabs(x_(0)) < 1e-4) {
    x_(0) = 1e-2;
    std::cout << "Init for px is too small\n";
  }

  // cap py if it is too small
  if (fabs(x_(1)) < 1e-4) {
    x_(1) = 1e-2;
    std::cout << "Init for py is too small\n";
  }
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

  P_aug_.topRightCorner(n_x_, n_x_) = P_;
  P_aug_(5, 5) = pow(std_a_, 2);
  P_aug_(6, 6) = pow(std_yawdd_, 2);
}

void UKF::GenerateSigmaPoints() {

  Xsig_aug_.col(0) = x_aug_;

  // create square root matrix
  MatrixXd sqP_aug = P_aug_.llt().matrixL();

  for (int i = 1; i < n_aug_; ++i) {
    Xsig_aug_.col(i) = x_aug_ + sqrt(lambda_ + n_x_) * sqP_aug.col(i);
    Xsig_aug_.col(i + n_aug_) = x_aug_ - sqrt(lambda_ + n_x_) * sqP_aug.col(i);
  }
}

void UKF::PredictSigmaPoints(double delta_t) {
  double vel, psi, psi_dot, nu_a, nu_psi;

  for (int k = 0; k < 2 * n_aug_ + 1; ++k) {
    vel = Xsig_aug_.col(k)(2);
    psi = Xsig_aug_.col(k)(3);
    psi_dot = Xsig_aug_.col(k)(4);
    nu_a = Xsig_aug_.col(k)(5);
    nu_psi = Xsig_aug_.col(k)(6);

    Xsig_pred_.col(k) = Xsig_aug_.col(k).head(5);
    // case for driving curves
    if (fabs(Xsig_aug_.col(k)(4)) >= 1e-3) {
      Xsig_pred_.col(k)(0) += vel / psi_dot * (sin(psi + psi_dot * delta_t) - sin(psi))
                              + 0.5 * pow(delta_t, 2) * cos(psi) * nu_a;
      Xsig_pred_.col(k)(1) += vel / psi_dot * (-cos(psi + psi_dot * delta_t) + cos(psi))
                              + 0.5 * pow(delta_t, 2) * sin(psi) * nu_a;
      Xsig_pred_.col(k)(2) += delta_t * nu_a;
      Xsig_pred_.col(k)(3) += psi_dot * delta_t + 0.5 * pow(delta_t, 2) * nu_psi;
      Xsig_pred_.col(k)(4) += delta_t * nu_psi;
    }
    // case for driving a straight line
    else if (fabs(Xsig_aug_.col(k)(4)) < 1e-3) {
      Xsig_pred_.col(k)(0) += vel * cos(psi) * delta_t
                              + 0.5 * pow(delta_t, 2) * cos(psi) * nu_a;
      Xsig_pred_.col(k)(1) += vel * sin(psi) * delta_t
                              + 0.5 * pow(delta_t, 2) * sin(psi) * nu_a;
      Xsig_pred_.col(k)(2) += delta_t * nu_a;
      Xsig_pred_.col(k)(3) += psi_dot * delta_t
                              + 0.5 * pow(delta_t, 2) * nu_psi;
      Xsig_pred_.col(k)(4) += delta_t * nu_psi;
    }
  }
}

void UKF::PredictMeanCovariance() {
  // calculate weights
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (int i = 1; i < 2 * n_aug_ + 1; ++i) {
    weights_(i) = 1 / (2 * (lambda_ + n_aug_));
  }
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    x_ += weights_(i) * Xsig_pred_.col(i);
  }

  // predict state and covariance
  VectorXd x_diff(n_x_);
  for (int j = 0; j < 2 * n_aug_ + 1; ++j) {
    x_diff = Xsig_pred_.col(j) - x_;
    // normalize angles
    while (x_diff(3) > M_PI) x_diff(3) -= 2 * M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2 * M_PI;

    P_ += weights_(j) * x_diff * x_diff.transpose();
  }



}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  // Predict Measurement
  //###############################################
  // map sigma points to measurement space
  MatrixXd Zsig(n_z_laser_, 2 * n_aug_ + 1);
  Zsig = Xsig_pred_.topLeftCorner(2, 2 * n_aug_ + 1);

  // map predictions to measurement space
  VectorXd z_pred(n_z_laser_);
  StateToMeasurementSpace(z_pred, Zsig);

  // map covariance to measurement space
  MatrixXd S(n_z_laser_, n_z_laser_);
  LidarSpaceCovariance(S, Zsig, z_pred);

  // Update State
  //###############################################
  // calculate Cross-correlation Matrix
  MatrixXd Tc(n_x_, n_z_laser_);
  CalculateCcMatrix(Tc, Zsig, z_pred);

  UpdateState(meas_package.raw_measurements_, Tc, S, z_pred, n_z_laser_, false);

}

void UKF::LidarSpaceCovariance(MatrixXd& S, MatrixXd& Zsig, VectorXd& z_pred) {
  MatrixXd R(2, 2);
  R(0,0) = pow(std_laspx_, 2);
  R(1,1) = pow(std_laspy_, 2);


  //calculate measurement covariance matrix S
  VectorXd z_diff(n_z_laser_);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    z_diff = Zsig.col(i) - z_pred;
    // normalizing angle
    while (z_diff(1) > M_PI) z_diff(1) -= 2 * M_PI;
    while (z_diff(1) < M_PI) z_diff(1) += 2 * M_PI;

    S += weights_(i) * z_diff * z_diff.transpose();
  }
  S += R;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
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
  SimgaPointsToRadarSpace(Zsig);

  // map predictions to measurement space
  VectorXd z_pred(n_z_radar_);
  StateToMeasurementSpace(z_pred, Zsig);

  MatrixXd S(n_z_radar_, n_z_radar_);
  RadarSpaceCovariance(S, Zsig, z_pred);

  // Update State
  //###############################################
  // calculate Cross-correlation Matrix
  MatrixXd Tc(n_x_, n_z_radar_);
  CalculateCcMatrix(Tc, Zsig, z_pred);

  UpdateState(meas_package.raw_measurements_, Tc, S, z_pred, n_z_radar_, true);

}


void UKF::SimgaPointsToRadarSpace(MatrixXd& Zsig) {
  // matrix for sigma points in measurement space

  double px, py, vel, psi, psi_dot;
  for (int j = 0; j < 2 * n_aug_ + 1; ++j) {
    px = Xsig_pred_.col(j)(0);
    py = Xsig_pred_.col(j)(1);
    vel = Xsig_pred_.col(j)(2);
    psi = Xsig_pred_.col(j)(3);
    psi_dot = Xsig_pred_.col(j)(4);

    Zsig.col(j)(0) = sqrt(pow(px, 2) + pow(py, 2));
    Zsig.col(j)(1) = atan2(py, px);
    Zsig.col(j)(2) = (px * cos(psi) * vel + py * sin(psi) * vel) / sqrt(pow(px, 2) + pow(py, 2));

    // TODO check if necessary
    // clip predictions too close to zero
    if (fabs(px) < 1e-4 or fabs(py) < 1e-4) {
      if (fabs(px) < 1e-4) {
        px = 1e-4;
      }

      if (fabs(py) < 1e-4) {
        py = 1e-4;
      }
      Zsig.col(j)(0) = sqrt(pow(px, 2) + pow(py, 2));
      Zsig.col(j)(1) = 0;
      Zsig.col(j)(2) = 0;
    }
  }
}

void UKF::StateToMeasurementSpace(VectorXd &z_pred, MatrixXd &Zsig) {
  //calculate mean predicted measurement
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    z_pred += weights_(i) * Zsig.col(i);
  }
}

void UKF::RadarSpaceCovariance(MatrixXd& S, MatrixXd& Zsig, VectorXd& z_pred) {
  MatrixXd R(3,3);
  R(0,0) = pow(std_radr_, 2);
  R(1,1) = pow(std_radphi_, 2);
  R(2,2) = pow(std_radrd_, 2);

  //calculate measurement covariance matrix S
  VectorXd z_diff(n_z_radar_);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    z_diff = Zsig.col(i) - z_pred;
    // normalizing angle
    while (z_diff(1) > M_PI) z_diff(1) -= 2 * M_PI;
    while (z_diff(1) < M_PI) z_diff(1) += 2 * M_PI;

    S += weights_(i) * z_diff * z_diff.transpose();

  }
  S += R;
}

void UKF::CalculateCcMatrix(MatrixXd& Tc, MatrixXd& Zsig, VectorXd& z_pred) {
  VectorXd x_diff(n_x_);
  VectorXd z_diff(n_x_);

  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    x_diff = Xsig_pred_.col(i) - x_;
    z_diff = Zsig.col(i) - z_pred;

    while (x_diff(3) > M_PI) x_diff(3) -= 2 * M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2 * M_PI;

    while (z_diff(1) > M_PI) z_diff(1) -= 2 * M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2 * M_PI;

    Tc += weights_(i) * x_diff * z_diff.transpose();
  }
}

void UKF::UpdateState(VectorXd& z, MatrixXd& Tc,
                      MatrixXd& S, VectorXd& z_pred,
                      int n_z, bool isRadar)
{
  //calculate Kalman gain K;
  MatrixXd K(n_x_, n_z);
  VectorXd z_diff(n_x_);

  K = Tc * S.inverse();

  //update state mean and covariance matrix
  z_diff = z - z_pred;

  if (isRadar) {
    while (z_diff(1) > M_PI) z_diff(1) -= 2 * M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2 * M_PI;
  }

  x_ += K * z_diff;
  P_ -= K * S * K.transpose();
}