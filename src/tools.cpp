#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */

  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;

  unsigned long n = estimations.size();
  for (int i=0; i < n; i++){
    rmse += (estimations[i] - ground_truth[i]).array().pow(2).matrix();
  }

  return (rmse / n).array().sqrt();
}