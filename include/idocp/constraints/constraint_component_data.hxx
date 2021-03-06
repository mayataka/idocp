#ifndef IDOCP_CONSTRAINT_COMPONENT_DATA_HXX_
#define IDOCP_CONSTRAINT_COMPONENT_DATA_HXX_

#include "idocp/constraints/constraint_component_data.hpp"

#include <stdexcept>
#include <iostream>


namespace idocp {

inline ConstraintComponentData::ConstraintComponentData(const int dimc)
  : slack(Eigen::VectorXd::Zero(dimc)),
    dual(Eigen::VectorXd::Zero(dimc)),
    residual(Eigen::VectorXd::Zero(dimc)),
    duality(Eigen::VectorXd::Zero(dimc)),
    dslack(Eigen::VectorXd::Zero(dimc)),
    ddual(Eigen::VectorXd::Zero(dimc)),
    r(),
    J(),
    dimc_(dimc) {
  try {
    if (dimc <= 0) {
      throw std::out_of_range("invalid argment: dimc must be positive!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
}


inline ConstraintComponentData::ConstraintComponentData()
  : slack(),
    dual(),
    residual(),
    duality(),
    dslack(),
    ddual(),
    r(),
    J(),
    dimc_(0) {
}


inline ConstraintComponentData::~ConstraintComponentData() {
}


inline int ConstraintComponentData::dimc() const {
  return dimc_;
}


inline bool ConstraintComponentData::checkDimensionalConsistency() const {
  if (slack.size() != dimc_) {
    return false;
  }
  if (dual.size() != dimc_) {
    return false;
  }
  if (residual.size() != dimc_) {
    return false;
  }
  if (duality.size() != dimc_) {
    return false;
  }
  if (dslack.size() != dimc_) {
    return false;
  }
  if (ddual.size() != dimc_) {
    return false;
  }
  return true;
}


inline bool ConstraintComponentData::isApprox(
    const ConstraintComponentData& other) const {
  if (!slack.isApprox(other.slack)) {
    return false;
  }
  if (!dual.isApprox(other.dual)) {
    return false;
  }
  if (!residual.isApprox(other.residual)) {
    return false;
  }
  if (!duality.isApprox(other.duality)) {
    return false;
  }
  if (!dslack.isApprox(other.dslack)) {
    return false;
  }
  if (!ddual.isApprox(other.ddual)) {
    return false;
  }
  return true;
}

} // namespace idocp

#endif // IDOCP_CONSTRAINT_COMPONENT_DATA_HXX_