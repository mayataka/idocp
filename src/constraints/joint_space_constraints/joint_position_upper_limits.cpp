#include "constraints/joint_space_constraints/joint_position_upper_limits.hpp"
#include "constraints/pdipm_func.hpp"

#include <cmath>
#include <assert.h>


namespace idocp {
namespace pdipm {

JointPositionUpperLimits::JointPositionUpperLimits(const Robot& robot, 
                                                   const double barrier)
  : dimq_(robot.dimq()),
    dimv_(robot.dimv()),
    dimc_(robot.upperJointPositionLimit().size()),
    barrier_(barrier),
    qmax_(robot.upperJointPositionLimit()),
    slack_(qmax_-Eigen::VectorXd::Constant(qmax_.size(), barrier)),
    dual_(Eigen::VectorXd::Constant(qmax_.size(), barrier)),
    residual_(Eigen::VectorXd::Zero(qmax_.size())),
    duality_(Eigen::VectorXd::Zero(qmax_.size())),
    dslack_(Eigen::VectorXd::Zero(qmax_.size())), 
    ddual_(Eigen::VectorXd::Zero(qmax_.size())) {
  assert(barrier_ > 0);
}


bool JointPositionUpperLimits::isFeasible(const Robot& robot, 
                                          const Eigen::VectorXd& q) {
  assert(q.size() == robot.dimq());
  for (int i=0; i<dimc_; ++i) {
    if (q.coeff(i) > qmax_.coeff(i)) {
      return false;
    }
  }
  return true;
}


void JointPositionUpperLimits::setSlackAndDual(const Robot& robot, 
                                               const double dtau, 
                                               const Eigen::VectorXd& q) {
  assert(dtau > 0);
  assert(q.size() == robot.dimq());
  slack_ = dtau * (qmax_-q);
  pdipmfunc::SetSlackAndDualPositive(dimc_, barrier_, slack_, dual_);
}


void JointPositionUpperLimits::condenseSlackAndDual(const Robot& robot, 
                                                    const double dtau, 
                                                    const Eigen::VectorXd& q,
                                                    Eigen::MatrixXd& Cqq, 
                                                    Eigen::VectorXd& Cq) {
  assert(dtau > 0);
  assert(Cqq.rows() == robot.dimv());
  assert(Cqq.cols() == robot.dimv());
  assert(Cq.size() == robot.dimv());
  for (int i=0; i<dimv_; ++i) {
    Cqq.coeffRef(i, i) += dtau * dtau * dual_.coeff(i) / slack_.coeff(i);
  }
  residual_ = dtau * (q-qmax_) + slack_;
  pdipmfunc::ComputeDualityResidual(barrier_, slack_, dual_, duality_);
  Cq.array() += dtau * (dual_.array()*residual_.array()-duality_.array()) / slack_.array();
}


void JointPositionUpperLimits::computeSlackAndDualDirection(
    const Robot& robot, const double dtau, const Eigen::VectorXd& dq) {
  assert(dtau > 0);
  assert(dq.size() == robot.dimv());
  dslack_ = - dtau * dq - residual_;
  pdipmfunc::ComputeDualDirection(dual_, slack_, dslack_, duality_, ddual_);
}


double JointPositionUpperLimits::maxSlackStepSize(const double margin_rate) {
  assert(margin_rate > 0);
  return pdipmfunc::FractionToBoundary(dimc_, margin_rate, slack_, dslack_);
}


double JointPositionUpperLimits::maxDualStepSize(const double margin_rate) {
  assert(margin_rate > 0);
  return pdipmfunc::FractionToBoundary(dimc_, margin_rate, dual_, ddual_);
}


void JointPositionUpperLimits::updateSlack(const double step_size) {
  assert(step_size > 0);
  slack_.noalias() += step_size * dslack_;
}


void JointPositionUpperLimits::updateDual(const double step_size) {
  assert(step_size > 0);
  dual_.noalias() += step_size * ddual_;
}


double JointPositionUpperLimits::costSlackBarrier() {
  return pdipmfunc::SlackBarrierCost(dimc_, barrier_, slack_);
}


double JointPositionUpperLimits::costSlackBarrier(const double step_size) {
  return pdipmfunc::SlackBarrierCost(dimc_, barrier_, slack_+step_size*dslack_);
}


void JointPositionUpperLimits::augmentDualResidual(const Robot& robot, 
                                                   const double dtau, 
                                                   Eigen::VectorXd& Cq) {
  assert(dtau > 0);
  assert(Cq.size() == robot.dimv());
  Cq.noalias() += dtau * dual_;
}


double JointPositionUpperLimits::residualL1Nrom(const Robot& robot, 
                                                const double dtau, 
                                                const Eigen::VectorXd& q) {
  assert(dtau > 0);
  assert(q.size() == robot.dimq());
  residual_ = dtau * (q-qmax_) + slack_;
  return residual_.lpNorm<1>();
}


double JointPositionUpperLimits::residualSquaredNrom(const Robot& robot, 
                                                     const double dtau, 
                                                     const Eigen::VectorXd& q) {
  assert(dtau > 0);
  assert(q.size() == robot.dimq());
  residual_ = dtau * (q-qmax_) + slack_;
  double error = 0;
  error += residual_.squaredNorm();
  pdipmfunc::ComputeDualityResidual(barrier_, slack_, dual_, residual_);
  error += residual_.squaredNorm();
  return error;
}

} // namespace pdipm
} // namespace idocp