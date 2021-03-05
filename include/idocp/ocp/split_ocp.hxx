#ifndef IDOCP_SPLIT_OCP_HXX_
#define IDOCP_SPLIT_OCP_HXX_

#include "idocp/ocp/split_ocp.hpp"

#include <cassert>

namespace idocp {

inline SplitOCP::SplitOCP(const Robot& robot, 
                          const std::shared_ptr<CostFunction>& cost,
                          const std::shared_ptr<Constraints>& constraints, 
                          const double dt) 
  : cost_(cost),
    cost_data_(cost->createCostFunctionData(robot)),
    constraints_(constraints),
    constraints_data_(constraints->createConstraintsData(robot, 0)),
    contact_dynamics_(robot, dt),
    has_floating_base_(robot.hasFloatingBase()),
    use_kinematics_(false) {
  if (cost_->useKinematics() || constraints_->useKinematics() 
                             || robot.maxPointContacts() > 0) {
    use_kinematics_ = true;
  }
}


inline SplitOCP::SplitOCP() 
  : cost_(),
    cost_data_(),
    constraints_(),
    constraints_data_(),
    contact_dynamics_(),
    has_floating_base_(false),
    use_kinematics_(false) {
}


inline SplitOCP::~SplitOCP() {
}


inline bool SplitOCP::isFeasible(Robot& robot, const SplitSolution& s) {
  return constraints_->isFeasible(robot, constraints_data_, s);
}


inline void SplitOCP::initConstraints(Robot& robot, 
                                      const ContactStatus& contact_status,
                                      const int time_step, SplitSolution& s) { 
  assert(time_step >= 0);
  contact_dynamics_.forwardDynamics(robot, contact_status, s);
  constraints_data_ = constraints_->createConstraintsData(robot, time_step);
  constraints_->setSlackAndDual(robot, constraints_data_, s);
}


template <typename SplitSolutionType>
inline void SplitOCP::linearizeOCP(Robot& robot, 
                                   const ContactStatus& contact_status,  
                                   const double t, const double dt, 
                                   const Eigen::VectorXd& q_prev, 
                                   SplitSolution& s, 
                                   const SplitSolutionType& s_next,
                                   SplitKKTMatrix& kkt_matrix, 
                                   SplitKKTResidual& kkt_residual) {
  assert(dt > 0);
  assert(q_prev.size() == robot.dimq());

  contact_dynamics_.forwardDynamics(robot, contact_status, s);

  kkt_matrix.setContactStatus(contact_status);
  kkt_residual.setContactStatus(contact_status);
  kkt_matrix.setImpulseStatus();
  kkt_residual.setImpulseStatus();
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  kkt_matrix.setZero();
  kkt_residual.setZero();
  cost_->computeStageCostDerivatives(robot, cost_data_, t, dt, s, kkt_residual);
  constraints_->augmentDualResidual(robot, constraints_data_, dt, s, kkt_residual);
  stateequation::linearizeForwardEuler(robot, dt, q_prev, s, s_next, 
                                       kkt_matrix, kkt_residual);
  stateequation::condenseForwardEuler(robot, dt, s, s_next.q, 
                                      kkt_matrix, kkt_residual);
  contact_dynamics_.linearizeContactDynamics(robot, contact_status, dt, s, 
                                             kkt_residual);
  cost_->computeStageCostHessian(robot, cost_data_, t, dt, s, kkt_matrix);
  constraints_->condenseSlackAndDual(robot, constraints_data_, dt, s, 
                                     kkt_matrix, kkt_residual);
  contact_dynamics_.condenseContactDynamics(robot, contact_status, dt, 
                                            kkt_matrix, kkt_residual, true);
}


inline void SplitOCP::computeCondensedPrimalDirection(Robot& robot, 
                                                      const double dt, 
                                                      SplitSolution& s, 
                                                      SplitDirection& d) {
  d.setContactStatusByDimension(s.dimf());
  contact_dynamics_.computeCondensedPrimalDirection(robot, d);
  constraints_->computeSlackAndDualDirection(robot, constraints_data_, s, d);
}


template <typename SplitDirectionType>
inline void SplitOCP::computeCondensedDualDirection(
    const Robot& robot, const double dt, const SplitKKTMatrix& kkt_matrix, 
    SplitKKTResidual& kkt_residual, const SplitDirectionType& d_next, 
    SplitDirection& d) {
  assert(dt > 0);
  contact_dynamics_.computeCondensedDualDirection(robot, dt, kkt_matrix,
                                                  kkt_residual, 
                                                  d_next.dgmm(), d);
  stateequation::correctCostateDirectionForwardEuler(robot, kkt_matrix, 
                                                     kkt_residual, d.dlmd());
}


inline double SplitOCP::maxPrimalStepSize() {
  return constraints_->maxSlackStepSize(constraints_data_);
}


inline double SplitOCP::maxDualStepSize() {
  return constraints_->maxDualStepSize(constraints_data_);
}


inline void SplitOCP::updatePrimal(const Robot& robot, 
                                   const double primal_step_size, 
                                   const SplitDirection& d, 
                                   SplitSolution& s) {
  assert(primal_step_size > 0);
  assert(primal_step_size <= 1);
  s.integrate(robot, primal_step_size, d);
  constraints_->updateSlack(constraints_data_, primal_step_size);
}


inline void SplitOCP::updateDual(const double dual_step_size) {
  assert(dual_step_size > 0);
  assert(dual_step_size <= 1);
  constraints_->updateDual(constraints_data_, dual_step_size);
}


template <typename SplitSolutionType>
inline void SplitOCP::computeKKTResidual(Robot& robot, 
                                         const ContactStatus& contact_status, 
                                         const double t, const double dt, 
                                         const Eigen::VectorXd& q_prev, 
                                         SplitSolution& s,
                                         const SplitSolutionType& s_next, 
                                         SplitKKTMatrix& kkt_matrix,
                                         SplitKKTResidual& kkt_residual) {
  assert(dt > 0);

  contact_dynamics_.forwardDynamics(robot, contact_status, s);

  kkt_matrix.setContactStatus(contact_status);
  kkt_residual.setContactStatus(contact_status);
  kkt_matrix.setImpulseStatus();
  kkt_residual.setImpulseStatus();
  kkt_residual.setZero();
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  cost_->computeStageCostDerivatives(robot, cost_data_, t, dt, s, kkt_residual);
  constraints_->computePrimalAndDualResidual(robot, constraints_data_, s);
  constraints_->augmentDualResidual(robot, constraints_data_, dt, s, kkt_residual);
  stateequation::linearizeForwardEuler(robot, dt, q_prev, s, s_next, 
                                       kkt_matrix, kkt_residual);
  contact_dynamics_.linearizeContactDynamics(robot, contact_status, dt, s, 
                                             kkt_residual);
}


inline double SplitOCP::squaredNormKKTResidual(
    const SplitKKTResidual& kkt_residual, const double dt) const {
  assert(dt > 0);
  double error = 0;
  error += kkt_residual.lx().squaredNorm();
  error += kkt_residual.la.squaredNorm();
  error += kkt_residual.lf().squaredNorm();
  if (has_floating_base_) {
    error += kkt_residual.lu_passive.squaredNorm();
  }
  error += kkt_residual.lu().squaredNorm();
  error += stateequation::squaredNormStateEuqationResidual(kkt_residual);
  error += contact_dynamics_.squaredNormContactDynamicsResidual(dt);
  error += dt * dt * constraints_->squaredNormPrimalAndDualResidual(constraints_data_);
  return error;
}


inline double SplitOCP::stageCost(Robot& robot, const double t,  
                                  const double dt, const SplitSolution& s, 
                                  const double primal_step_size) {
  assert(dt > 0);
  assert(primal_step_size >= 0);
  assert(primal_step_size <= 1);
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  double cost = 0;
  cost += cost_->computeStageCost(robot, cost_data_, t, dt, s);
  if (primal_step_size > 0) {
    cost += dt * constraints_->costSlackBarrier(constraints_data_, 
                                                primal_step_size);
  }
  else {
    cost += dt * constraints_->costSlackBarrier(constraints_data_);
  }
  return cost;
}


inline double SplitOCP::constraintViolation(Robot& robot, 
                                            const ContactStatus& contact_status, 
                                            const double t, const double dt, 
                                            const SplitSolution& s, 
                                            const Eigen::VectorXd& q_next, 
                                            const Eigen::VectorXd& v_next,
                                            SplitKKTResidual& kkt_residual) {
  assert(dt > 0);
  kkt_residual.setContactStatus(contact_status);
  kkt_residual.setImpulseStatus();
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  constraints_->computePrimalAndDualResidual(robot, constraints_data_, s);
  stateequation::computeForwardEulerResidual(robot, dt, s, q_next, v_next, 
                                              kkt_residual);
  contact_dynamics_.computeContactDynamicsResidual(robot, contact_status, s);
  double violation = 0;
  violation += stateequation::l1NormStateEuqationResidual(kkt_residual);
  violation += contact_dynamics_.l1NormContactDynamicsResidual(dt);
  violation += dt * constraints_->l1NormPrimalResidual(constraints_data_);
  return violation;
}


inline double SplitOCP::constraintViolation(Robot& robot, 
                                            const ContactStatus& contact_status, 
                                            const double t, const double dt, 
                                            const SplitSolution& s, 
                                            const Eigen::VectorXd& q_next, 
                                            const Eigen::VectorXd& v_next,
                                            SplitKKTResidual& kkt_residual,
                                            const ImpulseStatus& impulse_status, 
                                            const double dt_next) {
  assert(dt > 0);
  assert(dt_next > 0);
  kkt_residual.setContactStatus(contact_status);
  kkt_residual.setImpulseStatus(impulse_status);
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  constraints_->computePrimalAndDualResidual(robot, constraints_data_, s);
  stateequation::computeForwardEulerResidual(robot, dt, s, q_next, v_next, 
                                             kkt_residual);
  contact_dynamics_.computeContactDynamicsResidual(robot, contact_status, s);
  double violation = 0;
  violation += stateequation::l1NormStateEuqationResidual(kkt_residual);
  violation += contact_dynamics_.l1NormContactDynamicsResidual(dt);
  violation += dt * constraints_->l1NormPrimalResidual(constraints_data_);
  return violation;
}

} // namespace idocp

#endif // IDOCP_SPLIT_OCP_HXX_