#ifndef IDOCP_ROBOT_HXX_
#define IDOCP_ROBOT_HXX_

#include "idocp/robot/robot.hpp"

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/kinematics-derivatives.hpp"
#include "pinocchio/algorithm/frames.hpp"
#include "pinocchio/algorithm/frames-derivatives.hpp"
#include "pinocchio/algorithm/crba.hpp"
#include "pinocchio/algorithm/rnea.hpp"
#include "pinocchio/algorithm/rnea-derivatives.hpp"
#include "pinocchio/algorithm/aba.hpp"
#include "pinocchio/algorithm/cholesky.hpp"
#include "pinocchio/algorithm/compute-all-terms.hpp"
#include "pinocchio/algorithm/contact-dynamics.hpp"

#include <stdexcept>
#include <cassert>

namespace idocp {

template <typename TangentVectorType, typename ConfigVectorType>
inline void Robot::integrateConfiguration(
    const Eigen::MatrixBase<TangentVectorType>& v, 
    const double integration_length, 
    const Eigen::MatrixBase<ConfigVectorType>& q) const {
  assert(v.size() == dimv_);
  assert(q.size() == dimq_);
  if (floating_base_.hasFloatingBase()) {
    const Eigen::VectorXd q_tmp = q;
    pinocchio::integrate(model_, q_tmp, integration_length*v, 
                         const_cast<Eigen::MatrixBase<ConfigVectorType>&>(q));
  }
  else {
    (const_cast<Eigen::MatrixBase<ConfigVectorType>&>(q)).noalias() 
        += integration_length * v;
  }
}


template <typename ConfigVectorType1, typename TangentVectorType,  
          typename ConfigVectorType2>
inline void Robot::integrateConfiguration(
    const Eigen::MatrixBase<ConfigVectorType1>& q, 
    const Eigen::MatrixBase<TangentVectorType>& v, 
    const double integration_length, 
    const Eigen::MatrixBase<ConfigVectorType2>& q_integrated) const {
  assert(q.size() == dimq_);
  assert(v.size() == dimv_);
  assert(q_integrated.size() == dimq_);
  if (floating_base_.hasFloatingBase()) {
    pinocchio::integrate(
        model_, q, integration_length*v, 
        const_cast<Eigen::MatrixBase<ConfigVectorType2>&>(q_integrated));
  }
  else {
    (const_cast<Eigen::MatrixBase<ConfigVectorType2>&>(q_integrated)).noalias() 
        = q + integration_length * v;
  }
}


template <typename ConfigVectorType, typename TangentVectorType, 
          typename MatrixType>
inline void Robot::dIntegratedConfiguration(
    const Eigen::MatrixBase<ConfigVectorType>& q,
    const Eigen::MatrixBase<TangentVectorType>& v,
    const Eigen::MatrixBase<MatrixType>& dIntegrate_dq) const {
  assert(q.size() == dimq_);
  assert(v.size() == dimv_);
  assert(dIntegrate_dq.rows() == dimv_);
  assert(dIntegrate_dq.cols() == dimv_);
  pinocchio::dIntegrate(
      model_, q, v, const_cast<Eigen::MatrixBase<MatrixType>&>(dIntegrate_dq),
      pinocchio::ARG0);
}


template <typename ConfigVectorType, typename TangentVectorType, 
          typename MatrixType>
inline void Robot::dIntegratedVelocity(
    const Eigen::MatrixBase<ConfigVectorType>& q,
    const Eigen::MatrixBase<TangentVectorType>& v,
    const Eigen::MatrixBase<MatrixType>& dIntegrate_dv) const {
  assert(q.size() == dimq_);
  assert(v.size() == dimv_);
  assert(dIntegrate_dv.rows() == dimv_);
  assert(dIntegrate_dv.cols() == dimv_);
  pinocchio::dIntegrate(
      model_, q, v, const_cast<Eigen::MatrixBase<MatrixType>&>(dIntegrate_dv),
      pinocchio::ARG1);
}


template <typename ConfigVectorType1, typename ConfigVectorType2, 
          typename TangentVectorType>
inline void Robot::subtractConfiguration(
    const Eigen::MatrixBase<ConfigVectorType1>& q_plus, 
    const Eigen::MatrixBase<ConfigVectorType2>& q_minus,
    const Eigen::MatrixBase<TangentVectorType>& difference) const {
  assert(q_plus.size() == dimq_);
  assert(q_minus.size() == dimq_);
  assert(difference.size() == dimv_);
  if (floating_base_.hasFloatingBase()) {
    pinocchio::difference(
        model_, q_minus, q_plus, 
        const_cast<Eigen::MatrixBase<TangentVectorType>&>(difference));
  }
  else {
    const_cast<Eigen::MatrixBase<TangentVectorType>&>(difference) 
        = q_plus - q_minus;
  }
}


template <typename ConfigVectorType1, typename ConfigVectorType2, 
          typename MatrixType>
inline void Robot::dSubtractdConfigurationPlus(
    const Eigen::MatrixBase<ConfigVectorType1>& q_plus,
    const Eigen::MatrixBase<ConfigVectorType2>& q_minus,
    const Eigen::MatrixBase<MatrixType>& dSubtract_dqplus) const {
  assert(q_plus.size() == dimq_);
  assert(q_minus.size() == dimq_);
  assert(dSubtract_dqplus.rows() == dimv_);
  assert(dSubtract_dqplus.cols() == dimv_);
  pinocchio::dDifference(
      model_, q_minus, q_plus, 
      const_cast<Eigen::MatrixBase<MatrixType>&>(dSubtract_dqplus),
      pinocchio::ARG1);
}


template <typename ConfigVectorType1, typename ConfigVectorType2, 
          typename MatrixType>
inline void Robot::dSubtractdConfigurationMinus(
    const Eigen::MatrixBase<ConfigVectorType1>& q_plus,
    const Eigen::MatrixBase<ConfigVectorType2>& q_minus,
    const Eigen::MatrixBase<MatrixType>& dSubtract_dqminus) const {
  assert(q_plus.size() == dimq_);
  assert(q_minus.size() == dimq_);
  assert(dSubtract_dqminus.rows() == dimv_);
  assert(dSubtract_dqminus.cols() == dimv_);
  pinocchio::dDifference(
      model_, q_minus, q_plus, 
      const_cast<Eigen::MatrixBase<MatrixType>&>(dSubtract_dqminus),
      pinocchio::ARG0);
}


template <typename MatrixType1, typename MatrixType2>
inline void Robot::dSubtractdConfigurationInverse(
    const Eigen::MatrixBase<MatrixType1>& dSubtract_dq,
    const Eigen::MatrixBase<MatrixType2>& dSubtract_dq_inv) {
  const_cast<Eigen::MatrixBase<MatrixType2>&>(dSubtract_dq_inv).template topLeftCorner<3, 3>().noalias()
      = dSubtract_dq.template topLeftCorner<3, 3>().inverse();
  const_cast<Eigen::MatrixBase<MatrixType2>&>(dSubtract_dq_inv).template block<3, 3>(3, 3).noalias()
      = dSubtract_dq.template block<3, 3>(3, 3).inverse();
  mat_3d_.noalias() = dSubtract_dq.template block<3, 3>(0, 3) 
                        * dSubtract_dq_inv.template block<3, 3>(3, 3);
  const_cast<Eigen::MatrixBase<MatrixType2>&>(dSubtract_dq_inv).template block<3, 3>(0, 3).noalias()
      = - dSubtract_dq_inv.template topLeftCorner<3, 3>() * mat_3d_;
}


inline const Eigen::Vector3d& Robot::framePosition(const int frame_id) const {
  return data_.oMf[frame_id].translation();
}


inline const Eigen::Matrix3d& Robot::frameRotation(const int frame_id) const {
  return data_.oMf[frame_id].rotation();
}


inline const pinocchio::SE3& Robot::framePlacement(const int frame_id) const {
  return data_.oMf[frame_id];
}


template <typename MatrixType>
inline void Robot::getFrameJacobian(const int frame_id, 
                                    const Eigen::MatrixBase<MatrixType>& J) {
  assert(J.rows() == 6);
  assert(J.cols() == dimv_);
  pinocchio::getFrameJacobian(model_, data_, frame_id, pinocchio::LOCAL, 
                              const_cast<Eigen::MatrixBase<MatrixType>&>(J));
}


template <typename ConfigVectorType, typename TangentVectorType1, 
          typename TangentVectorType2>
inline void Robot::updateKinematics(
    const Eigen::MatrixBase<ConfigVectorType>& q, 
    const Eigen::MatrixBase<TangentVectorType1>& v, 
    const Eigen::MatrixBase<TangentVectorType2>& a) {
  assert(q.size() == dimq_);
  assert(v.size() == dimv_);
  assert(a.size() == dimv_);
  pinocchio::forwardKinematics(model_, data_, q, v, a);
  pinocchio::updateFramePlacements(model_, data_);
  pinocchio::computeForwardKinematicsDerivatives(model_, data_, q, v, a);
  pinocchio::jacobianCenterOfMass(model_, data_, false);
}


template <typename ConfigVectorType, typename TangentVectorType>
inline void Robot::updateKinematics(
    const Eigen::MatrixBase<ConfigVectorType>& q, 
    const Eigen::MatrixBase<TangentVectorType>& v) {
  assert(q.size() == dimq_);
  assert(v.size() == dimv_);
  pinocchio::forwardKinematics(model_, data_, q, v);
  pinocchio::updateFramePlacements(model_, data_);
  pinocchio::computeForwardKinematicsDerivatives(model_, data_, q, v, 
                                                 Eigen::VectorXd::Zero(dimv_));
  pinocchio::jacobianCenterOfMass(model_, data_, false);
}


template <typename ConfigVectorType>
inline void Robot::updateKinematics(
    const Eigen::MatrixBase<ConfigVectorType>& q) {
  assert(q.size() == dimq_);
  pinocchio::computeJointJacobians(model_, data_, q);
  pinocchio::framesForwardKinematics(model_, data_, q);
  pinocchio::jacobianCenterOfMass(model_, data_, false);
}


template <typename ConfigVectorType>
inline void Robot::updateFrameKinematics(
    const Eigen::MatrixBase<ConfigVectorType>& q) {
  assert(q.size() == dimq_);
  pinocchio::framesForwardKinematics(model_, data_, q);
}


template <typename VectorType>
inline void Robot::computeBaumgarteResidual(
    const ContactStatus& contact_status, const double time_step,
    const std::vector<Eigen::Vector3d>& contact_points,
    const Eigen::MatrixBase<VectorType>& baumgarte_residual) const {
  assert(time_step > 0);
  assert(contact_points.size() == point_contacts_.size());
  int num_active_contacts = 0;
  for (int i=0; i<point_contacts_.size(); ++i) {
    if (contact_status.isContactActive(i)) {
      point_contacts_[i].computeBaumgarteResidual(
          model_, data_, 
          baumgarte_weight_on_velocity_, baumgarte_weight_on_position_, 
          contact_points[i],
          (const_cast<Eigen::MatrixBase<VectorType>&>(baumgarte_residual))
              .template segment<3>(3*num_active_contacts));
      ++num_active_contacts;
    }
  }
}


template <typename MatrixType1, typename MatrixType2, typename MatrixType3>
inline void Robot::computeBaumgarteDerivatives(
    const ContactStatus& contact_status, const double time_step,
    const Eigen::MatrixBase<MatrixType1>& baumgarte_partial_dq, 
    const Eigen::MatrixBase<MatrixType2>& baumgarte_partial_dv, 
    const Eigen::MatrixBase<MatrixType3>& baumgarte_partial_da) {
  assert(baumgarte_partial_dq.cols() == dimv_);
  assert(baumgarte_partial_dv.cols() == dimv_);
  assert(baumgarte_partial_da.cols() == dimv_);
  int num_active_contacts = 0;
  for (int i=0; i<point_contacts_.size(); ++i) {
    if (contact_status.isContactActive(i)) {
      point_contacts_[i].computeBaumgarteDerivatives(
          model_, data_, 
          baumgarte_weight_on_velocity_, baumgarte_weight_on_position_,
          (const_cast<Eigen::MatrixBase<MatrixType1>&>(baumgarte_partial_dq))
              .block(3*num_active_contacts, 0, 3, dimv_),
          (const_cast<Eigen::MatrixBase<MatrixType2>&>(baumgarte_partial_dv))
              .block(3*num_active_contacts, 0, 3, dimv_),
          (const_cast<Eigen::MatrixBase<MatrixType3>&>(baumgarte_partial_da))
              .block(3*num_active_contacts, 0, 3, dimv_));
      ++num_active_contacts;
    }
  }
}


template <typename VectorType>
inline void Robot::computeImpulseVelocityResidual(
    const ImpulseStatus& impulse_status, 
    const Eigen::MatrixBase<VectorType>& velocity_residual) const {
  int num_active_impulse = 0;
  for (int i=0; i<point_contacts_.size(); ++i) {
    if (impulse_status.isImpulseActive(i)) {
      point_contacts_[i].computeContactVelocityResidual(
          model_, data_,
          (const_cast<Eigen::MatrixBase<VectorType>&>(velocity_residual))
              .template segment<3>(3*num_active_impulse));
      ++num_active_impulse;
    }
  }
}


template <typename MatrixType1, typename MatrixType2>
inline void Robot::computeImpulseVelocityDerivatives(
    const ImpulseStatus& impulse_status, 
    const Eigen::MatrixBase<MatrixType1>& velocity_partial_dq, 
    const Eigen::MatrixBase<MatrixType2>& velocity_partial_dv) {
  assert(velocity_partial_dq.cols() == dimv_);
  assert(velocity_partial_dv.cols() == dimv_);
  int num_active_impulse = 0;
  for (int i=0; i<point_contacts_.size(); ++i) {
    if (impulse_status.isImpulseActive(i)) {
      point_contacts_[i].computeContactVelocityDerivatives(
          model_, data_, 
          (const_cast<Eigen::MatrixBase<MatrixType1>&>(velocity_partial_dq))
              .block(3*num_active_impulse, 0, 3, dimv_),
          (const_cast<Eigen::MatrixBase<MatrixType2>&>(velocity_partial_dv))
              .block(3*num_active_impulse, 0, 3, dimv_));
      ++num_active_impulse;
    }
  }
}


template <typename VectorType>
inline void Robot::computeContactResidual(
    const ImpulseStatus& impulse_status, 
    const std::vector<Eigen::Vector3d>& contact_points,
    const Eigen::MatrixBase<VectorType>& contact_residual) const {
  assert(contact_points.size() == point_contacts_.size());
  int num_active_impulse = 0;
  for (int i=0; i<point_contacts_.size(); ++i) {
    if (impulse_status.isImpulseActive(i)) {
      point_contacts_[i].computeContactResidual(
          model_, data_, contact_points[i],
          (const_cast<Eigen::MatrixBase<VectorType>&>(contact_residual))
              .template segment<3>(3*num_active_impulse));
      ++num_active_impulse;
    }
  }
}


template <typename MatrixType>
inline void Robot::computeContactDerivative(
    const ImpulseStatus& impulse_status, 
    const Eigen::MatrixBase<MatrixType>& contact_partial_dq) {
  assert(contact_partial_dq.cols() == dimv_);
  int num_active_impulse = 0;
  for (int i=0; i<point_contacts_.size(); ++i) {
    if (impulse_status.isImpulseActive(i)) {
      point_contacts_[i].computeContactDerivative(
          model_, data_, 
          (const_cast<Eigen::MatrixBase<MatrixType>&>(contact_partial_dq))
              .block(3*num_active_impulse, 0, 3, dimv_));
      ++num_active_impulse;
    }
  }
}


inline void Robot::setFrictionCoefficient(
    const std::vector<double>& friction_coefficient) {
  for (int i=0; i<point_contacts_.size(); ++i) {
    point_contacts_[i].setFrictionCoefficient(friction_coefficient[i]);
  }
}


inline double Robot::frictionCoefficient(const int contact_index) const {
  assert(contact_index >= 0);
  assert(contact_index < point_contacts_.size());
  try {
    if (point_contacts_.empty()) {
      throw std::runtime_error(
          "invalid function call: robot has no point contacts!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  return point_contacts_[contact_index].frictionCoefficient();
}


inline void Robot::setRestitutionCoefficient(
    const std::vector<double>& restitution_coefficient) {
  for (int i=0; i<point_contacts_.size(); ++i) {
    point_contacts_[i].setRestitutionCoefficient(restitution_coefficient[i]);
  }
}


inline double Robot::restitutionCoefficient(const int contact_index) const {
  assert(contact_index >= 0);
  assert(contact_index < point_contacts_.size());
  try {
    if (point_contacts_.empty()) {
      throw std::runtime_error(
          "invalid function call: robot has no point contacts!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  return point_contacts_[contact_index].restitutionCoefficient();
}


inline void Robot::setContactForces(const ContactStatus& contact_status, 
                                    const std::vector<Eigen::Vector3d>& f) {
  assert(f.size() == point_contacts_.size());
  int num_active_contacts = 0;
  for (int i=0; i<point_contacts_.size(); ++i) {
    if (contact_status.isContactActive(i)) {
      point_contacts_[i].computeJointForceFromContactForce(f[i], fjoint_);
      ++num_active_contacts;
    }
    else {
      point_contacts_[i].computeJointForceFromContactForce(
          Eigen::Vector3d::Zero(), fjoint_);
    }
  }
}


inline void Robot::setImpulseForces(const ImpulseStatus& impulse_status, 
                                    const std::vector<Eigen::Vector3d>& f) {
  assert(f.size() == point_contacts_.size());
  int num_active_impulse = 0;
  for (int i=0; i<point_contacts_.size(); ++i) {
    if (impulse_status.isImpulseActive(i)) {
      point_contacts_[i].computeJointForceFromContactForce(f[i], fjoint_);
      ++num_active_impulse;
    }
    else {
      point_contacts_[i].computeJointForceFromContactForce(
          Eigen::Vector3d::Zero(), fjoint_);
    }
  }
}


template <typename ConfigVectorType, typename TangentVectorType1, 
          typename TangentVectorType2, typename TangentVectorType3>
inline void Robot::RNEA(const Eigen::MatrixBase<ConfigVectorType>& q, 
                        const Eigen::MatrixBase<TangentVectorType1>& v, 
                        const Eigen::MatrixBase<TangentVectorType2>& a, 
                        const Eigen::MatrixBase<TangentVectorType3>& tau) {
  assert(q.size() == dimq_);
  assert(v.size() == dimv_);
  assert(a.size() == dimv_);
  assert(tau.size() == dimv_);
  if (point_contacts_.empty()) {
    const_cast<Eigen::MatrixBase<TangentVectorType3>&>(tau)
        = pinocchio::rnea(model_, data_, q, v, a);
  }
  else {
    const_cast<Eigen::MatrixBase<TangentVectorType3>&>(tau)
        = pinocchio::rnea(model_, data_, q, v, a, fjoint_);
  }
}


template <typename ConfigVectorType, typename TangentVectorType1, 
          typename TangentVectorType2, typename MatrixType1, 
          typename MatrixType2, typename MatrixType3>
inline void Robot::RNEADerivatives(
    const Eigen::MatrixBase<ConfigVectorType>& q, 
    const Eigen::MatrixBase<TangentVectorType1>& v, 
    const Eigen::MatrixBase<TangentVectorType2>& a, 
    const Eigen::MatrixBase<MatrixType1>& dRNEA_partial_dq, 
    const Eigen::MatrixBase<MatrixType2>& dRNEA_partial_dv, 
    const Eigen::MatrixBase<MatrixType3>& dRNEA_partial_da) {
  assert(q.size() == dimq_);
  assert(v.size() == dimv_);
  assert(a.size() == dimv_);
  assert(dRNEA_partial_dq.cols() == dimv_);
  assert(dRNEA_partial_dq.rows() == dimv_);
  assert(dRNEA_partial_dv.cols() == dimv_);
  assert(dRNEA_partial_dv.rows() == dimv_);
  assert(dRNEA_partial_da.cols() == dimv_);
  assert(dRNEA_partial_da.rows() == dimv_);
  if (point_contacts_.empty()) {
      pinocchio::computeRNEADerivatives(
          model_, data_, q, v, a, 
          const_cast<Eigen::MatrixBase<MatrixType1>&>(dRNEA_partial_dq),
          const_cast<Eigen::MatrixBase<MatrixType2>&>(dRNEA_partial_dv),
          const_cast<Eigen::MatrixBase<MatrixType3>&>(dRNEA_partial_da));
  }
  else {
      pinocchio::computeRNEADerivatives(
          model_, data_, q, v, a, fjoint_,
          const_cast<Eigen::MatrixBase<MatrixType1>&>(dRNEA_partial_dq),
          const_cast<Eigen::MatrixBase<MatrixType2>&>(dRNEA_partial_dv),
          const_cast<Eigen::MatrixBase<MatrixType3>&>(dRNEA_partial_da));
  }
  (const_cast<Eigen::MatrixBase<MatrixType3>&>(dRNEA_partial_da)) 
      .template triangularView<Eigen::StrictlyLower>() 
      = (const_cast<Eigen::MatrixBase<MatrixType3>&>(dRNEA_partial_da)).transpose()
          .template triangularView<Eigen::StrictlyLower>();
}


template <typename ConfigVectorType, typename TangentVectorType1, 
          typename TangentVectorType2>
inline void Robot::RNEAImpulse(
    const Eigen::MatrixBase<ConfigVectorType>& q, 
    const Eigen::MatrixBase<TangentVectorType1>& dv,
    const Eigen::MatrixBase<TangentVectorType2>& res) {
  assert(q.size() == dimq_);
  assert(dv.size() == dimv_);
  assert(res.size() == dimv_);
  const_cast<Eigen::MatrixBase<TangentVectorType2>&>(res)
      = pinocchio::rnea(impulse_model_, impulse_data_, q, 
                        Eigen::VectorXd::Zero(dimv_),  dv, fjoint_);
}


template <typename ConfigVectorType, typename TangentVectorType, 
          typename MatrixType1, typename MatrixType2>
inline void Robot::RNEAImpulseDerivatives(
    const Eigen::MatrixBase<ConfigVectorType>& q, 
    const Eigen::MatrixBase<TangentVectorType>& dv, 
    const Eigen::MatrixBase<MatrixType1>& dRNEA_partial_dq, 
    const Eigen::MatrixBase<MatrixType2>& dRNEA_partial_ddv) {
  assert(q.size() == dimq_);
  assert(dv.size() == dimv_);
  assert(dRNEA_partial_dq.cols() == dimv_);
  assert(dRNEA_partial_dq.rows() == dimv_);
  assert(dRNEA_partial_ddv.cols() == dimv_);
  assert(dRNEA_partial_ddv.rows() == dimv_);
  pinocchio::computeRNEADerivatives(
      impulse_model_, impulse_data_, q, Eigen::VectorXd::Zero(dimv_), dv, 
      fjoint_, const_cast<Eigen::MatrixBase<MatrixType1>&>(dRNEA_partial_dq),
      dimpulse_dv_,
      const_cast<Eigen::MatrixBase<MatrixType2>&>(dRNEA_partial_ddv));
  (const_cast<Eigen::MatrixBase<MatrixType2>&>(dRNEA_partial_ddv)) 
      .template triangularView<Eigen::StrictlyLower>() 
      = (const_cast<Eigen::MatrixBase<MatrixType2>&>(dRNEA_partial_ddv)).transpose()
          .template triangularView<Eigen::StrictlyLower>();
}


template <typename MatrixType>
inline void Robot::dRNEAPartialdFext(
    const ContactStatus& contact_status,
    const Eigen::MatrixBase<MatrixType>& dRNEA_partial_dfext) {
  assert(dRNEA_partial_dfext.rows() == dimv_);
  int num_active_contacts = 0;
  for (int i=0; i<point_contacts_.size(); ++i) {
    if (contact_status.isContactActive(i)) {
      point_contacts_[i].getContactJacobian(
          model_, data_,  -1, 
          (const_cast<Eigen::MatrixBase<MatrixType>&>(dRNEA_partial_dfext))
              .block(0, 3*num_active_contacts, dimv_, 3), true);
      ++num_active_contacts;
    }
  }
}


template <typename ConfigVectorType, typename TangentVectorType>
inline void Robot::computeAllTerms(const Eigen::MatrixBase<ConfigVectorType>& q, 
                                   const Eigen::MatrixBase<TangentVectorType>& v) {
  pinocchio::computeAllTerms(model_, data_, q, v);
  pinocchio::framesForwardKinematics(model_, data_, q);
}


template <typename MatrixType>
inline void Robot::getContactJacobian(const ContactStatus& contact_status,
                                      const Eigen::MatrixBase<MatrixType>& J) {
  assert(J.cols() == dimv_);
  int num_active_contacts = 0;
  for (int i=0; i<point_contacts_.size(); ++i) {
    if (contact_status.isContactActive(i)) {
      point_contacts_[i].getContactJacobian(
          model_, data_, 
          (const_cast<Eigen::MatrixBase<MatrixType>&>(J))
              .block(3*num_active_contacts, 0, 3, dimv_), false);
      ++num_active_contacts;
    }
  }
}


template <typename MatrixType>
inline void Robot::getContactJacobian(const ImpulseStatus& impulse_status,
                                      const Eigen::MatrixBase<MatrixType>& J) {
  assert(J.cols() == dimv_);
  int num_active_contacts = 0;
  for (int i=0; i<point_contacts_.size(); ++i) {
    if (impulse_status.isImpulseActive(i)) {
      point_contacts_[i].getContactJacobian(
          model_, data_, 
          (const_cast<Eigen::MatrixBase<MatrixType>&>(J))
              .block(3*num_active_contacts, 0, 3, dimv_), false);
      ++num_active_contacts;
    }
  }
}


template <typename TangentVectorType, typename MatrixType, typename DriftVectorType, 
          typename ResultVectorType1, typename ResultVectorType2>
inline void Robot::forwardDynamics(const Eigen::MatrixBase<TangentVectorType>& tau, 
                                   const Eigen::MatrixBase<MatrixType>& J, 
                                   const Eigen::MatrixBase<DriftVectorType>& gmm,
                                   const Eigen::MatrixBase<ResultVectorType1>& a,
                                   const Eigen::MatrixBase<ResultVectorType2>& f) {
  const int dimf = f.size();

  // Compute the UDUt decomposition of data.M
  pinocchio::cholesky::decompose(model_, data_);

  // Compute the dynamic drift (control - nle)
  data_.torque_residual = tau - data_.nle;
  pinocchio::cholesky::solve(model_, data_, data_.torque_residual);

  data_.sDUiJt.leftCols(dimf) = J.transpose();
  // Compute U^-1 * J.T
  pinocchio::cholesky::Uiv(model_, data_, data_.sDUiJt);
  for (Eigen::DenseIndex k=0; k<model_.nv; ++k) {
    data_.sDUiJt.leftCols(dimf).row(k) /= sqrt(data_.D[k]);
  }

  data_.JMinvJt.topLeftCorner(dimf, dimf).noalias() 
      = data_.sDUiJt.leftCols(dimf).transpose() * data_.sDUiJt.leftCols(dimf);

  // data_.JMinvJt.topLeftCorner(dimf, dimf).diagonal().array() += inv_damping;
  data_.llt_JMinvJt.compute(data_.JMinvJt.topLeftCorner(dimf, dimf));

  // Compute the Lagrange Multipliers
  data_.lambda_c.head(dimf).noalias() = - J * data_.torque_residual;
  data_.lambda_c.head(dimf).noalias() -= gmm;
  data_.llt_JMinvJt.solveInPlace(data_.lambda_c.head(dimf));

  // Compute the joint acceleration
  data_.ddq.noalias() = J.transpose() * data_.lambda_c.head(dimf);
  pinocchio::cholesky::solve(model_, data_, data_.ddq);
  data_.ddq.noalias() += data_.torque_residual;

  (const_cast<Eigen::MatrixBase<ResultVectorType1>&>(a)) = data_.ddq;
  (const_cast<Eigen::MatrixBase<ResultVectorType2>&>(f)) = data_.lambda_c;
}


template <typename ConfigVectorType, typename TangentVectorType, typename MatrixType, 
          typename ResultVectorType1, typename ResultVectorType2>
inline void Robot::impulseDynamics(const Eigen::MatrixBase<ConfigVectorType>& q, 
                                   const Eigen::MatrixBase<TangentVectorType>& v, 
                                   const Eigen::MatrixBase<MatrixType>& J, 
                                   const Eigen::MatrixBase<ResultVectorType1>& dv,
                                   const Eigen::MatrixBase<ResultVectorType2>& lmd) {
  pinocchio::crba(model_, data_, q);

  const int dimf = lmd.size();
  // Compute the UDUt decomposition of data.M
  pinocchio::cholesky::decompose(model_, data_);

  data_.sDUiJt.leftCols(dimf) = J.transpose();
  // Compute U^-1 * J.T
  pinocchio::cholesky::Uiv(model_, data_, data_.sDUiJt.leftCols(dimf));
  for(int k=0; k<model_.nv; ++k) {
    data_.sDUiJt.leftCols(dimf).row(k) /= sqrt(data_.D[k]);
  }

  data_.JMinvJt.topLeftCorner(dimf, dimf).noalias() 
      = data_.sDUiJt.leftCols(dimf).transpose() * data_.sDUiJt.leftCols(dimf);

  data_.llt_JMinvJt.compute(data_.JMinvJt.topLeftCorner(dimf, dimf));

  // Compute the Lagrange Multipliers related to the contact impulses
  data_.impulse_c.noalias() = - J * v;
  data_.llt_JMinvJt.solveInPlace(data_.impulse_c.head(dimf));

  // Compute the joint velocity after impacts
  data_.dq_after.noalias() = J.transpose() * data_.impulse_c.head(dimf);
  pinocchio::cholesky::solve(model_, data_, data_.dq_after);
  data_.dq_after.noalias() += v;

  (const_cast<Eigen::MatrixBase<ResultVectorType1>&>(dv)) = data_.dq_after - v;
  (const_cast<Eigen::MatrixBase<ResultVectorType2>&>(lmd)) = data_.impulse_c;
}


template <typename MatrixType1, typename MatrixType2>
inline void Robot::computeMinv(const Eigen::MatrixBase<MatrixType1>& M, 
                               const Eigen::MatrixBase<MatrixType2>& Minv) {
  assert(M.rows() == dimv_);
  assert(M.cols() == dimv_);
  assert(Minv.rows() == dimv_);
  assert(Minv.cols() == dimv_);
  data_.M = M;
  pinocchio::cholesky::decompose(model_, data_);
  pinocchio::cholesky::computeMinv(
      model_, data_, const_cast<Eigen::MatrixBase<MatrixType2>&>(Minv));
}


template <typename MatrixType1, typename MatrixType2, typename MatrixType3>
inline void Robot::computeMJtJinv(
    const Eigen::MatrixBase<MatrixType1>& M, 
    const Eigen::MatrixBase<MatrixType2>& J, 
    const Eigen::MatrixBase<MatrixType3>& MJtJinv) {
  assert(M.rows() == dimv_);
  assert(M.cols() == dimv_);
  assert(J.rows() <= max_dimf_);
  assert(J.cols() == dimv_);
  assert(MJtJinv.rows() == M.rows()+J.rows());
  assert(MJtJinv.cols() == M.rows()+J.rows());
  const int dimf = J.rows();
  data_.M = M;
  pinocchio::cholesky::decompose(model_, data_);
  data_.sDUiJt.leftCols(dimf) = J.transpose();
  pinocchio::cholesky::Uiv(model_, data_, data_.sDUiJt.leftCols(dimf));
  for (Eigen::DenseIndex k=0; k<dimv_; ++k) {
    data_.sDUiJt.leftCols(dimf).row(k) /= std::sqrt(data_.D[k]);
  }
  data_.JMinvJt.topLeftCorner(dimf, dimf).noalias() 
      = data_.sDUiJt.leftCols(dimf).transpose() * data_.sDUiJt.leftCols(dimf);
  data_.llt_JMinvJt.compute(data_.JMinvJt.topLeftCorner(dimf, dimf));
  assert(data_.llt_JMinvJt.info() == Eigen::Success);
  Eigen::Block<MatrixType3> topLeft 
      = const_cast<Eigen::MatrixBase<MatrixType3>&>(MJtJinv).topLeftCorner(dimv_, dimv_);
  Eigen::Block<MatrixType3> topRight 
      = const_cast<Eigen::MatrixBase<MatrixType3>&>(MJtJinv).topRightCorner(dimv_, dimf);
  Eigen::Block<MatrixType3> bottomLeft 
      = const_cast<Eigen::MatrixBase<MatrixType3>&>(MJtJinv).bottomLeftCorner(dimf, dimv_);
  Eigen::Block<MatrixType3> bottomRight 
      = const_cast<Eigen::MatrixBase<MatrixType3>&>(MJtJinv).bottomRightCorner(dimf, dimf);
  bottomRight = - pinocchio::Data::MatrixXs::Identity(dimf, dimf);
  topLeft.setIdentity();
  data_.llt_JMinvJt.solveInPlace(bottomRight);
  pinocchio::cholesky::solve(model_, data_, topLeft);
  bottomLeft.noalias() = J * topLeft;
  topRight.noalias() = bottomLeft.transpose() * (-bottomRight);
  topLeft.noalias() -= topRight*bottomLeft;
  bottomLeft = topRight.transpose();
  assert(!MJtJinv.hasNaN());
}


inline Eigen::VectorXd Robot::generateFeasibleConfiguration() const {
  Eigen::VectorXd q_min = model_.lowerPositionLimit;
  Eigen::VectorXd q_max = model_.upperPositionLimit;
  if (floating_base_.hasFloatingBase()) {
    q_min.head(7) = - Eigen::VectorXd::Ones(7);
    q_max.head(7) = Eigen::VectorXd::Ones(7);
  }
  return pinocchio::randomConfiguration(model_, q_min, q_max);
}


template <typename ConfigVectorType>
inline void Robot::normalizeConfiguration(
    const Eigen::MatrixBase<ConfigVectorType>& q) const {
  assert(q.size() == dimq_);
  if (floating_base_.hasFloatingBase()) {
    if (q.template segment<4>(3).squaredNorm() 
          <= std::numeric_limits<double>::epsilon()) {
      (const_cast<Eigen::MatrixBase<ConfigVectorType>&> (q)).coeffRef(3) = 1;
    }
    pinocchio::normalize(model_, 
                         const_cast<Eigen::MatrixBase<ConfigVectorType>&>(q));
  }
}


inline Eigen::VectorXd Robot::jointEffortLimit() const {
  return joint_effort_limit_;
}


inline Eigen::VectorXd Robot::jointVelocityLimit() const {
  return joint_velocity_limit_;
}


inline Eigen::VectorXd Robot::lowerJointPositionLimit() const {
  return lower_joint_position_limit_;
}


inline Eigen::VectorXd Robot::upperJointPositionLimit() const {
  return upper_joint_position_limit_;
}


inline int Robot::dimq() const {
  return dimq_;
}


inline int Robot::dimv() const {
  return dimv_;
}


inline int Robot::dimu() const {
  return dimu_;
}


inline int Robot::max_dimf() const {
  return max_dimf_;
}


inline int Robot::dim_passive() const {
  return floating_base_.dim_passive();
}


inline bool Robot::hasFloatingBase() const {
  return floating_base_.hasFloatingBase();
}


inline int Robot::maxPointContacts() const {
  return point_contacts_.size();
}


inline void Robot::initializeJointLimits() {
  const int dim_joint = model_.nv - floating_base_.dim_passive();
  joint_effort_limit_.resize(dim_joint);
  joint_velocity_limit_.resize(dim_joint);
  lower_joint_position_limit_.resize(dim_joint);
  upper_joint_position_limit_.resize(dim_joint);
  joint_effort_limit_ = model_.effortLimit.tail(dim_joint);
  joint_velocity_limit_ = model_.velocityLimit.tail(dim_joint);
  lower_joint_position_limit_ = model_.lowerPositionLimit.tail(dim_joint);
  upper_joint_position_limit_ = model_.upperPositionLimit.tail(dim_joint);
}


inline std::vector<int> Robot::contactFramesIndices() const {
  std::vector<int> contact_frames_indices;
  for (const auto& e : point_contacts_) {
    contact_frames_indices.push_back(e.contact_frame_id());
  }
  return contact_frames_indices;
}


inline ContactStatus Robot::createContactStatus() const {
  return ContactStatus(point_contacts_.size());
}


inline ImpulseStatus Robot::createImpulseStatus() const {
  return ImpulseStatus(point_contacts_.size());
}


template <typename ContactStatusType>
inline void Robot::setContactPoints(ContactStatusType& contact_status) const {
  for (int i=0; i<point_contacts_.size(); ++i) {
    contact_status.setContactPoint(i, point_contacts_[i].contactPoint(data_));
  }
}


inline void Robot::getContactPoints(
    std::vector<Eigen::Vector3d>& contact_points) const {
  for (int i=0; i<point_contacts_.size(); ++i) {
    contact_points[i] = point_contacts_[i].contactPoint(data_);
  }
}


inline double Robot::totalWeight() const {
  return (- pinocchio::computeTotalMass(model_) * model_.gravity981.coeff(2));
}


template <typename VectorType>
inline void Robot::get_com(const Eigen::MatrixBase<VectorType>& com) const {
  const_cast<Eigen::MatrixBase<VectorType>&>(com) = data_.com[0];
}

inline const Eigen::Vector3d& Robot::com() const {
  return data_.com[0];
}

template <typename MatrixType>
inline void Robot::get_Jcom(const Eigen::MatrixBase<MatrixType>& J_com) const {
  const_cast<Eigen::MatrixBase<MatrixType>&>(J_com) = data_.Jcom;
}

} // namespace idocp

#endif // IDOCP_ROBOT_HXX_ 