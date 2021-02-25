#ifndef IDOCP_IMPULSE_DYNAMICS_FORWARD_EULER_DATA_HXX_
#define IDOCP_IMPULSE_DYNAMICS_FORWARD_EULER_DATA_HXX_

#include "idocp/impulse/impulse_dynamics_forward_euler_data.hpp"

namespace idocp {

inline ImpulseDynamicsForwardEulerData::ImpulseDynamicsForwardEulerData(
    const Robot& robot) 
  : dImDddv(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    dImDCdqv_full_(Eigen::MatrixXd::Zero(robot.dimv()+robot.max_dimf(), 
                                         2*robot.dimv())), 
    MJtJinv_full_(Eigen::MatrixXd::Zero(robot.dimv()+robot.max_dimf(), 
                                        robot.dimv()+robot.max_dimf())), 
    MJtJinv_dImDCdqv_full_(Eigen::MatrixXd::Zero(robot.dimv()+robot.max_dimf(), 
                                                 2*robot.dimv())), 
    Qdvfqv_full_(Eigen::MatrixXd::Zero(robot.dimv()+robot.max_dimf(), 
                                       2*robot.dimv())), 
    ImDC_full_(Eigen::VectorXd::Zero(robot.dimv()+robot.max_dimf())),
    MJtJinv_ImDC_full_(Eigen::VectorXd::Zero(robot.dimv()+robot.max_dimf())),
    ldvf_full_(Eigen::VectorXd::Zero(robot.dimv()+robot.max_dimf())),
    Pq_full_(Eigen::MatrixXd::Zero(robot.max_dimf(), robot.dimv())),
    P_full_(Eigen::VectorXd::Zero(robot.max_dimf())),
    dimv_(robot.dimv()),
    dimf_(0),
    dimvf_(robot.dimv()) {
}


inline ImpulseDynamicsForwardEulerData::ImpulseDynamicsForwardEulerData() 
  : dImDddv(),
    dImDCdqv_full_(), 
    MJtJinv_full_(), 
    MJtJinv_dImDCdqv_full_(), 
    Qdvfqv_full_(), 
    ImDC_full_(),
    MJtJinv_ImDC_full_(),
    ldvf_full_(),
    Pq_full_(),
    P_full_(),
    dimv_(0),
    dimf_(0),
    dimvf_(0) {
}


inline ImpulseDynamicsForwardEulerData::~ImpulseDynamicsForwardEulerData() {
}


inline void ImpulseDynamicsForwardEulerData::setImpulseStatus(
    const ImpulseStatus& impulse_status) {
  dimf_ = impulse_status.dimf();
  dimvf_ = dimv_ + dimf_;
}


inline Eigen::Block<Eigen::MatrixXd> 
ImpulseDynamicsForwardEulerData::dImDCdqv() {
  return dImDCdqv_full_.topLeftCorner(dimvf_, 2*dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseDynamicsForwardEulerData::dImDCdqv() const {
  return dImDCdqv_full_.topLeftCorner(dimvf_, 2*dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> 
ImpulseDynamicsForwardEulerData::dImDCdq() {
  return dImDCdqv_full_.topLeftCorner(dimvf_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseDynamicsForwardEulerData::dImDCdq() const {
  return dImDCdqv_full_.topLeftCorner(dimvf_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> 
ImpulseDynamicsForwardEulerData::dImDdq() {
  return dImDCdqv_full_.topLeftCorner(dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseDynamicsForwardEulerData::dImDdq() const {
  return dImDCdqv_full_.topLeftCorner(dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> 
ImpulseDynamicsForwardEulerData::dCdq() {
  return dImDCdqv_full_.block(dimv_, 0, dimf_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseDynamicsForwardEulerData::dCdq() const {
  return dImDCdqv_full_.block(dimv_, 0, dimf_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> 
ImpulseDynamicsForwardEulerData::dCdv() {
  return dImDCdqv_full_.block(dimv_, dimv_, dimf_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseDynamicsForwardEulerData::dCdv() const {
  return dImDCdqv_full_.block(dimv_, dimv_, dimf_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> 
ImpulseDynamicsForwardEulerData::MJtJinv() {
  return MJtJinv_full_.topLeftCorner(dimvf_, dimvf_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseDynamicsForwardEulerData::MJtJinv() const {
  return MJtJinv_full_.topLeftCorner(dimvf_, dimvf_);
}


inline Eigen::Block<Eigen::MatrixXd> 
ImpulseDynamicsForwardEulerData::MJtJinv_dImDCdqv() {
  return MJtJinv_dImDCdqv_full_.topLeftCorner(dimvf_, 2*dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseDynamicsForwardEulerData::MJtJinv_dImDCdqv() const {
  return MJtJinv_dImDCdqv_full_.topLeftCorner(dimvf_, 2*dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseDynamicsForwardEulerData::Qdvfqv() {
  return Qdvfqv_full_.topLeftCorner(dimvf_, 2*dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseDynamicsForwardEulerData::Qdvfqv() const {
  return Qdvfqv_full_.topLeftCorner(dimvf_, 2*dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
ImpulseDynamicsForwardEulerData::ImDC() {
  return ImDC_full_.head(dimvf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseDynamicsForwardEulerData::ImDC() const {
  return ImDC_full_.head(dimvf_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
ImpulseDynamicsForwardEulerData::ImD() {
  return ImDC_full_.head(dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseDynamicsForwardEulerData::ImD() const {
  return ImDC_full_.head(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
ImpulseDynamicsForwardEulerData::C() {
  return ImDC_full_.segment(dimv_, dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseDynamicsForwardEulerData::C() const {
  return ImDC_full_.segment(dimv_, dimf_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
ImpulseDynamicsForwardEulerData::MJtJinv_ImDC() {
  return MJtJinv_ImDC_full_.head(dimvf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseDynamicsForwardEulerData::MJtJinv_ImDC() const {
  return MJtJinv_ImDC_full_.head(dimvf_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
ImpulseDynamicsForwardEulerData::ldvf() {
  return ldvf_full_.head(dimvf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseDynamicsForwardEulerData::ldvf() const {
  return ldvf_full_.head(dimvf_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
ImpulseDynamicsForwardEulerData::ldv() {
  return ldvf_full_.head(dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseDynamicsForwardEulerData::ldv() const {
  return ldvf_full_.head(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
ImpulseDynamicsForwardEulerData::lf() {
  return ldvf_full_.segment(dimv_, dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseDynamicsForwardEulerData::lf() const {
  return ldvf_full_.segment(dimv_, dimf_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseDynamicsForwardEulerData::Pq() {
  return Pq_full_.topLeftCorner(dimf_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseDynamicsForwardEulerData::Pq() const {
  return Pq_full_.topLeftCorner(dimf_, dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
ImpulseDynamicsForwardEulerData::P() {
  return P_full_.head(dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseDynamicsForwardEulerData::P() const {
  return P_full_.head(dimf_);
}

} // namespace idocp 

#endif // IDOCP_IMPULSE_DYNAMICS_FORWARD_EULER_DATA_HXX_ 