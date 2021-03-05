#ifndef IDOCP_IMPULSE_DYNAMICS_FORWARD_EULER_HPP_
#define IDOCP_IMPULSE_DYNAMICS_FORWARD_EULER_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/impulse_status.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"
#include "idocp/impulse/impulse_split_kkt_residual.hpp"
#include "idocp/impulse/impulse_split_kkt_matrix.hpp"
#include "idocp/impulse/impulse_dynamics_forward_euler_data.hpp"


namespace idocp {

class ImpulseDynamicsForwardEuler {
public:
  ImpulseDynamicsForwardEuler(const Robot& robot, const double p=1.0e07);

  ImpulseDynamicsForwardEuler();

  ~ImpulseDynamicsForwardEuler();

  ImpulseDynamicsForwardEuler(const ImpulseDynamicsForwardEuler&) = default;

  ImpulseDynamicsForwardEuler& operator=(const ImpulseDynamicsForwardEuler&) 
      = default;
 
  ImpulseDynamicsForwardEuler(ImpulseDynamicsForwardEuler&&) noexcept = default;

  ImpulseDynamicsForwardEuler& operator=(ImpulseDynamicsForwardEuler&&) noexcept 
      = default;

  void impulseDynamics(Robot& robot, const ImpulseStatus& impulse_status, 
                       ImpulseSplitSolution& s);

  void impulseDynamicsDual(Robot& robot, const ImpulseStatus& impulse_status, 
                           const ImpulseSplitKKTResidual& kkt_residual,
                           ImpulseSplitSolution& s);

  void linearizeImpulseDynamics(Robot& robot, 
                                const ImpulseStatus& impulse_status, 
                                const ImpulseSplitSolution& s, 
                                ImpulseSplitKKTMatrix& kkt_matrix, 
                                ImpulseSplitKKTResidual& kkt_residual);

  static void linearizeInverseImpulseDynamics(
      Robot& robot, const ImpulseStatus& impulse_status, 
      const ImpulseSplitSolution& s, ImpulseDynamicsForwardEulerData& data);

  static void linearizeImpulseVelocityConstraint(
      Robot& robot, const ImpulseStatus& impulse_status, 
      ImpulseDynamicsForwardEulerData& data);

  void linearizeSwitchingConstraints(
      Robot& robot, const ImpulseStatus& impulse_status,
      ImpulseSplitKKTResidual& kkt_residual);

  void condenseImpulseDynamics(Robot& robot, 
                               const ImpulseStatus& impulse_status,
                               ImpulseSplitKKTMatrix& kkt_matrix, 
                               ImpulseSplitKKTResidual& kkt_residual);

  static void condensing(const Robot& robot, 
                         const ImpulseStatus& impulse_status,
                         ImpulseDynamicsForwardEulerData& data, 
                         ImpulseSplitKKTMatrix& kkt_matrix, 
                         ImpulseSplitKKTResidual& kkt_residual);

  void computeCondensedPrimalDirection(const Robot& robot, 
                                       ImpulseSplitDirection& d);

  template <typename VectorType>
  void computeCondensedDualDirection(const Robot& robot, 
                                     const ImpulseSplitKKTMatrix& kkt_matrix, 
                                     const ImpulseSplitKKTResidual& kkt_residual, 
                                     const Eigen::MatrixBase<VectorType>& dgmm,
                                     ImpulseSplitDirection& d);

  static void expansionPrimal(const Robot& robot, 
                              const ImpulseDynamicsForwardEulerData& data, 
                              ImpulseSplitDirection& d);

  template <typename VectorType>
  static void expansionDual(const Robot& robot, 
                            ImpulseDynamicsForwardEulerData& data,
                            const ImpulseSplitKKTMatrix& kkt_matrix, 
                            const ImpulseSplitKKTResidual& kkt_residual,
                            const Eigen::MatrixBase<VectorType>& dgmm,
                            ImpulseSplitDirection& d);

  void computeImpulseDynamicsResidual(Robot& robot, 
                                      const ImpulseStatus& impulse_status,
                                      const ImpulseSplitSolution& s, 
                                      ImpulseSplitKKTResidual& kkt_residual);

  void computeSwitchingConstraintsResidual(Robot& robot, 
                                           const ImpulseStatus& impulse_status);

  double l1NormImpulseDynamicsResidual(
      const ImpulseSplitKKTResidual& kkt_residual) const;

  double squaredNormImpulseDynamicsResidual(
      const ImpulseSplitKKTResidual& kkt_residual) const;

//   double l1NormSwitchingConstraintsResidual() const;

//   double squaredNormSwitchingConstraintsResidual() const;

  double penaltyCost() const;

  void set_penalty(const double penalty);

  double get_penalty() const;

private:
  ImpulseDynamicsForwardEulerData data_;
  double p_;

  void setImpulseStatus(const ImpulseStatus& impulse_status);

};

} // namespace idocp 

#include "idocp/impulse/impulse_dynamics_forward_euler.hxx"

#endif // IDOCP_IMPULSE_DYNAMICS_FORWARD_EULER_HPP_ 