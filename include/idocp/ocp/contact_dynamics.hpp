#ifndef IDOCP_CONTACT_DYNAMICS_HPP_
#define IDOCP_CONTACT_DYNAMICS_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/ocp/contact_dynamics_data.hpp"
#include "idocp/ocp/split_state_constraint_jacobian.hpp"

#include <limits>
#include <cmath>


namespace idocp {

class ContactDynamics {
public:
  ContactDynamics(const Robot& robot, const double baumgarte_time_step);

  ContactDynamics();

  ~ContactDynamics();

  ContactDynamics(const ContactDynamics&) = default;

  ContactDynamics& operator=(const ContactDynamics&) = default;
 
  ContactDynamics(ContactDynamics&&) noexcept = default;

  ContactDynamics& operator=(ContactDynamics&&) noexcept = default;

  void forwardDynamics(Robot& robot, const ContactStatus& contact_status, 
                       SplitSolution& s);

  void linearizeContactDynamics(Robot& robot, 
                                const ContactStatus& contact_status, 
                                const double dt, const SplitSolution& s, 
                                SplitKKTResidual& kkt_residual);

  static void linearizeInverseDynamics(Robot& robot, 
                                       const ContactStatus& contact_status,
                                       const SplitSolution& s, 
                                       ContactDynamicsData& data);

  static void linearizeContactConstraint(Robot& robot, 
                                         const ContactStatus& contact_status, 
                                         const double baumgarte_time_step, 
                                         ContactDynamicsData& data);

  void condenseContactDynamics(Robot& robot, 
                               const ContactStatus& contact_status, 
                               const double dt, SplitKKTMatrix& kkt_matrix, 
                               SplitKKTResidual& kkt_residual,
                               const bool is_forward_euler);

  void computeCondensedPrimalDirection(const Robot& robot, 
                                       SplitDirection& d) const;

  template <typename VectorType>
  void computeCondensedDualDirection(const Robot& robot, const double dt, 
                                     const SplitKKTMatrix& kkt_matrix, 
                                     const SplitKKTResidual& kkt_residual, 
                                     const Eigen::MatrixBase<VectorType>& dgmm,
                                     SplitDirection& d);

  void condenseSwitchingConstraint(SplitKKTResidual& kkt_residual, 
                                   SplitStateConstraintJacobian& jac) const;

  void computeContactDynamicsResidual(Robot& robot, 
                                      const ContactStatus& contact_status,
                                      const SplitSolution& s);

  double l1NormContactDynamicsResidual(const double dt) const;

  double squaredNormContactDynamicsResidual(const double dt) const;

private:
  ContactDynamicsData data_;
  bool has_floating_base_, has_active_contacts_;
  double baumgarte_time_step_;
  int dimv_, dimu_, dim_passive_;
  static constexpr int kDimFloatingBase = 6;

  void setContactStatus(const ContactStatus& contact_status);

};

} // namespace idocp 

#include "idocp/ocp/contact_dynamics.hxx"

#endif // IDOCP_CONTACT_DYNAMICS_HPP_ 