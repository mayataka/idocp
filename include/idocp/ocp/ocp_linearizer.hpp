#ifndef IDOCP_OCP_LINEARIZER_HPP_ 
#define IDOCP_OCP_LINEARIZER_HPP_

#include <vector>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/hybrid/hybrid_container.hpp"
#include "idocp/ocp/state_constraint_jacobian.hpp"
#include "idocp/hybrid/contact_sequence.hpp"
#include "idocp/hybrid/ocp_discretizer.hpp"


namespace idocp {

///
/// @class OCPLinearizer
/// @brief Linearize of the optimal control problem. 
///
class OCPLinearizer {
public:
  ///
  /// @brief Construct optimal control problem solver.
  /// @param[in] N Number of discretization of the horizon. Must be more than 1. 
  /// @param[in] max_num_impulse Maximum number of the impulse on the horizon. 
  /// Must be non-negative. 
  /// @param[in] nthreads Number of the threads in solving the optimal control 
  /// problem. Must be positive. 
  ///
  OCPLinearizer(const int N, const int max_num_impulse, const int nthreads);

  ///
  /// @brief Default constructor. 
  ///
  OCPLinearizer();

  ///
  /// @brief Destructor. 
  ///
  ~OCPLinearizer();

  ///
  /// @brief Default copy constructor. 
  ///
  OCPLinearizer(const OCPLinearizer&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  OCPLinearizer& operator=(const OCPLinearizer&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  OCPLinearizer(OCPLinearizer&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  OCPLinearizer& operator=(OCPLinearizer&&) noexcept = default;

  void initConstraints(OCP& ocp, std::vector<Robot>& robots,
                       const ContactSequence& contact_sequence, Solution& s) const;

  void linearizeOCP(OCP& ocp, std::vector<Robot>& robots,
                    const ContactSequence& contact_sequence,
                    const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                    Solution& s, KKTMatrix& kkt_matrix, 
                    KKTResidual& kkt_residual,
                    StateConstraintJacobian& jac) const;

  void computeKKTResidual(OCP& ocp, std::vector<Robot>& robots, 
                          const ContactSequence& contact_sequence,
                          const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                          Solution& s, KKTMatrix& kkt_matrix, 
                          KKTResidual& kkt_residual,
                          StateConstraintJacobian& jac) const;

  double KKTError(const OCP& ocp, const KKTResidual& kkt_residual);

  void printKKTError() const;

  void integrateSolution(OCP& ocp, const std::vector<Robot>& robots,
                         const KKTMatrix& kkt_matrix, KKTResidual& kkt_residual,
                         const double primal_step_size,
                         const double dual_step_size,
                         Direction& d, Solution& s) const;

  static const Eigen::VectorXd& q_prev(const OCPDiscretizer& ocp_discretizer, 
                                       const Eigen::VectorXd& q, 
                                       const Solution& s, const int time_stage);


private:
  template <typename Algorithm>
  void runParallel(OCP& ocp, std::vector<Robot>& robots,
                   const ContactSequence& contact_sequence,
                   const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                   Solution& s, KKTMatrix& kkt_matrix, 
                   KKTResidual& kkt_residual,
                   StateConstraintJacobian& jac) const;

  int max_num_impulse_, nthreads_;
  Eigen::VectorXd kkt_error_;
};

} // namespace idocp 

#include "idocp/ocp/ocp_linearizer.hxx"

#endif // IDOCP_OCP_LINEARIZER_HPP_