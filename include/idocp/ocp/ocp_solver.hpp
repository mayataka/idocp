#ifndef IDOCP_OCP_SOLVER_HPP_
#define IDOCP_OCP_SOLVER_HPP_

#include <vector>
#include <memory>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/hybrid/contact_sequence.hpp"
#include "idocp/hybrid/hybrid_container.hpp"
#include "idocp/ocp/state_constraint_jacobian.hpp"
#include "idocp/ocp/ocp_linearizer.hpp"
#include "idocp/ocp/riccati_recursion_solver.hpp"
#include "idocp/line_search/line_search.hpp"

namespace idocp {

///
/// @class OCPSolver
/// @brief Optimal control problem solver by Riccati recursion. 
///
class OCPSolver {
public:
  ///
  /// @brief Construct optimal control problem solver.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] cost Shared ptr to the cost function.
  /// @param[in] constraints Shared ptr to the constraints.
  /// @param[in] T Length of the horizon. Must be positive.
  /// @param[in] N Number of discretization of the horizon. Must be more than 1. 
  /// @param[in] nthreads Number of the threads in solving the optimal control 
  /// problem. Must be positive. Default is 1.
  ///
  OCPSolver(const Robot& robot, const std::shared_ptr<CostFunction>& cost,
            const std::shared_ptr<Constraints>& constraints, const double T, 
            const int N, const int max_num_impulse=0, const int nthreads=1);

  ///
  /// @brief Default constructor. 
  ///
  OCPSolver();

  ///
  /// @brief Destructor. 
  ///
  ~OCPSolver();

  ///
  /// @brief Default copy constructor. 
  ///
  OCPSolver(const OCPSolver&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  OCPSolver& operator=(const OCPSolver&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  OCPSolver(OCPSolver&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  OCPSolver& operator=(OCPSolver&&) noexcept = default;

  ///
  /// @brief Initializes the inequality constraints, i.e., set slack variables 
  /// and the Lagrange multipliers of inequality constraints. Based on the 
  /// current solution.
  ///
  void initConstraints(const double t);

  ///
  /// @brief Updates the solution by computing the primal-dual Newon direction.
  /// @param[in] t Initial time of the horizon. Current time in MPC. 
  /// @param[in] q Initial configuration. Size must be Robot::dimq().
  /// @param[in] v Initial velocity. Size must be Robot::dimv().
  /// @param[in] line_search If true, filter line search is enabled. If false
  /// filter line search is disabled. Default is false.
  ///
  void updateSolution(const double t, const Eigen::VectorXd& q, 
                      const Eigen::VectorXd& v, const bool line_search=false);

  ///
  /// @brief Get the const reference to the split solution of a time stage. 
  /// For example, you can get the const reference to the control input torques 
  /// at the initial stage via ocp.getSolution(0).u.
  /// @param[in] stage Time stage of interest. Must be more than 0 and less 
  /// than N.
  /// @return Const reference to the split solution of the specified time stage.
  ///
  const SplitSolution& getSolution(const int stage) const;

  ///
  /// @brief Gets the state-feedback gain for the control input torques.
  /// @param[in] stage Time stage of interest. Must be more than 0 and less 
  /// than N-1.
  /// @param[out] Kq Gain with respec to the configuration. Size must be 
  /// Robot::dimv() x Robot::dimv().
  /// @param[out] Kv Gain with respec to the velocity. Size must be
  /// Robot::dimv() x Robot::dimv().
  ///
  void getStateFeedbackGain(const int stage, Eigen::MatrixXd& Kq, 
                            Eigen::MatrixXd& Kv) const;

  void setSolution(const std::string& name, const Eigen::VectorXd& value);

  void setContactStatusUniformly(const ContactStatus& contact_status);

  void pushBackContactStatus(const ContactStatus& contact_status, 
                             const double switching_time);

  void setContactPoints(const int contact_phase, 
                        const std::vector<Eigen::Vector3d>& contact_points);

  ///
  /// @brief Pop back a contact status. 
  ///
  void popBackContactStatus();

  ///
  /// @brief Pop front a contact status. 
  ///
  void popFrontContactStatus();

  ///
  /// @brief Clear the line search filter. 
  ///
  void clearLineSearchFilter();

  ///
  /// @brief Computes the KKT residula of the optimal control problem. 
  /// @param[in] t Current time. 
  /// @param[in] q Initial configuration. Size must be Robot::dimq().
  /// @param[in] v Initial velocity. Size must be Robot::dimv().
  ///
  void computeKKTResidual(const double t, const Eigen::VectorXd& q, 
                          const Eigen::VectorXd& v);

  ///
  /// @brief Returns the squared KKT error norm by using previously computed 
  /// results computed by updateSolution(). The result is not exactly the 
  /// same as the squared KKT error norm of the original optimal control 
  /// problem. The result is the squared norm of the condensed residual. 
  /// However, this variables is sufficiently close to the original KKT error norm.
  /// @return The squared norm of the condensed KKT residual.
  ///
  double KKTError();

  ///
  /// @brief Return true if the current solution is feasible under the 
  /// inequality constraints. Return false if it is not feasible.
  /// @return true if the current solution is feasible under the inequality 
  /// constraints. false if it is not feasible.
  ///
  bool isCurrentSolutionFeasible();

  ///
  /// @brief Get the solution vector. This function is not suitable for 
  /// real-time application, e.g., MPC, since this function reconstructs the 
  /// solution vector object.
  /// @param[in] name Name of the printed variable. 
  ///
  std::vector<Eigen::VectorXd> getSolution(const std::string& name) const;

  ///
  /// @brief Shows the information of the optimal control problem.
  ///
  void showInfo() const;

  double cost();

private:
  std::vector<Robot> robots_;
  ContactSequence contact_sequence_;
  OCPLinearizer ocp_linearizer_;
  RiccatiRecursionSolver riccati_solver_;
  LineSearch line_search_;
  OCP ocp_;
  KKTMatrix kkt_matrix_;
  KKTResidual kkt_residual_;
  StateConstraintJacobian jac_;
  RiccatiFactorization riccati_factorization_;
  Solution s_;
  Direction d_;

  void discretizeSolution();

};

} // namespace idocp 


#endif // IDOCP_OCP_SOLVER_HPP_ 