#ifndef IDOCP_SPLIT_OCP_HPP_
#define IDOCP_SPLIT_OCP_HPP_

#include <memory>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/constraints/constraints_data.hpp"
#include "idocp/ocp/state_equation.hpp"
#include "idocp/ocp/contact_dynamics.hpp"


namespace idocp {

///
/// @class SplitOCP
/// @brief Split optimal control problem of a single stage. 
///
class SplitOCP {
public:
  ///
  /// @brief Construct a split optimal control problem.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] cost Shared ptr to the cost function.
  /// @param[in] constraints Shared ptr to the constraints.
  /// @param[in] dt Step size of the discretization of the horizon. Only used 
  /// to determine the weight parameters of the Baumgarte's stabilization method.
  ///
  SplitOCP(const Robot& robot, const std::shared_ptr<CostFunction>& cost,
           const std::shared_ptr<Constraints>& constraints, const double dt);

  ///
  /// @brief Default constructor.  
  ///
  SplitOCP();

  ///
  /// @brief Destructor. 
  ///
  ~SplitOCP();

  ///
  /// @brief Default copy constructor. 
  ///
  SplitOCP(const SplitOCP&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  SplitOCP& operator=(const SplitOCP&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  SplitOCP(SplitOCP&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  SplitOCP& operator=(SplitOCP&&) noexcept = default;

  ///
  /// @brief Check whether the solution is feasible under inequality constraints.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] s Split solution of this stage.
  ///
  bool isFeasible(Robot& robot, const SplitSolution& s);

  ///
  /// @brief Initialize the constraints, i.e., set slack and dual variables. 
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] time_step Time step of this stage.
  /// @param[in] s Split solution of this stage.
  ///
  void initConstraints(Robot& robot, const ContactStatus& contact_status, 
                       const int time_step, SplitSolution& s);

  ///
  /// @brief Linearize the OCP for Newton's method around the current solution, 
  /// i.e., computes the KKT residual and Hessian.
  /// @tparam SplitSolutionType Type of the split solution at the next stage.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] contact_status Contact status of robot at this stage. 
  /// @param[in] t Current time of this stage. 
  /// @param[in] dt Length of the discretization of the horizon.
  /// @param[in] q_prev Configuration of the previous stage.
  /// @param[in] s Split solution of this stage.
  /// @param[in] s_next Split solution of the next stage.
  /// @param[out] kkt_matrix KKT matrix of this stage.
  /// @param[out] kkt_residual KKT residual of this stage.
  ///
  template <typename SplitSolutionType>
  void linearizeOCP(Robot& robot, const ContactStatus& contact_status, 
                    const double t, const double dt, 
                    const Eigen::VectorXd& q_prev, const SplitSolution& s, 
                    const SplitSolutionType& s_next, SplitKKTMatrix& kkt_matrix,
                    SplitKKTResidual& kkt_residual);

  ///
  /// @brief Computes the Newton direction of the condensed primal variables  
  /// at this stage.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] dt Length of the discretization of the horizon.
  /// @param[in] s Split solution of this stage.
  /// @param[in, out] d Split direction of this stage.
  /// 
  void computeCondensedPrimalDirection(Robot& robot, const double dt, 
                                       const SplitSolution& s, 
                                       SplitDirection& d);

  ///
  /// @brief Computes the Newton direction of the condensed dual variables 
  /// at this stage.
  /// @tparam SplitDirectionType Type of the split direction at the next stage.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] dt Length of the discretization of the horizon.
  /// @param[in] kkt_matrix KKT matrix of this stage.
  /// @param[in] kkt_residual KKT residual of this stage.
  /// @param[in] d_next Split direction of the next stage.
  /// @param[in, out] d Split direction of this stage.
  /// 
  template <typename SplitDirectionType>
  void computeCondensedDualDirection(const Robot& robot, const double dt, 
                                     const SplitKKTMatrix& kkt_matrix, 
                                     SplitKKTResidual& kkt_residual,
                                     const SplitDirectionType& d_next,
                                     SplitDirection& d);

  ///
  /// @brief Returns maximum stap size of the primal variables that satisfies 
  /// the inequality constraints.
  /// @return Maximum stap size of the primal variables that satisfies 
  /// the inequality constraints.
  ///
  double maxPrimalStepSize();

  ///
  /// @brief Returns maximum stap size of the dual variables that satisfies 
  /// the inequality constraints.
  /// @return Maximum stap size of the dual variables that satisfies 
  /// the inequality constraints.
  ///
  double maxDualStepSize();

  ///
  /// @brief Updates primal variables of this stage.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] primal_step_size Primal step size of the OCP. 
  /// @param[in] d Split direction of this stage.
  /// @param[in, out] s Split solution of this stage.
  ///
  void updatePrimal(const Robot& robot, const double primal_step_size, 
                    const SplitDirection& d, SplitSolution& s);

  ///
  /// @brief Updates dual variables of the inequality constraints.
  /// @param[in] dual_step_size Dula step size of the OCP. 
  ///
  void updateDual(const double dual_step_size);

  ///
  /// @brief Computes the KKT residual of the OCP at this stage.
  /// @tparam SplitSolutionType Type of the split solution at the next stage.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] contact_status Contact status of robot at this stage. 
  /// @param[in] t Current time of this stage. 
  /// @param[in] dt Length of the discretization of the horizon.
  /// @param[in] s Split solution of this stage.
  /// @param[in] q_prev Configuration of the previous stage.
  /// @param[in] s Split solution of this stage.
  /// @param[in] s_next Split solution of the next stage.
  /// @param[in, out] kkt_matrix KKT matrix of this stage.
  /// @param[in, out] kkt_residual KKT residual of this stage.
  ///
  template <typename SplitSolutionType>
  void computeKKTResidual(Robot& robot, const ContactStatus& contact_status,
                          const double t, const double dt, 
                          const Eigen::VectorXd& q_prev, const SplitSolution& s, 
                          const SplitSolutionType& s_next,
                          SplitKKTMatrix& kkt_matrix, 
                          SplitKKTResidual& kkt_residual);

  ///
  /// @brief Returns the KKT residual of the OCP at this stage. Before calling 
  /// this function, SplitOCP::linearizeOCP() or SplitOCP::computeKKTResidual()
  /// must be called.
  /// @return The squared norm of the kKT residual.
  ///
  double squaredNormKKTResidual(const SplitKKTResidual& kkt_residual, 
                                const double dt) const;

  ///
  /// @brief Computes the stage cost of this stage for line search.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] t Current time of this stage. 
  /// @param[in] dt Length of the discretization of the horizon.
  /// @param[in] s Split solution of this stage.
  /// @param[in] primal_step_size Primal step size of the OCP. Default is 0.
  /// @return Stage cost of this stage.
  /// 
  double stageCost(Robot& robot, const double t, const double dt, 
                   const SplitSolution& s, const double primal_step_size=0);

  ///
  /// @brief Computes and returns the constraint violation of the OCP at this 
  /// stage for line search.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] contact_status Contact status of robot at this stage. 
  /// @param[in] t Current time of this stage. 
  /// @param[in] dt Length of the discretization of the horizon.
  /// @param[in] s Split solution of this stage.
  /// @param[in] q_next Configuration at the next stage.
  /// @param[in] v_next Generaized velocity at the next stage.
  /// @param[in, out] kkt_residual KKT residual of this stage.
  /// @return Constraint violation of this stage.
  ///
  double constraintViolation(Robot& robot, const ContactStatus& contact_status, 
                             const double t, const double dt, 
                             const SplitSolution& s, 
                             const Eigen::VectorXd& q_next,
                             const Eigen::VectorXd& v_next,
                             SplitKKTResidual& kkt_residual);

  ///
  /// @brief Computes and returns the constraint violation of the OCP at this 
  /// stage for line search.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] contact_status Contact status of robot at this stage. 
  /// @param[in] t Current time of this stage. 
  /// @param[in] dt Length of the discretization of the horizon.
  /// @param[in] s Split solution of this stage.
  /// @param[in] q_next Configuration at the next stage.
  /// @param[in] v_next Generaized velocity at the next stage.
  /// @param[in, out] kkt_residual KKT residual of this stage.
  /// @param[in] impulse_status Impulse status of the next impulse. 
  /// @param[in] dt_next Length of the discretization at the next stage. 
  /// @return Constraint violation of this stage.
  ///
  double constraintViolation(Robot& robot, const ContactStatus& contact_status, 
                             const double t, const double dt, 
                             const SplitSolution& s, 
                             const Eigen::VectorXd& q_next,
                             const Eigen::VectorXd& v_next,
                             SplitKKTResidual& kkt_residual,
                             const ImpulseStatus& impulse_status, 
                             const double dt_next);

private:
  std::shared_ptr<CostFunction> cost_;
  CostFunctionData cost_data_;
  std::shared_ptr<Constraints> constraints_;
  ConstraintsData constraints_data_;
  ContactDynamics contact_dynamics_;
  bool use_kinematics_, has_floating_base_;
  double stage_cost_, constraint_violation_;

};

} // namespace idocp

#include "idocp/ocp/split_ocp.hxx"

#endif // IDOCP_SPLIT_OCP_HPP_