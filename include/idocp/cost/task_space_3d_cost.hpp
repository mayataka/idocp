#ifndef IDOCP_TASK_SPACE_3D_COST_HPP_
#define IDOCP_TASK_SPACE_3D_COST_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function_component_base.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/kkt_matrix.hpp"


namespace idocp {

class TaskSpace3DCost final : public CostFunctionComponentBase {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  TaskSpace3DCost(const Robot& robot, const int frame_id);

  TaskSpace3DCost();

  ~TaskSpace3DCost();

  // Use defalut copy constructor.
  TaskSpace3DCost(const TaskSpace3DCost&) = default;

  // Use defalut copy operator.
  TaskSpace3DCost& operator=(const TaskSpace3DCost&) = default;

  // Use defalut move constructor.
  TaskSpace3DCost(TaskSpace3DCost&&) noexcept = default;

  // Use defalut copy operator.
  TaskSpace3DCost& operator=(TaskSpace3DCost&&) noexcept = default;

  void set_q_ref(const Eigen::Vector3d& q_ref);

  void set_q_weight(const Eigen::Vector3d& q_weight);

  void set_qf_weight(const Eigen::Vector3d& qf_weight);

  double l(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s) const override;

  double phi(const Robot& robot, CostFunctionData& data, const double t, 
             const SplitSolution& s) const override; 

  void lq(const Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const SplitSolution& s, 
          KKTResidual& kkt_residual) const override;

  void lv(const Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const SplitSolution& s, 
          KKTResidual& kkt_residual) const override {}

  void la(const Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const SplitSolution& s,
          KKTResidual& kkt_residual) const override {}

  void lf(const Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const SplitSolution& s, 
          KKTResidual& kkt_residual) const override {}

  void lqq(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           KKTMatrix& kkt_matrix) const override;

  void lvv(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           KKTMatrix& kkt_matrix) const override {}

  void laa(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           KKTMatrix& kkt_matrix) const override {}

  void lff(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           KKTMatrix& kkt_matrix) const override {}

  void phiq(const Robot& robot, CostFunctionData& data, const double t, 
            const SplitSolution& s, KKTResidual& kkt_residual) const override;

  void phiv(const Robot& robot, CostFunctionData& data, const double t, 
            const SplitSolution& s, KKTResidual& kkt_residual) const override {}

  void phiqq(const Robot& robot, CostFunctionData& data, const double t, 
             const SplitSolution& s, KKTMatrix& kkt_matrix) const override;

  void phivv(const Robot& robot, CostFunctionData& data, const double t, 
             const SplitSolution& s, KKTMatrix& kkt_matrix) const override {}

  void lu(const Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const Eigen::VectorXd& u, 
          Eigen::VectorXd& lu) const override {}

  void luu(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const Eigen::VectorXd& u, 
           Eigen::MatrixXd& Quu) const override {}

private:
  int frame_id_;
  Eigen::Vector3d q_ref_, q_weight_, qf_weight_;

};

} // namespace idocp


#endif // IDOCP_TASK_SPACE_3D_COST_HPP_