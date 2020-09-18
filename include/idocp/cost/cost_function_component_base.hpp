#ifndef IDOCP_COST_FUNCTION_COMPONENT_BASE_HPP_
#define IDOCP_COST_FUNCTION_COMPONENT_BASE_HPP_

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/kkt_matrix.hpp"


namespace idocp {

///
/// @class CostFunctionComponentBase
/// @brief Base class of components of cost function.
///
class CostFunctionComponentBase {
public:

  ///
  /// @brief Default constructor. 
  ///
  CostFunctionComponentBase() {}

  ///
  /// @brief Destructor. 
  ///
  virtual ~CostFunctionComponentBase() {}

  ///
  /// @brief Default copy constructor. 
  ///
  CostFunctionComponentBase(const CostFunctionComponentBase&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  CostFunctionComponentBase& operator=(const CostFunctionComponentBase&) 
      = default;

  ///
  /// @brief Default move constructor. 
  ///
  CostFunctionComponentBase(CostFunctionComponentBase&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  CostFunctionComponentBase& operator=(CostFunctionComponentBase&&) noexcept 
      = default;

  ///
  /// @brief Return true if the cost function component requres kinematics of 
  /// robot model. Return false if not.
  /// @return true if the cost function component requres kinematics of 
  /// robot model. false if not.
  ///
  virtual bool useKinematics() const = 0;

  virtual double l(Robot& robot, CostFunctionData& data, const double t, 
                   const double dtau, const SplitSolution& s) const = 0;

  virtual double phi(Robot& robot, CostFunctionData& data, const double t, 
                     const SplitSolution& s) const = 0;

  virtual void lq(Robot& robot, CostFunctionData& data, const double t, 
                  const double dtau, const SplitSolution& s, 
                  KKTResidual& kkt_residual) const = 0;

  virtual void lv(Robot& robot, CostFunctionData& data, const double t, 
                  const double dtau, const SplitSolution& s, 
                  KKTResidual& kkt_residual) const = 0;

  virtual void la(Robot& robot, CostFunctionData& data, const double t, 
                  const double dtau, const SplitSolution& s, 
                  KKTResidual& kkt_residual) const = 0;

  virtual void lf(Robot& robot, CostFunctionData& data, const double t, 
                  const double dtau, const SplitSolution& s, 
                  KKTResidual& kkt_residual) const = 0;

  virtual void lqq(Robot& robot, CostFunctionData& data, const double t, 
                   const double dtau, const SplitSolution& s, 
                   KKTMatrix& kkt_matrix) const = 0;

  virtual void lvv(Robot& robot, CostFunctionData& data, const double t, 
                   const double dtau, const SplitSolution& s, 
                   KKTMatrix& kkt_matrix) const = 0;

  virtual void laa(Robot& robot, CostFunctionData& data, const double t, 
                   const double dtau, const SplitSolution& s, 
                   KKTMatrix& kkt_matrix) const = 0;

  virtual void lff(Robot& robot, CostFunctionData& data, const double t, 
                   const double dtau, const SplitSolution& s, 
                   KKTMatrix& kkt_matrix) const = 0;

  virtual void phiq(Robot& robot, CostFunctionData& data,  
                    const double t, const SplitSolution& s, 
                    KKTResidual& kkt_residual) const = 0;

  virtual void phiv(Robot& robot, CostFunctionData& data,  
                    const double t, const SplitSolution& s, 
                    KKTResidual& kkt_residual) const = 0;

  virtual void phiqq(Robot& robot, CostFunctionData& data,  
                     const double t, const SplitSolution& s,
                     KKTMatrix& kkt_matrix) const = 0;

  virtual void phivv(Robot& robot, CostFunctionData& data,  
                     const double t, const SplitSolution& s,
                     KKTMatrix& kkt_matrix) const = 0;

  virtual void lu(Robot& robot, CostFunctionData& data, const double t, 
                  const double dtau, const Eigen::VectorXd& u, 
                  Eigen::VectorXd& lu) const = 0;

  virtual void luu(Robot& robot, CostFunctionData& data, const double t, 
                   const double dtau, const Eigen::VectorXd& u, 
                   Eigen::MatrixXd& Quu) const  = 0;

};

} // namespace idocp

#endif // IDOCP_COST_FUNCTION_COMPONENT_BASE_HPP_