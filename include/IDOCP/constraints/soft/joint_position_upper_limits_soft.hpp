#ifndef IDOCP_JOINT_POSITION_UPPER_LIMITS_SOFT_HPP_
#define IDOCP_JOINT_POSITION_UPPER_LIMITS_SOFT_HPP_

#include "Eigen/Core"

#include "robot/robot.hpp"


namespace idocp {
namespace softconstraints {

class JointPositionUpperLimits {
public:
  JointPositionUpperLimits(const Robot& robot, const double barrier);

  // Use default copy constructor.
  JointPositionUpperLimits(const JointPositionUpperLimits&) = default;

  // Use default copy operator.
  JointPositionUpperLimits& operator=(const JointPositionUpperLimits&) = default;

  void setSlackAndDual(const Robot& robot, const double dtau, 
                       const Eigen::VectorXd& q);

  void linearizeConstraint(const Robot& robot, const double dtau, 
                           const Eigen::VectorXd& q);

  void updateSlackAndDual(const Robot& robot, const double dtau, 
                          const Eigen::VectorXd& dq); 

  void condenseSlackAndDual(const Robot& robot, const double dtau, 
                            Eigen::MatrixXd& Cqq, Eigen::VectorXd& Cq) const;

  void augmentDualResidual(const Robot& robot, const double dtau, 
                           Eigen::VectorXd& Cq);

  double squaredConstraintsResidualNrom(const Robot& robot, const double dtau, 
                                        const Eigen::VectorXd& q);

private:
  unsigned int dimq_, dimv_, dimc_;
  double barrier_;
  Eigen::VectorXd qmax_, slack_, dual_, residual_, FB_residual_, dslack_, 
                  ddual_; 
};

} // namespace softconstraints
} // namespace idocp


#endif // IDOCP_JOINT_POSITION_UPPER_LIMITS_SOFT_HPP_