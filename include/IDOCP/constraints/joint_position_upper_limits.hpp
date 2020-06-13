#ifndef IDOCP_JOINT_POSITION_UPPER_LIMITS_HPP_
#define IDOCP_JOINT_POSITION_UPPER_LIMITS_HPP_

#include "Eigen/Core"

#include "robot/robot.hpp"


namespace idocp {

class JointPositionUpperLimits {
public:
  JointPositionUpperLimits(const Robot& robot);

  // Use default copy constructor.
  JointPositionUpperLimits(const JointPositionUpperLimits&) = default;

  // Use default copy operator.
  JointPositionUpperLimits& operator=(const JointPositionUpperLimits&) = default;

  void setSlackAndDual(const Robot& robot, const double barrier, 
                       const Eigen::VectorXd& q);

  void linearizeConstraint(const Robot& robot, const double barrier, 
                           const Eigen::VectorXd& q);

  void updateSlackAndDual(const Robot& robot, const Eigen::VectorXd& dq); 

  void condenseSlackAndDual(const Robot& robot, Eigen::MatrixXd& Cqq, 
                            Eigen::VectorXd& Cq) const;

  void augmentDualResidual(const Robot& robot, Eigen::VectorXd& Cq);

  double squaredConstraintsResidualNrom(const Robot& robot, 
                                        const Eigen::VectorXd& q);

private:
  unsigned int dimq_, dimv_, dimc_;
  Eigen::VectorXd qmax_, slack_, dual_, residual_, FB_residual_, dslack_, 
                  ddual_; 
};

} // namespace idocp


#endif // IDOCP_JOINT_POSITION_UPPER_LIMITS_HPP_