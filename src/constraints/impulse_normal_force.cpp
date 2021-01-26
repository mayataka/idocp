#include "idocp/constraints/impulse_normal_force.hpp"

#include <iostream>
#include <cassert>


namespace idocp {

ImpulseNormalForce::ImpulseNormalForce(const Robot& robot, const double barrier,
                                       const double fraction_to_boundary_rate)
  : ImpulseConstraintComponentBase(barrier, fraction_to_boundary_rate),
    dimc_(robot.maxPointContacts()),
    fraction_to_boundary_rate_(fraction_to_boundary_rate) {
}


ImpulseNormalForce::ImpulseNormalForce()
  : ImpulseConstraintComponentBase(),
    dimc_(0),
    fraction_to_boundary_rate_(0) {
}


ImpulseNormalForce::~ImpulseNormalForce() {
}


KinematicsLevel ImpulseNormalForce::kinematicsLevel() const {
  return KinematicsLevel::AccelerationLevel;
}


bool ImpulseNormalForce::isFeasible(Robot& robot, ConstraintComponentData& data, 
                                    const ImpulseSplitSolution& s) const {
  for (int i=0; i<dimc_; ++i) {
    if (s.isImpulseActive(i)) {
      if (s.f[i].coeff(2) < 0) {
        return false;
      }
    }
  }
  return true;
}


void ImpulseNormalForce::setSlackAndDual(
    Robot& robot, ConstraintComponentData& data, 
    const ImpulseSplitSolution& s) const {
  for (int i=0; i<dimc_; ++i) {
    data.slack.coeffRef(i) = s.f[i].coeff(2);
  }
  setSlackAndDualPositive(data);
}


void ImpulseNormalForce::augmentDualResidual(
    Robot& robot, ConstraintComponentData& data, 
    const ImpulseSplitSolution& s, ImpulseSplitKKTResidual& kkt_residual) const {
  int dimf_stack = 0;
  for (int i=0; i<dimc_; ++i) {
    if (s.isImpulseActive(i)) {
      kkt_residual.lf().coeffRef(dimf_stack+2) -= data.dual.coeff(i);
      dimf_stack += 3;
    }
  }
}


void ImpulseNormalForce::condenseSlackAndDual(
    Robot& robot, ConstraintComponentData& data, 
    const ImpulseSplitSolution& s, ImpulseSplitKKTMatrix& kkt_matrix, 
    ImpulseSplitKKTResidual& kkt_residual) const {
  int dimf_stack = 0;
  for (int i=0; i<dimc_; ++i) {
    if (s.isImpulseActive(i)) {
      kkt_matrix.Qff().coeffRef(dimf_stack+2, dimf_stack+2) 
          += data.dual.coeff(i) / data.slack.coeff(i);
      data.residual.coeffRef(i) = - s.f[i].coeff(2) + data.slack.coeff(i);
      data.duality.coeffRef(i) = computeDuality(data.slack.coeff(i), 
                                                data.dual.coeff(i));
      kkt_residual.lf().coeffRef(dimf_stack+2) 
          -= (data.dual.coeff(i)*data.residual.coeff(i)-data.duality.coeff(i)) 
             / data.slack.coeff(i);
      dimf_stack += 3;
    }
  }
}


void ImpulseNormalForce::computeSlackAndDualDirection(
    Robot& robot, ConstraintComponentData& data,  
    const ImpulseSplitSolution& s, const ImpulseSplitDirection& d) const {
  int dimf_stack = 0;
  for (int i=0; i<dimc_; ++i) {
    if (s.isImpulseActive(i)) {
      data.dslack.coeffRef(i) 
          = d.df().coeff(dimf_stack+2) - data.residual.coeff(i);
      data.ddual.coeffRef(i) = computeDualDirection(data.slack.coeff(i), 
                                                    data.dual.coeff(i), 
                                                    data.dslack.coeff(i), 
                                                    data.duality.coeff(i));
      dimf_stack += 3;
    }
    else {
      data.residual.coeffRef(i) = 0;
      data.duality.coeffRef(i)  = 0;
      data.slack.coeffRef(i)    = 1.0;
      data.dslack.coeffRef(i)   = fraction_to_boundary_rate_;
      data.dual.coeffRef(i)     = 1.0;
      data.ddual.coeffRef(i)    = fraction_to_boundary_rate_;
    }
  }
}


void ImpulseNormalForce::computePrimalAndDualResidual(
    Robot& robot, ConstraintComponentData& data, 
    const ImpulseSplitSolution& s) const {
  for (int i=0; i<dimc_; ++i) {
    if (s.isImpulseActive(i)) {
      data.residual.coeffRef(i) = - s.f[i].coeff(2) + data.slack.coeff(i);
      data.duality.coeffRef(i)  = computeDuality(data.slack.coeff(i), 
                                                 data.dual.coeff(i));
    }
    else {
      data.residual.coeffRef(i) = 0;
      data.duality.coeffRef(i)  = 0;
    }
  }
}


int ImpulseNormalForce::dimc() const {
  return dimc_;
}

} // namespace idocp