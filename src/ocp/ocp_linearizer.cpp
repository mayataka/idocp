#include "idocp/ocp/ocp_linearizer.hpp"

#include <omp.h>
#include <stdexcept>
#include <cassert>

namespace idocp{

OCPLinearizer::OCPLinearizer(const int N, const int max_num_impulse, 
                             const int nthreads) 
  : max_num_impulse_(max_num_impulse),
    nthreads_(nthreads),
    kkt_error_(Eigen::VectorXd::Zero(N+1+3*max_num_impulse)) {
  try {
    if (max_num_impulse < 0) {
      throw std::out_of_range("invalid value: max_num_impulse must be non-negative!");
    }
    if (nthreads <= 0) {
      throw std::out_of_range("invalid value: nthreads must be positive!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
}


OCPLinearizer::OCPLinearizer()
  : max_num_impulse_(0),
    nthreads_(0),
    kkt_error_() {
}


OCPLinearizer::~OCPLinearizer() {
}


void OCPLinearizer::initConstraints(OCP& ocp, std::vector<Robot>& robots, 
                                    const ContactSequence& contact_sequence, 
                                    Solution& s) const {
  const int N = ocp.discrete().N();
  const int N_impulse = ocp.discrete().N_impulse();
  const int N_lift = ocp.discrete().N_lift();
  const int N_all = N + 1 + 2*N_impulse + N_lift;
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<N_all; ++i) {
    if (i < N) {
      ocp[i].initConstraints(
          robots[omp_get_thread_num()], 
          contact_sequence.contactStatus(ocp.discrete().contactPhase(i)), 
          i, s[i]);
    }
    else if (i == N) {
      ocp.terminal.initConstraints(robots[omp_get_thread_num()], N, s[N]);
    }
    else if (i < N+1+N_impulse) {
      const int impulse_index  = i - (N+1);
      ocp.impulse[impulse_index].initConstraints(
          robots[omp_get_thread_num()], 
          contact_sequence.impulseStatus(impulse_index), 
          s.impulse[impulse_index]);
    }
    else if (i < N+1+2*N_impulse) {
      const int impulse_index  = i - (N+1+N_impulse);
      const int time_stage_after_impulse 
          = ocp.discrete().timeStageAfterImpulse(impulse_index);
      ocp.aux[impulse_index].initConstraints(
          robots[omp_get_thread_num()],  
          contact_sequence.contactStatus(
              ocp.discrete().contactPhaseAfterImpulse(impulse_index)), 
          0, s.aux[impulse_index]);
    }
    else {
      const int lift_index = i - (N+1+2*N_impulse);
      ocp.lift[lift_index].initConstraints(
          robots[omp_get_thread_num()], 
          contact_sequence.contactStatus(
              ocp.discrete().contactPhaseAfterLift(lift_index)), 
          0, s.lift[lift_index]);
    }
  }
}


void OCPLinearizer::linearizeOCP(OCP& ocp, std::vector<Robot>& robots, 
                                 const ContactSequence& contact_sequence, 
                                 const Eigen::VectorXd& q, 
                                 const Eigen::VectorXd& v, Solution& s, 
                                 KKTMatrix& kkt_matrix, 
                                 KKTResidual& kkt_residual,
                                 StateConstraintJacobian& jac) const {
  runParallel<internal::LinearizeOCP>(ocp, robots, contact_sequence, q, v, 
                                      s, kkt_matrix, kkt_residual, jac);
}


void OCPLinearizer::computeKKTResidual(OCP& ocp, std::vector<Robot>& robots, 
                                       const ContactSequence& contact_sequence, 
                                       const Eigen::VectorXd& q, 
                                       const Eigen::VectorXd& v, 
                                       Solution& s, KKTMatrix& kkt_matrix, 
                                       KKTResidual& kkt_residual,
                                       StateConstraintJacobian& jac) const {
  runParallel<internal::ComputeKKTResidual>(ocp, robots, contact_sequence, q, v,  
                                            s, kkt_matrix, kkt_residual, jac);
}


double OCPLinearizer::KKTError(const OCP& ocp, 
                               const KKTResidual& kkt_residual) {
  const int N = ocp.discrete().N();
  const int N_impulse = ocp.discrete().N_impulse();
  const int N_lift = ocp.discrete().N_lift();
  const int N_all = N + 1 + 2*N_impulse + N_lift;
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<N_all; ++i) {
    if (i < N) {
      kkt_error_.coeffRef(i) 
          = ocp[i].squaredNormKKTResidual(kkt_residual[i], ocp.discrete().dt(i));
    }
    else if (i == N) {
      kkt_error_.coeffRef(N) 
          = ocp.terminal.squaredNormKKTResidual(kkt_residual[N]);
    }
    else if (i < N+1+N_impulse) {
      const int impulse_index  = i - (N+1);
      kkt_error_.coeffRef(i) 
          = ocp.impulse[impulse_index].squaredNormKKTResidual(
              kkt_residual.impulse[impulse_index]);
    }
    else if (i < N+1+2*N_impulse) {
      const int impulse_index  = i - (N+1+N_impulse);
      kkt_error_.coeffRef(i) 
          = ocp.aux[impulse_index].squaredNormKKTResidual(
                kkt_residual.aux[impulse_index], 
                ocp.discrete().dt_aux(impulse_index));
    }
    else {
      const int lift_index = i - (N+1+2*N_impulse);
      kkt_error_.coeffRef(i) 
          = ocp.lift[lift_index].squaredNormKKTResidual(
              kkt_residual.lift[lift_index], ocp.discrete().dt_lift(lift_index));
    }
  }
  return std::sqrt(kkt_error_.head(N_all).sum());
}


void OCPLinearizer::integrateSolution(OCP& ocp, 
                                      const std::vector<Robot>& robots, 
                                      const KKTMatrix& kkt_matrix, 
                                      KKTResidual& kkt_residual, 
                                      const double primal_step_size, 
                                      const double dual_step_size, 
                                      Direction& d, Solution& s) const {
  assert(robots.size() == nthreads_);
  const int N = ocp.discrete().N();
  const int N_impulse = ocp.discrete().N_impulse();
  const int N_lift = ocp.discrete().N_lift();
  const int N_all = N + 1 + 2*N_impulse + N_lift;
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<N_all; ++i) {
    if (i < N) {
      if (ocp.discrete().isTimeStageBeforeImpulse(i)) {
        ocp[i].computeCondensedDualDirection(
            robots[omp_get_thread_num()], ocp.discrete().dt(i), 
            kkt_matrix[i], kkt_residual[i], 
            d.impulse[ocp.discrete().impulseIndexAfterTimeStage(i)], d[i]);
      }
      else if (ocp.discrete().isTimeStageBeforeLift(i)) {
        ocp[i].computeCondensedDualDirection(
            robots[omp_get_thread_num()], ocp.discrete().dt(i), 
            kkt_matrix[i], kkt_residual[i], 
            d.lift[ocp.discrete().liftIndexAfterTimeStage(i)], d[i]);
      }
      else {
        ocp[i].computeCondensedDualDirection(
            robots[omp_get_thread_num()], ocp.discrete().dt(i), 
            kkt_matrix[i], kkt_residual[i], d[i+1], d[i]);
      }
      ocp[i].updatePrimal(robots[omp_get_thread_num()], primal_step_size, 
                          d[i], s[i]);
      ocp[i].updateDual(dual_step_size);
    }
    else if (i == N) {
      ocp.terminal.computeCondensedDualDirection(robots[omp_get_thread_num()], 
                                                 kkt_matrix[N], 
                                                 kkt_residual[N], d[N]);
      ocp.terminal.updatePrimal(robots[omp_get_thread_num()], primal_step_size, 
                                d[N], s[N]);
      ocp.terminal.updateDual(dual_step_size);
    }
    else if (i < N+1+N_impulse) {
      const int impulse_index  = i - (N+1);
      ocp.impulse[impulse_index].computeCondensedDualDirection(
          robots[omp_get_thread_num()], kkt_matrix.impulse[impulse_index], 
          kkt_residual.impulse[impulse_index], d.aux[impulse_index], 
          d.impulse[impulse_index]);
      ocp.impulse[impulse_index].updatePrimal(robots[omp_get_thread_num()], 
                                              primal_step_size, 
                                              d.impulse[impulse_index], 
                                              s.impulse[impulse_index]);
      ocp.impulse[impulse_index].updateDual(dual_step_size);
    }
    else if (i < N+1+2*N_impulse) {
      const int impulse_index  = i - (N+1+N_impulse);
      ocp.aux[impulse_index].computeCondensedDualDirection(
          robots[omp_get_thread_num()], ocp.discrete().dt_aux(impulse_index),
          kkt_matrix.aux[impulse_index], kkt_residual.aux[impulse_index],
          d[ocp.discrete().timeStageAfterImpulse(impulse_index)], 
          d.aux[impulse_index]);
      ocp.aux[impulse_index].updatePrimal(robots[omp_get_thread_num()], 
                                          primal_step_size, 
                                          d.aux[impulse_index], 
                                          s.aux[impulse_index]);
      ocp.aux[impulse_index].updateDual(dual_step_size);
    }
    else {
      const int lift_index = i - (N+1+2*N_impulse);
      ocp.lift[lift_index].computeCondensedDualDirection(
          robots[omp_get_thread_num()], ocp.discrete().dt_lift(lift_index), 
          kkt_matrix.lift[lift_index], kkt_residual.lift[lift_index], 
          d[ocp.discrete().timeStageAfterLift(lift_index)], d.lift[lift_index]);
      ocp.lift[lift_index].updatePrimal(robots[omp_get_thread_num()], 
                                        primal_step_size, 
                                        d.lift[lift_index], s.lift[lift_index]);
      ocp.lift[lift_index].updateDual(dual_step_size);
    }
  }
}

} // namespace idocp