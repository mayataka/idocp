#ifndef IDOCP_OCP_LINEARIZER_HXX_ 
#define IDOCP_OCP_LINEARIZER_HXX_

#include "idocp/ocp/ocp_linearizer.hpp"

#include <omp.h>
#include <stdexcept>
#include <cassert>


namespace idocp {
namespace internal {

struct LinearizeOCP {
  template <typename SplitSolutionType>
  static inline void run(SplitOCP& split_ocp, Robot& robot, 
                         const ContactStatus& contact_status, const double t, 
                         const double dt, const Eigen::VectorXd& q_prev, 
                         SplitSolution& s, 
                         const SplitSolutionType& s_next, 
                         SplitKKTMatrix& kkt_matrix, 
                         SplitKKTResidual& kkt_residual) {
    split_ocp.linearizeOCP(robot, contact_status, t, dt, q_prev, s, s_next,
                           kkt_matrix, kkt_residual);
  }

  static inline void run(TerminalOCP& terminal_ocp, Robot& robot, 
                         const double t, const Eigen::VectorXd& q_prev, 
                         SplitSolution& s, SplitKKTMatrix& kkt_matrix, 
                         SplitKKTResidual& kkt_residual) {
    terminal_ocp.linearizeOCP(robot, t, q_prev, s, kkt_matrix, kkt_residual);
  }

  static inline void run(ImpulseSplitOCP& impulse_split_ocp, Robot& robot, 
                         const ImpulseStatus& impulse_status, const double t, 
                         const Eigen::VectorXd& q_prev, 
                         ImpulseSplitSolution& s, 
                         const SplitSolution& s_next, 
                         ImpulseSplitKKTMatrix& kkt_matrix, 
                         ImpulseSplitKKTResidual& kkt_residual) {
    impulse_split_ocp.linearizeOCP(robot, impulse_status, t, q_prev, s, s_next,
                                   kkt_matrix, kkt_residual);
  }
};


struct ComputeKKTResidual {
  template <typename SplitSolutionType>
  static inline void run(SplitOCP& split_ocp, Robot& robot, 
                         const ContactStatus& contact_status, const double t, 
                         const double dt, const Eigen::VectorXd& q_prev, 
                         SplitSolution& s, 
                         const SplitSolutionType& s_next, 
                         SplitKKTMatrix& kkt_matrix, 
                         SplitKKTResidual& kkt_residual) {
    split_ocp.computeKKTResidual(robot, contact_status, t, dt, q_prev, s, 
                                 s_next, kkt_matrix, kkt_residual);
  }

  static inline void run(TerminalOCP& terminal_ocp, Robot& robot,  
                         const double t, const Eigen::VectorXd& q_prev,
                         SplitSolution& s, SplitKKTMatrix& kkt_matrix, 
                         SplitKKTResidual& kkt_residual) {
    terminal_ocp.computeKKTResidual(robot, t, q_prev, s, kkt_matrix, kkt_residual);
  }

  static inline void run(ImpulseSplitOCP& impulse_split_ocp, Robot& robot, 
                         const ImpulseStatus& impulse_status, const double t, 
                         const Eigen::VectorXd& q_prev, 
                         ImpulseSplitSolution& s, 
                         const SplitSolution& s_next, 
                         ImpulseSplitKKTMatrix& kkt_matrix, 
                         ImpulseSplitKKTResidual& kkt_residual) {
    impulse_split_ocp.computeKKTResidual(robot, impulse_status, t, q_prev, s, 
                                         s_next, kkt_matrix, kkt_residual);
  }
};

} // namespace internal
} // namespace idocp


namespace idocp {

template <typename Algorithm>
inline void OCPLinearizer::runParallel(OCP& ocp, std::vector<Robot>& robots,
                                       const ContactSequence& contact_sequence,
                                       const Eigen::VectorXd& q, 
                                       const Eigen::VectorXd& v, 
                                       Solution& s, KKTMatrix& kkt_matrix,
                                       KKTResidual& kkt_residual,
                                       StateConstraintJacobian& jac) const {
  assert(robots.size() == nthreads_);
  assert(q.size() == robots[0].dimq());
  assert(v.size() == robots[0].dimv());
  const int N = ocp.discrete().N();
  const int N_impulse = ocp.discrete().N_impulse();
  const int N_lift = ocp.discrete().N_lift();
  const int N_all = N + 1 + 2*N_impulse + N_lift;
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<N_all; ++i) {
    if (i < N) {
      if (ocp.discrete().isTimeStageBeforeImpulse(i)) {
        Algorithm::run(
            ocp[i], robots[omp_get_thread_num()], 
            contact_sequence.contactStatus(ocp.discrete().contactPhase(i)), 
            ocp.discrete().t(i), ocp.discrete().dt(i), 
            q_prev(ocp.discrete(), q, s, i), 
            s[i], s.impulse[ocp.discrete().impulseIndexAfterTimeStage(i)], 
            kkt_matrix[i], kkt_residual[i]);
      }
      else if (ocp.discrete().isTimeStageBeforeLift(i)) {
        Algorithm::run(
            ocp[i], robots[omp_get_thread_num()], 
            contact_sequence.contactStatus(ocp.discrete().contactPhase(i)), 
            ocp.discrete().t(i), ocp.discrete().dt(i), 
            q_prev(ocp.discrete(), q, s, i), s[i], 
            s.lift[ocp.discrete().liftIndexAfterTimeStage(i)], 
            kkt_matrix[i], kkt_residual[i]);
      }
      else {
        Algorithm::run(
            ocp[i], robots[omp_get_thread_num()], 
            contact_sequence.contactStatus(ocp.discrete().contactPhase(i)), 
            ocp.discrete().t(i), ocp.discrete().dt(i), 
            q_prev(ocp.discrete(), q, s, i), 
            s[i], s[i+1], kkt_matrix[i], kkt_residual[i]);
      }
    }
    else if (i == N) {
      Algorithm::run(ocp.terminal, robots[omp_get_thread_num()], 
                     ocp.discrete().t(N), q_prev(ocp.discrete(), q, s, N), s[N], 
                     kkt_matrix[N], kkt_residual[N]);
    }
    else if (i < N+1+N_impulse) {
      const int impulse_index  = i - (N+1);
      const int time_stage_before_impulse 
          = ocp.discrete().timeStageBeforeImpulse(impulse_index);
      Algorithm::run(ocp.impulse[impulse_index], robots[omp_get_thread_num()], 
                     contact_sequence.impulseStatus(impulse_index), 
                     ocp.discrete().t_impulse(impulse_index), 
                     s[time_stage_before_impulse].q, s.impulse[impulse_index], 
                     s.aux[impulse_index], kkt_matrix.impulse[impulse_index], 
                     kkt_residual.impulse[impulse_index]);
    }
    else if (i < N+1+2*N_impulse) {
      const int impulse_index  = i - (N+1+N_impulse);
      const int time_stage_after_impulse 
          = ocp.discrete().timeStageAfterImpulse(impulse_index);
      Algorithm::run(
          ocp.aux[impulse_index], robots[omp_get_thread_num()], 
          contact_sequence.contactStatus(
              ocp.discrete().contactPhaseAfterImpulse(impulse_index)), 
          ocp.discrete().t_impulse(impulse_index), 
          ocp.discrete().dt_aux(impulse_index), s.impulse[impulse_index].q, 
          s.aux[impulse_index], s[time_stage_after_impulse], 
          kkt_matrix.aux[impulse_index], kkt_residual.aux[impulse_index]);
    }
    else {
      const int lift_index = i - (N+1+2*N_impulse);
      const int time_stage_after_lift
          = ocp.discrete().timeStageAfterLift(lift_index);
      Algorithm::run(
          ocp.lift[lift_index], robots[omp_get_thread_num()], 
          contact_sequence.contactStatus(
              ocp.discrete().contactPhaseAfterLift(lift_index)), 
          ocp.discrete().t_lift(lift_index), 
          ocp.discrete().dt_lift(lift_index), s[time_stage_after_lift-1].q, 
          s.lift[lift_index], s[time_stage_after_lift], 
          kkt_matrix.lift[lift_index], kkt_residual.lift[lift_index]);
    }
  }
}


inline const Eigen::VectorXd& OCPLinearizer::q_prev(
    const OCPDiscretizer& discretizer, const Eigen::VectorXd& q, 
    const Solution& s, const int time_stage) {
  assert(time_stage >= 0);
  assert(time_stage <= discretizer.N());
  if (time_stage == 0) {
    return q;
  }
  else if (discretizer.isTimeStageBeforeImpulse(time_stage-1)) {
    return s.aux[discretizer.impulseIndexAfterTimeStage(time_stage-1)].q;
  }
  else if (discretizer.isTimeStageBeforeLift(time_stage-1)) {
    return s.lift[discretizer.liftIndexAfterTimeStage(time_stage-1)].q;
  }
  else {
    return s[time_stage-1].q;
  }
}

} // namespace idocp 

#endif // IDOCP_OCP_LINEARIZER_HXX_ 