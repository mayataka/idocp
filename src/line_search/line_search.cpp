#include "idocp/line_search/line_search.hpp"

#include <stdexcept>
#include <cassert>


namespace idocp {

LineSearch::LineSearch(const Robot& robot, const int N, 
                       const int max_num_impulse, const int nthreads, 
                       const double step_size_reduction_rate,
                       const double min_step_size) 
  : filter_(),
    max_num_impulse_(max_num_impulse), 
    nthreads_(nthreads),
    step_size_reduction_rate_(step_size_reduction_rate), 
    min_step_size_(min_step_size),
    costs_(Eigen::VectorXd::Zero(N+1)), 
    costs_impulse_(Eigen::VectorXd::Zero(max_num_impulse)), 
    costs_aux_(Eigen::VectorXd::Zero(max_num_impulse)), 
    costs_lift_(Eigen::VectorXd::Zero(max_num_impulse)), 
    violations_(Eigen::VectorXd::Zero(N)), 
    violations_impulse_(Eigen::VectorXd::Zero(max_num_impulse)), 
    violations_aux_(Eigen::VectorXd::Zero(max_num_impulse)), 
    violations_lift_(Eigen::VectorXd::Zero(max_num_impulse)),
    s_try_(robot, N, max_num_impulse), 
    kkt_residual_(robot, N, max_num_impulse) {
}


LineSearch::LineSearch() 
  : filter_(),
    max_num_impulse_(0), 
    nthreads_(0),
    step_size_reduction_rate_(0), 
    min_step_size_(0),
    costs_(), 
    costs_impulse_(), 
    costs_aux_(), 
    costs_lift_(), 
    violations_(), 
    violations_impulse_(), 
    violations_aux_(), 
    violations_lift_(),
    s_try_(), 
    kkt_residual_() {
}


LineSearch::~LineSearch() {
}


void LineSearch::clearFilter() {
  filter_.clear();
}


bool LineSearch::isFilterEmpty() const {
  return filter_.isEmpty();
}


void LineSearch::computeCostAndViolation(
    OCP& ocp, std::vector<Robot>& robots, 
    const ContactSequence& contact_sequence, const Eigen::VectorXd& q, 
    const Eigen::VectorXd& v, const Solution& s, 
    const double primal_step_size) {
  assert(robots.size() == nthreads_);
  assert(q.size() == robots[0].dimq());
  assert(v.size() == robots[0].dimv());
  const int N = ocp.discrete().N();
  const int N_impulse = ocp.discrete().N_impulse();
  const int N_lift = ocp.discrete().N_lift();
  const int N_all = N + 1 + 2*N_impulse + N_lift;
  clearCosts();
  clearViolations();
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<N_all; ++i) {
    if (i < N) {
      costs_.coeffRef(i) = ocp[i].stageCost(robots[omp_get_thread_num()], 
                                            ocp.discrete().t(i), 
                                            ocp.discrete().dt(i), s[i], 
                                            primal_step_size);
      if (ocp.discrete().isTimeStageBeforeImpulse(i)) {
        violations_.coeffRef(i) = ocp[i].constraintViolation(
            robots[omp_get_thread_num()], 
            contact_sequence.contactStatus(ocp.discrete().contactPhase(i)),
            ocp.discrete().t(i), ocp.discrete().dt(i), s[i], 
            s.impulse[ocp.discrete().impulseIndexAfterTimeStage(i)].q,
            s.impulse[ocp.discrete().impulseIndexAfterTimeStage(i)].v,
            kkt_residual_[i]);
      }
      else if (ocp.discrete().isTimeStageBeforeLift(i)) {
        violations_.coeffRef(i) = ocp[i].constraintViolation(
            robots[omp_get_thread_num()], 
            contact_sequence.contactStatus(ocp.discrete().contactPhase(i)),
            ocp.discrete().t(i), ocp.discrete().dt(i), s[i], 
            s.lift[ocp.discrete().liftIndexAfterTimeStage(i)].q,
            s.lift[ocp.discrete().liftIndexAfterTimeStage(i)].v,
            kkt_residual_[i]);
      }
      if (ocp.discrete().isTimeStageBeforeImpulse(i+1)) {
        const int impulse_index  
            = ocp.discrete().impulseIndexAfterTimeStage(i+1);
        violations_.coeffRef(i) = ocp[i].constraintViolation(
            robots[omp_get_thread_num()], 
            contact_sequence.contactStatus(ocp.discrete().contactPhase(i)),
            ocp.discrete().t(i), ocp.discrete().dt(i), s[i], s[i+1].q, s[i+1].v, 
            kkt_residual_[i], contact_sequence.impulseStatus(impulse_index), 
            ocp.discrete().dt(i+1));
      }
      else {
        violations_.coeffRef(i) = ocp[i].constraintViolation(
            robots[omp_get_thread_num()], 
            contact_sequence.contactStatus(ocp.discrete().contactPhase(i)),
            ocp.discrete().t(i), ocp.discrete().dt(i), s[i], s[i+1].q, s[i+1].v, 
            kkt_residual_[i]);
      }
    }
    else if (i == N) {
      costs_.coeffRef(i) = ocp.terminal.terminalCost(robots[omp_get_thread_num()], 
                                                     ocp.discrete().t(i), s[i]);
    }
    else if (i < N+1+N_impulse) {
      const int impulse_index = i - (N+1);
      const int time_stage_before_impulse 
          = ocp.discrete().timeStageBeforeImpulse(impulse_index);
      violations_impulse_.coeffRef(impulse_index) 
          = ocp.impulse[impulse_index].constraintViolation(
              robots[omp_get_thread_num()], 
              contact_sequence.impulseStatus(impulse_index), 
              ocp.discrete().t_impulse(impulse_index), s.impulse[impulse_index],
              s.aux[impulse_index].q, s.aux[impulse_index].v,
              kkt_residual_.impulse[impulse_index]);
      costs_impulse_.coeffRef(impulse_index) 
          = ocp.impulse[impulse_index].stageCost(
              robots[omp_get_thread_num()], 
              ocp.discrete().t_impulse(impulse_index), s.impulse[impulse_index],
              primal_step_size);
    }
    else if (i < N+1+2*N_impulse) {
      const int impulse_index  = i - (N+1+N_impulse);
      const int time_stage_after_impulse 
          = ocp.discrete().timeStageAfterImpulse(impulse_index);
      costs_aux_.coeffRef(impulse_index) = ocp.aux[impulse_index].stageCost(
              robots[omp_get_thread_num()], 
              ocp.discrete().t_impulse(impulse_index),
              ocp.discrete().dt_aux(impulse_index), s.aux[impulse_index],
              primal_step_size);
      violations_aux_.coeffRef(impulse_index) 
          = ocp.aux[impulse_index].constraintViolation(
              robots[omp_get_thread_num()], 
              contact_sequence.contactStatus(
                  ocp.discrete().contactPhaseAfterImpulse(impulse_index)), 
              ocp.discrete().t_impulse(impulse_index), 
              ocp.discrete().dt_aux(impulse_index), s.aux[impulse_index],
              s[time_stage_after_impulse].q, s[time_stage_after_impulse].v,
              kkt_residual_.aux[impulse_index]);
    }
    else {
      const int lift_index = i - (N+1+2*N_impulse);
      const int time_stage_after_lift
          = ocp.discrete().timeStageAfterLift(lift_index);
      costs_lift_.coeffRef(lift_index) = ocp.lift[lift_index].stageCost(
              robots[omp_get_thread_num()], ocp.discrete().t_lift(lift_index), 
              ocp.discrete().dt_lift(lift_index), s.lift[lift_index], 
              primal_step_size);
      if (ocp.discrete().isTimeStageBeforeImpulse(time_stage_after_lift)) {
        const int impulse_index
            = ocp.discrete().impulseIndexAfterTimeStage(time_stage_after_lift);
        violations_lift_.coeffRef(lift_index) 
            = ocp.lift[lift_index].constraintViolation(
                robots[omp_get_thread_num()], 
                contact_sequence.contactStatus(
                    ocp.discrete().contactPhaseAfterLift(lift_index)), 
                ocp.discrete().t_lift(lift_index), 
                ocp.discrete().dt_lift(lift_index), s.lift[lift_index],
                s[time_stage_after_lift].q, s[time_stage_after_lift].v,
                kkt_residual_.lift[lift_index], 
                contact_sequence.impulseStatus(impulse_index), 
                ocp.discrete().dt(time_stage_after_lift));
      }
      else {
        violations_lift_.coeffRef(lift_index) 
            = ocp.lift[lift_index].constraintViolation(
                robots[omp_get_thread_num()], 
                contact_sequence.contactStatus(
                    ocp.discrete().contactPhaseAfterLift(lift_index)), 
                ocp.discrete().t_lift(lift_index), 
                ocp.discrete().dt_lift(lift_index), s.lift[lift_index],
                s[time_stage_after_lift].q, s[time_stage_after_lift].v,
                kkt_residual_.lift[lift_index]);
      }
    }
  }
}


void LineSearch::computeCostAndViolation(
    ParNMPC& parnmpc, std::vector<Robot>& robots, 
    const ContactSequence& contact_sequence, const Eigen::VectorXd& q, 
    const Eigen::VectorXd& v, const Solution& s, 
    const double primal_step_size) {
  assert(robots.size() == nthreads_);
  assert(q.size() == robots[0].dimq());
  assert(v.size() == robots[0].dimv());
  const int N = parnmpc.discrete().N();
  const int N_impulse = parnmpc.discrete().N_impulse();
  const int N_lift = parnmpc.discrete().N_lift();
  const int N_all = N + 2*N_impulse + N_lift;
  clearCosts();
  clearViolations();
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<N_all; ++i) {
    if (i < N-1) {
      costs_.coeffRef(i) = parnmpc[i].stageCost(robots[omp_get_thread_num()], 
                                                parnmpc.discrete().t(i), 
                                                parnmpc.discrete().dt(i), 
                                                s[i], primal_step_size);
      violations_.coeffRef(i) = parnmpc[i].constraintViolation(
          robots[omp_get_thread_num()], 
          contact_sequence.contactStatus(parnmpc.discrete().contactPhase(i)),
          parnmpc.discrete().t(i), parnmpc.discrete().dt(i), 
          q_prev(parnmpc.discrete(), q, s, i), 
          v_prev(parnmpc.discrete(), v, s, i), s[i], kkt_residual_[i]);
    }
    else if (i == N-1) {
      costs_.coeffRef(i) = parnmpc.terminal.stageCost(
          robots[omp_get_thread_num()], parnmpc.discrete().t(i), 
          parnmpc.discrete().dt(i), s[i], primal_step_size);
      violations_.coeffRef(i) = parnmpc.terminal.constraintViolation(
          robots[omp_get_thread_num()], 
          contact_sequence.contactStatus(parnmpc.discrete().contactPhase(i)),
          parnmpc.discrete().t(i), parnmpc.discrete().dt(i), 
          q_prev(parnmpc.discrete(), q, s, i), 
          v_prev(parnmpc.discrete(), v, s, i), s[i], kkt_residual_[i]);
    }
    else if (i < N+N_impulse) {
      const int impulse_index = i - N;
      const int time_stage_after_impulse 
          = parnmpc.discrete().timeStageAfterImpulse(impulse_index);
      costs_impulse_.coeffRef(impulse_index) 
          = parnmpc.impulse[impulse_index].stageCost(
              robots[omp_get_thread_num()], 
              parnmpc.discrete().t_impulse(impulse_index), 
              s.impulse[impulse_index], primal_step_size);
      violations_impulse_.coeffRef(impulse_index) 
          = parnmpc.impulse[impulse_index].constraintViolation(
              robots[omp_get_thread_num()], 
              contact_sequence.impulseStatus(impulse_index), 
              parnmpc.discrete().t_impulse(impulse_index), 
              s.aux[impulse_index].q, s.aux[impulse_index].v, 
              s.impulse[impulse_index], kkt_residual_.impulse[impulse_index]);
    }
    else if (i < N+2*N_impulse) {
      const int impulse_index  = i - (N+N_impulse);
      const int time_stage_before_impulse 
          = parnmpc.discrete().timeStageBeforeImpulse(impulse_index);
      costs_aux_.coeffRef(impulse_index) 
          = parnmpc.aux[impulse_index].stageCost(
              robots[omp_get_thread_num()], 
              parnmpc.discrete().t_impulse(impulse_index), 
              parnmpc.discrete().dt_aux(impulse_index), 
              s.aux[impulse_index], primal_step_size);
      if (time_stage_before_impulse >= 0) {
        violations_aux_.coeffRef(impulse_index) 
            = parnmpc.aux[impulse_index].constraintViolation(
                robots[omp_get_thread_num()], 
                contact_sequence.contactStatus(
                    parnmpc.discrete().contactPhaseBeforeImpulse(impulse_index)), 
                parnmpc.discrete().t_impulse(impulse_index), 
                parnmpc.discrete().dt_aux(impulse_index), 
                s[time_stage_before_impulse].q, s[time_stage_before_impulse].v, 
                s.aux[impulse_index], kkt_residual_.aux[impulse_index],
                contact_sequence.impulseStatus(impulse_index));
      }
      else {
        assert(time_stage_before_impulse == -1);
        violations_aux_.coeffRef(impulse_index) 
            = parnmpc.aux[impulse_index].constraintViolation(
                robots[omp_get_thread_num()], 
                contact_sequence.contactStatus(
                    parnmpc.discrete().contactPhaseBeforeImpulse(impulse_index)), 
                parnmpc.discrete().t_impulse(impulse_index), 
                parnmpc.discrete().dt_aux(impulse_index), 
                q, v, s.aux[impulse_index], kkt_residual_.aux[impulse_index]);
      }
    }
    else {
      const int lift_index = i - (N+2*N_impulse);
      const int time_stage_before_lift
          = parnmpc.discrete().timeStageBeforeLift(lift_index);
      costs_lift_.coeffRef(lift_index) = parnmpc.lift[lift_index].stageCost(
          robots[omp_get_thread_num()], parnmpc.discrete().t_lift(lift_index), 
          parnmpc.discrete().dt_lift(lift_index), s.lift[lift_index], 
          primal_step_size);
      if (time_stage_before_lift >= 0) {
        violations_lift_.coeffRef(lift_index) 
            = parnmpc.lift[lift_index].constraintViolation(
                robots[omp_get_thread_num()], 
                contact_sequence.contactStatus(
                    parnmpc.discrete().contactPhaseBeforeLift(lift_index)), 
                parnmpc.discrete().t_lift(lift_index), 
                parnmpc.discrete().dt_lift(lift_index), 
                s[time_stage_before_lift].q, s[time_stage_before_lift].v, 
                s.lift[lift_index], kkt_residual_.lift[lift_index]);
      }
      else {
        assert(time_stage_before_lift == -1);
        violations_lift_.coeffRef(lift_index) 
            = parnmpc.lift[lift_index].constraintViolation(
                robots[omp_get_thread_num()], 
                contact_sequence.contactStatus(
                    parnmpc.discrete().contactPhaseBeforeLift(lift_index)), 
                parnmpc.discrete().t_lift(lift_index), 
                parnmpc.discrete().dt_lift(lift_index), q, v, 
                s.lift[lift_index], kkt_residual_.lift[lift_index]);
      }
    }
  }
}


void LineSearch::computeSolution(const OCP& ocp, 
                                 const std::vector<Robot>& robots, 
                                 const Solution& s, const Direction& d, 
                                 const double step_size) {
  assert(robots.size() == nthreads_);
  assert(step_size > 0);
  assert(step_size <= 1);
  const int N = ocp.discrete().N();
  const int N_impulse = ocp.discrete().N_impulse();
  const int N_lift = ocp.discrete().N_lift();
  const int N_all = N + 1 + 2*N_impulse + N_lift;
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<N_all; ++i) {
    if (i <= N) {
      computeSolution(robots[omp_get_thread_num()], s[i], d[i], step_size, 
                      s_try_[i]);
    }
    else if (i < N+1+N_impulse) {
      const int impulse_index = i - (N+1);
      computeSolution(robots[omp_get_thread_num()], s.impulse[impulse_index], 
                      d.impulse[impulse_index], step_size, 
                      s_try_.impulse[impulse_index]);
    }
    else if (i < N+1+2*N_impulse) {
      const int impulse_index  = i - (N+1+N_impulse);
      computeSolution(robots[omp_get_thread_num()], s.aux[impulse_index], 
                      d.aux[impulse_index], step_size, 
                      s_try_.aux[impulse_index]);
    }
    else {
      const int lift_index = i - (N+1+2*N_impulse);
      computeSolution(robots[omp_get_thread_num()], s.lift[lift_index], 
                      d.lift[lift_index], step_size, s_try_.lift[lift_index]);
    }
  }
}


void LineSearch::computeSolution(const ParNMPC& parnmpc, 
                                 const std::vector<Robot>& robots, 
                                 const Solution& s, const Direction& d, 
                                 const double step_size) {
  assert(robots.size() == nthreads_);
  assert(step_size > 0);
  assert(step_size <= 1);
  const int N = parnmpc.discrete().N();
  const int N_impulse = parnmpc.discrete().N_impulse();
  const int N_lift = parnmpc.discrete().N_lift();
  const int N_all = N + 2*N_impulse + N_lift;
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<N_all; ++i) {
    if (i < N) {
      computeSolution(robots[omp_get_thread_num()], s[i], d[i], step_size, 
                      s_try_[i]);
    }
    else if (i < N+N_impulse) {
      const int impulse_index = i - N;
      computeSolution(robots[omp_get_thread_num()], s.impulse[impulse_index], 
                      d.impulse[impulse_index], step_size, 
                      s_try_.impulse[impulse_index]);
    }
    else if (i < N+2*N_impulse) {
      const int impulse_index  = i - (N+N_impulse);
      computeSolution(robots[omp_get_thread_num()], s.aux[impulse_index], 
                      d.aux[impulse_index], step_size, 
                      s_try_.aux[impulse_index]);
    }
    else {
      const int lift_index = i - (N+2*N_impulse);
      computeSolution(robots[omp_get_thread_num()], s.lift[lift_index], 
                      d.lift[lift_index], step_size, s_try_.lift[lift_index]);
    }
  }
}

} // namespace idocp