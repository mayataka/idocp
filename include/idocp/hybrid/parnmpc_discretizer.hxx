#ifndef IDOCP_PARNMPC_DISCRETIZER_HXX_
#define IDOCP_PARNMPC_DISCRETIZER_HXX_

#include "idocp/hybrid/parnmpc_discretizer.hpp"

#include <cmath>
#include <cassert>
#include <algorithm>

namespace idocp {

inline ParNMPCDiscretizer::ParNMPCDiscretizer(const double T, const int N, 
                                              const int max_events) 
  : T_(T),
    dt_ideal_(T/N), 
    max_dt_(dt_ideal_-min_dt_),
    N_(N),
    N_ideal_(N),
    N_impulse_(0),
    N_lift_(0),
    max_events_(max_events),
    contact_phase_index_from_time_stage_(N+1, 0), 
    impulse_index_before_time_stage_(N, -1), 
    lift_index_before_time_stage_(N, -1), 
    time_stage_after_impulse_(max_events+1, -1), 
    time_stage_after_lift_(max_events+1, -1),
    is_time_stage_after_impulse_(N, false),
    is_time_stage_after_lift_(N, false),
    t_(N, 0),
    t_impulse_(max_events+1, 0),
    t_lift_(max_events+1, 0),
    dt_(N, (double)(T/N)),
    dt_aux_(max_events+1, 0),
    dt_lift_(max_events+1, 0) {
}


inline ParNMPCDiscretizer::ParNMPCDiscretizer()
  : T_(0),
    dt_ideal_(0), 
    max_dt_(0),
    N_(0),
    N_ideal_(0),
    N_impulse_(0),
    N_lift_(0),
    max_events_(0),
    contact_phase_index_from_time_stage_(), 
    impulse_index_before_time_stage_(), 
    lift_index_before_time_stage_(), 
    time_stage_after_impulse_(), 
    time_stage_after_lift_(),
    is_time_stage_after_impulse_(),
    is_time_stage_after_lift_(),
    t_(),
    t_impulse_(),
    t_lift_(),
    dt_(),
    dt_aux_(),
    dt_lift_() {
}


inline ParNMPCDiscretizer::~ParNMPCDiscretizer() {
}


inline void ParNMPCDiscretizer::discretizeOCP(
    const ContactSequence& contact_sequence, const double t) {
  countDiscreteEvents(contact_sequence, t);
  countTimeSteps(t);
  countTimeStages();
  countContactPhase();
  assert(isWellDefined());
}


inline int ParNMPCDiscretizer::N() const {
  return N_;
}


inline int ParNMPCDiscretizer::N_impulse() const {
  return N_impulse_;
}


inline int ParNMPCDiscretizer::N_lift() const {
  return N_lift_;
}


inline int ParNMPCDiscretizer::N_all() const {
  return (N()+2*N_impulse()+N_lift());
}


inline int ParNMPCDiscretizer::N_ideal() const {
  return N_ideal_;
}


inline int ParNMPCDiscretizer::contactPhase(const int time_stage) const {
  assert(time_stage >= -1);
  assert(time_stage <= N());
  return contact_phase_index_from_time_stage_[time_stage+1];
}


inline int ParNMPCDiscretizer::contactPhaseBeforeImpulse(
    const int impulse_index) const {
  return contactPhase(timeStageBeforeImpulse(impulse_index));
}


inline int ParNMPCDiscretizer::contactPhaseBeforeLift(
    const int lift_index) const {
  return contactPhase(timeStageBeforeLift(lift_index));
}


inline int ParNMPCDiscretizer::impulseIndexBeforeTimeStage(
    const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage < N());
  return impulse_index_before_time_stage_[time_stage];
}


inline int ParNMPCDiscretizer::impulseIndexAfterTimeStage(
    const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage < N());
  return impulseIndexBeforeTimeStage(time_stage+1);
}


inline int ParNMPCDiscretizer::liftIndexBeforeTimeStage(
    const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage < N());
  return lift_index_before_time_stage_[time_stage];
}


inline int ParNMPCDiscretizer::liftIndexAfterTimeStage(
    const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage < N());
  return liftIndexBeforeTimeStage(time_stage+1);
}


inline int ParNMPCDiscretizer::timeStageBeforeImpulse(
    const int impulse_index) const {
  assert(impulse_index >= 0);
  assert(impulse_index < N_impulse());
  return timeStageAfterImpulse(impulse_index) - 1;
}


inline int ParNMPCDiscretizer::timeStageAfterImpulse(
    const int impulse_index) const {
  return time_stage_after_impulse_[impulse_index];
}


inline int ParNMPCDiscretizer::timeStageBeforeLift(const int lift_index) const {
  assert(lift_index >= 0);
  assert(lift_index < N_lift());
  return timeStageAfterLift(lift_index) - 1;
}


inline int ParNMPCDiscretizer::timeStageAfterLift(const int lift_index) const {
  return time_stage_after_lift_[lift_index];
}


inline bool ParNMPCDiscretizer::isTimeStageBeforeImpulse(
    const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage < N());
  return isTimeStageAfterImpulse(time_stage+1);
}


inline bool ParNMPCDiscretizer::isTimeStageAfterImpulse(
    const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage < N());
  return is_time_stage_after_impulse_[time_stage];
}


inline bool ParNMPCDiscretizer::isTimeStageBeforeLift(
    const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage < N());
  return isTimeStageAfterLift(time_stage+1);
}


inline bool ParNMPCDiscretizer::isTimeStageAfterLift(
    const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage < N());
  return is_time_stage_after_lift_[time_stage];
}


inline double ParNMPCDiscretizer::t(const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage < N());
  return t_[time_stage];
}


inline double ParNMPCDiscretizer::t_impulse(const int impulse_index) const {
  assert(impulse_index >= 0);
  assert(impulse_index < N_impulse());
  return t_impulse_[impulse_index];
}


inline double ParNMPCDiscretizer::t_lift(const int lift_index) const {
  assert(lift_index >= 0);
  assert(lift_index < N_lift());
  return t_lift_[lift_index];
}


inline double ParNMPCDiscretizer::dt(const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage <= N());
  return dt_[time_stage];
}


inline double ParNMPCDiscretizer::dt_aux(const int impulse_index) const {
  assert(impulse_index >= 0);
  assert(impulse_index < N_impulse());
  return dt_aux_[impulse_index];
}


inline double ParNMPCDiscretizer::dt_lift(const int lift_index) const {
  assert(lift_index >= 0);
  assert(lift_index < N_lift());
  return dt_lift_[lift_index];
}


inline bool ParNMPCDiscretizer::isWellDefined() const {
  for (int i=0; i<N(); ++i) {
    if (isTimeStageAfterImpulse(i) && isTimeStageAfterLift(i)) {
      return false;
    }
  }
  for (int i=0; i<N()-1; ++i) {
    if (isTimeStageBeforeImpulse(i) && isTimeStageBeforeImpulse(i+1)) {
      return false;
    }
  }
  return true;
}


inline void ParNMPCDiscretizer::countDiscreteEvents(
    const ContactSequence& contact_sequence, const double t) {
  N_impulse_ = contact_sequence.numImpulseEvents();
  assert(N_impulse_ <= max_events_);
  for (int i=0; i<N_impulse_; ++i) {
    t_impulse_[i] = contact_sequence.impulseTime(i);
    time_stage_after_impulse_[i] = std::floor((t_impulse_[i]-t)/dt_ideal_);
  }
  N_lift_ = contact_sequence.numLiftEvents();
  assert(N_lift_ <= max_events_);
  for (int i=0; i<N_lift_; ++i) {
    t_lift_[i] = contact_sequence.liftTime(i);
    time_stage_after_lift_[i] = std::floor((t_lift_[i]-t)/dt_ideal_);
  }
  N_ = N_ideal_;
  assert(isWellDefined());
}


inline void ParNMPCDiscretizer::countTimeSteps(const double t) {
  int impulse_index = 0;
  int lift_index = 0;
  int num_events_on_grid = 0;
  for (int i=0; i<N_ideal_; ++i) {
    const int stage = i - num_events_on_grid;
    if (i == time_stage_after_impulse_[impulse_index]) {
      dt_[stage] = (i+1) * dt_ideal_ + t - t_impulse_[impulse_index];
      assert(dt_[stage] >= -min_dt_);
      assert(dt_[stage] <= dt_ideal_+min_dt_);
      if (dt_[stage] <= min_dt_) {
        time_stage_after_impulse_[impulse_index] = i + 1;
        t_[stage] = t + (i+1) * dt_ideal_;
      }
      else if (dt_[stage] >= max_dt_) {
        time_stage_after_impulse_[impulse_index] = stage - 1;
        dt_aux_[impulse_index] = dt_ideal_;
        t_[stage] = t + i * dt_ideal_;
        ++num_events_on_grid;
        ++impulse_index;
      }
      else {
        time_stage_after_impulse_[impulse_index] = stage;
        dt_aux_[impulse_index] = dt_ideal_ - dt_[stage];
        t_[stage] = t + (i+1) * dt_ideal_;
        ++impulse_index;
      }
    }
    else if (i == time_stage_after_lift_[lift_index]) {
      dt_[stage] = (i+1) * dt_ideal_ + t - t_lift_[lift_index];
      assert(dt_[stage] >= -min_dt_);
      assert(dt_[stage] <= dt_ideal_+min_dt_);
      if (dt_[stage] <= min_dt_) {
        time_stage_after_lift_[lift_index] = i + 1;
        t_[stage] = t + (i+1) * dt_ideal_;
      }
      else if (dt_[stage] >= max_dt_) {
        time_stage_after_lift_[lift_index] = stage - 1;
        dt_lift_[lift_index] = dt_ideal_;
        t_[stage] = t + i * dt_ideal_;
        ++num_events_on_grid;
        ++lift_index;
      }
      else {
        time_stage_after_lift_[lift_index] = stage;
        dt_lift_[lift_index] = dt_ideal_ - dt_[stage];
        t_[stage] = t + (i+1) * dt_ideal_;
        ++lift_index;
      }
    }
    else {
      dt_[stage] = dt_ideal_;
      t_[stage] = t + (i+1) * dt_ideal_;
    }
  }
  N_ = N_ideal_ - num_events_on_grid;
  t_[N_-1] = t + T_;
}


inline void ParNMPCDiscretizer::countTimeStages() {
  int impulse_index = 0;
  int lift_index = 0;
  for (int i=0; i<N(); ++i) {
    if (impulse_index < N_impulse()) {
      if (i == timeStageAfterImpulse(impulse_index)) {
        is_time_stage_after_impulse_[i] = true;
        impulse_index_before_time_stage_[i] = impulse_index;
        ++impulse_index;
      }
      else {
        is_time_stage_after_impulse_[i] = false;
        impulse_index_before_time_stage_[i] = -1;
      }
    }
    else {
      is_time_stage_after_impulse_[i] = false;
      impulse_index_before_time_stage_[i] = -1;
    }
    if (lift_index < N_lift()) {
      if (i == timeStageAfterLift(lift_index)) {
        is_time_stage_after_lift_[i] = true;
        lift_index_before_time_stage_[i] = lift_index;
        ++lift_index;
      }
      else {
        is_time_stage_after_lift_[i] = false;
        lift_index_before_time_stage_[i] = -1;
      }
    }
    else {
      is_time_stage_after_lift_[i] = false;
      lift_index_before_time_stage_[i] = -1;
    }
  }
}


inline void ParNMPCDiscretizer::countContactPhase() {
  contact_phase_index_from_time_stage_[0] = 0;
  int num_events = 0;
  for (int i=0; i<N(); ++i) {
    if (isTimeStageAfterImpulse(i) || isTimeStageAfterLift(i)) {
      ++num_events; 
    }
    contact_phase_index_from_time_stage_[i+1] = num_events;
  }
}

} // namespace idocp

#endif // IDOCP_PARNMPC_DISCRETIZER_HXX_ 