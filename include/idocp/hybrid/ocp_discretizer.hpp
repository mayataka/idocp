#ifndef IDOCP_OCP_DISCRETIZER_HPP_
#define IDOCP_OCP_DISCRETIZER_HPP_

#include "idocp/hybrid/contact_sequence.hpp"

#include <vector>
#include <limits>
#include <cmath>

namespace idocp {

class OCPDiscretizer {
public:
  OCPDiscretizer(const double T, const int N, const int max_events);

  ///
  /// @brief Default constructor. 
  ///
  OCPDiscretizer();

  ///
  /// @brief Destructor. 
  ///
  ~OCPDiscretizer();

  ///
  /// @brief Default copy constructor. 
  ///
  OCPDiscretizer(const OCPDiscretizer&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  OCPDiscretizer& operator=(const OCPDiscretizer&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  OCPDiscretizer(OCPDiscretizer&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  OCPDiscretizer& operator=(OCPDiscretizer&&) noexcept = default;

  void discretizeOCP(const ContactSequence& contact_sequence, const double t);

  int N() const;

  int N_impulse() const;

  int N_lift() const; 

  int N_all() const;

  int N_ideal() const;

  int contactPhase(const int time_stage) const;

  int contactPhaseAfterImpulse(const int impulse_index) const;

  int contactPhaseAfterLift(const int lift_index) const;

  int impulseIndexAfterTimeStage(const int time_stage) const;

  int liftIndexAfterTimeStage(const int time_stage) const;

  int timeStageBeforeImpulse(const int impulse_index) const;

  int timeStageAfterImpulse(const int impulse_index) const;

  int timeStageBeforeLift(const int lift_index) const;

  int timeStageAfterLift(const int lift_index) const;

  bool isTimeStageBeforeImpulse(const int time_stage) const;

  bool isTimeStageAfterImpulse(const int time_stage) const;

  bool isTimeStageBeforeLift(const int time_stage) const;

  bool isTimeStageAfterLift(const int time_stage) const;

  double t(const int time_stage) const;

  double t_impulse(const int impulse_index) const;

  double t_lift(const int lift_index) const;

  double dt(const int time_stage) const;

  double dt_aux(const int impulse_index) const;

  double dt_lift(const int lift_index) const;

  bool isWellDefined() const;

  void showInfo() const;

  void showInfo(const ContactSequence& contact_sequence) const;

private:
  double T_, dt_ideal_, max_dt_;
  int N_, N_ideal_, N_impulse_, N_lift_, max_events_;
  std::vector<int> contact_phase_index_from_time_stage_, 
                   impulse_index_after_time_stage_, 
                   lift_index_after_time_stage_, time_stage_before_impulse_, 
                   time_stage_before_lift_;
  std::vector<bool> is_time_stage_before_impulse_, is_time_stage_before_lift_;
  std::vector<double> t_, t_impulse_, t_lift_, dt_, dt_aux_, dt_lift_;

  static constexpr double min_dt_ 
      = std::sqrt(std::numeric_limits<double>::epsilon());

  void countDiscreteEvents(const ContactSequence& contact_sequence, 
                           const double t);

  void countTimeSteps(const double t);

  void countTimeStages();

  void countContactPhase();

};

} // namespace idocp

#include "idocp/hybrid/ocp_discretizer.hxx"

#endif // IDOCP_OCP_DISCRETIZER_HPP_ 