#include <iostream>
#include <string>
#include <memory>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/ocp_solver.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/cost/trotting_configuration_space_cost.hpp"
#include "idocp/cost/contact_force_cost.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/constraints/joint_position_lower_limit.hpp"
#include "idocp/constraints/joint_position_upper_limit.hpp"
#include "idocp/constraints/joint_velocity_lower_limit.hpp"
#include "idocp/constraints/joint_velocity_upper_limit.hpp"
#include "idocp/constraints/joint_torques_lower_limit.hpp"
#include "idocp/constraints/joint_torques_upper_limit.hpp"
#include "idocp/constraints/linearized_friction_cone.hpp"
#include "idocp/constraints/linearized_impulse_friction_cone.hpp"

#include "idocp/utils/ocp_benchmarker.hpp"

#ifdef ENABLE_VIEWER
#include "idocp/utils/trajectory_viewer.hpp"
#endif 


int main(int argc, char *argv[]) {
  std::vector<int> contact_frames = {14, 24, 34, 44}; // LF, LH, RF, RH
  const std::string path_to_urdf = "../anymal_b_simple_description/urdf/anymal.urdf";
  idocp::Robot robot(path_to_urdf, contact_frames);
  const double baum = 1.55 / 30;
  const double alpha = 2 / baum;
  const double beta = 1 / (baum*baum);
  robot.setBaumgarteWeights(alpha, beta);

  const double step_length = 0.15;
  const double t_start = 0.5;
  const double t_period = 0.5;

  auto cost = std::make_shared<idocp::CostFunction>();
  Eigen::VectorXd q_standing(Eigen::VectorXd::Zero(robot.dimq()));
  q_standing << 0, 0, 0.4792, 0, 0, 0, 1, 
                -0.1,  0.7, -1.0, 
                -0.1, -0.7,  1.0, 
                 0.1,  0.7, -1.0, 
                 0.1, -0.7,  1.0;
  Eigen::VectorXd q_weight(Eigen::VectorXd::Zero(robot.dimv()));
  q_weight << 10, 10, 10, 10, 10, 10, 
              10, 10, 10,
              10, 10, 10,
              10, 10, 10,
              10, 10, 10;
  Eigen::VectorXd v_weight(Eigen::VectorXd::Zero(robot.dimv()));
  v_weight << 1, 1, 1, 1, 1, 1, 
              0.1, 0.1, 0.1,
              0.1, 0.1, 0.1,
              0.1, 0.1, 0.1,
              0.1, 0.1, 0.1;
  Eigen::VectorXd a_weight(Eigen::VectorXd::Zero(robot.dimv()));
  a_weight << 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 
              0.01, 0.01, 0.01,
              0.01, 0.01, 0.01,
              0.01, 0.01, 0.01,
              0.01, 0.01, 0.01;

  idocp::TrottingSwingAngles swing_angles;
  swing_angles.front_swing_knee   = 1.7; 
  swing_angles.hip_swing_knee     = 1.7;
  auto config_cost = std::make_shared<idocp::TrottingConfigurationSpaceCost>(robot);
  config_cost->set_ref(t_start, t_period, q_standing, step_length, swing_angles);
  config_cost->set_q_weight(q_weight);
  config_cost->set_qf_weight(q_weight);
  config_cost->set_qi_weight(q_weight);
  config_cost->set_v_weight(v_weight);
  config_cost->set_vf_weight(v_weight);
  config_cost->set_vi_weight(v_weight);
  config_cost->set_a_weight(a_weight);
  config_cost->set_dvi_weight(a_weight);
  cost->push_back(config_cost);

  auto contact_cost = std::make_shared<idocp::ContactForceCost>(robot);
  std::vector<Eigen::Vector3d> f_weight;
  for (int i=0; i<contact_frames.size(); ++i) {
    Eigen::Vector3d fw; 
    fw << 0.001, 0.001, 0.001;
    f_weight.push_back(fw);
  }
  contact_cost->set_f_weight(f_weight);
  contact_cost->set_fi_weight(f_weight);
  contact_cost->set_f_ref(robot);
  cost->push_back(contact_cost);

  auto constraints           = std::make_shared<idocp::Constraints>();
  auto joint_position_lower  = std::make_shared<idocp::JointPositionLowerLimit>(robot);
  auto joint_position_upper  = std::make_shared<idocp::JointPositionUpperLimit>(robot);
  auto joint_velocity_lower  = std::make_shared<idocp::JointVelocityLowerLimit>(robot);
  auto joint_velocity_upper  = std::make_shared<idocp::JointVelocityUpperLimit>(robot);
  auto joint_torques_lower   = std::make_shared<idocp::JointTorquesLowerLimit>(robot);
  auto joint_torques_upper   = std::make_shared<idocp::JointTorquesUpperLimit>(robot);
  const double mu = 0.7;
  auto friction_cone         = std::make_shared<idocp::LinearizedFrictionCone>(robot, mu);
  auto impulse_friction_cone = std::make_shared<idocp::LinearizedImpulseFrictionCone>(robot, mu);
  constraints->push_back(joint_position_lower);
  constraints->push_back(joint_position_upper);
  constraints->push_back(joint_velocity_lower);
  constraints->push_back(joint_velocity_upper);
  constraints->push_back(joint_torques_lower);
  constraints->push_back(joint_torques_upper);
  constraints->push_back(friction_cone);
  constraints->push_back(impulse_friction_cone);

  // 11 steps
  // const double T = 6.05; // t_start + max_num_impulse_phase * t_period + 0.05;
  // const int N = 120;
  // const int max_num_impulse_phase = 11;

  // // 10 steps
  // const double T = 5.55; // t_start + max_num_impulse_phase * t_period + 0.05;
  // const int N = 110;
  // const int max_num_impulse_phase = 10;

  // // 8 steps
  // const double T = 4.55; // t_start + max_num_impulse_phase * t_period + 0.05;
  // const int N = 90;
  // const int max_num_impulse_phase = 8;

  // // 6 steps
  // const double T = 3.55; // t_start + max_num_impulse_phase * t_period + 0.05;
  // const int N = 70;
  // const int max_num_impulse_phase = 6;

  // // 4 steps
  // const double T = 2.55; // t_start + max_num_impulse_phase * t_period + 0.05;
  // const int N = 50;
  // const int max_num_impulse_phase = 4;

  // 2 steps
  const double T = 1.55; // t_start + max_num_impulse_phase * t_period + 0.05;
  const int N = 30;
  const int max_num_impulse_phase = 2;

  const int nthreads = 4;
  const double t = 0;
  idocp::OCPSolver ocp_solver(robot, cost, constraints, T, N, max_num_impulse_phase+1, nthreads);

  robot.updateFrameKinematics(q_standing);
  std::vector<Eigen::Vector3d> contact_points(robot.maxPointContacts(), Eigen::Vector3d::Zero());
  robot.getContactPoints(contact_points);
  auto contact_status_initial = robot.createContactStatus();
  contact_status_initial.activateContacts({0, 1, 2, 3});
  contact_status_initial.setContactPoints(contact_points);
  ocp_solver.setContactStatusUniformly(contact_status_initial);

  auto contact_status_even = robot.createContactStatus();
  contact_status_even.activateContacts({1, 2});
  contact_status_even.setContactPoints(contact_points);
  ocp_solver.pushBackContactStatus(contact_status_even, t_start);

  auto contact_status_odd = robot.createContactStatus();
  contact_points[0].coeffRef(0) += 0.5 * step_length;
  contact_points[3].coeffRef(0) += 0.5 * step_length;
  contact_status_odd.activateContacts({0, 3});
  contact_status_odd.setContactPoints(contact_points);
  ocp_solver.pushBackContactStatus(contact_status_odd, t_start+t_period);

  for (int i=2; i<=max_num_impulse_phase; ++i) {
    if (i % 2 == 0) {
      contact_points[1].coeffRef(0) += step_length;
      contact_points[2].coeffRef(0) += step_length;
      contact_status_even.setContactPoints(contact_points);
      ocp_solver.pushBackContactStatus(contact_status_even, t_start+i*t_period);
    }
    else {
      contact_points[0].coeffRef(0) += step_length;
      contact_points[3].coeffRef(0) += step_length;
      contact_status_odd.setContactPoints(contact_points);
      ocp_solver.pushBackContactStatus(contact_status_odd, t_start+i*t_period);
    }
  }

  Eigen::VectorXd q(Eigen::VectorXd::Zero(robot.dimq()));
  q << 0, 0, 0.4792, 0, 0, 0, 1, 
       -0.1,  0.7, -1.0, // LF
       -0.1, -0.7,  1.0, // LH
        0.1,  0.7, -1.0, // RF
        0.1, -0.7,  1.0; // RH
  Eigen::VectorXd v(Eigen::VectorXd::Zero(robot.dimv()));

  ocp_solver.setSolution("q", q);
  ocp_solver.setSolution("v", v);
  Eigen::Vector3d f_init;
  f_init << 0, 0, 0.25*robot.totalWeight();
  ocp_solver.setSolution("f", f_init);

  ocp_solver.initConstraints(t);

  ocp_solver.setPenalty(1.0e02);


  const bool line_search = false;
  idocp::ocpbenchmarker::Convergence(ocp_solver, t, q, v, 25, line_search);

#ifdef ENABLE_VIEWER
  if (argc != 2) {
    std::cout << "Invalid argment!" << std::endl;
    std::cout << "Package serach path must be specified as the second argment!" << std::endl;
  }
  else {
    const std::string pkg_search_path = argv[1];
    idocp::TrajectoryViewer viewer(pkg_search_path, path_to_urdf);
    const double dt = T/N;
    viewer.display(ocp_solver.getSolution("q"), dt);
  }
#endif 

  return 0;
}