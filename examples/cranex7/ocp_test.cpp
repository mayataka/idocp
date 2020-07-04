#include <iostream>
#include <string>
#include <chrono>

#include "Eigen/Core"

#include "ocp/ocp.hpp"
#include "robot/robot.hpp"
#include "manipulator/cost_function.hpp"
#include "manipulator/constraints.hpp"


int main() {
  srand((unsigned int) time(0));
  const std::string urdf_file_name = "../crane_x7_description/urdf/crane_x7.urdf";
  idocp::Robot robot(urdf_file_name);
  robot.setJointDamping(Eigen::VectorXd::Constant(7, 1.0e-06));
  idocp::manipulator::CostFunction cost(robot);
  idocp::manipulator::Constraints constraints(robot);
  Eigen::VectorXd q_ref = 1.5 * Eigen::VectorXd::Random(robot.dimq());
  q_ref.fill(1);
  cost.set_q_ref(q_ref);
  const double T = 1;
  const unsigned int N = 50;
  const unsigned int num_proc = 4;
  idocp::OCP ocp_(robot, &cost, &constraints, T, N, num_proc);
  const double t = 0;
  Eigen::VectorXd q = Eigen::VectorXd::Zero(robot.dimq());
  Eigen::VectorXd v = Eigen::VectorXd::Zero(robot.dimv());
  // ocp_.setStateTrajectory(q, v);
  // std::cout << ocp_.KKTError(t, q, v) << std::endl;
  const int num_iteration = 20;
  std::chrono::system_clock::time_point start_clock, end_clock;
  start_clock = std::chrono::system_clock::now();
  for (int i=0; i<num_iteration; ++i) {
    ocp_.solveSQP(t, q, v, false);
    std::cout << ocp_.KKTError(t, q, v) << std::endl;
  }
  end_clock = std::chrono::system_clock::now();
  // ocp_.printSolution();
  std::cout << "q_ref = " << q_ref.transpose() << std::endl;
  std::cout << "total CPU time: " << 1e-03 * std::chrono::duration_cast<std::chrono::microseconds>(end_clock-start_clock).count() << "[ms]" << std::endl;
  std::cout << "CPU time per update: " << 1e-03 * std::chrono::duration_cast<std::chrono::microseconds>(end_clock-start_clock).count() / num_iteration << "[ms]" << std::endl;
  Eigen::VectorXd u;
  u.resize(robot.dimv());
  ocp_.getInitialControlInput(u);
  std::cout << "u = " << u.transpose() << std::endl;

  return 0;
}