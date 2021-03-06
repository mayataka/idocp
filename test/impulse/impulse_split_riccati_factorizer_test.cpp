#include <string>
#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/impulse/impulse_split_kkt_matrix.hpp"
#include "idocp/impulse/impulse_split_kkt_residual.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/split_riccati_factorization.hpp"
#include "idocp/impulse/impulse_backward_riccati_recursion_factorizer.hpp"
#include "idocp/impulse/impulse_split_riccati_factorizer.hpp"


namespace idocp {

class ImpulseSplitRiccatiFactorizerTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
    fixed_base_robot = Robot(fixed_base_urdf);
    floating_base_robot = Robot(floating_base_urdf);
  }

  virtual void TearDown() {
  }

  static ImpulseSplitKKTMatrix createKKTMatrix(const Robot& robot);
  static ImpulseSplitKKTResidual createKKTResidual(const Robot& robot);
  static SplitRiccatiFactorization createRiccatiFactorization(const Robot& robot);

  static void testBackwardRecursion(const Robot& robot);
  static void testForwardRecursion(const Robot& robot);

  std::string fixed_base_urdf, floating_base_urdf;
  Robot fixed_base_robot, floating_base_robot;
};


ImpulseSplitKKTMatrix ImpulseSplitRiccatiFactorizerTest::createKKTMatrix(const Robot& robot) {
  const int dimv = robot.dimv();
  Eigen::MatrixXd seed = Eigen::MatrixXd::Random(dimv, dimv);
  ImpulseSplitKKTMatrix kkt_matrix(robot);
  kkt_matrix.Qqq() = seed * seed.transpose();
  kkt_matrix.Qqv().setRandom();
  kkt_matrix.Qvq() = kkt_matrix.Qqv().transpose();
  seed = Eigen::MatrixXd::Random(dimv, dimv);
  kkt_matrix.Qvv() = seed * seed.transpose();
  if (robot.hasFloatingBase()) {
    kkt_matrix.Fqq().setIdentity();
    kkt_matrix.Fqq().topLeftCorner(robot.dim_passive(), robot.dim_passive()).setRandom();
  }
  kkt_matrix.Fqv().setZero();
  kkt_matrix.Fvq().setRandom();
  kkt_matrix.Fvv().setRandom();
  return kkt_matrix;
}


ImpulseSplitKKTResidual ImpulseSplitRiccatiFactorizerTest::createKKTResidual(const Robot& robot) {
  ImpulseSplitKKTResidual kkt_residual(robot);
  kkt_residual.lx().setRandom();
  kkt_residual.Fx().setRandom();
  return kkt_residual;
}


SplitRiccatiFactorization ImpulseSplitRiccatiFactorizerTest::createRiccatiFactorization(const Robot& robot) {
  SplitRiccatiFactorization riccati(robot);
  const int dimv = robot.dimv();
  Eigen::MatrixXd seed = Eigen::MatrixXd::Random(dimv, dimv);
  riccati.Pqq = seed * seed.transpose();
  riccati.Pqv = Eigen::MatrixXd::Random(dimv, dimv);
  riccati.Pvq = riccati.Pqv.transpose();
  seed = Eigen::MatrixXd::Random(dimv, dimv);
  riccati.Pvv = seed * seed.transpose();
  riccati.sq.setRandom();
  riccati.sv.setRandom();
  return riccati; 
}


void ImpulseSplitRiccatiFactorizerTest::testBackwardRecursion(const Robot& robot) {
  const int dimv = robot.dimv();
  const SplitRiccatiFactorization riccati_next = createRiccatiFactorization(robot);
  ImpulseSplitKKTMatrix kkt_matrix = createKKTMatrix(robot);
  ImpulseSplitKKTResidual kkt_residual = createKKTResidual(robot);
  ImpulseSplitKKTMatrix kkt_matrix_ref = kkt_matrix;
  ImpulseSplitKKTResidual kkt_residual_ref = kkt_residual;
  ImpulseSplitRiccatiFactorizer factorizer(robot);
  ImpulseBackwardRiccatiRecursionFactorizer backward_recursion_ref(robot);
  SplitRiccatiFactorization riccati = createRiccatiFactorization(robot);
  SplitRiccatiFactorization riccati_ref = riccati;
  factorizer.backwardRiccatiRecursion(riccati_next, kkt_matrix, kkt_residual, riccati);
  backward_recursion_ref.factorizeKKTMatrix(riccati_next, kkt_matrix_ref);
  backward_recursion_ref.factorizeRiccatiFactorization(riccati_next, kkt_matrix_ref, kkt_residual_ref, riccati_ref);
  EXPECT_TRUE(riccati.Pqq.isApprox(riccati_ref.Pqq));
  EXPECT_TRUE(riccati.Pqv.isApprox(riccati_ref.Pqv));
  EXPECT_TRUE(riccati.Pvq.isApprox(riccati_ref.Pvq));
  EXPECT_TRUE(riccati.Pvv.isApprox(riccati_ref.Pvv));
  EXPECT_TRUE(riccati.sq.isApprox(riccati_ref.sq));
  EXPECT_TRUE(riccati.sv.isApprox(riccati_ref.sv));
  EXPECT_TRUE(riccati.Pqq.isApprox(riccati.Pqq.transpose()));
  EXPECT_TRUE(riccati.Pvv.isApprox(riccati.Pvv.transpose()));
  EXPECT_TRUE(riccati.Pvq.isApprox(riccati.Pqv.transpose()));
  EXPECT_TRUE(kkt_matrix.Qxx().isApprox(kkt_matrix.Qxx().transpose()));
}


void ImpulseSplitRiccatiFactorizerTest::testForwardRecursion(const Robot& robot) {
  const int dimv = robot.dimv();
  SplitRiccatiFactorization riccati_next = createRiccatiFactorization(robot);
  SplitRiccatiFactorization riccati_next_ref = riccati_next;
  ImpulseSplitKKTMatrix kkt_matrix = createKKTMatrix(robot);
  ImpulseSplitKKTResidual kkt_residual = createKKTResidual(robot);
  ImpulseSplitRiccatiFactorizer factorizer(robot);
  SplitRiccatiFactorization riccati = createRiccatiFactorization(robot);
  SplitRiccatiFactorization riccati_ref = riccati;
  factorizer.backwardRiccatiRecursion(riccati_next, kkt_matrix, kkt_residual, riccati);
  ImpulseSplitKKTMatrix kkt_matrix_ref = kkt_matrix;
  ImpulseSplitKKTResidual kkt_residual_ref = kkt_residual;
  ImpulseSplitDirection d(robot);
  d.setRandom();
  SplitDirection d_next(robot);
  d_next.setRandom();
  SplitDirection d_next_ref = d_next;
  factorizer.forwardRiccatiRecursion(kkt_matrix, kkt_residual, d, d_next);
  if (!robot.hasFloatingBase()) {
    kkt_matrix_ref.Fqq().setIdentity();
  }
  d_next_ref.dx() = kkt_matrix_ref.Fxx() * d.dx() + kkt_residual_ref.Fx();
  EXPECT_TRUE(d_next.isApprox(d_next_ref));
  ImpulseSplitDirection d_ref = d;
  factorizer.computeCostateDirection(riccati, d);
  d_ref.dlmd() = riccati.Pqq * d.dq() + riccati.Pqv * d.dv() - riccati.sq;
  d_ref.dgmm() = riccati.Pvq * d.dq() + riccati.Pvv * d.dv() - riccati.sv;
  EXPECT_TRUE(d.isApprox(d_ref));
}


TEST_F(ImpulseSplitRiccatiFactorizerTest, fixedBase) {
  testBackwardRecursion(fixed_base_robot);
  testForwardRecursion(fixed_base_robot);
}


TEST_F(ImpulseSplitRiccatiFactorizerTest, floating_base) {
  testBackwardRecursion(floating_base_robot);
  testForwardRecursion(floating_base_robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}