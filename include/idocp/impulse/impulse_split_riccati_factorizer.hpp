#ifndef IDOCP_IMPULSE_SPLIT_RICCATI_FACTORIZER_HPP_
#define IDOCP_IMPULSE_SPLIT_RICCATI_FACTORIZER_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_riccati_factorization.hpp"
#include "idocp/impulse/impulse_split_kkt_matrix.hpp"
#include "idocp/impulse/impulse_split_kkt_residual.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/impulse/impulse_backward_riccati_recursion_factorizer.hpp"
#include "idocp/impulse/impulse_forward_riccati_recursion_factorizer.hpp"


namespace idocp {

///
/// @class ImpulseSplitRiccatiFactorizer
/// @brief Riccati factorizer for ImpulseSplitOCP.
///
class ImpulseSplitRiccatiFactorizer {
public:
  ///
  /// @brief Construct factorizer.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  ///
  ImpulseSplitRiccatiFactorizer(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  ImpulseSplitRiccatiFactorizer();

  ///
  /// @brief Destructor. 
  ///
  ~ImpulseSplitRiccatiFactorizer();
 
  ///
  /// @brief Default copy constructor. 
  ///
  ImpulseSplitRiccatiFactorizer(
      const ImpulseSplitRiccatiFactorizer&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  ImpulseSplitRiccatiFactorizer& operator=(
      const ImpulseSplitRiccatiFactorizer&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  ImpulseSplitRiccatiFactorizer(
      ImpulseSplitRiccatiFactorizer&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  ImpulseSplitRiccatiFactorizer& operator=(
      ImpulseSplitRiccatiFactorizer&&) noexcept = default;

  ///
  /// @brief Performs the backward Riccati recursion. 
  /// @param[in] riccati_next Riccati factorization at the next time stage. 
  /// @param[in, out] kkt_matrix KKT matrix at the current impulse stage. 
  /// @param[in, out] kkt_residual KKT residual at the current impulse stage. 
  /// @param[out] riccati Riccati factorization at the current impulse stage. 
  ///
  void backwardRiccatiRecursion(const SplitRiccatiFactorization& riccati_next, 
                                ImpulseSplitKKTMatrix& kkt_matrix, 
                                ImpulseSplitKKTResidual& kkt_residual, 
                                SplitRiccatiFactorization& riccati);

  ///
  /// @brief Performs the serial part of the forward Riccati recursion due to 
  /// pure-state equality constraints. 
  /// @param[in] riccati Riccati factorization at the current impulse stage. 
  /// @param[in] kkt_matrix KKT matrix at the current impulse stage. 
  /// @param[in] kkt_residual KKT residual at the current impulse stage. 
  /// @param[out] riccati_next Riccati factorization at the next time stage. 
  /// @param[in] exist_state_constraint If true, the factorization for state
  /// constraints are also performed. 
  ///
  void forwardStateConstraintFactorization(
      const SplitRiccatiFactorization& riccati, 
      const ImpulseSplitKKTMatrix& kkt_matrix, 
      const ImpulseSplitKKTResidual& kkt_residual, 
      SplitRiccatiFactorization& riccati_next, 
      const bool exist_state_constraint);

  ///
  /// @brief Performs the backward factorization of matrices related to the 
  /// pure-state equality constraints. 
  /// @param[in] kkt_matrix KKT matrix at the current impulse stage. 
  /// @param[in] T_next A factorization at the next time stage. 
  /// @param[out] T A factorization at the current impulse stage. 
  ///
  template <typename MatrixType1, typename MatrixType2>
  void backwardStateConstraintFactorization(
      const Eigen::MatrixBase<MatrixType1>& T_next,  
      const ImpulseSplitKKTMatrix& kkt_matrix, 
      const Eigen::MatrixBase<MatrixType2>& T) const;

  ///
  /// @brief Performs forward Riccati recursion and computes state direction. 
  /// @param[in] kkt_matrix KKT matrix at the current time stage. 
  /// @param[in] kkt_residual KKT residual at the current time stage. 
  /// @param[in] d Split direction at the current time stage. 
  /// @param[out] d_next Split direction at the next time stage. 
  ///
  void forwardRiccatiRecursion(const ImpulseSplitKKTMatrix& kkt_matrix, 
                               const ImpulseSplitKKTResidual& kkt_residual,
                               const ImpulseSplitDirection& d, 
                               SplitDirection& d_next) const;

  ///
  /// @brief Computes the Newton direction of the costate vector. 
  /// @param[in] riccati Riccati factorization at the current impulse stage. 
  /// @param[out] d Split direction of the current impulse stage. 
  /// @param[in] exist_state_constraint If true, the factorization for state
  /// constraints are also performed.
  ///
  ///
  static void computeCostateDirection(const SplitRiccatiFactorization& riccati, 
                                      ImpulseSplitDirection& d,
                                      const bool exist_state_constraint);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  bool has_floating_base_;
  int dimv_;
  ImpulseBackwardRiccatiRecursionFactorizer backward_recursion_;
  ImpulseForwardRiccatiRecursionFactorizer forward_recursion_;

};

} // namespace idocp

#include "idocp/impulse/impulse_split_riccati_factorizer.hxx"

#endif // IDOCP_IMPULSE_SPLIT_RICCATI_FACTORIZER_HPP_ 