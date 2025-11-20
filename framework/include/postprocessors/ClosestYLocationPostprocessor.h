#pragma once

#include "GeneralPostprocessor.h"
#include "libmesh/system.h"
#include <vector>
#include <utility>

/**
 * Postprocessor that finds the smallest y-location above a specified wall
 * where the interface criterion |c| < interface_tol is met.
 *
 * Performance notes:
 *  - We precompute candidate (dof_id, y) pairs once in initialSetup()
 *    to avoid re-walking the mesh each execute().
 *  - We still check DOF ownership in execute() to be safe under repartitioning.
 *  - Per-thread minima are joined in threadJoin(), then an MPI min in finalize().
 */
class ClosestYLocationPostprocessor : public GeneralPostprocessor
{
public:
  static InputParameters validParams();

  ClosestYLocationPostprocessor(const InputParameters & parameters);

  // GeneralPostprocessor interface
  virtual void initialSetup() override;
  virtual void initialize() override;
  virtual void execute() override;
  virtual void threadJoin(const UserObject & other_uo) override;
  virtual void finalize() override;
  virtual Real getValue() const override;

protected:
  /// Variable name
  const VariableName _c_var;

  /// Interface tolerance
  const Real _interface_tol;

  /// Wall y-location
  const Real _wall_y;

  /// Running minimum (per thread, then per rank, then global)
  Real _closest_y;

  /// System for variable c
  libMesh::System & _c_sys;

  /// Precomputed candidate list: (dof_id, y)
  std::vector<std::pair<libMesh::dof_id_type, Real>> _candidates;

  /// Flag to indicate the candidate list has been built
  bool _candidates_built = false;
};
