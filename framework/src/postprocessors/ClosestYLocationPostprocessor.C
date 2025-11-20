#include "ClosestYLocationPostprocessor.h"

#include "MooseMesh.h"
#include "SubProblem.h"

#include "libmesh/mesh_base.h"
#include "libmesh/node.h"
#include "libmesh/numeric_vector.h"

#include <algorithm>
#include <cmath>
#include <limits>

registerMooseObject("MooseApp", ClosestYLocationPostprocessor);

InputParameters
ClosestYLocationPostprocessor::validParams()
{
  InputParameters params = GeneralPostprocessor::validParams();
  params.addRequiredParam<VariableName>("c", "The phase-field or order parameter variable name.");
  params.addParam<Real>("interface_tol", 1e-3, "Tolerance for identifying |c| ≈ 0.");
  params.addParam<Real>("wall_y", 0.0, "Wall y-location; nodes with y <= wall_y are ignored.");
  return params;
}

ClosestYLocationPostprocessor::ClosestYLocationPostprocessor(const InputParameters & parameters)
  : GeneralPostprocessor(parameters),
    _c_var(getParam<VariableName>("c")),
    _interface_tol(getParam<Real>("interface_tol")),
    _wall_y(getParam<Real>("wall_y")),
    _closest_y(std::numeric_limits<Real>::max()),
    _c_sys(_subproblem.getSystem(_c_var))
{
}

void
ClosestYLocationPostprocessor::initialSetup()
{
  // Build the candidate DOF list once.
  // If your mesh can change/partition at runtime, add logic to rebuild when needed.
  const libMesh::MeshBase & mesh = _subproblem.mesh().getMesh();

  const unsigned int sys_num = _c_sys.number();
  const unsigned int var_num = _c_sys.variable_number(_c_var);

  _candidates.clear();
  _candidates.reserve(mesh.n_local_nodes()); // upper bound; we'll push_back conditionally

  for (const auto & node : mesh.local_node_ptr_range())
  {
    if (!node || !node->active())
      continue;

    if (!node->has_dofs(sys_num))
      continue;

    const Real y = (*node)(1);
    if (y <= _wall_y)
      continue;

    // One nodal DOF for CG variables
    const libMesh::dof_id_type dof = node->dof_number(sys_num, var_num, 0);

    // We store the pair regardless of ownership; ownership can change if repartitioning occurs.
    _candidates.emplace_back(dof, y);
  }

  _candidates_built = true;
}

void
ClosestYLocationPostprocessor::initialize()
{
  _closest_y = std::numeric_limits<Real>::max();
}

void
ClosestYLocationPostprocessor::execute()
{
  // If, for any reason, initialSetup() didn't run, build now.
  if (!_candidates_built)
    initialSetup();

  const libMesh::NumericVector<libMesh::Number> & c_sol = *_c_sys.current_local_solution;

  // Fast path: loop only over prefiltered candidates
  for (const auto & pr : _candidates)
  {
    const libMesh::dof_id_type dof = pr.first;
    const Real y = pr.second;

    // Only owned DOFs: avoid invalid access and remote gets
    if (dof < c_sol.first_local_index() || dof >= c_sol.last_local_index())
      continue;

    const Real c_val = c_sol.el(dof);

    if (std::abs(c_val) < _interface_tol && y < _closest_y)
      _closest_y = y;
  }
}

void
ClosestYLocationPostprocessor::threadJoin(const UserObject & other_uo)
{
  const auto & other = static_cast<const ClosestYLocationPostprocessor &>(other_uo);
  _closest_y = std::min(_closest_y, other._closest_y);
}

void
ClosestYLocationPostprocessor::finalize()
{
  // In-place MPI min reduction
  _communicator.min(_closest_y);
}

Real
ClosestYLocationPostprocessor::getValue() const
{
  if (_closest_y == std::numeric_limits<Real>::max())
  {
    mooseWarning("No interface node found with |c| < ",
                 _interface_tol,
                 " and y > ",
                 _wall_y);
    return std::numeric_limits<Real>::quiet_NaN();
  }
  return _closest_y;
}
