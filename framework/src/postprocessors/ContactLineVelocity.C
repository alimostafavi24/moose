#include "ContactLineVelocity.h"
#include "MooseMesh.h"
#include "SubProblem.h"
#include "libmesh/system.h"
#include "libmesh/numeric_vector.h"

registerMooseObject("MooseApp", ContactLineVelocity);

InputParameters
ContactLineVelocity::validParams()
{
  InputParameters params = GeneralPostprocessor::validParams();
  params.addRequiredParam<VariableName>("c", "The order parameter (e.g., phase field)");
  params.addRequiredParam<VariableName>("variable", "Variable to extract (e.g., vel_x)");
  params.addParam<Real>("interface_tol", 1e-3, "Tolerance for detecting c=0 interface");
  params.addParam<Real>("wall_y", 0.0, "Wall y-location (points below are ignored)");
  return params;
}

ContactLineVelocity::ContactLineVelocity(const InputParameters & parameters)
  : GeneralPostprocessor(parameters),
    _c_var(getParam<VariableName>("c")),
    _target_var(getParam<VariableName>("variable")),
    _interface_tol(getParam<Real>("interface_tol")),
    _wall_y(getParam<Real>("wall_y")),
    _value_at_contact_line(std::numeric_limits<Real>::quiet_NaN()),
    _best_point(Point(std::numeric_limits<Real>::max(), std::numeric_limits<Real>::max(), 0)),
    _c_sys(_subproblem.getSystem(_c_var)),
    _target_sys(_subproblem.getSystem(_target_var))
{
}

void
ContactLineVelocity::initialize()
{
  _value_at_contact_line = std::numeric_limits<Real>::quiet_NaN();
  _best_point = Point(std::numeric_limits<Real>::max(), std::numeric_limits<Real>::max(), 0);
}

void
ContactLineVelocity::execute()
{
  const MeshBase & mesh = _subproblem.mesh().getMesh();
  const NumericVector<Number> & c_sol = *_c_sys.current_local_solution;
  const NumericVector<Number> & var_sol = *_target_sys.current_local_solution;

  const unsigned int c_sys_num = _c_sys.number();
  const unsigned int c_var_num = _c_sys.variable_number(_c_var);

  const unsigned int var_sys_num = _target_sys.number();
  const unsigned int var_var_num = _target_sys.variable_number(_target_var);

  const auto c_start = c_sol.first_local_index();
  const auto c_end   = c_sol.last_local_index();

  const auto var_start = var_sol.first_local_index();
  const auto var_end   = var_sol.last_local_index();

  Real min_y = std::numeric_limits<Real>::max();

  for (const auto & node : mesh.local_node_ptr_range())
  {
    if (!node || !node->active())
      continue;

    if (!node->has_dofs(c_sys_num) || !node->has_dofs(var_sys_num))
      continue;

    const Point & pt = *node;
    const Real y = pt(1);

    if (y <= _wall_y)
      continue;

    dof_id_type c_dof = node->dof_number(c_sys_num, c_var_num, 0);
    if (c_dof < c_start || c_dof >= c_end)
      continue;

    Real c_val = c_sol.el(c_dof);

    if (std::abs(c_val) < _interface_tol && y < min_y)
    {
      dof_id_type var_dof = node->dof_number(var_sys_num, var_var_num, 0);
      if (var_dof < var_start || var_dof >= var_end)
        continue;

      _value_at_contact_line = var_sol.el(var_dof);
      _best_point = pt;
      min_y = y;
    }
  }

  if (std::isnan(_value_at_contact_line))
    mooseWarning("No interface node found with |c| < ", _interface_tol, " and y > ", _wall_y);
  else
    mooseInfo("Contact line velocity evaluated at point (x, y) = (",
              _best_point(0), ", ", _best_point(1), ")");
}

Real
ContactLineVelocity::getValue() const
{
  return _value_at_contact_line;
}
