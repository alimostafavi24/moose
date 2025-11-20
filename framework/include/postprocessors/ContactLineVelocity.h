#pragma once

#include "GeneralPostprocessor.h"

/**
 * ContactLineVelocity:
 * Finds the closest point above a wall where |c| < interface_tol,
 * and returns the value of a user-specified variable at that location.
 */
class ContactLineVelocity : public GeneralPostprocessor
{
public:
  static InputParameters validParams();

  ContactLineVelocity(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual void finalize() override {}
  virtual Real getValue() const override;

protected:
  const VariableName _c_var;
  const VariableName _target_var;
  const Real _interface_tol;
  const Real _wall_y;

  Real _value_at_contact_line;
  Point _best_point;

  const libMesh::System & _c_sys;
  const libMesh::System & _target_sys;
};
