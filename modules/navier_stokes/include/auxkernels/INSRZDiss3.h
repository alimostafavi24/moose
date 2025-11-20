#pragma once

#include "AuxKernel.h"

// Forward Declarations

/**
 * Computes h_min / |u|
 */
class INSRZDiss3: public AuxKernel
{
public:
  static InputParameters validParams();

  INSRZDiss3(const InputParameters & parameters);

  virtual ~INSRZDiss3() {}

protected:
  virtual Real computeValue();

  // Velocity gradients
  const VariableValue & _vel_x;
  const VariableGradient & _grad_vel_x;
  const VariableSecond & _second_vel_x;

  const VariableValue & _vel_y;
  const VariableGradient & _grad_vel_y;
  const VariableSecond & _second_vel_y;
  const MaterialProperty<Real> & _mu;
};
