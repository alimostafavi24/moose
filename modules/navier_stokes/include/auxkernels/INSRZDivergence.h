#pragma once

#include "AuxKernel.h"

// Forward Declarations

/**
 * Computes h_min / |u|
 */
class INSRZDivergence : public AuxKernel
{
public:
  static InputParameters validParams();

  INSRZDivergence(const InputParameters & parameters);

  virtual ~INSRZDivergence() {}

protected:
  virtual Real computeValue();

  // Velocity gradients
  const VariableValue & _vel_x;
  const VariableGradient & _grad_vel_x;
  const VariableGradient & _grad_vel_y;

};
