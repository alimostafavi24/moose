#pragma once

#include "AuxKernel.h"

// Forward Declarations

/**
 * Computes h_min / |u|
 */
class NormalVectorPF: public AuxKernel
{
public:
  static InputParameters validParams();

  NormalVectorPF(const InputParameters & parameters);

  virtual ~NormalVectorPF() {}

protected:
  virtual Real computeValue();

  // Velocity gradients
  const VariableValue & _c;
  const VariableGradient & _grad_c;

  unsigned int _component;
};
