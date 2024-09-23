//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "IntegratedBC.h"
#include "TimeKernel.h"

/**
 * A FluxBC which is consistent with the boundary terms arising from
 * the Diffusion Kernel. The residual contribution is:
 *
 * \f$ F(u) = - \int_{\Gamma} \nabla u * \hat n * \phi d\Gamma \f$
 *
 * This class is essentially identical to the DiffusionFluxBC, but it
 * is not a part of the FluxBC hierarchy. It does not actually impose
 * any boundary condition, instead it computes the residual
 * contribution due to the boundary term arising from integration by
 * parts of the Diffusion Kernel.
 */
class WettingBCEnergy : public IntegratedBC
{
public:
  static InputParameters validParams();

  WettingBCEnergy(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
 // virtual Real computeQpOffDiagJacobian(unsigned jvar);
  //const VariableValue & _u;
  //const VariableGradient & _grad_u;
  const VariableValue & _u_dot;
  const VariableValue & _du_dot_du;
  /// Time derivative of u
 // const VariableValue & _u_dot;

  /// Derivative of u_dot with respect to u
  //const VariableValue & _du_dot_du;

  const VariableValue & _vel_x;
  const VariableValue & _vel_y;

  const Real _cos;
  const Real _sigma;
  const Real _gamma;
  const Real _eta;

    // Variable numberings
  
  //unsigned _c_var;  
 // unsigned _vel_x_vel_var_number;
 // unsigned _vel_y_vel_var_number;
//  unsigned _w_vel_var_number;

};
