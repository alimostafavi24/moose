//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "AdvectionKernel222Multi.h"
#include "MooseMesh.h"

registerMooseObject("NavierStokesApp", AdvectionKernel222Multi);

InputParameters
AdvectionKernel222Multi::validParams()
{
  InputParameters params = Kernel::validParams();

  params.addClassDescription("Adds advection term to the CH equation.");
  // Coupled variables
  params.addRequiredCoupledVar("u", "x-velocity");
  params.addCoupledVar("v", "y-velocity"); // only required in 2D and 3D
  params.addCoupledVar("w", "z-velocity"); // only required in 3D
  params.addRequiredCoupledVar("c", "grad c"); //c is the order parameter

  return params;
}

AdvectionKernel222Multi::AdvectionKernel222Multi(const InputParameters & parameters)
  : Kernel(parameters),

    // Coupled variables
    _grad_u_vel(coupledGradient("u")),
    _grad_v_vel(coupledGradient("v")),
    _grad_w_vel(coupledGradient("w")),
    _c(coupledValue("c")),

    // Variable numberings
     _c_var(coupled("c"))
{
  
}

Real
AdvectionKernel222Multi::computeQpResidual()
{

  Real convective_part =_c[_qp] * (_grad_u_vel[_qp](0)+_grad_v_vel[_qp](1) + +_grad_w_vel[_qp](2)) *
                         _test[_i][_qp];



  return convective_part;
}

Real
AdvectionKernel222Multi::computeQpJacobian()
{
 return 0;
}



Real
AdvectionKernel222Multi::computeQpOffDiagJacobian(unsigned jvar)
{
 if (jvar == _c_var)
  {
    Real convective_part = _phi[_j][_qp]  * (_grad_u_vel[_qp](0)+_grad_v_vel[_qp](1)) *
                         _test[_i][_qp];;
     return convective_part;
  }
  else
    return 0;
}