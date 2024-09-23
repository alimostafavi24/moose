//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "AdvectionKernel2PF.h"
#include "MooseMesh.h"

registerMooseObject("NavierStokesApp", AdvectionKernel2PF);

InputParameters
AdvectionKernel2PF::validParams()
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

AdvectionKernel2PF::AdvectionKernel2PF(const InputParameters & parameters)
  : Kernel(parameters),


    // Coupled variables
    _u_vel(coupledValue("u")),
    _v_vel(_mesh.dimension() >= 2 ? coupledValue("v") : _zero),
    _w_vel(_mesh.dimension() == 3 ? coupledValue("w") : _zero),
    _grad_c(coupledGradient("c")),

    // Variable numberings
     _c_var(coupled("c"))
  {
}

Real
AdvectionKernel2PF::computeQpResidual()
{

  Real convective_part =(_u_vel[_qp] * _grad_c[_qp](0) + _v_vel[_qp] * _grad_c[_qp](1) +
                          _w_vel[_qp] * _grad_c[_qp](2)) *
                         _test[_i][_qp];



  return convective_part;
}

Real
AdvectionKernel2PF::computeQpJacobian()
{
 return 0;
}



Real
AdvectionKernel2PF::computeQpOffDiagJacobian(unsigned jvar)
{
 if (jvar == _c_var)
  {
    Real convective_part = (_u_vel[_qp] * _grad_phi[_j][_qp](0) + _v_vel[_qp] * _grad_phi[_j][_qp](1) + _w_vel[_qp] * _grad_phi[_j][_qp](2)) * _test[_i][_qp];
     return convective_part;
  }
  else
    return 0;
}