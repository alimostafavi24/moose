//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "AdvectionKernel2.h"
#include "MooseMesh.h"

registerMooseObject("NavierStokesApp", AdvectionKernel2);

InputParameters
AdvectionKernel2::validParams()
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

AdvectionKernel2::AdvectionKernel2(const InputParameters & parameters)
  : Kernel(parameters),

    // Coupled variables
    _u_vel(coupledValue("u")),
    _v_vel(_mesh.dimension() >= 2 ? coupledValue("v") : _zero),
    _w_vel(_mesh.dimension() == 3 ? coupledValue("w") : _zero),
    _grad_c(coupledGradient("c")),

    // Variable numberings
     _c_var(coupled("c")),
    _u_vel_var_number(coupled("u")),
    _v_vel_var_number(_mesh.dimension() >= 2 ? coupled("v") : libMesh::invalid_uint),
    _w_vel_var_number(_mesh.dimension() == 3 ? coupled("w") : libMesh::invalid_uint)
{
  
}

Real
AdvectionKernel2::computeQpResidual()
{

  Real convective_part =(_u_vel[_qp] * _grad_c[_qp](0) + _v_vel[_qp] * _grad_c[_qp](1) +
                          _w_vel[_qp] * _grad_c[_qp](2)) *
                         _test[_i][_qp];



  return convective_part;
}

Real
AdvectionKernel2::computeQpJacobian()
{
 return 0;
}



Real
AdvectionKernel2::computeQpOffDiagJacobian(unsigned jvar)
{
  if (jvar == _u_vel_var_number)
  {
    Real convective_part =  _phi[_j][_qp] * _grad_c[_qp](0) * _test[_i][_qp];
    return convective_part;
  }

  else if (jvar == _v_vel_var_number)
  {
    Real convective_part = _phi[_j][_qp] * _grad_c[_qp](1) * _test[_i][_qp];
    return convective_part;
  }

  else if (jvar == _w_vel_var_number)
  {
    Real convective_part =  _phi[_j][_qp] * _grad_c[_qp](2) * _test[_i][_qp];
    return convective_part;
  }
  else if (jvar == _c_var)
  {
    Real convective_part = (_u_vel[_qp] * _grad_phi[_j][_qp](0) + _v_vel[_qp] * _grad_phi[_j][_qp](1) + _w_vel[_qp] * _grad_phi[_j][_qp](2)) * _test[_i][_qp];
     return convective_part;
  }
  else
    return 0;
}