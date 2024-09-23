//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "INSRZDivergence.h"
#include "MooseMesh.h"
#include "NS.h"

registerMooseObject("NavierStokesApp", INSRZDivergence);

InputParameters
INSRZDivergence::validParams()
{
  InputParameters params = AuxKernel::validParams();

  params.addClassDescription("This class computes the stress component based on "
                             "pressure and velocity for incompressible Navier-Stokes");
  params.addCoupledVar("vel_x", "vel_x");
  params.addCoupledVar("vel_y", "vel_y");

  return params;
}

INSRZDivergence::INSRZDivergence(const InputParameters & parameters)
  : AuxKernel(parameters),
    _vel_x(coupledValue("vel_x")),
    _grad_vel_x(coupledGradient("vel_x")),
    _grad_vel_y(coupledGradient("vel_y"))

{
}

Real
INSRZDivergence::computeValue()
{
  //const auto r = _q_point[_qp](0);

  return _grad_vel_x[_qp](0)+_grad_vel_y[_qp](1)+_vel_x[_qp]/_q_point[_qp](0);
}
