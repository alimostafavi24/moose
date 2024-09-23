//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "INSRZDiss.h"
#include "MooseMesh.h"
#include "NS.h"

registerMooseObject("NavierStokesApp", INSRZDiss);

InputParameters
INSRZDiss::validParams()
{
  InputParameters params = AuxKernel::validParams();

  params.addClassDescription("This class computes the stress component based on "
                             "pressure and velocity for incompressible Navier-Stokes");
  params.addCoupledVar("vel_x", "vel_x");
  params.addCoupledVar("vel_y", "vel_y");
  params.addParam<MaterialPropertyName>("mu_name", "mu", "The name of the mu");


  return params;
}

INSRZDiss::INSRZDiss(const InputParameters & parameters)
  : AuxKernel(parameters),
    _vel_x(coupledValue("vel_x")),
    _grad_vel_x(coupledGradient("vel_x")),
    _second_vel_x(coupledSecond("vel_x")),
    _vel_y(coupledValue("vel_y")),
    _grad_vel_y(coupledGradient("vel_y")), 
    _second_vel_y(coupledSecond("vel_y")),
    _mu(getMaterialProperty<Real>("mu_name"))
{
}

Real
INSRZDiss::computeValue()
{
  //const auto r = _q_point[_qp](0);
  //return _mu[_qp] * ( 2 * ( _grad_vel_x[_qp](0)*_grad_vel_x[_qp](0) + _vel_x[_qp]*_vel_x[_qp] + _grad_vel_y[_qp](1)*_grad_vel_y[_qp](1))
  //                + (_grad_vel_x[_qp](1)+_grad_vel_y[_qp](0)) * (_grad_vel_x[_qp](1)+_grad_vel_y[_qp](0))
 //                 - 2/3 *(_grad_vel_x[_qp](0)+_grad_vel_y[_qp](1)+_vel_x[_qp]/_q_point[_qp](0))*(_grad_vel_x[_qp](0)+_grad_vel_y[_qp](1)+_vel_x[_qp]/_q_point[_qp](0)) );
  return +1 * _mu[_qp] * 
  ((_vel_x[_qp] * (_second_vel_x[_qp](0,0) + _second_vel_x[_qp](1,1) + (1/_q_point[_qp](0)*_grad_vel_x[_qp](0)) - (_vel_x[_qp] / (_q_point[_qp](0)*_q_point[_qp](0)) ) ) +
  _vel_y[_qp] * (_second_vel_y[_qp](0,0) + _second_vel_y[_qp](1,1) + (1/_q_point[_qp](0)*_grad_vel_y[_qp](0)) )));
}
