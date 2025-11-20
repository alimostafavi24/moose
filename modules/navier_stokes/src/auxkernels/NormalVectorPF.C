//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "NormalVectorPF.h"
#include "MooseMesh.h"
#include "NS.h"

registerMooseObject("NavierStokesApp", NormalVectorPF);

InputParameters
NormalVectorPF::validParams()
{
  InputParameters params = AuxKernel::validParams();

  params.addClassDescription("This class computes Normal Vector");
  params.addCoupledVar("c", "c");
  params.addRequiredParam<unsigned int>("component", "number of component (0 = x, 1 = y, 2 = z)");


  return params;
}

NormalVectorPF::NormalVectorPF(const InputParameters & parameters)
  : AuxKernel(parameters),
    _c(coupledValue("c")),
    _grad_c(coupledGradient("c")),
  _component(getParam<unsigned int>("component")) 
{
}

Real
NormalVectorPF::computeValue()
{

  return +1 * _grad_c[_qp](_component) / (_grad_c[_qp].norm() + std::numeric_limits<Real>::epsilon());
}
