

//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "Electric.h"
#include "NS.h"
#include "MooseVariableFE.h"
#include "SystemBase.h"

registerMooseObject("NavierStokesApp", Electric);

InputParameters
Electric::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Adds the Columb force to the NS equation.");  
  params.addRequiredParam<unsigned int>("component", "number of component (0 = x, 1 = y, 2 = z)");
  params.addRequiredParam<MaterialPropertyName>("permittivity", "permittivity");
  params.addRequiredParam<MaterialPropertyName>("dpermittivitydc", "dpermittivitydc");
  params.addRequiredParam<MaterialPropertyName>("d2permittivitydc2", "d2permittivitydc2");
  params.addRequiredCoupledVar("phit", "Potential");
  params.addRequiredCoupledVar("rhoe", "Charge density");
  params.addRequiredCoupledVar("c", "c");


  return params;
}

Electric::Electric(const InputParameters & parameters)
  : Kernel(parameters),
    _permittivity(getMaterialProperty<Real>("permittivity")),
    _dpermittivitydc(getMaterialProperty<Real>("dpermittivitydc")),
    _d2permittivitydc2(getMaterialProperty<Real>("d2permittivitydc2")),

    _grad_phit(coupledGradient("phit")),
    _rhoe(coupledValue("rhoe")),
    _grad_c(coupledGradient("c")),

    _component(getParam<unsigned int>("component")),
    _phit_var(coupled("phit")),
    _rhoe_var(coupled("rhoe")),
    _c_var(coupled("c"))
{
}

Real
Electric::computeQpResidual()
{

 Real Electric_1 = (-0.5 * _grad_phit[_qp] * _grad_phit[_qp] * _dpermittivitydc[_qp] * _grad_c[_qp](_component));

 Real Electric_2 = (-1.0 * _rhoe[_qp] * _grad_phit[_qp](_component));

 return (Electric_1 + Electric_2) * -1 * _test[_i][_qp];

}

Real
Electric::computeQpJacobian()
{

  return 0.;
}

Real
Electric::computeQpOffDiagJacobian(unsigned jvar)
{


  if (jvar == _phit_var)
  {

 Real Electric_1 = (-1.0 * _grad_phit[_qp] * _grad_phi[_j][_qp] * _dpermittivitydc[_qp] * _grad_c[_qp](_component));

 Real Electric_2 = (-1.0 * _rhoe[_qp] * _grad_phi[_j][_qp](_component));

 return (Electric_1 + Electric_2) * -1 * _test[_i][_qp];

   }

     if (jvar == _rhoe_var)
  {
 Real Electric_1 = 0;

 Real Electric_2 = (-1.0 * _phi[_j][_qp] * _grad_phit[_qp](_component));

 return (Electric_1 + Electric_2) * -1 * _test[_i][_qp];
   }


        if (jvar == _c_var)
  {

 Real Electric_1 = (-0.5 * _grad_phit[_qp] * _grad_phit[_qp]) *
  (_d2permittivitydc2[_qp] * _grad_c[_qp](_component) + _dpermittivitydc[_qp] * _grad_phi[_j][_qp](_component));

 Real Electric_2 = 0;

 return (Electric_1 + Electric_2) * -1 * _test[_i][_qp];

   }
   
  else
  return 0.;
}
