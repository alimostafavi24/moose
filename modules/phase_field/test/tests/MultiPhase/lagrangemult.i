[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 14
  ny = 10
  nz = 0
  xmin = 10
  xmax = 40
  ymin = 15
  ymax = 35
  elem_type = QUAD4
[]

[Variables]
  [./c]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = SmoothCircleIC
      x1 = 25.0
      y1 = 25.0
      radius = 6.0
      invalue = 0.9
      outvalue = 0.1
      int_width = 3.0
    [../]
  [../]
  [./w]
    order = FIRST
    family = LAGRANGE
  [../]

  [./eta1]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = SmoothCircleIC
      x1 = 30.0
      y1 = 25.0
      radius = 4.0
      invalue = 0.9
      outvalue = 0.1
      int_width = 2.0
    [../]
  [../]
  [./eta2]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.5
  [../]

  [./lambda]
    order = FIRST
    family = LAGRANGE
    initial_condition = 1.0
  [../]
[]

[Kernels]
  [./deta1dt]
    type = TimeDerivative
    variable = eta1
  [../]
  [./ACBulk1]
    type = AllenCahn
    variable = eta1
    coupled_variables = 'c eta2'
    f_name = F
  [../]
  [./ACInterface1]
    type = ACInterface
    variable = eta1
    kappa_name = kappa_eta
  [../]
  [./lagrange1]
    type = SwitchingFunctionConstraintEta
    variable = eta1
    h_name   = h1
    lambda = lambda
  [../]

  [./deta2dt]
    type = TimeDerivative
    variable = eta2
  [../]
  [./ACBulk2]
    type = AllenCahn
    variable = eta2
    coupled_variables = 'c eta1'
    f_name = F
  [../]
  [./ACInterface2]
    type = ACInterface
    variable = eta2
    kappa_name = kappa_eta
  [../]
  [./lagrange2]
    type = SwitchingFunctionConstraintEta
    variable = eta2
    h_name   = h2
    lambda = lambda
  [../]

  [./lagrange]
    type = SwitchingFunctionConstraintLagrange
    variable = lambda
    etas    = 'eta1 eta2'
    h_names = 'h1   h2'
    epsilon = 0
  [../]

  [./c_res]
    type = SplitCHParsed
    variable = c
    f_name = F
    kappa_name = kappa_c
    w = w
    coupled_variables = 'eta1 eta2'
  [../]
  [./w_res]
    type = SplitCHWRes
    variable = w
    mob_name = M
  [../]
  [./time1]
    type = CoupledTimeDerivative
    variable = w
    v = c
  [../]
[]

[BCs]
  [./Periodic]
    [./All]
      auto_direction = 'x y'
    [../]
  [../]
[]

[Materials]
  [./consts]
    type = GenericConstantMaterial
    prop_names  = 'L kappa_eta'
    prop_values = '1 1        '
  [../]
  [./consts2]
    type = GenericConstantMaterial
    prop_names  = 'M kappa_c'
    prop_values = '1 1'
  [../]

  [./switching1]
    type = SwitchingFunctionMaterial
    function_name = h1
    eta = eta1
    h_order = SIMPLE
    outputs = exodus
  [../]
  [./switching2]
    type = SwitchingFunctionMaterial
    function_name = h2
    eta = eta2
    h_order = SIMPLE
    outputs = exodus
  [../]

  [./barrier]
    type = MultiBarrierFunctionMaterial
    etas = 'eta1 eta2'
  [../]

  [./free_energy_A]
    type = DerivativeParsedMaterial
    property_name = Fa
    coupled_variables = 'c'
    expression = '(c-0.1)^2'
    derivative_order = 2
    enable_jit = true
  [../]
  [./free_energy_B]
    type = DerivativeParsedMaterial
    property_name = Fb
    coupled_variables = 'c'
    expression = '(c-0.9)^2'
    derivative_order = 2
    enable_jit = true
  [../]

  [./free_energy]
    type = DerivativeMultiPhaseMaterial
    f_name = F
    fi_names = 'Fa   Fb'
    hi_names = 'h1   h2'
    etas     = 'eta1 eta2'
    coupled_variables = 'c'
    derivative_order = 2
  [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  scheme = 'bdf2'

  solve_type = 'PJFNK'
  #petsc_options = '-snes_ksp -snes_ksp_ew'
  #petsc_options = '-ksp_monitor_snes_lg-snes_ksp_ew'
  #petsc_options_iname = '-ksp_gmres_restart'
  #petsc_options_value = '1000              '

  l_max_its = 15
  l_tol = 1.0e-6

  nl_max_its = 50
  nl_rel_tol = 1.0e-8
  nl_abs_tol = 1.0e-10

  start_time = 0.0
  num_steps = 1
  dt = 0.01
  dtmin = 0.01
[]

[Debug]
  # show_var_residual_norms = true
[]

[Outputs]
  execute_on = 'timestep_end'
  exodus = true
[]
