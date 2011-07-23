!> \file
!> \author Chris Bradley
!> \brief This module handles all problem wide constants.
!>
!> \section LICENSE
!>
!> Version: MPL 1.1/GPL 2.0/LGPL 2.1
!>
!> The contents of this file are subject to the Mozilla Public License
!> Version 1.1 (the "License"); you may not use this file except in
!> compliance with the License. You may obtain a copy of the License at
!> http://www.mozilla.org/MPL/
!>
!> Software distributed under the License is distributed on an "AS IS"
!> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
!> License for the specific language governing rights and limitations
!> under the License.
!>
!> The Original Code is OpenCMISS
!>
!> The Initial Developer of the Original Code is University of Auckland,
!> Auckland, New Zealand, the University of Oxford, Oxford, United
!> Kingdom and King's College, London, United Kingdom. Portions created
!> by the University of Auckland, the University of Oxford and King's
!> College, London are Copyright (C) 2007-2010 by the University of
!> Auckland, the University of Oxford and King's College, London.
!> All Rights Reserved.
!>
!> Contributor(s):
!>
!> Alternatively, the contents of this file may be used under the terms of
!> either the GNU General Public License Version 2 or later (the "GPL"), or
!> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
!> in which case the provisions of the GPL or the LGPL are applicable instead
!> of those above. If you wish to allow use of your version of this file only
!> under the terms of either the GPL or the LGPL, and not to allow others to
!> use your version of this file under the terms of the MPL, indicate your
!> decision by deleting the provisions above and replace them with the notice
!> and other provisions required by the GPL or the LGPL. If you do not delete
!> the provisions above, a recipient may use your version of this file under
!> the terms of any one of the MPL, the GPL or the LGPL.
!>

!> This module handles all problem wide constants.
MODULE PROBLEM_CONSTANTS
  
  USE KINDS

  IMPLICIT NONE

  !Module parameters

  !Problem Classes
  INTEGER(INTG), PARAMETER :: PROBLEM_NO_CLASS=0  
  INTEGER(INTG), PARAMETER :: PROBLEM_ELASTICITY_CLASS=1
  INTEGER(INTG), PARAMETER :: PROBLEM_FLUID_MECHANICS_CLASS=2
  INTEGER(INTG), PARAMETER :: PROBLEM_ELECTROMAGNETICS_CLASS=3
  INTEGER(INTG), PARAMETER :: PROBLEM_CLASSICAL_FIELD_CLASS=4  
  INTEGER(INTG), PARAMETER :: PROBLEM_BIOELECTRICS_CLASS=5
  INTEGER(INTG), PARAMETER :: PROBLEM_MODAL_CLASS=6
  INTEGER(INTG), PARAMETER :: PROBLEM_FITTING_CLASS=7
  INTEGER(INTG), PARAMETER :: PROBLEM_OPTIMISATION_CLASS=8
  INTEGER(INTG), PARAMETER :: PROBLEM_MULTI_PHYSICS_CLASS=9

  !Problem types
  INTEGER(INTG), PARAMETER :: PROBLEM_NO_TYPE=0
  !Elasticity class
  INTEGER(INTG), PARAMETER :: PROBLEM_LINEAR_ELASTICITY_TYPE=1
  INTEGER(INTG), PARAMETER :: PROBLEM_FINITE_ELASTICITY_TYPE=2
  !Fluid mechanics class
  INTEGER(INTG), PARAMETER :: PROBLEM_STOKES_EQUATION_TYPE=1
  INTEGER(INTG), PARAMETER :: PROBLEM_NAVIER_STOKES_EQUATION_TYPE=2
  INTEGER(INTG), PARAMETER :: PROBLEM_DARCY_EQUATION_TYPE=3
  INTEGER(INTG), PARAMETER :: PROBLEM_DARCY_PRESSURE_EQUATION_TYPE=4
  INTEGER(INTG), PARAMETER :: PROBLEM_POISEUILLE_EQUATION_TYPE=5
  INTEGER(INTG), PARAMETER :: PROBLEM_BURGERS_EQUATION_TYPE=6
  !Electromagnetics class
  INTEGER(INTG), PARAMETER :: PROBLEM_ELECTROSTATIC_TYPE=1
  INTEGER(INTG), PARAMETER :: PROBLEM_MAGNETOSTATIC_TYPE=2
  INTEGER(INTG), PARAMETER :: PROBLEM_MAXWELLS_EQUATIONS_TYPE=3
  !Classical field class
  INTEGER(INTG), PARAMETER :: PROBLEM_LAPLACE_EQUATION_TYPE=1
  INTEGER(INTG), PARAMETER :: PROBLEM_POISSON_EQUATION_TYPE=2
  INTEGER(INTG), PARAMETER :: PROBLEM_HELMHOLTZ_EQUATION_TYPE=3
  INTEGER(INTG), PARAMETER :: PROBLEM_WAVE_EQUATION_TYPE=4
  INTEGER(INTG), PARAMETER :: PROBLEM_DIFFUSION_EQUATION_TYPE=5
  INTEGER(INTG), PARAMETER :: PROBLEM_ADVECTION_DIFFUSION_EQUATION_TYPE=6
  INTEGER(INTG), PARAMETER :: PROBLEM_REACTION_DIFFUSION_EQUATION_TYPE=7
  INTEGER(INTG), PARAMETER :: PROBLEM_BIHARMONIC_EQUATION_TYPE=8
  INTEGER(INTG), PARAMETER :: PROBLEM_DATA_FITTING_TYPE=9
  INTEGER(INTG), PARAMETER :: PROBLEM_HJ_EQUATION_TYPE=10

  !Bioelectrics class
  INTEGER(INTG), PARAMETER :: PROBLEM_MONODOMAIN_EQUATION_TYPE=1
  INTEGER(INTG), PARAMETER :: PROBLEM_BIDOMAIN_EQUATION_TYPE=2
  INTEGER(INTG), PARAMETER :: PROBLEM_MONODOMAIN_STRANG_SPLITTING_EQUATION_TYPE=3
  !Modal class
  INTEGER(INTG), PARAMETER :: PROBLEM_LINEAR_ELASTIC_MODAL_TYPE=1
  !Multi physics class
  INTEGER(INTG), PARAMETER :: PROBLEM_FINITE_ELASTICITY_DARCY_TYPE=1
  INTEGER(INTG), PARAMETER :: PROBLEM_FINITE_ELASTICITY_STOKES_TYPE=2
  INTEGER(INTG), PARAMETER :: PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_TYPE=3
  INTEGER(INTG), PARAMETER :: PROBLEM_DIFFUSION_DIFFUSION_TYPE=4
  INTEGER(INTG), PARAMETER :: PROBLEM_DIFFUSION_ADVECTION_DIFFUSION_TYPE=5
  INTEGER(INTG), PARAMETER :: PROBLEM_MULTI_COMPARTMENT_TRANSPORT_TYPE=6 !<Problem type for the multi-compartment coupled transport, comprising either/or/both advection-diffusion & diffusion 
  INTEGER(INTG), PARAMETER :: PROBLEM_FINITE_ELASTICITY_FLUID_PRESSURE_TYPE=7
  INTEGER(INTG), PARAMETER :: PROBLEM_BIOELECTRIC_FINITE_ELASTICITY_TYPE=8
  !Problem subtypes
  INTEGER(INTG), PARAMETER :: PROBLEM_NO_SUBTYPE=0
  !Elasticity class
  !  Linear elasticity - uses PROBLEM_NO_SUBTYPE=0
  !  Finite elasticity
  INTEGER(INTG), PARAMETER :: PROBLEM_QUASISTATIC_FINITE_ELASTICITY_SUBTYPE=1
  INTEGER(INTG), PARAMETER :: PROBLEM_FINITE_ELASTICITY_CELLML_SUBTYPE=2

  !Fluid mechanics class
  !  Stokes equations
  INTEGER(INTG), PARAMETER :: PROBLEM_STATIC_STOKES_SUBTYPE=1
  INTEGER(INTG), PARAMETER :: PROBLEM_LAPLACE_STOKES_SUBTYPE=2
  INTEGER(INTG), PARAMETER :: PROBLEM_TRANSIENT_STOKES_SUBTYPE=3
  INTEGER(INTG), PARAMETER :: PROBLEM_OPTIMISED_STOKES_SUBTYPE=4
  INTEGER(INTG), PARAMETER :: PROBLEM_ALE_STOKES_SUBTYPE=5
  INTEGER(INTG), PARAMETER :: PROBLEM_PGM_STOKES_SUBTYPE=6
!  INTEGER(INTG), PARAMETER :: PROBLEM_HJ_SUBTYPE=6
  !  Navier-Stokes equations
  INTEGER(INTG), PARAMETER :: PROBLEM_STATIC_NAVIER_STOKES_SUBTYPE=1
  INTEGER(INTG), PARAMETER :: PROBLEM_LAPLACE_NAVIER_STOKES_SUBTYPE=2
  INTEGER(INTG), PARAMETER :: PROBLEM_TRANSIENT_NAVIER_STOKES_SUBTYPE=3
  INTEGER(INTG), PARAMETER :: PROBLEM_OPTIMISED_NAVIER_STOKES_SUBTYPE=4
  INTEGER(INTG), PARAMETER :: PROBLEM_ALE_NAVIER_STOKES_SUBTYPE=5
  INTEGER(INTG), PARAMETER :: PROBLEM_PGM_NAVIER_STOKES_SUBTYPE=6
  INTEGER(INTG), PARAMETER :: PROBLEM_1DTRANSIENT_NAVIER_STOKES_SUBTYPE=8
  INTEGER(INTG), PARAMETER :: PROBLEM_QUASISTATIC_NAVIER_STOKES_SUBTYPE=7
  !  Darcy equation
  INTEGER(INTG), PARAMETER :: PROBLEM_STANDARD_DARCY_SUBTYPE=1
  INTEGER(INTG), PARAMETER :: PROBLEM_QUASISTATIC_DARCY_SUBTYPE=2
  INTEGER(INTG), PARAMETER :: PROBLEM_ALE_DARCY_SUBTYPE=3
  INTEGER(INTG), PARAMETER :: PROBLEM_TRANSIENT_DARCY_SUBTYPE=4
  INTEGER(INTG), PARAMETER :: PROBLEM_PGM_DARCY_SUBTYPE=5
  INTEGER(INTG), PARAMETER :: PROBLEM_PGM_TRANSIENT_DARCY_SUBTYPE=6
  !  Poiseuille equation
  INTEGER(INTG), PARAMETER :: PROBLEM_STATIC_POISEUILLE_SUBTYPE=1
  INTEGER(INTG), PARAMETER :: PROBLEM_DYNAMIC_POISEUILLE_SUBTYPE=2
  !  Burgers equation
  INTEGER(INTG), PARAMETER :: PROBLEM_STATIC_BURGERS_SUBTYPE=1
  INTEGER(INTG), PARAMETER :: PROBLEM_DYNAMIC_BURGERS_SUBTYPE=2
  !Electromagnetics class
  !Classical field class
  !  Laplace equation
!!TODO: We don't really have two problem types here? Maybe a different type with nonlinear boundary conditions???
  INTEGER(INTG), PARAMETER :: PROBLEM_STANDARD_LAPLACE_SUBTYPE=1
  INTEGER(INTG), PARAMETER :: PROBLEM_GENERALISED_LAPLACE_SUBTYPE=2
  INTEGER(INTG), PARAMETER :: PROBLEM_MOVING_MESH_LAPLACE_SUBTYPE=3
  !  Hamilton-Jacobi equation
  INTEGER(INTG), PARAMETER :: PROBLEM_STANDARD_HJ_SUBTYPE=1
  INTEGER(INTG), PARAMETER :: PROBLEM_GENERALISED_HJ_SUBTYPE=2
  INTEGER(INTG), PARAMETER :: PROBLEM_MOVING_MESH_HJ_SUBTYPE=3
  !  Poisson equation
  INTEGER(INTG), PARAMETER :: PROBLEM_LINEAR_SOURCE_POISSON_SUBTYPE=1
  INTEGER(INTG), PARAMETER :: PROBLEM_NONLINEAR_SOURCE_POISSON_SUBTYPE=2
  INTEGER(INTG), PARAMETER :: PROBLEM_LINEAR_PRESSURE_POISSON_SUBTYPE=3
  INTEGER(INTG), PARAMETER :: PROBLEM_NONLINEAR_PRESSURE_POISSON_SUBTYPE=4
  INTEGER(INTG), PARAMETER :: PROBLEM_ALE_PRESSURE_POISSON_SUBTYPE=5
  INTEGER(INTG), PARAMETER :: PROBLEM_FITTED_PRESSURE_POISSON_SUBTYPE=6
  !  Helmholtz equation
  INTEGER(INTG), PARAMETER :: PROBLEM_STANDARD_HELMHOLTZ_SUBTYPE=2
  INTEGER(INTG), PARAMETER :: PROBLEM_GENERALISED_HELMHOLTZ_SUBTYPE=3
  !  Wave equation
  !  Diffusion equation
  INTEGER(INTG), PARAMETER :: PROBLEM_NO_SOURCE_DIFFUSION_SUBTYPE=1
  INTEGER(INTG), PARAMETER :: PROBLEM_LINEAR_SOURCE_DIFFUSION_SUBTYPE=2
  INTEGER(INTG), PARAMETER :: PROBLEM_NONLINEAR_SOURCE_DIFFUSION_SUBTYPE=3
  INTEGER(INTG), PARAMETER :: PROBLEM_NO_SOURCE_ALE_DIFFUSION_SUBTYPE=4
  INTEGER(INTG), PARAMETER :: PROBLEM_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE=5
  INTEGER(INTG), PARAMETER :: PROBLEM_NONLINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE=6
  !  Reaction-diffusion equation
  INTEGER(INTG), PARAMETER :: PROBLEM_CELLML_REAC_INTEG_REAC_DIFF_STRANG_SPLIT_SUBTYPE=1
  INTEGER(INTG), PARAMETER :: PROBLEM_CELLML_REAC_EVAL_REAC_DIFF_NO_SPLIT_SUBTYPE=2
  INTEGER(INTG), PARAMETER :: PROBLEM_CONSTANT_REAC_DIFF_NO_SPLIT_SUBTYPE=3
  !  Advection-diffusion equation
  INTEGER(INTG), PARAMETER :: PROBLEM_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE=1
  INTEGER(INTG), PARAMETER :: PROBLEM_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE=2
  INTEGER(INTG), PARAMETER :: PROBLEM_NONLINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE=3
  INTEGER(INTG), PARAMETER :: PROBLEM_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE=4
  INTEGER(INTG), PARAMETER :: PROBLEM_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE=5
  INTEGER(INTG), PARAMETER :: PROBLEM_NONLINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE=6
  !Subtypes for steady-state advection-diffusion equation
  INTEGER(INTG), PARAMETER :: PROBLEM_NO_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE=7
  INTEGER(INTG), PARAMETER :: PROBLEM_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE=8
  INTEGER(INTG), PARAMETER :: PROBLEM_NONLINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE=9
  

  !Bioelectric class
  !  Monodomain equation
  INTEGER(INTG), PARAMETER :: PROBLEM_MONODOMAIN_GUDUNOV_SPLIT_SUBTYPE=1
  INTEGER(INTG), PARAMETER :: PROBLEM_MONODOMAIN_STRANG_SPLIT_SUBTYPE=2   
  !  Bidomain equation
  INTEGER(INTG), PARAMETER :: PROBLEM_BIDOMAIN_GUDUNOV_SPLIT_SUBTYPE=1
  INTEGER(INTG), PARAMETER :: PROBLEM_BIDOMAIN_STRANG_SPLIT_SUBTYPE=2

  !Monodomain Strang splitting type : specific cell models
  INTEGER(INTG), PARAMETER :: PROBLEM_MONODOMAIN_BUENOOROVIO_SUBTYPE=1
  INTEGER(INTG), PARAMETER :: PROBLEM_MONODOMAIN_TENTUSSCHER06_SUBTYPE=2

  !Modal class
  !  Galerkin projection
  INTEGER(INTG), PARAMETER :: PROBLEM_STANDARD_DATA_FITTING_SUBTYPE=1
  INTEGER(INTG), PARAMETER :: PROBLEM_GENERALISED_DATA_FITTING_SUBTYPE=2
  INTEGER(INTG), PARAMETER :: PROBLEM_MAT_PROPERTIES_DATA_FITTING_SUBTYPE=3
  INTEGER(INTG), PARAMETER :: PROBLEM_VECTOR_DATA_FITTING_SUBTYPE=4
  INTEGER(INTG), PARAMETER :: PROBLEM_DIV_FREE_VECTOR_DATA_FITTING_SUBTYPE=5
  INTEGER(INTG), PARAMETER :: PROBLEM_VECTOR_DATA_PRE_FITTING_SUBTYPE=6
  INTEGER(INTG), PARAMETER :: PROBLEM_DIV_FREE_VECTOR_DATA_PRE_FITTING_SUBTYPE=7
  !  Multi physics (subtype numbers must be different from Darcy ones)
  INTEGER(INTG), PARAMETER :: PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE=101
  INTEGER(INTG), PARAMETER :: PROBLEM_PGM_ELASTICITY_DARCY_SUBTYPE=102
  INTEGER(INTG), PARAMETER :: PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE=103
  INTEGER(INTG), PARAMETER :: PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE=104
  INTEGER(INTG), PARAMETER :: PROBLEM_COUPLED_SOURCE_DIFFUSION_DIFFUSION_SUBTYPE=111
  INTEGER(INTG), PARAMETER :: PROBLEM_COUPLED_SOURCE_ALE_DIFFUSION_DIFFUSION_SUBTYPE=112
  INTEGER(INTG), PARAMETER :: PROBLEM_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE=121
  INTEGER(INTG), PARAMETER :: PROBLEM_COUPLED_SOURCE_ALE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE=122
  INTEGER(INTG), PARAMETER :: PROBLEM_STANDARD_MULTI_COMPARTMENT_TRANSPORT_SUBTYPE=131
  INTEGER(INTG), PARAMETER :: PROBLEM_STANDARD_ALE_MULTI_COMPARTMENT_TRANSPORT_SUBTYPE=132    
  INTEGER(INTG), PARAMETER :: PROBLEM_STANDARD_ELASTICITY_FLUID_PRESSURE_SUBTYPE=133
  INTEGER(INTG), PARAMETER :: PROBLEM_GUDUNOV_MONODOMAIN_SIMPLE_ELASTICITY_SUBTYPE=141

  !> \addtogroup PROBLEM_CONSTANTS_SetupTypes PROBLEM_CONSTANTS::SetupTypes
  !> \brief Setup type parameters
  !> \see PROBLEM_CONSTANTS
  !>@{
  INTEGER(INTG), PARAMETER :: PROBLEM_SETUP_INITIAL_TYPE=1 !<Initial setup for a problem. \see PROBLEM_CONSTANTS_SetupTypes,PROBLEM_CONSTANTS
  INTEGER(INTG), PARAMETER :: PROBLEM_SETUP_CONTROL_TYPE=2 !<Solver setup for a problem. \see PROBLEM_CONSTANTS_SetupTypes,PROBLEM_ROU
  INTEGER(INTG), PARAMETER :: PROBLEM_SETUP_SOLVERS_TYPE=3 !<Solver setup for a problem. \see PROBLEM_CONSTANTS_SetupTypes,PROBLEM_CONSTANTS
  INTEGER(INTG), PARAMETER :: PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE=4 !<Solver equations setup for a problem. \see PROBLEM_CONSTANTS_SetupTypes,PROBLEM_CONSTANTS
  INTEGER(INTG), PARAMETER :: PROBLEM_SETUP_CELLML_EQUATIONS_TYPE=5 !<CellML equations setup for a problem. \see PROBLEM_CONSTANTS_SetupTypes,PROBLEM_CONSTANTS
  !>@}
  
  !> \addtogroup PROBLEM_CONSTANTS_SetupActionTypes PROBLEM_CONSTANTS::SetupActionTypes
  !> \brief Setup action type parameters
  !> \see PROBLEM_CONSTANTS
  !>@{
  INTEGER(INTG), PARAMETER :: PROBLEM_SETUP_START_ACTION=1 !<Start setup action. \see PROBLEM_CONSTANTS_SetupActionTypes,CONSTANTS_ROUTINES
  INTEGER(INTG), PARAMETER :: PROBLEM_SETUP_FINISH_ACTION=2 !<Finish setup action. \see PROBLEM_CONSTANTS_SetupActionTypes,CONSTANTS_ROUTINES
  !>@}

  !> \addtogroup PROBLEM_CONSTANTS_ControlLoopTypes PROBLEM_CONSTANTS::ControlLoopTypes
  !> \brief Control loop type parameters
  !> \see PROBLEM_CONSTANTS
  !>@{
  INTEGER(INTG), PARAMETER :: PROBLEM_CONTROL_SIMPLE_TYPE=1 !<Simple, one iteration control loop. \see PROBLEM_CONSTANTS_ControlLoopTypes,PROBLEM_CONSTANTS
  INTEGER(INTG), PARAMETER :: PROBLEM_CONTROL_FIXED_LOOP_TYPE=2 !<Fixed iteration control loop. \see PROBLEM_CONSTANTS_ControlLoopTypes,PROBLEM_CONSTANTS
  INTEGER(INTG), PARAMETER :: PROBLEM_CONTROL_TIME_LOOP_TYPE=3 !<Time control loop. \see PROBLEM_CONSTANTS_ControlLoopTypes,PROBLEM_CONSTANTS
  INTEGER(INTG), PARAMETER :: PROBLEM_CONTROL_WHILE_LOOP_TYPE=4 !<While control loop. \see PROBLEM_CONSTANTS_ControlLoopTypes,PROBLEM_CONSTANTS
  INTEGER(INTG), PARAMETER :: PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE=5 !<Load increment control loop. \see PROBLEM_CONSTANTS_ControlLoopTypes,PROBLEM_CONSTANTS
  !>@}
   
  !> \addtogroup PROBLEM_CONSTANTS_LinearityTypes PROBLEM_CONSTANTS::LinearityTypes
  !> \brief Setup type parameters
  !> \see PROBLEM_CONSTANTS
  !>@{
  INTEGER(INTG), PARAMETER :: PROBLEM_SOLVER_LINEAR=1 !<Linear problem. \see PROBLEM_CONSTANTS_LinearityTypes,PROBLEM_CONSTANTS
  INTEGER(INTG), PARAMETER :: PROBLEM_SOLVER_NONLINEAR=2 !<Nonlinear problem. \see PROBLEM_CONSTANTS_LinearityTypes,PROBLEM_CONSTANTS
  !>@}
  
  
  !> \addtogroup PROBLEM_CONSTANTS_EquationsLinearityTypes PROBLEM_CONSTANTS::EquationsLinearityTypes
  !> \brief The solver equations linearity types 
  !> \see PROBLEM_CONSTANTS
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_EQUATIONS_LINEAR=1 !<Solver equations are linear \see PROBLEM_CONSTANTS_EquationLinearityTypes,PROBLEM_CONSTANTS
  INTEGER(INTG), PARAMETER :: SOLVER_EQUATIONS_NONLINEAR=2 !<Solver equations are nonlinear \see PROBLEM_CONSTANTS_EquationLinearityTypes,PROBLEM_CONSTANTS
  !>@}

  !> \addtogroup PROBLEM_CONSTANTS_EquationsTimeDependenceTypes PROBLEM_CONSTANTS::EquationsTimeDependenceTypes
  !> \brief The solver equations time dependence types 
  !> \see PROBLEM_CONSTANTS
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_EQUATIONS_STATIC=1 !<Solver equations are static \see PROBLEM_CONSTANTS_EquationTimeDependenceTypes,PROBLEM_CONSTANTS
  INTEGER(INTG), PARAMETER :: SOLVER_EQUATIONS_QUASISTATIC=2 !<Solver equations are quasistatic \see PROBLEM_CONSTANTS_EquationTimeDependenceTypes,PROBLEM_CONSTANTS
  INTEGER(INTG), PARAMETER :: SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC=3 !<Solver equations are first order dynamic \see PROBLEM_CONSTANTS_EquationTimeDependenceTypes,PROBLEM_CONSTANTS
  INTEGER(INTG), PARAMETER :: SOLVER_EQUATIONS_SECOND_ORDER_DYNAMIC=4 !<Solver equations are second order dynamic \see PROBLEM_CONSTANTS_EquationTimeDependenceTypes,PROBLEM_CONSTANTS
  !>@}

   
END MODULE PROBLEM_CONSTANTS
