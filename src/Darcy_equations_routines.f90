!> \file
!> \author Christian Michler
!> \brief This module handles all Darcy equations routines.
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
!> Contributor(s): Chris Bradley
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

!>This module handles all Darcy equations routines.

MODULE DARCY_EQUATIONS_ROUTINES

  USE BaseRoutines
  USE BASIS_ROUTINES
  USE BOUNDARY_CONDITIONS_ROUTINES
  USE Constants
  USE CONTROL_LOOP_ROUTINES
  USE ControlLoopAccessRoutines
  USE ComputationEnvironment
  USE COORDINATE_ROUTINES
  USE DISTRIBUTED_MATRIX_VECTOR
  USE DOMAIN_MAPPINGS
  USE EquationsRoutines
  USE EquationsAccessRoutines
  USE EquationsMappingRoutines
  USE EquationsMatricesRoutines
  USE EQUATIONS_SET_CONSTANTS
  USE EquationsSetAccessRoutines
  USE FIELD_ROUTINES
  USE FieldAccessRoutines
  USE FINITE_ELASTICITY_ROUTINES
  USE FLUID_MECHANICS_IO_ROUTINES
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE Kinds
  USE Maths
  USE MATRIX_VECTOR
  USE MESH_ROUTINES
  USE NODE_ROUTINES
  USE PROBLEM_CONSTANTS
  USE Strings
  USE SOLVER_ROUTINES
  USE SolverAccessRoutines
  USE Timer
  USE Types

#include "macros.h"


  IMPLICIT NONE

  PUBLIC DARCY_EQUATION_EQUATIONS_SET_SETUP
  PUBLIC Darcy_EquationsSetSpecificationSet
  PUBLIC Darcy_EquationsSetSolutionMethodSet
  PUBLIC Darcy_BoundaryConditionsAnalyticCalculate

  PUBLIC DARCY_EQUATION_PROBLEM_SETUP
  PUBLIC Darcy_ProblemSpecificationSet

  PUBLIC DARCY_EQUATION_FINITE_ELEMENT_CALCULATE

  PUBLIC DARCY_EQUATION_PRE_SOLVE
  PUBLIC DARCY_EQUATION_POST_SOLVE
  PUBLIC DARCY_EQUATION_POST_SOLVE_OUTPUT_DATA

  PUBLIC DARCY_CONTROL_TIME_LOOP_PRE_LOOP

  PUBLIC Darcy_PreSolveStorePreviousIterate

  PUBLIC DARCY_EQUATION_MONITOR_CONVERGENCE

  INTEGER(INTG) :: SOLVER_NUMBER_SOLID,SOLVER_NUMBER_MAT_PROPERTIES,SOLVER_NUMBER_DARCY
  INTEGER(INTG) :: SOLVER_INDEX_SOLID,SOLVER_INDEX_MAT_PROPERTIES,SOLVER_INDEX_DARCY

  REAL(DP) :: RESIDUAL_NORM_0

  LOGICAL :: idebug1, idebug2, idebug3


CONTAINS

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a Darcy equation type of a fluid mechanics equations set class.
  SUBROUTINE Darcy_EquationsSetSolutionMethodSet(EQUATIONS_SET,SOLUTION_METHOD,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equations_SET !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: SOLUTION_METHOD !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Darcy_EquationsSetSolutionMethodSet",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a Darcy type equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
      CASE(EQUATIONS_SET_STANDARD_DARCY_SUBTYPE, EQUATIONS_SET_QUASISTATIC_DARCY_SUBTYPE, EQUATIONS_SET_ALE_DARCY_SUBTYPE, &
        & EQUATIONS_SET_TRANSIENT_DARCY_SUBTYPE,EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE, &
        & EQUATIONS_SET_TRANSIENT_ALE_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
        & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE,EQUATIONS_SET_MULTI_COMPARTMENT_DARCY_SUBTYPE, &
        & EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
        SELECT CASE(SOLUTION_METHOD)
        CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
          EQUATIONS_SET%SOLUTION_METHOD=EQUATIONS_SET_FEM_SOLUTION_METHOD
        CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE DEFAULT
          localError="The specified solution method of "//TRIM(NumberToVString(SOLUTION_METHOD,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="Equations set subtype of "//TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
          & " is not valid for a Darcy equation type of a fluid mechanics equations set class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF

    EXITS("Darcy_EquationsSetSolutionMethodSet")
    RETURN
999 ERRORSEXITS("Darcy_EquationsSetSolutionMethodSet",err,error)
    RETURN 1
  END SUBROUTINE Darcy_EquationsSetSolutionMethodSet

  !
  !================================================================================================================================
  !

  !>Sets up the Darcy equation.
  SUBROUTINE DARCY_EQUATION_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equations_SET !<A pointer to the equations set to setup
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: equations_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: GEOMETRIC_SCALING_TYPE, GEOMETRIC_MESH_COMPONENT,NUMBER_OF_COMPONENTS,NUMBER_OF_DARCY_COMPONENTS
    INTEGER(INTG) :: imy_matrix,Ncompartments
    TYPE(DECOMPOSITION_TYPE), POINTER :: GEOMETRIC_DECOMPOSITION
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EQUATIONS_SET_MATERIALS_TYPE), POINTER :: equations_MATERIALS
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EQUATIONS_SET_EQUATIONS_SET_FIELD_TYPE), POINTER :: equations_EQUATIONS_SET_FIELD
    TYPE(FIELD_TYPE), POINTER :: equations_SET_FIELD_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(EQUATIONS_SET_SOURCE_TYPE), POINTER :: equations_SOURCE
    TYPE(VARYING_STRING) :: localError
    INTEGER:: DEPENDENT_FIELD_NUMBER_OF_VARIABLES, DEPENDENT_FIELD_NUMBER_OF_COMPONENTS
    INTEGER:: DEPENDENT_FIELD_ELASTICITY_NUMBER_OF_COMPONENTS, DEPENDENT_FIELD_DARCY_NUMBER_OF_COMPONENTS
    INTEGER:: INDEPENDENT_FIELD_NUMBER_OF_VARIABLES, INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS
    INTEGER:: NUMBER_OF_DIMENSIONS, GEOMETRIC_COMPONENT_NUMBER
    INTEGER:: MATERIAL_FIELD_NUMBER_OF_VARIABLES, MATERIAL_FIELD_NUMBER_OF_COMPONENTS
    INTEGER:: MESH_COMPONENT,MATERIAL_FIELD_NUMBER_OF_U_VAR_COMPONENTS,MATERIAL_FIELD_NUMBER_OF_V_VAR_COMPONENTS, &
            & MATERIAL_FIELD_NUMBER_OF_U1_VAR_COMPONENTS
    INTEGER:: i,component_idx

    INTEGER(INTG) :: num_var,num_var_count
    INTEGER(INTG) :: equations_SET_FIELD_NUMBER_OF_VARIABLES,EQUATIONS_SET_FIELD_NUMBER_OF_COMPONENTS,NUMBER_OF_SOURCE_COMPONENTS
    INTEGER(INTG), POINTER :: equations_SET_FIELD_DATA(:)
    INTEGER(INTG), ALLOCATABLE :: VARIABLE_TYPES(:),VARIABLE_U_TYPES(:),COUPLING_MATRIX_STORAGE_TYPE(:), &
      & COUPLING_MATRIX_STRUCTURE_TYPE(:)

    ENTERS("DARCY_EQUATION_EQUATIONS_SET_SETUP",err,error,*999)

    NULLIFY(EQUATIONS)
    NULLIFY(vectorMapping)
    NULLIFY(vectorMatrices)
    NULLIFY(GEOMETRIC_DECOMPOSITION)
    NULLIFY(EQUATIONS_EQUATIONS_SET_FIELD)
    NULLIFY(EQUATIONS_SET_FIELD_FIELD)
    NULLIFY(EQUATIONS_SET_FIELD_DATA)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a Darcy type equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
      CASE(EQUATIONS_SET_STANDARD_DARCY_SUBTYPE, EQUATIONS_SET_QUASISTATIC_DARCY_SUBTYPE, EQUATIONS_SET_ALE_DARCY_SUBTYPE, &
        & EQUATIONS_SET_TRANSIENT_DARCY_SUBTYPE, EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE, &
        & EQUATIONS_SET_TRANSIENT_ALE_DARCY_SUBTYPE, EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
        & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE,EQUATIONS_SET_MULTI_COMPARTMENT_DARCY_SUBTYPE, &
        & EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
        SELECT CASE(EQUATIONS_SET_SETUP%SETUP_TYPE)

        !-----------------------------------------------------------------
        ! s o l u t i o n   m e t h o d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
          SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
          CASE(EQUATIONS_SET_STANDARD_DARCY_SUBTYPE, EQUATIONS_SET_QUASISTATIC_DARCY_SUBTYPE, EQUATIONS_SET_ALE_DARCY_SUBTYPE, &
            & EQUATIONS_SET_TRANSIENT_DARCY_SUBTYPE,EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE, &
            & EQUATIONS_SET_TRANSIENT_ALE_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
            & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE)
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
             CALL Darcy_EquationsSetSolutionMethodSet(EQUATIONS_SET,EQUATIONS_SET_FEM_SOLUTION_METHOD,err,error,*999)
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
               !do nothing
            CASE DEFAULT
              localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
                & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
                & " is invalid for a standard or quasistatic Darcy equation."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE(EQUATIONS_SET_MULTI_COMPARTMENT_DARCY_SUBTYPE,EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              CALL Darcy_EquationsSetSolutionMethodSet(EQUATIONS_SET,EQUATIONS_SET_FEM_SOLUTION_METHOD,err,error,*999)

              EQUATIONS_SET_FIELD_NUMBER_OF_VARIABLES = 1
              EQUATIONS_SET_FIELD_NUMBER_OF_COMPONENTS = 2
              EQUATIONS_EQUATIONS_SET_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD
              IF(EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_AUTO_CREATED) THEN
                !Create the auto created equations set field
                CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION, &
                  & EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,err,error,*999)
                EQUATIONS_SET_FIELD_FIELD=>EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD
                CALL FIELD_LABEL_SET(EQUATIONS_SET_FIELD_FIELD,"Equations Set Field",err,error,*999)
                CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET_FIELD_FIELD,FIELD_GENERAL_TYPE,&
                  & err,error,*999)
                CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET_FIELD_FIELD,&
                  & FIELD_INDEPENDENT_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_SET(EQUATIONS_SET_FIELD_FIELD, &
                  & EQUATIONS_SET_FIELD_NUMBER_OF_VARIABLES,err,error,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET_FIELD_FIELD,&
                  & [FIELD_U_VARIABLE_TYPE],err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET_FIELD_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET_FIELD_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_INTG_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET_FIELD_FIELD,&
                  & FIELD_U_VARIABLE_TYPE,EQUATIONS_SET_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
              ELSE
                !Check the user specified field
                CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,EQUATIONS_SET_FIELD_NUMBER_OF_VARIABLES, &
                  & err,error,*999)
                CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE],err,error,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_INTG_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                  & EQUATIONS_SET_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
              ENDIF
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              IF(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_AUTO_CREATED) THEN
                CALL FIELD_CREATE_FINISH(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,err,error,*999)
                CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,&
                  & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, 1, 1_INTG, ERR, ERROR, *999)
                CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,&
                  & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, 2, 1_INTG, ERR, ERROR, *999)
              ENDIF
!!TODO: Check valid setup
            CASE DEFAULT
              localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
                & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
                & " is invalid for a standard or quasistatic Darcy equation."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          END SELECT

        !-----------------------------------------------------------------
        ! g e o m e t r y   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
          SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
          CASE(EQUATIONS_SET_STANDARD_DARCY_SUBTYPE, EQUATIONS_SET_QUASISTATIC_DARCY_SUBTYPE, &
            & EQUATIONS_SET_TRANSIENT_DARCY_SUBTYPE)
          !Do nothing
          CASE(EQUATIONS_SET_ALE_DARCY_SUBTYPE,EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE, &
            & EQUATIONS_SET_TRANSIENT_ALE_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
            & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE,EQUATIONS_SET_MULTI_COMPARTMENT_DARCY_SUBTYPE, &
            & EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)

            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)

              FIELD_VARIABLE=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%ptr

              CALL Field_ParameterSetEnsureCreated(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD, FIELD_U_VARIABLE_TYPE, &
                & FIELD_INITIAL_VALUES_SET_TYPE, ERR, ERROR, *999)

              CALL Field_ParameterSetEnsureCreated(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD, FIELD_U_VARIABLE_TYPE, &
                & FIELD_PREVIOUS_VALUES_SET_TYPE, ERR, ERROR, *999)

              CALL Field_ParameterSetEnsureCreated(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD, FIELD_U_VARIABLE_TYPE, &
                & FIELD_MESH_DISPLACEMENT_SET_TYPE, ERR, ERROR, *999)

              CALL Field_ParameterSetEnsureCreated(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD, FIELD_U_VARIABLE_TYPE, &
                & FIELD_MESH_VELOCITY_SET_TYPE, ERR, ERROR, *999)

              CALL Field_ParameterSetEnsureCreated(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD, FIELD_U_VARIABLE_TYPE, &
                & FIELD_NEGATIVE_MESH_VELOCITY_SET_TYPE, ERR, ERROR, *999)

              IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_MULTI_COMPARTMENT_DARCY_SUBTYPE .OR. &
                 EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE) THEN
                !Create the equations set field for multi-compartment Darcy
                EQUATIONS_SET_FIELD_NUMBER_OF_COMPONENTS = 2

                EQUATIONS_EQUATIONS_SET_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD
                EQUATIONS_SET_FIELD_FIELD=>EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD

                IF(EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_AUTO_CREATED) THEN
                  CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,err,error,*999)
                  CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET_FIELD_FIELD,&
                    & GEOMETRIC_DECOMPOSITION,err,error,*999)
                  CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET_FIELD_FIELD,&
                    & EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & 1,GEOMETRIC_COMPONENT_NUMBER,err,error,*999)
                  DO component_idx = 1, EQUATIONS_SET_FIELD_NUMBER_OF_COMPONENTS
                    CALL FIELD_COMPONENT_MESH_COMPONENT_SET_AND_LOCK(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,component_idx,GEOMETRIC_COMPONENT_NUMBER,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,component_idx,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                  END DO

                  !Default the field scaling to that of the geometric field
                  CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                  CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,GEOMETRIC_SCALING_TYPE, &
                    & err,error,*999)
                ELSE
                  !Do nothing
                ENDIF
              ENDIF
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              ! do nothing
            CASE DEFAULT
              localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
                & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
                & " is invalid for a linear diffusion equation."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          END SELECT

        !-----------------------------------------------------------------
        ! d e p e n d e n t   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)
          SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
          CASE(EQUATIONS_SET_STANDARD_DARCY_SUBTYPE, EQUATIONS_SET_QUASISTATIC_DARCY_SUBTYPE, EQUATIONS_SET_ALE_DARCY_SUBTYPE, &
            & EQUATIONS_SET_TRANSIENT_DARCY_SUBTYPE,EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE, &
            & EQUATIONS_SET_TRANSIENT_ALE_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
            & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE, &
            & EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
                !Create the auto created dependent field
                CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION, &
                  & EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,err,error,*999)
                CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_GEOMETRIC_GENERAL_TYPE,err,error,*999)
                CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DEPENDENT_TYPE,err,error,*999)

                CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,err,error,*999)
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,GEOMETRIC_DECOMPOSITION, &
                  & err,error,*999)
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,err,error,*999)

                DEPENDENT_FIELD_NUMBER_OF_VARIABLES = 2  ! U and the normal component of its flux
                CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & DEPENDENT_FIELD_NUMBER_OF_VARIABLES,err,error,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,[FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DELUDELN_VARIABLE_TYPE],err,error,*999)
                CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,"U",err,error,*999)
                CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,"del U/del n", &
                    & err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)

                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD, FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS, ERR, ERROR, *999)
                DEPENDENT_FIELD_NUMBER_OF_COMPONENTS = NUMBER_OF_DIMENSIONS + 1
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, FIELD_U_VARIABLE_TYPE, &
                  & DEPENDENT_FIELD_NUMBER_OF_COMPONENTS, ERR, ERROR, *999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_DELUDELN_VARIABLE_TYPE, DEPENDENT_FIELD_NUMBER_OF_COMPONENTS, ERR, ERROR, *999)
                !Default to the geometric interpolation setup
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                  & GEOMETRIC_MESH_COMPONENT,err,error,*999)
                DO i=1,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS
                  IF( i < DEPENDENT_FIELD_NUMBER_OF_COMPONENTS ) THEN
                    !Set velocity mesh component (default to the geometric one)
                    MESH_COMPONENT = GEOMETRIC_MESH_COMPONENT
                    CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & i, MESH_COMPONENT,err,error,*999)
                    CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                      & i, MESH_COMPONENT,err,error,*999)
                  ELSE
                    !Set pressure mesh component (default to the geometric one)
                    MESH_COMPONENT = GEOMETRIC_MESH_COMPONENT
                    CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & i, MESH_COMPONENT,err,error,*999)
                    CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                      & i, MESH_COMPONENT,err,error,*999)
                  ENDIF
                ENDDO

                SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  DO i = 1, DEPENDENT_FIELD_NUMBER_OF_COMPONENTS
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,i,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_DELUDELN_VARIABLE_TYPE,i,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  ENDDO
                  !Default the scaling to the geometric field scaling
                  CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                  CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE DEFAULT
                  localError="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                    & " is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              ELSE
                !-----------------------------------
                ! DEPENDENT_FIELD: not AUTO_CREATED
                !-----------------------------------
                SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
                CASE(EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE, &
                  & EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
                  & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE)
                  !-----------------------------------------------------------------------
                  ! Check the shared dependent field set up in finite elasticity routines
                  !-----------------------------------------------------------------------
                  CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GEOMETRIC_GENERAL_TYPE,err,error,*999)
                  CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DEPENDENT_TYPE,err,error,*999)
                  CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,4,err,error,*999)
                  CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE, &
                    & FIELD_V_VARIABLE_TYPE,FIELD_DELVDELN_VARIABLE_TYPE],err,error,*999)
                  CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                    & err,error,*999)
                  CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                    & err,error,*999)
                  CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                    & err,error,*999)
                  CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELVDELN_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                    & err,error,*999)
                  CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                  CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                  CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                  CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELVDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD, FIELD_U_VARIABLE_TYPE, &
                    & NUMBER_OF_DIMENSIONS, ERR, ERROR, *999)

                  SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
                  CASE(EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE)  !compressible elasticity
                    DEPENDENT_FIELD_ELASTICITY_NUMBER_OF_COMPONENTS = NUMBER_OF_DIMENSIONS
                    DEPENDENT_FIELD_DARCY_NUMBER_OF_COMPONENTS = NUMBER_OF_DIMENSIONS + 2  !(u,v,w,p,m)
                  CASE(EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE)
                    DEPENDENT_FIELD_ELASTICITY_NUMBER_OF_COMPONENTS = NUMBER_OF_DIMENSIONS + 1
                    DEPENDENT_FIELD_DARCY_NUMBER_OF_COMPONENTS = NUMBER_OF_DIMENSIONS + 1
                  CASE(EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE)
                    DEPENDENT_FIELD_ELASTICITY_NUMBER_OF_COMPONENTS = NUMBER_OF_DIMENSIONS + 1 !(u1,u2,u3,p)
                    DEPENDENT_FIELD_DARCY_NUMBER_OF_COMPONENTS = NUMBER_OF_DIMENSIONS + 1 !(u,v,w,m)
                  END SELECT

                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                    & DEPENDENT_FIELD_ELASTICITY_NUMBER_OF_COMPONENTS,err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                    & DEPENDENT_FIELD_ELASTICITY_NUMBER_OF_COMPONENTS,err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE, &
                    & DEPENDENT_FIELD_DARCY_NUMBER_OF_COMPONENTS,err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELVDELN_VARIABLE_TYPE, &
                    & DEPENDENT_FIELD_DARCY_NUMBER_OF_COMPONENTS,err,error,*999)

                  SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                  CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1, &
                      & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                      !Mind that elastic hydrostatic pressure might be interpolated element-wise
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                      & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,1, &
                      & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELVDELN_VARIABLE_TYPE,1, &
                      & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                    CALL FlagError("Not implemented.",err,error,*999)
                  CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                    CALL FlagError("Not implemented.",err,error,*999)
                  CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                    CALL FlagError("Not implemented.",err,error,*999)
                  CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                    CALL FlagError("Not implemented.",err,error,*999)
                  CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                    CALL FlagError("Not implemented.",err,error,*999)
                  CASE DEFAULT
                    localError="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                      & " is invalid."
                    CALL FlagError(localError,err,error,*999)
                  END SELECT
                CASE(EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
                  !-----------------------------------------------------------------------
                  ! Check the shared dependent field set up in finite elasticity routines
                  ! Must have 2+2*Ncompartments number of variable types
                  !-----------------------------------------------------------------------
                  CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GEOMETRIC_GENERAL_TYPE,err,error,*999)
                  CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DEPENDENT_TYPE,err,error,*999)
                  !Get the number of Darcy compartments from the equations set field
                    EQUATIONS_SET_FIELD_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD
                    CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET_FIELD_FIELD,FIELD_U_VARIABLE_TYPE, &
                     & FIELD_VALUES_SET_TYPE,EQUATIONS_SET_FIELD_DATA,err,error,*999)
                    Ncompartments=EQUATIONS_SET_FIELD_DATA(2)
                  CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,(2+2*Ncompartments),err,error,*999)
                  ALLOCATE(VARIABLE_TYPES(2*Ncompartments+2))
                  DO num_var=1,Ncompartments+1
                    VARIABLE_TYPES(2*num_var-1)=FIELD_U_VARIABLE_TYPE+(FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(num_var-1))
                    VARIABLE_TYPES(2*num_var)=FIELD_DELUDELN_VARIABLE_TYPE+(FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(num_var-1))
                  ENDDO
                  CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,VARIABLE_TYPES,err,error,*999)

                  CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & NUMBER_OF_DIMENSIONS,err,error,*999)
                  NUMBER_OF_COMPONENTS=NUMBER_OF_DIMENSIONS+1
                  NUMBER_OF_DARCY_COMPONENTS=NUMBER_OF_DIMENSIONS+1

                  DO num_var=1,2*Ncompartments+2
                    CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,VARIABLE_TYPES(num_var),FIELD_VECTOR_DIMENSION_TYPE, &
                      & err,error,*999)
                    CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,VARIABLE_TYPES(num_var),FIELD_DP_TYPE,err,error,*999)
                    CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,VARIABLE_TYPES(num_var),NUMBER_OF_COMPONENTS, &
                      & err,error,*999)
                  ENDDO

                  SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                  CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                    !Elasticity:
                   DO component_idx=1,NUMBER_OF_DIMENSIONS
                     CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,component_idx, &
                       & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                     CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,component_idx,&
                       & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                   ENDDO !component_idx
                   !If solid hydrostatic pressure is driving Darcy flow, check that pressure uses node based interpolation
                   CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_COMPONENTS,&
                      & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                   CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                      & NUMBER_OF_COMPONENTS,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  DO num_var=3,2*Ncompartments+2
                    !Darcy:
                    DO component_idx=1,NUMBER_OF_DARCY_COMPONENTS
                      CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,VARIABLE_TYPES(num_var),component_idx, &
                        & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    ENDDO !component_idx
                  ENDDO
                  CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                    CALL FlagError("Not implemented.",err,error,*999)
                  CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                    CALL FlagError("Not implemented.",err,error,*999)
                  CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                    CALL FlagError("Not implemented.",err,error,*999)
                  CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                    CALL FlagError("Not implemented.",err,error,*999)
                  CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                    CALL FlagError("Not implemented.",err,error,*999)
                  CASE DEFAULT
                    localError="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                      & " is invalid."
                    CALL FlagError(localError,err,error,*999)
                  END SELECT
                CASE(EQUATIONS_SET_MULTI_COMPARTMENT_DARCY_SUBTYPE)
                  !Check the field created by Darcy routines for the multi-compartment model
                  CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                  CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DEPENDENT_TYPE,err,error,*999)
                  EQUATIONS_SET_FIELD_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD
                  CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET_FIELD_FIELD,FIELD_U_VARIABLE_TYPE, &
                   & FIELD_VALUES_SET_TYPE,EQUATIONS_SET_FIELD_DATA,err,error,*999)
                  Ncompartments=EQUATIONS_SET_FIELD_DATA(2)
                  CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,2*Ncompartments,err,error,*999)
                  !Create & populate array storing all of the relevant variable types against which to check the field variables
                  ALLOCATE(VARIABLE_TYPES(2*Ncompartments))
                  DO num_var=1,Ncompartments
                    VARIABLE_TYPES(2*num_var-1)=FIELD_U_VARIABLE_TYPE+(FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(num_var-1))
                    VARIABLE_TYPES(2*num_var)=FIELD_DELUDELN_VARIABLE_TYPE+(FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(num_var-1))
                  ENDDO
                  CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,VARIABLE_TYPES,err,error,*999)

                  CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & NUMBER_OF_DIMENSIONS,err,error,*999)
                  DO num_var=1,2*Ncompartments
                    CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,VARIABLE_TYPES(num_var), &
                       & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                    CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,VARIABLE_TYPES(num_var),FIELD_DP_TYPE,err,error,*999)
                    CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,VARIABLE_TYPES(num_var), &
                       & NUMBER_OF_DIMENSIONS+1,err,error,*999)
                  ENDDO

                  SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                  CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                      component_idx=1
                    DO num_var=1,2*Ncompartments
                     DO component_idx=1,NUMBER_OF_DIMENSIONS+1
                      CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,VARIABLE_TYPES(num_var),component_idx, &
                        & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                      !NOTE-pressure might use element based interpolation - need to account for this
                     ENDDO
                    ENDDO
                  CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                    CALL FlagError("Not implemented.",err,error,*999)
                  CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                    CALL FlagError("Not implemented.",err,error,*999)
                  CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                    CALL FlagError("Not implemented.",err,error,*999)
                  CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                    CALL FlagError("Not implemented.",err,error,*999)
                  CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                    CALL FlagError("Not implemented.",err,error,*999)
                  CASE DEFAULT
                    localError="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                      & " is invalid."
                    CALL FlagError(localError,err,error,*999)
                  END SELECT

                CASE DEFAULT
                  !--------------------------------
                  ! Check the user specified field
                  !--------------------------------
                  CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                  CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DEPENDENT_TYPE,err,error,*999)
                  CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,2,err,error,*999)
                  CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE],&
                    & err,error,*999)
                  CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                    & err,error,*999)
                  CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                    & err,error,*999)
                  CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                  CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD, FIELD_U_VARIABLE_TYPE, &
                    & NUMBER_OF_DIMENSIONS, ERR, ERROR, *999)
                  DEPENDENT_FIELD_NUMBER_OF_COMPONENTS = NUMBER_OF_DIMENSIONS + 1
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                    & DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                    & DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)

                  SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                  CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1, &
                      & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                      & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                    CALL FlagError("Not implemented.",err,error,*999)
                  CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                    CALL FlagError("Not implemented.",err,error,*999)
                  CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                    CALL FlagError("Not implemented.",err,error,*999)
                  CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                    CALL FlagError("Not implemented.",err,error,*999)
                  CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                    CALL FlagError("Not implemented.",err,error,*999)
                  CASE DEFAULT
                    localError="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                      & " is invalid."
                    CALL FlagError(localError,err,error,*999)
                  END SELECT
                END SELECT ! on (EQUATIONS_SET%SPECIFICATION(3))
              ENDIF ! on (EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED)
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
                CALL FIELD_CREATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,err,error,*999)
                CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                   & FIELD_INITIAL_VALUES_SET_TYPE,err,error,*999)
              ENDIF
              IF(EQUATIONS_SET%SPECIFICATION(3)/=EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)THEN
              !Actually, only needed for PGM (for elasticity_Darcy defined in elasticity V var):
                CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                   & FIELD_RELATIVE_VELOCITY_SET_TYPE,err,error,*999)
              ENDIF
            CASE DEFAULT
              localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
                & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
                & " is invalid for a standard, quasistatic or ALE Darcy equation"
              CALL FlagError(localError,err,error,*999)
            END SELECT
          END SELECT

        !-----------------------------------------------------------------
        ! I N d e p e n d e n t   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_INDEPENDENT_TYPE)
          SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
          CASE(EQUATIONS_SET_STANDARD_DARCY_SUBTYPE, EQUATIONS_SET_QUASISTATIC_DARCY_SUBTYPE, EQUATIONS_SET_ALE_DARCY_SUBTYPE, &
            & EQUATIONS_SET_TRANSIENT_DARCY_SUBTYPE,EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE, &
            & EQUATIONS_SET_TRANSIENT_ALE_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
            & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE, &
            & EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
            !\todo: revise: do they all need an independent field ?
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              IF(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
                !Create the auto created INdependent field
                CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION, &
                  & EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,err,error,*999)
                CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_GENERAL_TYPE,err,error,*999)

                CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_INDEPENDENT_TYPE,&
                  & err,error,*999)

                CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,err,error,*999)
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,GEOMETRIC_DECOMPOSITION, &
                  & err,error,*999)
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                  & EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,err,error,*999)

                INDEPENDENT_FIELD_NUMBER_OF_VARIABLES = 2  ! U and the normal component of its flux
                CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                  & INDEPENDENT_FIELD_NUMBER_OF_VARIABLES,err,error,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,[FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DELUDELN_VARIABLE_TYPE],err,error,*999)
                CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,"Independent U", &
                    & err,error,*999)
                CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                    & "Independent del U/del n",err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)

                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD, FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS, ERR, ERROR, *999)
                INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS = NUMBER_OF_DIMENSIONS !+ 1
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, FIELD_U_VARIABLE_TYPE, &
                  & INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS, ERR, ERROR, *999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                  & FIELD_DELUDELN_VARIABLE_TYPE, INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS, ERR, ERROR, *999)
                !Default to the geometric interpolation setup
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                  & GEOMETRIC_MESH_COMPONENT,err,error,*999)
                DO i=1,INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS
                  IF( i < INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS ) THEN
                    !Set velocity mesh component (default to the geometric one)
                    MESH_COMPONENT = GEOMETRIC_MESH_COMPONENT
                    CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & i, MESH_COMPONENT,err,error,*999)
                    CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,&
                      & FIELD_DELUDELN_VARIABLE_TYPE, i, MESH_COMPONENT,err,error,*999)
                  ELSE
                    !Set pressure mesh component (default to the geometric one)
                    MESH_COMPONENT = GEOMETRIC_MESH_COMPONENT
                    CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & i, MESH_COMPONENT,err,error,*999)
                    CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,&
                      & FIELD_DELUDELN_VARIABLE_TYPE, i, MESH_COMPONENT,err,error,*999)
                  ENDIF
                ENDDO

                SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  DO i = 1, INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,i,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                      & FIELD_DELUDELN_VARIABLE_TYPE,i,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  ENDDO
                  !Default the scaling to the geometric field scaling
                  CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                  CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE DEFAULT
                  localError="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                    & " is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              ELSE
                !Check the user specified field
                CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,err,error,*999)

                CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)

                CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,2,err,error,*999)
                CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE], &
                  & err,error,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD, FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS, ERR, ERROR, *999)
                INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS = NUMBER_OF_DIMENSIONS !+ 1
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                  & INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)

                SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE DEFAULT
                  localError="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                    & " is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              ENDIF
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              IF(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
                CALL FIELD_CREATE_FINISH(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,err,error,*999)
              ENDIF
            CASE DEFAULT
              localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
                & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
                & " is invalid for a standard, quasistatic or ALE Darcy equation"
              CALL FlagError(localError,err,error,*999)
            END SELECT
          END SELECT

        !-----------------------------------------------------------------
        !   m a t e r i a l   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
          SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
          CASE(EQUATIONS_SET_STANDARD_DARCY_SUBTYPE, EQUATIONS_SET_QUASISTATIC_DARCY_SUBTYPE, EQUATIONS_SET_ALE_DARCY_SUBTYPE, &
            & EQUATIONS_SET_TRANSIENT_DARCY_SUBTYPE,EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE, &
            & EQUATIONS_SET_TRANSIENT_ALE_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
            & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE, EQUATIONS_SET_MULTI_COMPARTMENT_DARCY_SUBTYPE)
            MATERIAL_FIELD_NUMBER_OF_VARIABLES = 1
            SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
            CASE(EQUATIONS_SET_STANDARD_DARCY_SUBTYPE,EQUATIONS_SET_QUASISTATIC_DARCY_SUBTYPE,EQUATIONS_SET_ALE_DARCY_SUBTYPE)
              !Porosity + scalar permeability/viscosity
              MATERIAL_FIELD_NUMBER_OF_COMPONENTS = 2
            CASE DEFAULT
              !Porosity + symmetric permeability/viscosity tensor
              MATERIAL_FIELD_NUMBER_OF_COMPONENTS = 7
            END SELECT
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
              IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
                IF(EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
                  !Create the auto created materials field
                  CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_MATERIALS% &
                    & MATERIALS_FIELD,err,error,*999)
                  CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_MATERIAL_TYPE,err,error,*999)
                  CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
                  CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,err,error,*999)
                  CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,GEOMETRIC_DECOMPOSITION, &
                    & err,error,*999)
                  CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,EQUATIONS_SET%GEOMETRY% &
                    & GEOMETRIC_FIELD,err,error,*999)
                  CALL FIELD_NUMBER_OF_VARIABLES_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD, &
                    & MATERIAL_FIELD_NUMBER_OF_VARIABLES,err,error,*999)
                  CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,[FIELD_U_VARIABLE_TYPE], &
                    & err,error,*999)
                  CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE,"Material", &
                      & err,error,*999)
                  CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                  CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_DP_TYPE,err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & MATERIAL_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & 1,GEOMETRIC_COMPONENT_NUMBER,err,error,*999)

                  !Auto-created / default is node_based_interpolation: that's an expensive default ...
                  !Maybe default should be constant; node_based should be requested by the user \todo
                  DO i = 1, MATERIAL_FIELD_NUMBER_OF_COMPONENTS
                    CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & i,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & i,GEOMETRIC_COMPONENT_NUMBER,err,error,*999)
                  END DO

                  !Default the field scaling to that of the geometric field
                  CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                  CALL FIELD_SCALING_TYPE_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                ELSE
                  !Check the user specified field
                  CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_MATERIAL_TYPE,err,error,*999)
                  CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
                  CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,MATERIAL_FIELD_NUMBER_OF_VARIABLES,err,error,*999)
                  CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE],err,error,*999)
                  CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                    & err,error,*999)
                  CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                    & MATERIAL_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
                ENDIF
              ELSE
                CALL FlagError("Equations set materials is not associated.",err,error,*999)
              ENDIF
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
              IF( ASSOCIATED(EQUATIONS_MATERIALS) ) THEN
                IF( EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED ) THEN
                  CALL FIELD_CREATE_FINISH(EQUATIONS_MATERIALS%MATERIALS_FIELD,err,error,*999)
                  !Set the default values for the materials field
                  DO i=1,MATERIAL_FIELD_NUMBER_OF_COMPONENTS
                    CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE, i, 1.0_DP, ERR, ERROR, *999)
                  ENDDO
                ENDIF
              ELSE
                CALL FlagError("Equations set materials is not associated.",err,error,*999)
              ENDIF
            CASE DEFAULT
              localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
                & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
                & " is invalid for a standard, quasistatic or ALE Darcy equation."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE(EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
          !Materials field needs two extra variable types
          !The V variable type stores the Darcy coupling coefficients that govern flux between compartments
          !The U1 variable type stores the parameters for the constitutive laws that determine the partial pressure in each compartment
          !For a first attempt at this, it will be assumed that the functional form of this law is the same for each compartment, with only the paramenters varying (default will be three components)
            EQUATIONS_SET_FIELD_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD
            CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET_FIELD_FIELD,FIELD_U_VARIABLE_TYPE, &
              & FIELD_VALUES_SET_TYPE,EQUATIONS_SET_FIELD_DATA,err,error,*999)
            Ncompartments=EQUATIONS_SET_FIELD_DATA(2)
            MATERIAL_FIELD_NUMBER_OF_VARIABLES = 3
            MATERIAL_FIELD_NUMBER_OF_U_VAR_COMPONENTS = 2
            MATERIAL_FIELD_NUMBER_OF_V_VAR_COMPONENTS = Ncompartments
            MATERIAL_FIELD_NUMBER_OF_U1_VAR_COMPONENTS = 3
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
              IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
                IF(EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
                  !Create the auto created materials field
                  CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_MATERIALS% &
                    & MATERIALS_FIELD,err,error,*999)
                  CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_MATERIAL_TYPE,err,error,*999)
                  CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
                  CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,err,error,*999)
                  CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,GEOMETRIC_DECOMPOSITION, &
                    & err,error,*999)
                  CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,EQUATIONS_SET%GEOMETRY% &
                    & GEOMETRIC_FIELD,err,error,*999)
                  CALL FIELD_NUMBER_OF_VARIABLES_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD, &
                    & MATERIAL_FIELD_NUMBER_OF_VARIABLES,err,error,*999)
                  CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,[FIELD_U_VARIABLE_TYPE, &
                    & FIELD_V_VARIABLE_TYPE,FIELD_U1_VARIABLE_TYPE],err,error,*999)
                  CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                  CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_DP_TYPE,err,error,*999)
                  CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                    & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                  CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                    & FIELD_DP_TYPE,err,error,*999)
                  CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U1_VARIABLE_TYPE, &
                    & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                  CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U1_VARIABLE_TYPE, &
                    & FIELD_DP_TYPE,err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & MATERIAL_FIELD_NUMBER_OF_U_VAR_COMPONENTS,err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                    & MATERIAL_FIELD_NUMBER_OF_V_VAR_COMPONENTS,err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U1_VARIABLE_TYPE, &
                    & MATERIAL_FIELD_NUMBER_OF_U1_VAR_COMPONENTS,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & 1,GEOMETRIC_COMPONENT_NUMBER,err,error,*999)

                  !Auto-created / default is node_based_interpolation: that's an expensive default ...
                  !Maybe default should be constant; node_based should be requested by the user \todo
                  DO i = 1, MATERIAL_FIELD_NUMBER_OF_U_VAR_COMPONENTS
                    CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & i,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & i,GEOMETRIC_COMPONENT_NUMBER,err,error,*999)
                  END DO
                  DO i = 1, MATERIAL_FIELD_NUMBER_OF_V_VAR_COMPONENTS
                    CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                      & i,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                      & i,GEOMETRIC_COMPONENT_NUMBER,err,error,*999)
                  END DO
                  DO i = 1, MATERIAL_FIELD_NUMBER_OF_U1_VAR_COMPONENTS
                    CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U1_VARIABLE_TYPE, &
                      & i,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U1_VARIABLE_TYPE, &
                      & i,GEOMETRIC_COMPONENT_NUMBER,err,error,*999)
                  END DO

                  !Default the field scaling to that of the geometric field
                  CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                  CALL FIELD_SCALING_TYPE_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                ELSE
                  !Check the user specified field
                  CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_MATERIAL_TYPE,err,error,*999)
                  CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
                  CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,MATERIAL_FIELD_NUMBER_OF_VARIABLES,err,error,*999)
                  CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE, &
                    & FIELD_V_VARIABLE_TYPE,FIELD_U1_VARIABLE_TYPE],err,error,*999)
                  CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                    & err,error,*999)
                  CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                    & err,error,*999)
                  CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U1_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                    & err,error,*999)
                  CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                  CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                  CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U1_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                     & MATERIAL_FIELD_NUMBER_OF_U_VAR_COMPONENTS,err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE, &
                     & MATERIAL_FIELD_NUMBER_OF_V_VAR_COMPONENTS,err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U1_VARIABLE_TYPE, &
                     & MATERIAL_FIELD_NUMBER_OF_U1_VAR_COMPONENTS,err,error,*999)
                ENDIF
              ELSE
                CALL FlagError("Equations set materials is not associated.",err,error,*999)
              ENDIF
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
              IF( ASSOCIATED(EQUATIONS_MATERIALS) ) THEN
                IF( EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED ) THEN
                  CALL FIELD_CREATE_FINISH(EQUATIONS_MATERIALS%MATERIALS_FIELD,err,error,*999)
                  !Set the default values for the materials field
                  DO i=1,MATERIAL_FIELD_NUMBER_OF_U_VAR_COMPONENTS
                    CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE, i, 1.0_DP, ERR, ERROR, *999)
                  ENDDO
                  DO i=1,MATERIAL_FIELD_NUMBER_OF_V_VAR_COMPONENTS
                    CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE, i, 0.0_DP, ERR, ERROR, *999)
                  ENDDO
                  DO i=1,MATERIAL_FIELD_NUMBER_OF_U1_VAR_COMPONENTS
                    CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U1_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE, i, 0.0_DP, ERR, ERROR, *999)
                  ENDDO
                ENDIF
              ELSE
                CALL FlagError("Equations set materials is not associated.",err,error,*999)
              ENDIF
            CASE DEFAULT
              localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
                & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
                & " is invalid for a standard, quasistatic or ALE Darcy equation."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          END SELECT

        !-----------------------------------------------------------------
        !   a n a l y t i c   f i e l d
        !-----------------------------------------------------------------

        CASE(EQUATIONS_SET_SETUP_ANALYTIC_TYPE)

          SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
          CASE(EQUATIONS_SET_STANDARD_DARCY_SUBTYPE)
              SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
                !Set start action
                CASE(EQUATIONS_SET_SETUP_START_ACTION)
                  IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
                    IF(ASSOCIATED(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD)) THEN
                      IF(ASSOCIATED(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD)) THEN
                        CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                          & NUMBER_OF_DIMENSIONS,err,error,*999)
                        SELECT CASE(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE)
                          CASE(EQUATIONS_SET_DARCY_EQUATION_TWO_DIM_1)
                            !Set analytic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_DARCY_EQUATION_TWO_DIM_1
                          CASE(EQUATIONS_SET_DARCY_EQUATION_TWO_DIM_2)
                            !Set analytic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_DARCY_EQUATION_TWO_DIM_2
                          CASE(EQUATIONS_SET_DARCY_EQUATION_TWO_DIM_3)
                            !Set analytic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_DARCY_EQUATION_TWO_DIM_3
                          CASE(EQUATIONS_SET_DARCY_EQUATION_THREE_DIM_1)
                            !Set analytic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_DARCY_EQUATION_THREE_DIM_1
                          CASE(EQUATIONS_SET_DARCY_EQUATION_THREE_DIM_2)
                            !Set analytic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_DARCY_EQUATION_THREE_DIM_2
                          CASE(EQUATIONS_SET_DARCY_EQUATION_THREE_DIM_3)
                            !Set analytic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_DARCY_EQUATION_THREE_DIM_3
                          CASE(EQUATIONS_SET_INCOMP_ELAST_DARCY_ANALYTIC_DARCY)
                            !Set analytic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_INCOMP_ELAST_DARCY_ANALYTIC_DARCY
                          CASE DEFAULT
                            localError="The specified analytic function type of "// &
                              & TRIM(NumberToVString(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",err,error))// &
                              & " is invalid for an analytic Darcy problem."
                            CALL FlagError(localError,err,error,*999)
                        END SELECT
                      ELSE
                        CALL FlagError("Equations set geometric field is not associated.",err,error,*999)
                      ENDIF
                    ELSE
                      CALL FlagError("Equations set dependent field is not associated.",err,error,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Equations set dependent field has not been finished.",err,error,*999)
                  ENDIF
                CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
                  IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
                    IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC%ANALYTIC_FIELD)) THEN
                      IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FIELD_AUTO_CREATED) THEN
                        !--- Why finish the dependent field and not the analytic one ???
                        CALL FIELD_CREATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,err,error,*999)
                      ENDIF
                    ENDIF
                  ELSE
                    CALL FlagError("Equations set analytic is not associated.",err,error,*999)
                  ENDIF
                CASE DEFAULT
                  localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
                    & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
                    & " is invalid for an analytic Darcy problem."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
          CASE(EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE)
              SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
                !Set start action
                CASE(EQUATIONS_SET_SETUP_START_ACTION)
                  IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
                    IF(ASSOCIATED(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD)) THEN
                      IF(ASSOCIATED(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD)) THEN
                        CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                          & NUMBER_OF_DIMENSIONS,err,error,*999)
                        !Initialise analytic parameter which stores value of time to zero - need to update this somewhere in a pre_solve routine
                        EQUATIONS_SET%ANALYTIC%ANALYTIC_USER_PARAMS(1)=0.0_DP
                        SELECT CASE(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE)
                        CASE(EQUATIONS_SET_INCOMP_ELAST_DARCY_ANALYTIC_DARCY)
                          !Set analytic function type
                          EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_INCOMP_ELAST_DARCY_ANALYTIC_DARCY
                        CASE DEFAULT
                          localError="The specified analytic function type of "// &
                            & TRIM(NumberToVString(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",err,error))// &
                            & " is invalid for an analytic Darcy problem."
                          CALL FlagError(localError,err,error,*999)
                        END SELECT
                      ELSE
                        CALL FlagError("Equations set geometric field is not associated.",err,error,*999)
                      ENDIF
                    ELSE
                      CALL FlagError("Equations set dependent field is not associated.",err,error,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Equations set dependent field has not been finished.",err,error,*999)
                  ENDIF
                CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
                  IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
                    IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC%ANALYTIC_FIELD)) THEN
                      IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FIELD_AUTO_CREATED) THEN
                        !--- Why finish the dependent field and not the analytic one ???
                        CALL FIELD_CREATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,err,error,*999)
                      ENDIF
                    ENDIF
                  ELSE
                    CALL FlagError("Equations set analytic is not associated.",err,error,*999)
                  ENDIF
                CASE DEFAULT
                  localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
                    & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
                    & " is invalid for an analytic Darcy problem."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
          CASE DEFAULT
            localError="The equation set subtype of "//TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a Darcy equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT

        !-----------------------------------------------------------------
        !   s o u r c e   t y p e   -   include gravity at some point
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_SOURCE_TYPE)
          SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
          CASE(EQUATIONS_SET_STANDARD_DARCY_SUBTYPE, EQUATIONS_SET_QUASISTATIC_DARCY_SUBTYPE, EQUATIONS_SET_ALE_DARCY_SUBTYPE, &
            & EQUATIONS_SET_TRANSIENT_DARCY_SUBTYPE,EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE, &
            & EQUATIONS_SET_TRANSIENT_ALE_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
            & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE, &
            & EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              EQUATIONS_SOURCE=>EQUATIONS_SET%SOURCE
              IF(ASSOCIATED(EQUATIONS_SOURCE)) THEN
                IF(EQUATIONS_SOURCE%SOURCE_FIELD_AUTO_CREATED) THEN
                  CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_SOURCE% &
                    & SOURCE_FIELD,err,error,*999)
                  CALL FIELD_LABEL_SET(EQUATIONS_SOURCE%SOURCE_FIELD,"Source Field",err,error,*999)
                  CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SOURCE%SOURCE_FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                  CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SOURCE%SOURCE_FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
                  CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,err,error,*999)
                  CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SOURCE%SOURCE_FIELD,GEOMETRIC_DECOMPOSITION, &
                    & err,error,*999)
                  CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SOURCE%SOURCE_FIELD,EQUATIONS_SET%GEOMETRY% &
                    & GEOMETRIC_FIELD,err,error,*999)
                  CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SOURCE%SOURCE_FIELD,1,err,error,*999)
                  CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SOURCE%SOURCE_FIELD,[FIELD_U_VARIABLE_TYPE], &
                    & err,error,*999)
                  CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE,"Source", &
                      & err,error,*999)
                  CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                  CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_DP_TYPE,err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD, FIELD_U_VARIABLE_TYPE, &
                    & NUMBER_OF_DIMENSIONS, ERR, ERROR, *999)
                    NUMBER_OF_SOURCE_COMPONENTS = NUMBER_OF_DIMENSIONS + 1
                  CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & NUMBER_OF_SOURCE_COMPONENTS,err,error,*999)

                  !Default the source components to the geometric interpolation setup with nodal interpolation
                  IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE .OR. &
                    & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE) THEN
                    !nodal / mesh based
                    DO component_idx=1,NUMBER_OF_DIMENSIONS !NUMBER_OF_SOURCE_COMPONENTS
                      CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                        & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                      CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                        & component_idx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                      CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                       & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                    ENDDO !component_idx
                    !Set source component 'NUMBER_OF_DIMENSIONS + 1' according to GEOMETRIC_MESH_COMPONENT 'NUMBER_OF_DIMENSIONS'
                    CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & NUMBER_OF_DIMENSIONS + 1,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & NUMBER_OF_DIMENSIONS + 1,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  ENDIF
                  !Default the field scaling to that of the geometric field
                  CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                  CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SOURCE%SOURCE_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                ELSE
                  !Check the user specified field
                  CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                  CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
                  CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,1,err,error,*999)
                  CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE],err,error,*999)
                  CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                    & err,error,*999)
                  CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD, FIELD_U_VARIABLE_TYPE, &
                    & NUMBER_OF_DIMENSIONS, ERR, ERROR, *999)
                    NUMBER_OF_SOURCE_COMPONENTS = NUMBER_OF_DIMENSIONS + 1
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                    & NUMBER_OF_SOURCE_COMPONENTS,err,error,*999)
                ENDIF
              ELSE
                CALL FlagError("Equations set source is not associated.",err,error,*999)
              ENDIF
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              EQUATIONS_SOURCE=>EQUATIONS_SET%SOURCE
              IF(ASSOCIATED(EQUATIONS_SOURCE)) THEN
                IF(EQUATIONS_SOURCE%SOURCE_FIELD_AUTO_CREATED) THEN
                  !Finish creating the source field
                  CALL FIELD_CREATE_FINISH(EQUATIONS_SOURCE%SOURCE_FIELD,err,error,*999)
                  !Set the default values for the source field
                  CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & NUMBER_OF_DIMENSIONS,err,error,*999)
                  IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE .OR. &
                    & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE) THEN
                    NUMBER_OF_SOURCE_COMPONENTS = NUMBER_OF_DIMENSIONS + 1
                  ELSE
                    NUMBER_OF_SOURCE_COMPONENTS=0
                  ENDIF
                  !Now set the source values to 0.0
                  DO component_idx=1,NUMBER_OF_SOURCE_COMPONENTS
                    CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,component_idx,0.0_DP,err,error,*999)
                  ENDDO !component_idx
                ENDIF
              ELSE
                CALL FlagError("Equations set source is not associated.",err,error,*999)
              ENDIF
            CASE DEFAULT
              localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
                & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
                & " is invalid for a standard, quasistatic or ALE Darcy equation."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          END SELECT

        !-----------------------------------------------------------------
        !   e q u a t i o n s   t y p e
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
          SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
          !-----------------------------------------------------------------
          !   s t a t i c
          !-----------------------------------------------------------------
          CASE(EQUATIONS_SET_STANDARD_DARCY_SUBTYPE,EQUATIONS_SET_MULTI_COMPARTMENT_DARCY_SUBTYPE)
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
              IF( ASSOCIATED(EQUATIONS_MATERIALS) ) THEN
                IF( EQUATIONS_MATERIALS%MATERIALS_FINISHED ) THEN
                  CALL Equations_CreateStart(EQUATIONS_SET,equations,err,error,*999)
                  CALL Equations_LinearityTypeSet(equations,EQUATIONS_LINEAR,err,error,*999)
                  CALL Equations_TimeDependenceTypeSet(equations,EQUATIONS_STATIC,err,error,*999)
                ELSE
                  CALL FlagError("Equations set materials has not been finished.",err,error,*999)
                ENDIF
              ELSE
                CALL FlagError("Equations set materials is not associated.",err,error,*999)
              ENDIF
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
                CASE(EQUATIONS_SET_MULTI_COMPARTMENT_DARCY_SUBTYPE)
                    !!!!!THE FOLLOWING IF STATEMENT IS ILLUSTRATIVE ONLY - need to implement the equation set field thing, and make a generalised case statement
                  CALL EquationsSet_EquationsGet(EQUATIONS_SET,equations,err,error,*999)
                  CALL Equations_CreateFinish(equations,err,error,*999)
                  NULLIFY(vectorEquations)
                  CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
                  !Create the equations mapping.
                  CALL EquationsMapping_VectorCreateStart(vectorEquations,FIELD_DELUDELN_VARIABLE_TYPE,vectorMapping,err,error,*999)
                  CALL EquationsMapping_LinearMatricesNumberSet(vectorMapping,1,err,error,*999)
                  EQUATIONS_SET_FIELD_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD
                  CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET_FIELD_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,EQUATIONS_SET_FIELD_DATA,err,error,*999)
                  imy_matrix = EQUATIONS_SET_FIELD_DATA(1)
                  Ncompartments = EQUATIONS_SET_FIELD_DATA(2)
                  ALLOCATE(VARIABLE_TYPES(2*Ncompartments))
                  DO num_var=1,Ncompartments
                    VARIABLE_TYPES(2*num_var-1)=FIELD_U_VARIABLE_TYPE+(FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(num_var-1))
                    VARIABLE_TYPES(2*num_var)=FIELD_DELUDELN_VARIABLE_TYPE+(FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(num_var-1))
                  ENDDO
                  CALL EquationsMapping_LinearMatricesVariableTypesSet(vectorMapping,[VARIABLE_TYPES(2*imy_matrix-1)], &
                    & err,error,*999)
                  CALL EquationsMapping_RHSVariableTypeSet(vectorMapping,VARIABLE_TYPES(2*imy_matrix),err,error,*999)
                  CALL EquationsMapping_VectorCreateFinish(vectorMapping,err,error,*999)
                  !Create the equations matrices
                  CALL EquationsMatrices_VectorCreateStart(vectorEquations,vectorMatrices,err,error,*999)
                  SELECT CASE(equations%sparsityType)
                  CASE(EQUATIONS_MATRICES_FULL_MATRICES)
                    CALL EquationsMatrices_LinearStorageTypeSet(vectorMatrices,[MATRIX_BLOCK_STORAGE_TYPE], &
                      & err,error,*999)
                  CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                    CALL EquationsMatrices_LinearStorageTypeSet(vectorMatrices,[MATRIX_COMPRESSED_ROW_STORAGE_TYPE], &
                      & err,error,*999)
                    CALL EquationsMatrices_LinearStructureTypeSet(vectorMatrices,[EQUATIONS_MATRIX_FEM_STRUCTURE], &
                      & err,error,*999)
                  CASE DEFAULT
                    localError="The equations matrices sparsity type of "// &
                      & TRIM(NumberToVString(equations%sparsityType,"*",err,error))//" is invalid."
                    CALL FlagError(localError,err,error,*999)
                  END SELECT
                  CALL EquationsMatrices_VectorCreateFinish(vectorMatrices,err,error,*999)
                CASE DEFAULT
                  !Finish the equations creation
                  CALL EquationsSet_EquationsGet(EQUATIONS_SET,equations,err,error,*999)
                  CALL Equations_CreateFinish(equations,err,error,*999)
                  !Create the equations mapping.
                  CALL EquationsMapping_VectorCreateStart(vectorEquations,FIELD_DELUDELN_VARIABLE_TYPE,vectorMapping,err,error,*999)
                  CALL EquationsMapping_LinearMatricesNumberSet(vectorMapping,1,err,error,*999)
                  CALL EquationsMapping_LinearMatricesVariableTypesSet(vectorMapping,[FIELD_U_VARIABLE_TYPE], &
                    & err,error,*999)
                  CALL EquationsMapping_RHSVariableTypeSet(vectorMapping,FIELD_DELUDELN_VARIABLE_TYPE,err,error,*999)
                  CALL EquationsMapping_VectorCreateFinish(vectorMapping,err,error,*999)
                  !Create the equations matrices
                  CALL EquationsMatrices_VectorCreateStart(vectorEquations,vectorMatrices,err,error,*999)
                  SELECT CASE(equations%sparsityType)
                  CASE(EQUATIONS_MATRICES_FULL_MATRICES)
                    CALL EquationsMatrices_LinearStorageTypeSet(vectorMatrices,[MATRIX_BLOCK_STORAGE_TYPE], &
                      & err,error,*999)
                  CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                    CALL EquationsMatrices_LinearStorageTypeSet(vectorMatrices,[MATRIX_COMPRESSED_ROW_STORAGE_TYPE], &
                      & err,error,*999)
                    CALL EquationsMatrices_LinearStructureTypeSet(vectorMatrices,[EQUATIONS_MATRIX_FEM_STRUCTURE], &
                      & err,error,*999)
                  CASE DEFAULT
                    localError="The equations matrices sparsity type of "// &
                      & TRIM(NumberToVString(equations%sparsityType,"*",err,error))//" is invalid."
                    CALL FlagError(localError,err,error,*999)
                  END SELECT
                  CALL EquationsMatrices_VectorCreateFinish(vectorMatrices,err,error,*999)
                END SELECT
              CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                CALL FlagError("Not implemented.",err,error,*999)
              CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                CALL FlagError("Not implemented.",err,error,*999)
              CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                CALL FlagError("Not implemented.",err,error,*999)
              CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                CALL FlagError("Not implemented.",err,error,*999)
              CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                CALL FlagError("Not implemented.",err,error,*999)
              CASE DEFAULT
                  localError="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                  & " is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT

            CASE DEFAULT
              localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
                & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
                & " is invalid for a standard Darcy equation."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          !-----------------------------------------------------------------
          !   q u a s i s t a t i c   and    A L E
          !-----------------------------------------------------------------
          CASE(EQUATIONS_SET_QUASISTATIC_DARCY_SUBTYPE, EQUATIONS_SET_ALE_DARCY_SUBTYPE, &
            & EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE)
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
              IF( ASSOCIATED(EQUATIONS_MATERIALS) ) THEN
                IF( EQUATIONS_MATERIALS%MATERIALS_FINISHED ) THEN
                  CALL Equations_CreateStart(EQUATIONS_SET,equations,err,error,*999)
                  CALL Equations_LinearityTypeSet(equations,EQUATIONS_LINEAR,err,error,*999)
                  CALL Equations_TimeDependenceTypeSet(equations,EQUATIONS_QUASISTATIC,err,error,*999)
                ELSE
                  CALL FlagError("Equations set materials has not been finished.",err,error,*999)
                ENDIF
              ELSE
                CALL FlagError("Equations set materials is not associated.",err,error,*999)
              ENDIF
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                !Finish the equations creation
                CALL EquationsSet_EquationsGet(EQUATIONS_SET,equations,err,error,*999)
                CALL Equations_CreateFinish(equations,err,error,*999)
                NULLIFY(vectorEquations)
                CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
                !Create the equations mapping.
                CALL EquationsMapping_VectorCreateStart(vectorEquations,FIELD_DELUDELN_VARIABLE_TYPE,vectorMapping,err,error,*999)
                CALL EquationsMapping_LinearMatricesNumberSet(vectorMapping,1,err,error,*999)
                SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
                CASE(EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE)
                  CALL EquationsMapping_LinearMatricesVariableTypesSet(vectorMapping,[FIELD_V_VARIABLE_TYPE], &
                    & err,error,*999)
                  CALL EquationsMapping_RHSVariableTypeSet(vectorMapping,FIELD_DELVDELN_VARIABLE_TYPE,err,error,*999)
                CASE DEFAULT
                  CALL EquationsMapping_LinearMatricesVariableTypesSet(vectorMapping,[FIELD_U_VARIABLE_TYPE], &
                    & err,error,*999)
                  CALL EquationsMapping_RHSVariableTypeSet(vectorMapping,FIELD_DELUDELN_VARIABLE_TYPE,err,error,*999)
                END SELECT
                CALL EquationsMapping_VectorCreateFinish(vectorMapping,err,error,*999)
                !Create the equations matrices
                CALL EquationsMatrices_VectorCreateStart(vectorEquations,vectorMatrices,err,error,*999)
                SELECT CASE(equations%sparsityType)
                CASE(EQUATIONS_MATRICES_FULL_MATRICES)
                  CALL EquationsMatrices_LinearStorageTypeSet(vectorMatrices,[MATRIX_BLOCK_STORAGE_TYPE], &
                    & err,error,*999)
                CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                  CALL EquationsMatrices_LinearStorageTypeSet(vectorMatrices,[MATRIX_COMPRESSED_ROW_STORAGE_TYPE], &
                    & err,error,*999)
                  CALL EquationsMatrices_LinearStructureTypeSet(vectorMatrices,[EQUATIONS_MATRIX_FEM_STRUCTURE], &
                    & err,error,*999)
                CASE DEFAULT
                  localError="The equations matrices sparsity type of "// &
                    & TRIM(NumberToVString(equations%sparsityType,"*",err,error))//" is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
                CALL EquationsMatrices_VectorCreateFinish(vectorMatrices,err,error,*999)
              CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                CALL FlagError("Not implemented.",err,error,*999)
              CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                CALL FlagError("Not implemented.",err,error,*999)
              CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                CALL FlagError("Not implemented.",err,error,*999)
              CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                CALL FlagError("Not implemented.",err,error,*999)
              CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                CALL FlagError("Not implemented.",err,error,*999)
              CASE DEFAULT
                localError="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                  & " is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            CASE DEFAULT
              localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
                & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
                & " is invalid for a quasistatic Darcy equation."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          !-----------------------------------------------------------------
          !   d y n a m i c
          !-----------------------------------------------------------------
          CASE(EQUATIONS_SET_TRANSIENT_DARCY_SUBTYPE,EQUATIONS_SET_TRANSIENT_ALE_DARCY_SUBTYPE, &
              & EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE,EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE, &
              & EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
              CASE(EQUATIONS_SET_SETUP_START_ACTION)
                EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
                IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
                  IF(EQUATIONS_MATERIALS%MATERIALS_FINISHED) THEN
                    CALL Equations_CreateStart(EQUATIONS_SET,equations,err,error,*999)
                    CALL Equations_LinearityTypeSet(equations,EQUATIONS_LINEAR,err,error,*999)
                    CALL Equations_TimeDependenceTypeSet(equations,EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
                  ELSE
                    CALL FlagError("Equations set materials has not been finished.",err,error,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Equations materials is not associated.",err,error,*999)
                ENDIF
              CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
                SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                  CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                    !Finish the equations creation
                    CALL EquationsSet_EquationsGet(EQUATIONS_SET,equations,err,error,*999)
                    CALL Equations_CreateFinish(equations,err,error,*999)
                    NULLIFY(vectorEquations)
                    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
                    !Create the equations mapping.
                    CALL EquationsMapping_VectorCreateStart(vectorEquations,FIELD_DELUDELN_VARIABLE_TYPE,vectorMapping, &
                      & err,error,*999)
                    IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE .OR. &
                      & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE) THEN
                      CALL EquationsMapping_LinearMatricesNumberSet(vectorMapping,0,err,error,*999)
                    ENDIF
                    CALL EquationsMapping_DynamicMatricesSet(vectorMapping,.TRUE.,.TRUE.,err,error,*999)
                    SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
                    CASE(EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
                        & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE)
                      CALL EquationsMapping_DynamicVariableTypeSet(vectorMapping,FIELD_V_VARIABLE_TYPE,err,error,*999)
                      CALL EquationsMapping_RHSVariableTypeSet(vectorMapping,FIELD_DELVDELN_VARIABLE_TYPE, &
                        & err,error,*999)
                      IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE) THEN
                        CALL EquationsMapping_SourceVariableTypeSet(vectorMapping,FIELD_U_VARIABLE_TYPE,err,error,*999)
                      ENDIF
                    CASE(EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
                      EQUATIONS_SET_FIELD_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD
                      CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET_FIELD_FIELD,FIELD_U_VARIABLE_TYPE, &
                         & FIELD_VALUES_SET_TYPE,EQUATIONS_SET_FIELD_DATA,err,error,*999)
                      imy_matrix = EQUATIONS_SET_FIELD_DATA(1)
                      Ncompartments = EQUATIONS_SET_FIELD_DATA(2)
                      CALL EquationsMapping_LinearMatricesNumberSet(vectorMapping,Ncompartments-1,err,error,*999)
                      ALLOCATE(VARIABLE_TYPES(2*Ncompartments+2))
                      ALLOCATE(VARIABLE_U_TYPES(Ncompartments-1))
                      DO num_var=1,Ncompartments+1
                        VARIABLE_TYPES(2*num_var-1)=FIELD_U_VARIABLE_TYPE+(FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(num_var-1))
                        VARIABLE_TYPES(2*num_var)=FIELD_DELUDELN_VARIABLE_TYPE+(FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(num_var-1))
                      ENDDO
                      num_var_count=0
                      DO num_var=2,Ncompartments+1
                        IF((num_var-1)/=imy_matrix)THEN
                          num_var_count=num_var_count+1
                          VARIABLE_U_TYPES(num_var_count)=VARIABLE_TYPES(2*num_var-1)
                        ENDIF
                      ENDDO
                      CALL EquationsMapping_DynamicVariableTypeSet(vectorMapping,VARIABLE_TYPES(2*imy_matrix+1), &
                         & err,error,*999)
                      CALL EquationsMapping_LinearMatricesVariableTypesSet(vectorMapping,VARIABLE_U_TYPES,err,error,*999)
                      CALL EquationsMapping_RHSVariableTypeSet(vectorMapping,VARIABLE_TYPES(2*imy_matrix+2),err,error,*999)
                      CALL EquationsMapping_SourceVariableTypeSet(vectorMapping,FIELD_U_VARIABLE_TYPE,err,error,*999)
                    CASE DEFAULT
                    CALL EquationsMapping_DynamicVariableTypeSet(vectorMapping,FIELD_U_VARIABLE_TYPE,err,error,*999)
                    CALL EquationsMapping_RHSVariableTypeSet(vectorMapping,FIELD_DELUDELN_VARIABLE_TYPE, &
                      & err,error,*999)
                    END SELECT
                    CALL EquationsMapping_VectorCreateFinish(vectorMapping,err,error,*999)
                    !Create the equations matrices
                    CALL EquationsMatrices_VectorCreateStart(vectorEquations,vectorMatrices,err,error,*999)
                    !Set up matrix storage and structure
                    IF(equations%lumpingType==EQUATIONS_LUMPED_MATRICES) THEN
                      !Set up lumping
                      CALL EquationsMatrices_DynamicLumpingTypeSet(vectorMatrices, &
                        & [EQUATIONS_MATRIX_UNLUMPED,EQUATIONS_MATRIX_LUMPED],err,error,*999)
                      CALL EquationsMatrices_DynamicStorageTypeSet(vectorMatrices, &
                        & [DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE] &
                        & ,err,error,*999)
                      CALL EquationsMatrices_DynamicStructureTypeSet(vectorMatrices, &
                        & [EQUATIONS_MATRIX_FEM_STRUCTURE,EQUATIONS_MATRIX_DIAGONAL_STRUCTURE],err,error,*999)
                    ELSE
                      SELECT CASE(equations%sparsityType)
                        CASE(EQUATIONS_MATRICES_FULL_MATRICES)
                        CALL EquationsMatrices_LinearStorageTypeSet(vectorMatrices, &
                            & [DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE,DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE],err,error,*999)
                        CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                          CALL EquationsMatrices_DynamicStorageTypeSet(vectorMatrices, &
                            & [DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE, &
                            & DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE],err,error,*999)
                          CALL EquationsMatrices_DynamicStructureTypeSet(vectorMatrices, &
                            & [EQUATIONS_MATRIX_FEM_STRUCTURE,EQUATIONS_MATRIX_FEM_STRUCTURE],err,error,*999)
                            IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)THEN
                              ALLOCATE(COUPLING_MATRIX_STORAGE_TYPE(Ncompartments-1))
                              ALLOCATE(COUPLING_MATRIX_STRUCTURE_TYPE(Ncompartments-1))
                              DO num_var=1,Ncompartments-1
                               COUPLING_MATRIX_STORAGE_TYPE(num_var)=DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE
                               COUPLING_MATRIX_STRUCTURE_TYPE(num_var)=EQUATIONS_MATRIX_FEM_STRUCTURE
                              ENDDO
                              CALL EquationsMatrices_LinearStorageTypeSet(vectorMatrices,COUPLING_MATRIX_STORAGE_TYPE, &
                                & err,error,*999)
                              CALL EquationsMatrices_LinearStructureTypeSet(vectorMatrices,COUPLING_MATRIX_STRUCTURE_TYPE, &
                                & err,error,*999)
                            ENDIF
                        CASE DEFAULT
                          localError="The equations matrices sparsity type of "// &
                            & TRIM(NumberToVString(equations%sparsityType,"*",err,error))//" is invalid."
                          CALL FlagError(localError,err,error,*999)
                      END SELECT
                    ENDIF
                    CALL EquationsMatrices_VectorCreateFinish(vectorMatrices,err,error,*999)
                  CASE DEFAULT
                    localError="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*", &
                      & err,error))//" is invalid."
                    CALL FlagError(localError,err,error,*999)
                END SELECT
              CASE DEFAULT
              localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
                  & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
                  & " is invalid for a Darcy equation."
                CALL FlagError(localError,err,error,*999)
            END SELECT
          !-----------------------------------------------------------------
          !   D e f a u l t
          !-----------------------------------------------------------------
          CASE DEFAULT
            localError="The equation set subtype of "//TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
              & " for a setup of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a Darcy equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT

        !-----------------------------------------------------------------
        !   c a s e   d e f a u l t
        !-----------------------------------------------------------------
        CASE DEFAULT
          localError="The setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
            & " is invalid for a standard, quasistatic, ALE or dynamic Darcy equation."
          CALL FlagError(localError,err,error,*999)

        END SELECT
      CASE DEFAULT
        localError="The equations set subtype of "//TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
          & " does not equal a standard, quasistatic, ALE or dynamic Darcy equation subtype."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF

    EXITS("DARCY_EQUATION_EQUATIONS_SET_SETUP")
    RETURN
999 ERRORSEXITS("DARCY_EQUATION_EQUATIONS_SET_SETUP",err,error)
    RETURN 1
  END SUBROUTINE DARCY_EQUATION_EQUATIONS_SET_SETUP

  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matrices and RHS for a Darcy equation finite element equations set.
  SUBROUTINE DARCY_EQUATION_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ELEMENT_NUMBER,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equations_SET !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    INTEGER(INTG) :: FIELD_VAR_TYPE,ng,mh,mhs,mi,ms,nh,nhs,ni,ns,idxdim,num_var_count,idx_tensor
    INTEGER(INTG) :: my_compartment,Ncompartments,imatrix
    INTEGER(INTG) :: component_idx,xi_idx,derivative_idx
    INTEGER(INTG) :: MESH_COMPONENT_NUMBER, global_element_idx
    INTEGER(INTG) :: MESH_COMPONENT_1, MESH_COMPONENT_2
    INTEGER(INTG) :: NDOFS, NUMBER_OF_VEL_PRESS_COMPONENTS
    INTEGER(INTG) :: FIELD_VAR_TYPES(99)
    INTEGER(INTG) :: equations_SET_SUBTYPE
    INTEGER(INTG), POINTER :: equations_SET_FIELD_DATA(:)
    REAL(DP) :: RWG,SUM,PGMSI(3),PGNSI(3),PGM,PGN,COUPLING_PARAM
    TYPE(BASIS_TYPE), POINTER :: DEPENDENT_BASIS,GEOMETRIC_BASIS
    TYPE(BASIS_TYPE), POINTER :: DEPENDENT_BASIS_1, DEPENDENT_BASIS_2
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatrixType), POINTER :: stiffnessMatrix, dampingMatrix
    TYPE(FIELD_TYPE), POINTER :: dependentField,geometricField,materialsField,EQUATIONS_SET_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(MESH_ELEMENT_TYPE), POINTER :: MESH_ELEMENT
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector
    TYPE(FIELD_TYPE), POINTER :: sourceField
    TYPE(FieldVariablePtrType) :: FIELD_VARIABLES(99)
    TYPE(EquationsMatrixPtrType) :: COUPLING_MATRICES(99)
    REAL(DP), ALLOCATABLE :: PRESSURE_COEFF(:),PRESSURE(:),GRAD_PRESSURE(:,:)

    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: QUADRATURE_SCHEME
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: QUADRATURE_SCHEME_1, QUADRATURE_SCHEME_2
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: GEOMETRIC_INTERPOLATED_POINT,MATERIALS_INTERPOLATED_POINT
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: REFERENCE_GEOMETRIC_INTERPOLATED_POINT
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: ELASTICITY_DEPENDENT_INTERPOLATION_PARAMETERS
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: ELASTICITY_DEPENDENT_INTERPOLATED_POINT
    TYPE(VARYING_STRING) :: localError

    REAL(DP):: SOURCE,INTER_COMP_SOURCE,INTER_COMP_PERM_1,INTER_COMP_PERM_2
    REAL(DP):: BETA_PARAM, P_SINK_PARAM

    REAL(DP):: PERM_OVER_VIS_PARAM, DARCY_RHO_0_F
    REAL(DP):: PERM_TENSOR_OVER_VIS(3,3), VIS_OVER_PERM_TENSOR(3,3), Jmat
    REAL(DP):: X(3), ARG(3), L, FACT
    REAL(DP):: LM_PRESSURE,GRAD_LM_PRESSURE(3)

    REAL(DP):: DXDY(3,3), DXDXI(3,3), DYDXI(3,3), DXIDY(3,3)
    REAL(DP):: Jxy, Jyxi

    REAL(DP):: Mfact, bfact, p0fact

    REAL(DP):: ffact !f(Jxy) of the INRIA model
    REAL(DP):: dfdJfact !dfdJfact = f'(Jxy) of the INRIA model

    REAL(DP):: C_PARAM

    LOGICAL :: STABILIZED


    !--- Parameter settings concerning the Finite Element implementation
    STABILIZED = .TRUE.
    DARCY%LENGTH = 10.0_DP
    L = DARCY%LENGTH

    !--- testcase: default
    DARCY%TESTCASE = 0
    DARCY%ANALYTIC = .FALSE.

    ENTERS("DARCY_EQUATION_FINITE_ELEMENT_CALCULATE",err,error,*999)

    !Parameters settings for coupled elasticity Darcy INRIA model:
    CALL GET_DARCY_FINITE_ELASTICITY_PARAMETERS(DARCY_RHO_0_F,Mfact,bfact,p0fact,err,error,*999)

    NULLIFY(DEPENDENT_BASIS,GEOMETRIC_BASIS)
    NULLIFY(DEPENDENT_BASIS_1, DEPENDENT_BASIS_2)
    NULLIFY(EQUATIONS)
    NULLIFY(vectorMapping)
    NULLIFY(linearMapping)
    NULLIFY(dynamicMapping)
    NULLIFY(vectorMatrices)
    NULLIFY(linearMatrices)
    NULLIFY(dynamicMatrices)
    NULLIFY(rhsVector)
    NULLIFY(stiffnessMatrix, dampingMatrix)
    NULLIFY(dependentField,geometricField,materialsField)
    NULLIFY(FIELD_VARIABLE)
    NULLIFY(QUADRATURE_SCHEME)
    NULLIFY(QUADRATURE_SCHEME_1, QUADRATURE_SCHEME_2)
    NULLIFY(GEOMETRIC_INTERPOLATED_POINT,MATERIALS_INTERPOLATED_POINT)
    NULLIFY(REFERENCE_GEOMETRIC_INTERPOLATED_POINT)
    NULLIFY(ELASTICITY_DEPENDENT_INTERPOLATION_PARAMETERS)
    NULLIFY(ELASTICITY_DEPENDENT_INTERPOLATED_POINT)
    NULLIFY(DECOMPOSITION,MESH_ELEMENT)
    NULLIFY(BOUNDARY_CONDITIONS,BOUNDARY_CONDITIONS_VARIABLE)
    NULLIFY(sourceVector,sourceField)
    NULLIFY(EQUATIONS_SET_FIELD_DATA)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        NULLIFY(vectorEquations)
        CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
        IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
          CALL FlagError("Equations set specification is not allocated.",err,error,*999)
        ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
          CALL FlagError("Equations set specification must have three entries for a Darcy type equations set.", &
            & err,error,*999)
        END IF
        EQUATIONS_SET_SUBTYPE=EQUATIONS_SET%SPECIFICATION(3)
        SELECT CASE(EQUATIONS_SET_SUBTYPE)
        CASE(EQUATIONS_SET_STANDARD_DARCY_SUBTYPE, EQUATIONS_SET_QUASISTATIC_DARCY_SUBTYPE, EQUATIONS_SET_ALE_DARCY_SUBTYPE, &
          & EQUATIONS_SET_TRANSIENT_DARCY_SUBTYPE,EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE, &
          & EQUATIONS_SET_TRANSIENT_ALE_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE,EQUATIONS_SET_MULTI_COMPARTMENT_DARCY_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
!!TODO: move these and scale factor adjustment out once generalised Darcy is put in.
          !Store all these in equations matrices/somewhere else?????
          dependentField=>equations%interpolation%dependentField
          geometricField=>equations%interpolation%geometricField
          materialsField=>equations%interpolation%materialsField
          IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE .OR. &
            & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE) THEN
            sourceField=>equations%interpolation%sourceField
          END IF

          vectorMatrices=>vectorEquations%vectorMatrices
          rhsVector=>vectorMatrices%rhsVector
          vectorMapping=>vectorEquations%vectorMapping

          IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE .OR. &
            & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE) THEN
            sourceVector=>vectorMatrices%sourceVector
            sourceVector%elementVector%vector = 0.0_DP
          END IF

          SELECT CASE(EQUATIONS_SET_SUBTYPE)
          CASE(EQUATIONS_SET_STANDARD_DARCY_SUBTYPE, EQUATIONS_SET_QUASISTATIC_DARCY_SUBTYPE, EQUATIONS_SET_ALE_DARCY_SUBTYPE, &
            EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE)
            linearMatrices=>vectorMatrices%linearMatrices
            stiffnessMatrix=>linearMatrices%matrices(1)%ptr

            linearMapping=>vectorMapping%linearMapping
            FIELD_VARIABLE=>linearMapping%equationsMatrixToVarMaps(1)%VARIABLE
            FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE

            stiffnessMatrix%elementMatrix%matrix=0.0_DP

          CASE(EQUATIONS_SET_TRANSIENT_DARCY_SUBTYPE,EQUATIONS_SET_TRANSIENT_ALE_DARCY_SUBTYPE, &
            & EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
            & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE)
            dynamicMatrices=>vectorMatrices%dynamicMatrices
            stiffnessMatrix=>dynamicMatrices%matrices(1)%ptr
            dampingMatrix=>dynamicMatrices%matrices(2)%ptr

            dynamicMapping=>vectorMapping%dynamicMapping
            FIELD_VARIABLE=>dynamicMapping%equationsMatrixToVarMaps(1)%VARIABLE
            FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE

            stiffnessMatrix%elementMatrix%matrix=0.0_DP
            dampingMatrix%elementMatrix%matrix=0.0_DP

            !Stuff used to check if this element is on the mesh boundary
            DECOMPOSITION => dependentField%DECOMPOSITION
            MESH_COMPONENT_NUMBER = DECOMPOSITION%MESH_COMPONENT_NUMBER
            global_element_idx = DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%ptr%MAPPINGS%ELEMENTS% &
              & LOCAL_TO_GLOBAL_MAP(ELEMENT_NUMBER)
            MESH_ELEMENT => DECOMPOSITION%MESH%TOPOLOGY(MESH_COMPONENT_NUMBER)%ptr%ELEMENTS%ELEMENTS(global_element_idx)

          CASE(EQUATIONS_SET_MULTI_COMPARTMENT_DARCY_SUBTYPE)
            EQUATIONS_SET_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD
            CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET_FIELD,FIELD_U_VARIABLE_TYPE, &
              & FIELD_VALUES_SET_TYPE,EQUATIONS_SET_FIELD_DATA,err,error,*999)

            my_compartment = EQUATIONS_SET_FIELD_DATA(1)
            Ncompartments  = EQUATIONS_SET_FIELD_DATA(2)
            linearMatrices=>vectorMatrices%linearMatrices
            stiffnessMatrix=>linearMatrices%matrices(1)%ptr

            linearMapping=>vectorMapping%linearMapping
            FIELD_VARIABLE=>linearMapping%equationsMatrixToVarMaps(1)%VARIABLE
            FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE

            stiffnessMatrix%elementMatrix%matrix=0.0_DP

          CASE(EQUATIONS_SET_ELASTICITY_MULTI_COMPARTMENT_DARCY_INRIA_SUBTYPE)

            EQUATIONS_SET_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD
            CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET_FIELD,FIELD_U_VARIABLE_TYPE, &
              & FIELD_VALUES_SET_TYPE,EQUATIONS_SET_FIELD_DATA,err,error,*999)

            my_compartment = EQUATIONS_SET_FIELD_DATA(1)
            Ncompartments  = EQUATIONS_SET_FIELD_DATA(2)

            !if Ncompartments>99 then flag error

            linearMatrices=>vectorMatrices%linearMatrices
            linearMapping=>vectorMapping%linearMapping

!             DO imatrix = 1,Ncompartments
!               COUPLING_MATRICES(imatrix)%ptr=>linearMatrices%matrices(imatrix)%ptr
!               FIELD_VARIABLES(imatrix)%ptr=>linearMapping%equationsMatrixToVarMaps(imatrix)%VARIABLE
!               FIELD_VAR_TYPES(imatrix)=FIELD_VARIABLES(imatrix)%ptr%VARIABLE_TYPE
!               COUPLING_MATRICES(imatrix)%ptr%elementMatrix%matrix=0.0_DP
!             END DO

            dynamicMatrices=>vectorMatrices%dynamicMatrices
            stiffnessMatrix=>dynamicMatrices%matrices(1)%ptr
            dampingMatrix=>dynamicMatrices%matrices(2)%ptr



            dynamicMapping=>vectorMapping%dynamicMapping
            FIELD_VARIABLE=>dynamicMapping%equationsMatrixToVarMaps(1)%VARIABLE
            FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE

            stiffnessMatrix%elementMatrix%matrix=0.0_DP
            dampingMatrix%elementMatrix%matrix=0.0_DP
          CASE(EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
            EQUATIONS_SET_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD
            CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET_FIELD,FIELD_U_VARIABLE_TYPE, &
              & FIELD_VALUES_SET_TYPE,EQUATIONS_SET_FIELD_DATA,err,error,*999)

            my_compartment = EQUATIONS_SET_FIELD_DATA(1)
            Ncompartments  = EQUATIONS_SET_FIELD_DATA(2)
            !These linear matrices are actually only required if we are coupling the momentum terms too
            !If it is just a mass coupling, then all of the additional terms are placed in the RHS of the mass-increase equation
            linearMatrices=>vectorMatrices%linearMatrices
            linearMapping=>vectorMapping%linearMapping

            num_var_count=0
            DO imatrix = 1,Ncompartments
             IF(imatrix/=my_compartment)THEN
              num_var_count=num_var_count+1
              COUPLING_MATRICES(num_var_count)%ptr=>linearMatrices%matrices(num_var_count)%ptr
              FIELD_VARIABLES(num_var_count)%ptr=>linearMapping%equationsMatrixToVarMaps(num_var_count)%VARIABLE
              FIELD_VAR_TYPES(num_var_count)=FIELD_VARIABLES(num_var_count)%ptr%VARIABLE_TYPE
              COUPLING_MATRICES(num_var_count)%ptr%elementMatrix%matrix=0.0_DP
             ENDIF
            END DO

            dynamicMatrices=>vectorMatrices%dynamicMatrices
            stiffnessMatrix=>dynamicMatrices%matrices(1)%ptr
            dampingMatrix=>dynamicMatrices%matrices(2)%ptr

            dynamicMapping=>vectorMapping%dynamicMapping
            FIELD_VARIABLE=>dynamicMapping%equationsMatrixToVarMaps(1)%VARIABLE
            FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE

            stiffnessMatrix%elementMatrix%matrix=0.0_DP
            dampingMatrix%elementMatrix%matrix=0.0_DP

              ALLOCATE(PRESSURE_COEFF(Ncompartments))

              ALLOCATE(PRESSURE(Ncompartments))
              ALLOCATE(GRAD_PRESSURE(3,Ncompartments))
              PRESSURE = 0.0_DP
              GRAD_PRESSURE = 0.0_DP
              PRESSURE_COEFF(1)=0.25_DP
              PRESSURE_COEFF(2)=0.25_DP
              PRESSURE_COEFF(3)=0.25_DP
              PRESSURE_COEFF(4)=0.25_DP
          END SELECT

          !\ToDo: DEPENDENT_BASIS, DEPENDENT_BASIS_1, DEPENDENT_BASIS_2 - consistency !!!

          GEOMETRIC_BASIS=>geometricField%DECOMPOSITION%DOMAIN(geometricField%DECOMPOSITION%MESH_COMPONENT_NUMBER)%ptr% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          DEPENDENT_BASIS=>dependentField%DECOMPOSITION%DOMAIN(dependentField%DECOMPOSITION%MESH_COMPONENT_NUMBER)%ptr% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS

          QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr

          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
            & geometricInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
            & materialsInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
          IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE) THEN
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
              & materialsInterpParameters(FIELD_V_VARIABLE_TYPE)%ptr,err,error,*999)
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
              & materialsInterpParameters(FIELD_U1_VARIABLE_TYPE)%ptr,err,error,*999)
          ENDIF

          IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE .OR. &
            & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE) THEN
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
              & sourceInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
          END IF

          IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE .OR. &
           & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE) THEN
            ELASTICITY_DEPENDENT_INTERPOLATION_PARAMETERS=>equations%interpolation% &
              & dependentInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER, &
              & ELASTICITY_DEPENDENT_INTERPOLATION_PARAMETERS,err,error,*999)
            ELASTICITY_DEPENDENT_INTERPOLATED_POINT=>equations%interpolation% &
              & dependentInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr
          ENDIF


          SELECT CASE(EQUATIONS_SET_SUBTYPE)
          CASE(EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE,EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE, &
             & EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
            NUMBER_OF_VEL_PRESS_COMPONENTS = FIELD_VARIABLE%NUMBER_OF_COMPONENTS - 1  !last component: mass increase
          CASE DEFAULT
            NUMBER_OF_VEL_PRESS_COMPONENTS = FIELD_VARIABLE%NUMBER_OF_COMPONENTS
          END SELECT

          !---------------------------------------------------------------------------------------------------------
          !Invoke penalty term to enforce impermeable BC
          !  should only be executed if THIS element lies on the surface
          !  (within the routine we check whether the element nodes have actually been set impermeable)
          SELECT CASE(EQUATIONS_SET_SUBTYPE)
          CASE(EQUATIONS_SET_TRANSIENT_DARCY_SUBTYPE,EQUATIONS_SET_TRANSIENT_ALE_DARCY_SUBTYPE, &
            & EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
            & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE)
            IF( MESH_ELEMENT%BOUNDARY_ELEMENT ) THEN
              CALL DARCY_EQUATION_IMPERMEABLE_BC_VIA_PENALTY(EQUATIONS_SET,ELEMENT_NUMBER,err,error,*999)
            ENDIF
          END SELECT
          !---------------------------------------------------------------------------------------------------------

          !--- Loop over gauss points
          !    Given that also materials field is interpolated, ensure sufficient number of Gauss points !!!
          DO ng=1,QUADRATURE_SCHEME%NUMBER_OF_GAUSS


            IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE.OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE) THEN
              !------------------------------------------------------------------------------
              !--- begin: Compute the Jacobian of the mapping

              !--- Interpolation of Reference Geometry
              CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_INITIAL_VALUES_SET_TYPE,ELEMENT_NUMBER, &
                & equations%interpolation%geometricInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
              REFERENCE_GEOMETRIC_INTERPOLATED_POINT => equations%interpolation%geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr
              CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng, &
                & REFERENCE_GEOMETRIC_INTERPOLATED_POINT,err,error,*999)
              !--- Retrieve local map DYDXI
              DO component_idx=1,DEPENDENT_BASIS%NUMBER_OF_XI
                DO xi_idx=1,DEPENDENT_BASIS%NUMBER_OF_XI
                  derivative_idx=PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xi_idx) !2,4,7
                  DYDXI(component_idx,xi_idx)=REFERENCE_GEOMETRIC_INTERPOLATED_POINT%VALUES(component_idx,derivative_idx) !dy/dxi (y = referential)
                ENDDO
              ENDDO

              !--- Interpolation of (actual) Geometry and Metrics
              CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER, &
                & equations%interpolation%geometricInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
              GEOMETRIC_INTERPOLATED_POINT => equations%interpolation%geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr
              CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng, &
                & GEOMETRIC_INTERPOLATED_POINT,err,error,*999)
              CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI, &
                & equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
              !--- Retrieve local map DXDXI
              DO component_idx=1,DEPENDENT_BASIS%NUMBER_OF_XI
                DO xi_idx=1,DEPENDENT_BASIS%NUMBER_OF_XI
                  derivative_idx=PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xi_idx) !2,4,7
                  DXDXI(component_idx,xi_idx)=GEOMETRIC_INTERPOLATED_POINT%VALUES(component_idx,derivative_idx) !dx/dxi
                ENDDO
              ENDDO

              !--- Compute deformation gradient tensor DXDY and its Jacobian Jxy
              CALL Invert(DYDXI,DXIDY,Jyxi,err,error,*999) !dy/dxi -> dxi/dy
              CALL MatrixProduct(DXDXI,DXIDY,DXDY,err,error,*999) !dx/dxi * dxi/dy = dx/dy (deformation gradient tensor, F)
              CALL Determinant(DXDY,Jxy,err,error,*999)

              IF( ABS(Jxy) < 1.0E-10_DP ) THEN
                localError="DARCY_EQUATION_FINITE_ELEMENT_CALCULATE: Jacobian Jxy is smaller than 1.0E-10_DP."
                CALL FlagError(localError,err,error,*999)
              END IF

              !ffact = f(Jxy) of the INRIA model, dfdJfact is not relevant here
              CALL EVALUATE_CHAPELLE_FUNCTION(Jxy,ffact,dfdJfact,err,error,*999)

              !--- end: Compute the Jacobian of the mapping
              !------------------------------------------------------------------------------
            END IF

            !--- Interpolate geometric and mesh velocity field (if applicable)
            GEOMETRIC_INTERPOLATED_POINT => equations%interpolation%geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng, &
              & GEOMETRIC_INTERPOLATED_POINT,err,error,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI, &
              & equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)

            !--- Calculate 'geometricInterpParameters' from 'FIELD_VALUES_SET_TYPE'
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER, &
              & equations%interpolation%geometricInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng, &
              & GEOMETRIC_INTERPOLATED_POINT,err,error,*999)


!             !--- Material Settings ---!
!             !*** If material is variable, need to account for this in deriving the variational statement ***!


            !--- Interpolate materials field
            !Get the Darcy permeability
            MATERIALS_INTERPOLATED_POINT => equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng, &
              & MATERIALS_INTERPOLATED_POINT,err,error,*999)
            IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE) THEN
            !Get the intercompartmental permeabilities
              CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
                & materialsInterpPoint(FIELD_V_VARIABLE_TYPE)%ptr,err,error,*999)
            !Get the material parameters for the constitutive law for each Darcy compartment (for determining the partial pressures)
              CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
                & materialsInterpPoint(FIELD_U1_VARIABLE_TYPE)%ptr,err,error,*999)
            ENDIF

            IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE .OR. &
              & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE) THEN
              CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
                & sourceInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            END IF

            SELECT CASE(EQUATIONS_SET_SUBTYPE)
            CASE(EQUATIONS_SET_STANDARD_DARCY_SUBTYPE,EQUATIONS_SET_QUASISTATIC_DARCY_SUBTYPE,EQUATIONS_SET_ALE_DARCY_SUBTYPE)
              !scalar permeability/viscosity
              PERM_TENSOR_OVER_VIS=0.0_DP
              PERM_TENSOR_OVER_VIS(1,1) = MATERIALS_INTERPOLATED_POINT%VALUES(2,NO_PART_DERIV)
              PERM_TENSOR_OVER_VIS(2,2) = MATERIALS_INTERPOLATED_POINT%VALUES(2,NO_PART_DERIV)
              PERM_TENSOR_OVER_VIS(3,3) = MATERIALS_INTERPOLATED_POINT%VALUES(2,NO_PART_DERIV)
              !Multiply by porosity
              PERM_TENSOR_OVER_VIS=PERM_TENSOR_OVER_VIS*MATERIALS_INTERPOLATED_POINT%VALUES(1,NO_PART_DERIV)
            CASE DEFAULT
              !symmetric permeability/viscosity tensor
              PERM_TENSOR_OVER_VIS(1,1) = MATERIALS_INTERPOLATED_POINT%VALUES(2,NO_PART_DERIV)
              PERM_TENSOR_OVER_VIS(1,2) = MATERIALS_INTERPOLATED_POINT%VALUES(3,NO_PART_DERIV)
              PERM_TENSOR_OVER_VIS(1,3) = MATERIALS_INTERPOLATED_POINT%VALUES(4,NO_PART_DERIV)
              PERM_TENSOR_OVER_VIS(2,2) = MATERIALS_INTERPOLATED_POINT%VALUES(5,NO_PART_DERIV)
              PERM_TENSOR_OVER_VIS(2,3) = MATERIALS_INTERPOLATED_POINT%VALUES(6,NO_PART_DERIV)
              PERM_TENSOR_OVER_VIS(3,3) = MATERIALS_INTERPOLATED_POINT%VALUES(7,NO_PART_DERIV)

              PERM_TENSOR_OVER_VIS(2,1) = PERM_TENSOR_OVER_VIS(1,2)
              PERM_TENSOR_OVER_VIS(3,1) = PERM_TENSOR_OVER_VIS(1,3)
              PERM_TENSOR_OVER_VIS(3,2) = PERM_TENSOR_OVER_VIS(2,3)
            END SELECT

            IF(DIAGNOSTICS3) THEN
              IF(idebug2) THEN
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"MATERIALS_INTERPOLATED_POINT%VALUES(1,NO_PART_DERIV) = ", &
                  & MATERIALS_INTERPOLATED_POINT%VALUES(1,NO_PART_DERIV),err,error,*999)
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"MATERIALS_INTERPOLATED_POINT%VALUES(2,NO_PART_DERIV) = ", &
                  & MATERIALS_INTERPOLATED_POINT%VALUES(2,NO_PART_DERIV),err,error,*999)
                CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"",err,error,*999)
                idebug2 = .FALSE.
              ENDIF
            ENDIF

            CALL Determinant(PERM_TENSOR_OVER_VIS,Jmat,err,error,*999)
            IF(Jmat>ZERO_TOLERANCE) THEN
              CALL INVERT(PERM_TENSOR_OVER_VIS,VIS_OVER_PERM_TENSOR,Jmat,err,error,*999)
            ELSE
              VIS_OVER_PERM_TENSOR = 0.0_DP
              DO idx_tensor=1,3
                VIS_OVER_PERM_TENSOR(idx_tensor,idx_tensor) = 1.0e10_DP
              END DO
!               CALL WRITE_STRING(GENERAL_OUTPUT_TYPE, &
!                 & "WARNING: Jmat<ZERO_TOLERANCE - Thus setting VIS_OVER_PERM_TENSOR(i,i) = 1.0e10_DP",err,error,*999)
            END IF


            !Two parameters that are used only for TESTCASE==3: VenousCompartment problem: Exclude this, too specific ???
            BETA_PARAM   = - DARCY%PERM_OVER_VIS * (2.0_DP * PI / DARCY%LENGTH) * (2.0_DP * PI / DARCY%LENGTH)
            P_SINK_PARAM = DARCY%P_SINK


            IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE) THEN
              CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng, &
                & ELASTICITY_DEPENDENT_INTERPOLATED_POINT,err,error,*999)
              !Mind the sign !!!
              !The minus sign derives from the convention of using "+ P * Jznu * AZU(i,j)"
              ! in the constitutive law in FINITE_ELASTICITY_GAUSS_CAUCHY_TENSOR
              LM_PRESSURE = -ELASTICITY_DEPENDENT_INTERPOLATED_POINT%VALUES(4,NO_PART_DERIV)
              DO xi_idx=1,DEPENDENT_BASIS%NUMBER_OF_XI
                derivative_idx=PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xi_idx) !2,4,7
                !gradient wrt. element coordinates xi
                GRAD_LM_PRESSURE(xi_idx) = -ELASTICITY_DEPENDENT_INTERPOLATED_POINT%VALUES(4,derivative_idx)
              ENDDO
            ENDIF

            !For multi-compartment model - determine pressure from partial derivative of constitutive law
            IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE) THEN
              write(*,*) 'NEED CONSTITUTIVE LAWS HERE!!!! THE FOLLOWING IS PLACEHOLDER ONLY!'
              !BEGIN PLACEHOLDER
              CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng, &
                & ELASTICITY_DEPENDENT_INTERPOLATED_POINT,err,error,*999)
              !Mind the sign !!!
              LM_PRESSURE = -ELASTICITY_DEPENDENT_INTERPOLATED_POINT%VALUES(4,NO_PART_DERIV)
              DO xi_idx=1,DEPENDENT_BASIS%NUMBER_OF_XI
                derivative_idx=PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xi_idx) !2,4,7
                !gradient wrt. element coordinates xi
                GRAD_LM_PRESSURE(xi_idx) = -ELASTICITY_DEPENDENT_INTERPOLATED_POINT%VALUES(4,derivative_idx)
              ENDDO
              !loop over compartments to determine the pressure in each one - this could be quite inefficient, as it will be calculated several times over
              !unless calculate the pressures in a pre-solve and store them in extra components/variables of the dependent field
              !these pressures should really be known immediately after the finite elasticity solve and not determined here
              !END PLACEHOLDER
              !The following pressure_coeff matrix is just for testing purposes and ultimately will be replaced with functions and materials field parameters (for present, sum of
              !coefficients should be 1).

              DO imatrix=1,Ncompartments

                PRESSURE(imatrix) =  PRESSURE_COEFF(imatrix)*LM_PRESSURE
                DO xi_idx=1,DEPENDENT_BASIS%NUMBER_OF_XI
                  GRAD_PRESSURE(xi_idx,imatrix) = PRESSURE_COEFF(imatrix)*GRAD_LM_PRESSURE(xi_idx)
                ENDDO
              ENDDO
            ENDIF


!             RWG = equations%interpolation%geometricInterpPointMetrics%JACOBIAN * QUADRATURE_SCHEME%GAUSS_WEIGHTS(ng)

            !Loop over element rows
            mhs=0
            DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS

              MESH_COMPONENT_1 = FIELD_VARIABLE%COMPONENTS(mh)%MESH_COMPONENT_NUMBER
              DEPENDENT_BASIS_1 => dependentField%DECOMPOSITION%DOMAIN(MESH_COMPONENT_1)%ptr% &
                & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
              QUADRATURE_SCHEME_1 => DEPENDENT_BASIS_1%QUADRATURE% &
                & QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
              RWG = equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%JACOBIAN * &
                & QUADRATURE_SCHEME_1%GAUSS_WEIGHTS(ng)

              DO ms=1,DEPENDENT_BASIS_1%NUMBER_OF_ELEMENT_PARAMETERS
                mhs=mhs+1

                !===================================================================================================================
                !stiffnessMatrix
                IF(stiffnessMatrix%updateMatrix) THEN

                  !Loop over element columns
                  nhs=0
                  DO nh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS

                    MESH_COMPONENT_2 = FIELD_VARIABLE%COMPONENTS(nh)%MESH_COMPONENT_NUMBER
                    DEPENDENT_BASIS_2 => dependentField%DECOMPOSITION%DOMAIN(MESH_COMPONENT_2)%ptr% &
                      & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
                    !--- We cannot use two different quadrature schemes here !!!
                    QUADRATURE_SCHEME_2 => DEPENDENT_BASIS_2%QUADRATURE% &
                     & QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
                    !RWG = equations%interpolation%geometricInterpPointMetrics%JACOBIAN * &
                    !  & QUADRATURE_SCHEME_2%GAUSS_WEIGHTS(ng)

                    DO ns=1,DEPENDENT_BASIS_2%NUMBER_OF_ELEMENT_PARAMETERS
                      nhs=nhs+1

                      SELECT CASE(EQUATIONS_SET_SUBTYPE)
                      !====================================================================================================
                      !  i n c o m p r e s s i b l e   e l a s t i c i t y   d r i v e n   D a r c y   :   M A T R I C E S
                      CASE(EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE)
                        !-------------------------------------------------------------------------------------------------------------
                        !velocity test function, velocity trial function
                        IF(mh==nh.AND.nh<FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN

                          SUM = 0.0_DP

                          PGM=QUADRATURE_SCHEME_1%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                          PGN=QUADRATURE_SCHEME_2%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)

                          SUM = SUM + VIS_OVER_PERM_TENSOR( mh, nh ) * PGM * PGN
                          !MIND: double check the matrix index order: (mh, nh) or (nh, mh)
                          !within this conditional: mh==nh anyway

                          stiffnessMatrix%elementMatrix%matrix(mhs,nhs) = stiffnessMatrix%elementMatrix%matrix(mhs,nhs) + &
                            & SUM * RWG
                        !-------------------------------------------------------------------------------------------------------------
                        !mass-increase test function, velocity trial function
                        ELSE IF(mh==FIELD_VARIABLE%NUMBER_OF_COMPONENTS.AND.nh<FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN

                          SUM = 0.0_DP

                          PGM=QUADRATURE_SCHEME_1%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                          PGN=QUADRATURE_SCHEME_2%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)

                          DO mi=1,DEPENDENT_BASIS_1%NUMBER_OF_XI
                            PGMSI(mi)=QUADRATURE_SCHEME_1%GAUSS_BASIS_FNS(ms,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(mi),ng)
                          ENDDO !mi

                          DO ni=1,DEPENDENT_BASIS_2%NUMBER_OF_XI
                            PGNSI(ni)=QUADRATURE_SCHEME_2%GAUSS_BASIS_FNS(ns,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)
                          ENDDO !ni

                          DO ni=1,DEPENDENT_BASIS_2%NUMBER_OF_XI
                            SUM = SUM + PGM * PGNSI(ni) * &
                              & equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%DXI_DX(ni,nh)
                          ENDDO !ni

                          stiffnessMatrix%elementMatrix%matrix(mhs,nhs) = stiffnessMatrix%elementMatrix%matrix(mhs,nhs) + &
                            & SUM * RWG
                        ENDIF
                        !=====================================================================================================================
                        !dampingMatrix
                        IF(dampingMatrix%updateMatrix) THEN
                          !MASS-INCREASE test function, mass-increase trial function
                          IF(mh==nh.AND.nh==FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
                            PGM=QUADRATURE_SCHEME_1%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                            PGN=QUADRATURE_SCHEME_2%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)

                            SUM = 0.0_DP

                            !To integrate the mass-increase term in the reference configuration, we divide by Jxy.
                            SUM = PGM * PGN / (Jxy * DARCY_RHO_0_F)

                            dampingMatrix%elementMatrix%matrix(mhs,nhs) = dampingMatrix%elementMatrix%matrix(mhs,nhs) + &
                              & SUM * RWG
                          END IF

!                           !Try out adding the inertia term ...
!                           IF(mh==nh.AND.mh<FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
!                             PGM=QUADRATURE_SCHEME_1%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
!                             PGN=QUADRATURE_SCHEME_2%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)
!
!                             SUM = 0.0_DP
!
!                             SUM = PGM*PGN*DARCY_RHO_0_F
!
!                             dampingMatrix%elementMatrix%matrix(mhs,nhs) = dampingMatrix%elementMatrix%matrix(mhs,nhs) + &
!                               & SUM * RWG
!                           END IF

                        END IF

                      ! matrices for multi-compartment poroelastic equations
                      CASE(EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
                        !velocity test function, velocity trial function
                        IF(mh==nh.AND.nh<FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN

                          SUM = 0.0_DP

                          PGM=QUADRATURE_SCHEME_1%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                          PGN=QUADRATURE_SCHEME_2%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)

                          SUM = SUM + VIS_OVER_PERM_TENSOR( mh, nh ) * PGM * PGN
                          !MIND: double check the matrix index order: (mh, nh) or (nh, mh)
                          !within this conditional: mh==nh anyway

                          stiffnessMatrix%elementMatrix%matrix(mhs,nhs) = stiffnessMatrix%elementMatrix%matrix(mhs,nhs) + &
                            & SUM * RWG
                        !-------------------------------------------------------------------------------------------------------------
                        !mass-increase test function, velocity trial function
                        ELSE IF(mh==FIELD_VARIABLE%NUMBER_OF_COMPONENTS.AND.nh<FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN

                          SUM = 0.0_DP

                          PGM=QUADRATURE_SCHEME_1%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                          PGN=QUADRATURE_SCHEME_2%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)

                          DO mi=1,DEPENDENT_BASIS_1%NUMBER_OF_XI
                            PGMSI(mi)=QUADRATURE_SCHEME_1%GAUSS_BASIS_FNS(ms,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(mi),ng)
                          ENDDO !mi

                          DO ni=1,DEPENDENT_BASIS_2%NUMBER_OF_XI
                            PGNSI(ni)=QUADRATURE_SCHEME_2%GAUSS_BASIS_FNS(ns,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)
                          ENDDO !ni

                          DO ni=1,DEPENDENT_BASIS_2%NUMBER_OF_XI
                            SUM = SUM + PGM * PGNSI(ni) * &
                              & equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%DXI_DX(ni,nh)
                          ENDDO !ni

                          stiffnessMatrix%elementMatrix%matrix(mhs,nhs) = stiffnessMatrix%elementMatrix%matrix(mhs,nhs) + &
                            & SUM * RWG
                        ENDIF
                        !=====================================================================================================================
                        !dampingMatrix
                        IF(dampingMatrix%updateMatrix) THEN
                          !MASS-INCREASE test function, mass-increase trial function
                          IF(mh==nh.AND.nh==FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
                            PGM=QUADRATURE_SCHEME_1%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                            PGN=QUADRATURE_SCHEME_2%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)

                            SUM = 0.0_DP

                            !To integrate the mass-increase term in the reference configuration, we divide by Jxy.
                            SUM = PGM * PGN / (Jxy * DARCY_RHO_0_F)

                            dampingMatrix%elementMatrix%matrix(mhs,nhs) = dampingMatrix%elementMatrix%matrix(mhs,nhs) + &
                              & SUM * RWG
                          END IF

!                           !Try out adding the inertia term ...
!                           IF(mh==nh.AND.mh<FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
!                             PGM=QUADRATURE_SCHEME_1%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
!                             PGN=QUADRATURE_SCHEME_2%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)
!
!                             SUM = 0.0_DP
!
!                             SUM = PGM*PGN*DARCY_RHO_0_F
!
!                             dampingMatrix%elementMatrix%matrix(mhs,nhs) = dampingMatrix%elementMatrix%matrix(mhs,nhs) + &
!                               & SUM * RWG
!                           END IF

                        END IF


                      !=================================================================================
                      !    d e f a u l t   :   M A T R I C E S
                      CASE DEFAULT
                        !-------------------------------------------------------------------------------------------------------------
                        !velocity test function, velocity trial function
                        IF(mh==nh.AND.nh<NUMBER_OF_VEL_PRESS_COMPONENTS) THEN

                          SUM = 0.0_DP

                          PGM=QUADRATURE_SCHEME_1%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                          PGN=QUADRATURE_SCHEME_2%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)

                          SUM = SUM + VIS_OVER_PERM_TENSOR( mh, nh ) * PGM * PGN
                          !MIND: double check the matrix index order: (mh, nh) or (nh, mh)
                          !within this conditional: mh==nh anyway

                          IF( STABILIZED ) THEN
                            SUM = SUM - 0.5_DP * VIS_OVER_PERM_TENSOR( mh, nh ) * PGM * PGN
                          END IF

                          stiffnessMatrix%elementMatrix%matrix(mhs,nhs) = stiffnessMatrix%elementMatrix%matrix(mhs,nhs) + &
                            & SUM * RWG

                        !-------------------------------------------------------------------------------------------------------------
                        !velocity test function, pressure trial function
                        ELSE IF(mh<NUMBER_OF_VEL_PRESS_COMPONENTS.AND.nh==NUMBER_OF_VEL_PRESS_COMPONENTS) THEN

                          SUM = 0.0_DP

                          PGM=QUADRATURE_SCHEME_1%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                          PGN=QUADRATURE_SCHEME_2%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)

                          DO mi=1,DEPENDENT_BASIS_1%NUMBER_OF_XI
                            PGMSI(mi)=QUADRATURE_SCHEME_1%GAUSS_BASIS_FNS(ms,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(mi),ng)
                          ENDDO !mi

                          DO ni=1,DEPENDENT_BASIS_2%NUMBER_OF_XI
                            PGNSI(ni)=QUADRATURE_SCHEME_2%GAUSS_BASIS_FNS(ns,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)
                          ENDDO !ni

                          DO mi=1,DEPENDENT_BASIS_1%NUMBER_OF_XI
                            SUM = SUM - PGMSI(mi) * PGN * &
                              & equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%DXI_DX(mi,mh)
                          ENDDO !mi

                          IF( STABILIZED ) THEN
                            DO ni=1,DEPENDENT_BASIS_2%NUMBER_OF_XI
                              SUM = SUM - 0.5_DP * PGM * PGNSI(ni) * &
                                & equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%DXI_DX(ni,mh)
                            ENDDO !ni
                          END IF

                          stiffnessMatrix%elementMatrix%matrix(mhs,nhs) = stiffnessMatrix%elementMatrix%matrix(mhs,nhs) + &
                            & SUM * RWG

                        !-------------------------------------------------------------------------------------------------------------
                        !pressure test function, velocity trial function
                        ELSE IF(mh==NUMBER_OF_VEL_PRESS_COMPONENTS.AND.nh<NUMBER_OF_VEL_PRESS_COMPONENTS) THEN

                          SUM = 0.0_DP

                          PGM=QUADRATURE_SCHEME_1%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                          PGN=QUADRATURE_SCHEME_2%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)

                          DO mi=1,DEPENDENT_BASIS_1%NUMBER_OF_XI
                            PGMSI(mi)=QUADRATURE_SCHEME_1%GAUSS_BASIS_FNS(ms,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(mi),ng)
                          ENDDO !mi

                          DO ni=1,DEPENDENT_BASIS_2%NUMBER_OF_XI
                            PGNSI(ni)=QUADRATURE_SCHEME_2%GAUSS_BASIS_FNS(ns,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)
                          ENDDO !ni

                          DO ni=1,DEPENDENT_BASIS_2%NUMBER_OF_XI
                            SUM = SUM + PGM * PGNSI(ni) * &
                              & equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%DXI_DX(ni,nh)
                          ENDDO !ni

                          IF( STABILIZED ) THEN
                            DO mi=1,DEPENDENT_BASIS_1%NUMBER_OF_XI
                              SUM = SUM + 0.5_DP * PGMSI(mi) * PGN * &
                                & equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%DXI_DX(mi,nh)
                            ENDDO !mi
                          END IF

                          stiffnessMatrix%elementMatrix%matrix(mhs,nhs) = stiffnessMatrix%elementMatrix%matrix(mhs,nhs) + &
                            & SUM * RWG

                        !-------------------------------------------------------------------------------------------------------------
                        !pressure test function, pressure trial function
                        ELSE IF(mh==nh.AND.nh==NUMBER_OF_VEL_PRESS_COMPONENTS) THEN

                          SUM = 0.0_DP

                          DO mi=1,DEPENDENT_BASIS_1%NUMBER_OF_XI
                            PGMSI(mi)=QUADRATURE_SCHEME_1%GAUSS_BASIS_FNS(ms,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(mi),ng)
                          ENDDO !mi

                          DO ni=1,DEPENDENT_BASIS_2%NUMBER_OF_XI
                            PGNSI(ni)=QUADRATURE_SCHEME_2%GAUSS_BASIS_FNS(ns,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)
                          ENDDO !ni

                          IF( STABILIZED ) THEN
                            DO idxdim =1,DEPENDENT_BASIS_1%NUMBER_OF_XI !number space dimension equiv. number of xi
                              DO mi=1,DEPENDENT_BASIS_1%NUMBER_OF_XI
                                DO ni=1,DEPENDENT_BASIS_2%NUMBER_OF_XI
                                  SUM = SUM + 0.5_DP * PERM_TENSOR_OVER_VIS( idxdim, idxdim ) * PGMSI(mi) * PGNSI(ni) * &
                                    & equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr% &
                                    & DXI_DX(mi,idxdim) * equations%interpolation% &
                                    & geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%DXI_DX(ni,idxdim)
                                ENDDO !ni
                              ENDDO !mi
                            ENDDO !idxdim
                          END IF

                          IF( DARCY%TESTCASE == 3 ) THEN
                            !This forms part of the pressure-dependent source term,
                            !thus it enters the LHS
                            PGM=QUADRATURE_SCHEME_1%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                            PGN=QUADRATURE_SCHEME_2%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)

                            SUM = SUM + BETA_PARAM * PGM * PGN
                          END IF

                          stiffnessMatrix%elementMatrix%matrix(mhs,nhs) = stiffnessMatrix%elementMatrix%matrix(mhs,nhs) + &
                            & SUM * RWG

                        !-------------------------------------------------------------------------------------------------------------
                        !For the INRIA model, and: mass-increase test function, pressure trial function
                        ELSE IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE.AND. &
                          & mh==FIELD_VARIABLE%NUMBER_OF_COMPONENTS.AND.nh==NUMBER_OF_VEL_PRESS_COMPONENTS) THEN

                          SUM = 0.0_DP

                          PGM=QUADRATURE_SCHEME_1%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                          PGN=QUADRATURE_SCHEME_2%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)

                          SUM = SUM - PGM * PGN / (Mfact * ffact)

                          stiffnessMatrix%elementMatrix%matrix(mhs,nhs) = stiffnessMatrix%elementMatrix%matrix(mhs,nhs) + &
                            & SUM * RWG

                        !-------------------------------------------------------------------------------------------------------------
                        !For the INRIA model, and: mass-increase test function, mass-increase trial function
                        ELSE IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE.AND. &
                          & mh==nh.AND.nh==FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN

                          SUM = 0.0_DP

                          PGM=QUADRATURE_SCHEME_1%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                          PGN=QUADRATURE_SCHEME_2%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)

                          SUM = SUM + PGM * PGN / DARCY_RHO_0_F

                          stiffnessMatrix%elementMatrix%matrix(mhs,nhs) = stiffnessMatrix%elementMatrix%matrix(mhs,nhs) + &
                            & SUM * RWG

                        !-------------------------------------------------------------------------------------------------------------
                        ELSE

                          stiffnessMatrix%elementMatrix%matrix(mhs,nhs) = 0.0_DP

                        ENDIF

                        !=====================================================================================================================
                        ! dampingMatrix

                        IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_TRANSIENT_DARCY_SUBTYPE) THEN
                          IF(dampingMatrix%updateMatrix) THEN
                            IF(mh==nh.AND.mh<NUMBER_OF_VEL_PRESS_COMPONENTS) THEN
                              PGM=QUADRATURE_SCHEME_1%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                              PGN=QUADRATURE_SCHEME_2%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)

                              SUM = 0.0_DP

                              SUM = PGM*PGN*DARCY_RHO_0_F

                              dampingMatrix%elementMatrix%matrix(mhs,nhs) = dampingMatrix%elementMatrix%matrix(mhs,nhs) + &
                                & SUM * RWG
                            END IF
                          END IF
                        ELSE IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE) THEN
                          IF(dampingMatrix%updateMatrix) THEN
                            !pressure test function, mass-increase trial function
                            IF(mh==NUMBER_OF_VEL_PRESS_COMPONENTS.AND.nh==FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
                              PGM=QUADRATURE_SCHEME_1%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                              PGN=QUADRATURE_SCHEME_2%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)

                              SUM = 0.0_DP

                              SUM = PGM * PGN / (Jxy * DARCY_RHO_0_F)

                              dampingMatrix%elementMatrix%matrix(mhs,nhs) = dampingMatrix%elementMatrix%matrix(mhs,nhs) + &
                                & SUM * RWG
                            END IF
                          END IF
                        END IF

                      END SELECT
                      !   e n d   s e l e c t   EQUATIONS_SET_SUBTYPE
                      !=================================================================================

                    ENDDO !ns
                  ENDDO !nh
                ENDIF
                !===================================================================================================================
                !rhsVector
                IF(rhsVector%updateVector) THEN

                  SELECT CASE(EQUATIONS_SET_SUBTYPE)
                  !==========================================================================================
                  !  i n c o m p r e s s i b l e   e l a s t i c i t y   d r i v e n   D a r c y   :   R H S
                  CASE(EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE)

                    !-----------------------------------------------------------------------------------------------------------------
                    !velocity test function
                    IF( mh<FIELD_VARIABLE%NUMBER_OF_COMPONENTS ) THEN

                      SUM = 0.0_DP

                      PGM = QUADRATURE_SCHEME_1%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)

                      !Term arising from the pressure / Lagrange Multiplier of elasticity (given):
                      DO mi=1,DEPENDENT_BASIS_1%NUMBER_OF_XI
                        SUM = SUM - PGM * GRAD_LM_PRESSURE(mi) * &
                          & equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%DXI_DX(mi,mh)
                      ENDDO !mi

                      rhsVector%elementVector%vector(mhs) = rhsVector%elementVector%vector(mhs) + SUM * RWG

                    !-----------------------------------------------------------------------------------------------------------------
                    !mass-increase test function
                    ELSE IF( mh==FIELD_VARIABLE%NUMBER_OF_COMPONENTS ) THEN

                      SUM = 0.0_DP

                      PGM=QUADRATURE_SCHEME_1%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)

                      ! + possible SOURCE AND SINK TERMS
                      SOURCE = 0.0_DP

                      SUM = SUM + PGM * SOURCE

                      rhsVector%elementVector%vector(mhs) = rhsVector%elementVector%vector(mhs) + SUM * RWG

                    ELSE

                      rhsVector%elementVector%vector(mhs) = 0.0_DP

                    END IF
                    !-------------------------------------------------------------------------------------------------------------
                  CASE(EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
                    !-----------------------------------------------------------------------------------------------------------------
                    !velocity test function
                    IF( mh<FIELD_VARIABLE%NUMBER_OF_COMPONENTS ) THEN

                      SUM = 0.0_DP

                      PGM = QUADRATURE_SCHEME_1%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)

                      !Term arising from the pressure / Lagrange Multiplier of elasticity (given):
                      !TO DO- need to read different grad p depending on the compartment of interest
                      DO mi=1,DEPENDENT_BASIS_1%NUMBER_OF_XI
!                         SUM = SUM - PGM * GRAD_LM_PRESSURE(mi) * &
!                           & equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%DXI_DX(mi,mh)
                        !this is the pressure gradient for the appropriate compartment
                        SUM = SUM - PGM * GRAD_PRESSURE(mi,my_compartment) * &
                          & equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%DXI_DX(mi,mh)
                      ENDDO !mi

                      rhsVector%elementVector%vector(mhs) = rhsVector%elementVector%vector(mhs) + SUM * RWG

                    !-----------------------------------------------------------------------------------------------------------------
                    !mass-increase test function
                    ELSE IF( mh==FIELD_VARIABLE%NUMBER_OF_COMPONENTS ) THEN

                      SUM = 0.0_DP

                      PGM=QUADRATURE_SCHEME_1%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)

                      ! n o   s o u r c e
                      !source terms need to be converted to use source field & vector
                      SOURCE = 0.0_DP


                     !Add in the source/sink terms due to the pressure difference between compartments
                     DO imatrix=1,Ncompartments
                       IF(imatrix/=my_compartment) THEN
                       !Interpolate the coupling material parameter from the V variable type of the materials field
                        INTER_COMP_PERM_1=equations%interpolation%materialsInterpPoint(FIELD_V_VARIABLE_TYPE)%ptr% &
                          & VALUES(my_compartment,NO_PART_DERIV)
                        INTER_COMP_PERM_2=equations%interpolation%materialsInterpPoint(FIELD_V_VARIABLE_TYPE)%ptr% &
                          & VALUES(imatrix,NO_PART_DERIV)
                       !Source term is coefficient*(p(my_compartment) - p(imatrix))
                       INTER_COMP_SOURCE=-INTER_COMP_PERM_1*PRESSURE(my_compartment) + INTER_COMP_PERM_2*PRESSURE(imatrix)
                       ENDIF
                     ENDDO


                      SUM = SUM + PGM * SOURCE + PGM * INTER_COMP_SOURCE

                      rhsVector%elementVector%vector(mhs) = rhsVector%elementVector%vector(mhs) + SUM * RWG

                    ELSE

                      rhsVector%elementVector%vector(mhs) = 0.0_DP

                    END IF
                  !=================================================================================
                  !    d e f a u l t   :   R H S
                  CASE DEFAULT
                    !-----------------------------------------------------------------------------------------------------------------
                    !velocity test function
                    IF( mh<NUMBER_OF_VEL_PRESS_COMPONENTS ) THEN

                      SUM = 0.0_DP

                      rhsVector%elementVector%vector(mhs) = rhsVector%elementVector%vector(mhs) + SUM * RWG

                    !-----------------------------------------------------------------------------------------------------------------
                    !pressure test function
                    ELSE IF( mh==NUMBER_OF_VEL_PRESS_COMPONENTS ) THEN

                      SUM = 0.0_DP

                      PGM=QUADRATURE_SCHEME_1%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)

                      ! n o   s o u r c e
                      SOURCE = 0.0_DP

                      SUM = SUM + PGM * SOURCE

                      rhsVector%elementVector%vector(mhs) = rhsVector%elementVector%vector(mhs) + SUM * RWG

                    !-------------------------------------------------------------------------------------------------------------
                    !For the INRIA model, and: mass-increase test function
                    ELSE IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE.AND. &
                      & mh==FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN

                      SUM = 0.0_DP

                      PGM=QUADRATURE_SCHEME_1%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)

                      SUM = SUM - PGM * bfact * (1.0_DP - Jxy)

                      SUM = SUM - PGM * p0fact / (Mfact * ffact)

                      rhsVector%elementVector%vector(mhs) = rhsVector%elementVector%vector(mhs) + SUM * RWG

                    ELSE

                      rhsVector%elementVector%vector(mhs) = 0.0_DP

                    END IF
                    !-------------------------------------------------------------------------------------------------------------
                  END SELECT
                  !   e n d   s e l e c t   EQUATIONS_SET_SUBTYPE
                  !=================================================================================

                END IF

                IF(ASSOCIATED(sourceVector)) THEN
                  IF(sourceVector%updateVector) THEN
                    IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE .OR. &
                      & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE) THEN

                        C_PARAM=equations%interpolation%sourceInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(mh, NO_PART_DERIV)

                        !IF(ABS(C_PARAM)>1.0E-08) WRITE(*,*)'C_PARAM = ',C_PARAM

                        SUM = 0.0_DP
                        PGM=QUADRATURE_SCHEME_1%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                        SUM = SUM + PGM * C_PARAM
                        sourceVector%elementVector%vector(mhs) = sourceVector%elementVector%vector(mhs) + SUM * RWG
                    ENDIF
                  END IF
                END IF
              ENDDO !ms
            ENDDO !mh

            IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE) THEN
            !Calculate the momentum coupling matrices

              !Loop over element rows
              mhs=0
              DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS !field_variable is the variable associated with the equations set under consideration

                MESH_COMPONENT_1 = FIELD_VARIABLE%COMPONENTS(mh)%MESH_COMPONENT_NUMBER
                DEPENDENT_BASIS_1 => dependentField%DECOMPOSITION%DOMAIN(MESH_COMPONENT_1)%ptr% &
                  & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
                QUADRATURE_SCHEME_1 => DEPENDENT_BASIS_1%QUADRATURE% &
                  & QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
                RWG = equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%JACOBIAN * &
                  & QUADRATURE_SCHEME_1%GAUSS_WEIGHTS(ng)

                DO ms=1,DEPENDENT_BASIS_1%NUMBER_OF_ELEMENT_PARAMETERS
                  mhs=mhs+1

                  num_var_count=0
                  DO imatrix = 1,Ncompartments
                  IF(imatrix/=my_compartment)THEN
                    num_var_count=num_var_count+1

!need to test for the case where imatrix==mycompartment
!the coupling terms then needs to be added into the stiffness matrix
                    IF(COUPLING_MATRICES(num_var_count)%ptr%updateMatrix) THEN

                      !Loop over element columns
                      nhs=0
                      DO nh=1,FIELD_VARIABLES(num_var_count)%ptr%NUMBER_OF_COMPONENTS

                        MESH_COMPONENT_2 = FIELD_VARIABLE%COMPONENTS(nh)%MESH_COMPONENT_NUMBER
                        DEPENDENT_BASIS_2 => dependentField%DECOMPOSITION%DOMAIN(MESH_COMPONENT_2)%ptr% &
                          & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
                        !--- We cannot use two different quadrature schemes here !!!
                        QUADRATURE_SCHEME_2 => DEPENDENT_BASIS_2%QUADRATURE% &
                         & QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
                        !RWG = equations%interpolation%geometricInterpPointMetrics%JACOBIAN * &
                        !  & QUADRATURE_SCHEME_2%GAUSS_WEIGHTS(ng)

                        DO ns=1,DEPENDENT_BASIS_2%NUMBER_OF_ELEMENT_PARAMETERS
                          nhs=nhs+1

!                           !-------------------------------------------------------------------------------------------------------------
!                           !concentration test function, concentration trial function
!                           !For now, this is only a dummy implementation - this still has to be properly set up.
!                           IF(mh==nh.AND.nh<NUMBER_OF_VEL_PRESS_COMPONENTS) THEN ! don't need this for diffusion equation

!                             SUM = 0.0_DP

                            PGM=QUADRATURE_SCHEME_1%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                            PGN=QUADRATURE_SCHEME_2%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)

                            !Get the coupling coefficients
                              COUPLING_PARAM=equations%interpolation%materialsInterpPoint(FIELD_V_VARIABLE_TYPE)%ptr% &
                                & VALUES(imatrix,NO_PART_DERIV)

!                              SUM = SUM + COUPLING_PARAM * PGM * PGN

                             COUPLING_MATRICES(num_var_count)%ptr%elementMatrix%matrix(mhs,nhs) = &
                               & COUPLING_MATRICES(num_var_count)%ptr%elementMatrix%matrix(mhs,nhs) + &
                               & COUPLING_PARAM * PGM * PGN * RWG
!                           ENDIF

                        ENDDO !ns
                      ENDDO !nh
                    ENDIF
                   ENDIF
                  ENDDO !imatrix
                ENDDO !ms
              ENDDO !mh
            ENDIF


            !-----------------------------------------------------------------------------------------------------------------------------------
            ! RIGHT HAND SIDE FOR ANALYTIC SOLUTION
            !-----------------------------------------------------------------------------------------------------------------------------------

            IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
              IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_DARCY_EQUATION_TWO_DIM_1.OR. &
                & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_DARCY_EQUATION_TWO_DIM_2.OR. &
                & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_DARCY_EQUATION_TWO_DIM_3.OR. &
                & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_DARCY_EQUATION_THREE_DIM_1.OR. &
                & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_DARCY_EQUATION_THREE_DIM_2.OR. &
                & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_DARCY_EQUATION_THREE_DIM_3) THEN

                mhs=0
                DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                  MESH_COMPONENT_1=FIELD_VARIABLE%COMPONENTS(mh)%MESH_COMPONENT_NUMBER
                  DEPENDENT_BASIS_1=>dependentField%DECOMPOSITION%DOMAIN(MESH_COMPONENT_1)%ptr% &
                    & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
                  QUADRATURE_SCHEME_1=>DEPENDENT_BASIS_1%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
                  RWG=equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%JACOBIAN* &
                    & QUADRATURE_SCHEME_1%GAUSS_WEIGHTS(ng)
                  DO ms=1,DEPENDENT_BASIS_1%NUMBER_OF_ELEMENT_PARAMETERS
                    mhs=mhs+1
                    PGM=QUADRATURE_SCHEME_1%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                    !note mh value derivative
                    SUM=0.0_DP

                    X(1) = equations%interpolation%geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(1,1)
                    X(2) = equations%interpolation%geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(2,1)
                    IF(DEPENDENT_BASIS_1%NUMBER_OF_XI==3) THEN
                      X(3) = equations%interpolation%geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(3,1)
                    END IF
                    IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_DARCY_EQUATION_TWO_DIM_1) THEN
                      SUM=0.0_DP
                    ELSE IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_DARCY_EQUATION_TWO_DIM_2) THEN
                      IF(mh==3) THEN
                        FACT   = PERM_OVER_VIS_PARAM / L
                        ARG(1) = X(1) / L
                        ARG(2) = X(2) / L
                        SOURCE = -2.0_DP / L * FACT * EXP( ARG(1) ) * EXP( ARG(2) )
                        SUM = PGM * SOURCE
                      ELSE
                        SUM = 0.0_DP
                      ENDIF
                    ELSE IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_DARCY_EQUATION_TWO_DIM_3) THEN
                      IF(mh==3) THEN
                        FACT   = 2.0_DP * PI * PERM_OVER_VIS_PARAM / L
                        ARG(1) = 2.0_DP * PI * X(1) / L
                        ARG(2) = 2.0_DP * PI * X(2) / L
                        SOURCE = +2.0_DP * (2.0_DP * PI / L) * FACT * SIN( ARG(1) ) * SIN( ARG(2) )
                        SUM = PGM * SOURCE
                      ELSE
                        SUM = 0.0_DP
                      ENDIF
                    ELSE IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_DARCY_EQUATION_THREE_DIM_1) THEN
                      SUM=0.0_DP
                    ELSE IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_DARCY_EQUATION_THREE_DIM_2) THEN
                      IF(mh==4) THEN
                        FACT   = PERM_OVER_VIS_PARAM / L
                        ARG(1) = X(1) / L
                        ARG(2) = X(2) / L
                        ARG(3) = X(3) / L
                        SOURCE = -3.0_DP / L * FACT * EXP( ARG(1) ) * EXP( ARG(2) ) * EXP( ARG(3) )
                        SUM = PGM * SOURCE
                      ELSE
                        SUM = 0.0_DP
                      ENDIF
                    ELSE IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_DARCY_EQUATION_THREE_DIM_3) THEN
                      IF(mh==4) THEN
                        FACT   = 2.0_DP * PI * PERM_OVER_VIS_PARAM / L
                        ARG(1) = 2.0_DP * PI * X(1) / L
                        ARG(2) = 2.0_DP * PI * X(2) / L
                        ARG(3) = 2.0_DP * PI * X(3) / L
                        SOURCE = +3.0_DP * ( 2.0_DP * PI / L ) * FACT * SIN( ARG(1) ) * SIN( ARG(2) ) * SIN( ARG(3) )
                        SUM = PGM * SOURCE
                      ELSE
                        SUM = 0.0_DP
                      END IF
                    ENDIF

                    !Calculate RHS VECTOR
                    rhsVector%elementVector%vector(mhs)=rhsVector%elementVector%vector(mhs)+SUM*RWG
                  ENDDO !ms
                ENDDO !mh
              ELSE
                rhsVector%elementVector%vector(mhs)=0.0_DP
              ENDIF
            ENDIF

            ! end: RIGHT HAND SIDE FOR ANALYTIC SOLUTION
            !-----------------------------------------------------------------------------------------------------------------------------------

!             !===================================================================================================================
!             !COUPLING_MATRICES
!             SELECT CASE(EQUATIONS_SET_SUBTYPE)
!             CASE(EQUATIONS_SET_MULTI_COMPARTMENT_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_MULTI_COMPARTMENT_DARCY_INRIA_SUBTYPE)
!
!               !Create FIELD_VARIABLES type, COUPLING_MATRICES type
!
!               !Loop over element rows
!               mhs=0
!               DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
!
!                 MESH_COMPONENT_1 = FIELD_VARIABLE%COMPONENTS(mh)%MESH_COMPONENT_NUMBER
!                 DEPENDENT_BASIS_1 => dependentField%DECOMPOSITION%DOMAIN(MESH_COMPONENT_1)%ptr% &
!                   & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
!                 QUADRATURE_SCHEME_1 => DEPENDENT_BASIS_1%QUADRATURE% &
!                   & QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
!                 RWG = equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%JACOBIAN * &
!                   & QUADRATURE_SCHEME_1%GAUSS_WEIGHTS(ng)
!
!                 DO ms=1,DEPENDENT_BASIS_1%NUMBER_OF_ELEMENT_PARAMETERS
!                   mhs=mhs+1
!
!                   DO imatrix=1,Ncompartments
!
!                     IF(COUPLING_MATRICES(imatrix)%ptr%updateMatrix) THEN
!
!                       !Loop over element columns
!                       nhs=0
! !                       DO nh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
!                       DO nh=1,FIELD_VARIABLES(imatrix)%ptr%NUMBER_OF_COMPONENTS
!
!                         MESH_COMPONENT_2 = FIELD_VARIABLE%COMPONENTS(nh)%MESH_COMPONENT_NUMBER
!                         DEPENDENT_BASIS_2 => dependentField%DECOMPOSITION%DOMAIN(MESH_COMPONENT_2)%ptr% &
!                           & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
!                         !--- We cannot use two different quadrature schemes here !!!
!                         QUADRATURE_SCHEME_2 => DEPENDENT_BASIS_2%QUADRATURE% &
!                          & QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
!                         !RWG = equations%interpolation%geometricInterpPointMetrics%JACOBIAN * &
!                         !  & QUADRATURE_SCHEME_2%GAUSS_WEIGHTS(ng)
!
!                         DO ns=1,DEPENDENT_BASIS_2%NUMBER_OF_ELEMENT_PARAMETERS
!                           nhs=nhs+1
!
!                           !-------------------------------------------------------------------------------------------------------------
!                           !velocity test function, velocity trial function
!                           !For now, this is only a dummy implementation - this still has to be properly set up.
!                           IF(mh==nh.AND.nh<NUMBER_OF_VEL_PRESS_COMPONENTS) THEN
!
!                             SUM = 0.0_DP
!
!                             PGM=QUADRATURE_SCHEME_1%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
!                             PGN=QUADRATURE_SCHEME_2%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)
!
!                             SUM = SUM + VIS_OVER_PERM_TENSOR( mh, nh ) * PGM * PGN
!
!                             COUPLING_MATRICES(imatrix)%ptr%elementMatrix%matrix(mhs,nhs) = &
!                               & COUPLING_MATRICES(imatrix)%ptr%elementMatrix%matrix(mhs,nhs) + SUM * RWG
!                           ENDIF
!
!                         ENDDO !ns
!                       ENDDO !nh
!                     ENDIF
!                   ENDDO !imatrix
!                 ENDDO !ms
!               ENDDO !mh
!             CASE DEFAULT
!               !Do nothing
!             END SELECT
          ENDDO !ng

          IF(rhsVector%updateVector) THEN
            ! Integrate pressure over faces, and add to RHS vector
            CALL Darcy_FiniteElementFaceIntegrate(EQUATIONS_SET,ELEMENT_NUMBER,FIELD_VARIABLE,err,error,*999)
          ENDIF

          ! CHECK STIFFNESS MATRIX WITH CMHEART
          IF(DIAGNOSTICS5) THEN
            IF( ELEMENT_NUMBER == 1 ) THEN
              NDOFS = 0
              DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                MESH_COMPONENT_1 = FIELD_VARIABLE%COMPONENTS(mh)%MESH_COMPONENT_NUMBER
                DEPENDENT_BASIS_1 => dependentField%DECOMPOSITION%DOMAIN(MESH_COMPONENT_1)%ptr% &
                  & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
                NDOFS = NDOFS + DEPENDENT_BASIS_1%NUMBER_OF_ELEMENT_PARAMETERS
              END DO

              CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Element Matrix for element number 1 (Darcy):",err,error,*999)
              DO mhs=1,NDOFS
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"row number = ",mhs,err,error,*999)
                CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NDOFS,NDOFS,NDOFS,&
                  & stiffnessMatrix%elementMatrix%matrix(mhs,:), &
                  & '("",4(X,E13.6))','4(4(X,E13.6))',err,error,*999)
                CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE," ",err,error,*999)
              END DO
            END IF
          END IF

          !Scale factor adjustment
          IF(dependentField%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
            CALL Field_InterpolationParametersScaleFactorsElementGet(ELEMENT_NUMBER,equations%interpolation% &
              & dependentInterpParameters(FIELD_VAR_TYPE)%ptr,err,error,*999)
            mhs=0
            DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
              !Loop over element rows
              MESH_COMPONENT_1=FIELD_VARIABLE%COMPONENTS(mh)%MESH_COMPONENT_NUMBER
              DEPENDENT_BASIS_1=>dependentField%DECOMPOSITION%DOMAIN(MESH_COMPONENT_1)%ptr% &
                & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
              DO ms=1,DEPENDENT_BASIS_1%NUMBER_OF_ELEMENT_PARAMETERS
                mhs=mhs+1
                nhs=0
                IF(ASSOCIATED(stiffnessMatrix).AND.ASSOCIATED(dampingMatrix)) THEN
                  IF(stiffnessMatrix%updateMatrix.OR.dampingMatrix%updateMatrix) THEN
                    !Loop over element columns
                    DO nh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                      MESH_COMPONENT_2=FIELD_VARIABLE%COMPONENTS(nh)%MESH_COMPONENT_NUMBER
                      DEPENDENT_BASIS_2=>dependentField%DECOMPOSITION%DOMAIN(MESH_COMPONENT_2)%ptr% &
                        & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
                      DO ns=1,DEPENDENT_BASIS_2%NUMBER_OF_ELEMENT_PARAMETERS
                        nhs=nhs+1
                        IF(stiffnessMatrix%updateMatrix)THEN
                          stiffnessMatrix%elementMatrix%matrix(mhs,nhs)=stiffnessMatrix%elementMatrix%matrix(mhs,nhs)* &
                            & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%SCALE_FACTORS(ms,mh)* &
                            & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%SCALE_FACTORS(ns,nh)
                        END IF
                        IF(dampingMatrix%updateMatrix)THEN
                          dampingMatrix%elementMatrix%matrix(mhs,nhs)=dampingMatrix%elementMatrix%matrix(mhs,nhs)* &
                            & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%SCALE_FACTORS(ms,mh)* &
                            & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%SCALE_FACTORS(ns,nh)
                        END IF
                      ENDDO !ns
                    ENDDO !nh
                  ENDIF
                ENDIF
                IF(ASSOCIATED(rhsVector)) THEN
                  IF(rhsVector%updateVector) rhsVector%elementVector%vector(mhs)=rhsVector%elementVector%vector(mhs)* &
                    & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%SCALE_FACTORS(ms,mh)
                ENDIF
                IF(ASSOCIATED(sourceVector)) THEN
                  IF(sourceVector%updateVector) sourceVector%elementVector%vector(mhs)= &
                     & sourceVector%elementVector%vector(mhs)* &
                     & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%SCALE_FACTORS(ms,mh)
                ENDIF
              ENDDO !ms
            ENDDO !mh
          ENDIF

         ! RESTORE ALL POINTERS CALL PARAMATER_SET_FIELD_DATA_RESTORE

        CASE DEFAULT
          localError="Equations set subtype "//TRIM(NumberToVString(EQUATIONS_SET_SUBTYPE,"*",err,error))// &
            & " is not valid for a Darcy equation type of a fluid mechanics equations set class."
          CALL FlagError(localError,err,error,*999)
        END SELECT

      ELSE
        CALL FlagError("Equations set equations is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF

  EXITS("DARCY_EQUATION_FINITE_ELEMENT_CALCULATE")
    RETURN
999 ERRORSEXITS("DARCY_EQUATION_FINITE_ELEMENT_CALCULATE",err,error)
    RETURN 1
  END SUBROUTINE DARCY_EQUATION_FINITE_ELEMENT_CALCULATE

  !
  !================================================================================================================================
  !

  !>Calculates the face integration term of the finite element formulation for Darcy's equation,
  !>required for pressure boundary conditions.
  SUBROUTINE Darcy_FiniteElementFaceIntegrate(equationsSet,elementNumber,dependentVariable,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<The equations set to calculate the RHS term for
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to calculat the RHS term for
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: dependentVariable
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: decomposition
    TYPE(DECOMPOSITION_ELEMENT_TYPE), POINTER :: decompElement
    TYPE(BASIS_TYPE), POINTER :: dependentBasis
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsMatricesVectorType), POINTER :: equationsMatrices
    TYPE(DECOMPOSITION_FACE_TYPE), POINTER :: face
    TYPE(BASIS_TYPE), POINTER :: faceBasis
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: dependentInterpolatedPoint
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: dependentInterpolationParameters
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: faceQuadratureScheme
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: geometricInterpolatedPoint
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: geometricInterpolationParameters
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_TYPE), POINTER :: pointMetrics
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    INTEGER(INTG) :: faceIdx, faceNumber
    INTEGER(INTG) :: componentIdx, gaussIdx
    INTEGER(INTG) :: elementBaseDofIdx, faceNodeIdx, elementNodeIdx
    INTEGER(INTG) :: faceNodeDerivativeIdx, meshComponentNumber, nodeDerivativeIdx, parameterIdx
    INTEGER(INTG) :: faceParameterIdx, elementDofIdx, normalComponentIdx
    REAL(DP) :: gaussWeight, normalProjection, pressureGauss

    ENTERS("Darcy_FiniteElementFaceIntegrate",err,error,*999)

    NULLIFY(decomposition)
    NULLIFY(decompElement)
    NULLIFY(dependentBasis)
    NULLIFY(equations)
    NULLIFY(equationsMatrices)
    NULLIFY(face)
    NULLIFY(faceBasis)
    NULLIFY(faceQuadratureScheme)
    NULLIFY(dependentInterpolatedPoint)
    NULLIFY(dependentInterpolationParameters)
    NULLIFY(geometricInterpolatedPoint)
    NULLIFY(geometricInterpolationParameters)
    NULLIFY(rhsVector)

    IF(ASSOCIATED(equationsSet)) THEN
      equations=>equationsSet%EQUATIONS
      IF(ASSOCIATED(equations)) THEN
        NULLIFY(vectorEquations)
        CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
        equationsMatrices=>vectorEquations%vectorMatrices
        IF(ASSOCIATED(equationsMatrices)) THEN
          rhsVector=>equationsMatrices%rhsVector
        END IF
      ELSE
        CALL FlagError("Equations set equations is not associated.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    END IF

    IF(.NOT.ALLOCATED(equationsSet%specification)) THEN
      CALL FlagError("Equations set specification is not allocated.",err,error,*999)
    ELSE IF(SIZE(equationsSet%specification,1)/=3) THEN
      CALL FlagError("Equations set specification must have three entries for a Darcy type equations set.", &
        & err,error,*999)
    END IF
    SELECT CASE(equationsSet%specification(3))
    CASE(EQUATIONS_SET_STANDARD_DARCY_SUBTYPE, &
        & EQUATIONS_SET_QUASISTATIC_DARCY_SUBTYPE, &
        & EQUATIONS_SET_ALE_DARCY_SUBTYPE, &
        & EQUATIONS_SET_TRANSIENT_DARCY_SUBTYPE)

      !Get the mesh decomposition and basis for this element
      decomposition=>dependentVariable%FIELD%DECOMPOSITION
      !These RHS terms are associated with the equations for the three velocity components,
      !rather than the pressure term
      meshComponentNumber=dependentVariable%COMPONENTS(1)%MESH_COMPONENT_NUMBER
      dependentBasis=>decomposition%DOMAIN(meshComponentNumber)%ptr%TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%BASIS
      decompElement=>DECOMPOSITION%TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)

      !Only add RHS terms if the face geometric parameters are calculated
      IF(decomposition%CALCULATE_FACES) THEN
        !Get interpolation parameters and point for Darcy pressure
        dependentInterpolationParameters=>equations%interpolation%dependentInterpParameters(dependentVariable%VARIABLE_TYPE)%ptr
        dependentInterpolatedPoint=>equations%interpolation%dependentInterpPoint(dependentVariable%VARIABLE_TYPE)%ptr

        DO faceIdx=1,dependentBasis%NUMBER_OF_LOCAL_FACES
          !Get the face normal and quadrature information
          IF(ALLOCATED(decompElement%ELEMENT_FACES)) THEN
            faceNumber=decompElement%ELEMENT_FACES(faceIdx)
          ELSE
            CALL FlagError("Decomposition element faces is not allocated.",err,error,*999)
          END IF
          face=>decomposition%TOPOLOGY%FACES%FACES(faceNumber)
          !This speeds things up but is also important, as non-boundary faces have an XI_DIRECTION that might
          !correspond to the other element.
          IF(.NOT.(face%BOUNDARY_FACE)) CYCLE
          CALL FIELD_INTERPOLATION_PARAMETERS_FACE_GET(FIELD_VALUES_SET_TYPE,faceNumber,dependentInterpolationParameters, &
            & err,error,*999)
          normalComponentIdx=ABS(face%XI_DIRECTION)
          faceBasis=>decomposition%DOMAIN(meshComponentNumber)%ptr%TOPOLOGY%FACES%FACES(faceNumber)%BASIS
          faceQuadratureScheme=>faceBasis%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr

          DO gaussIdx=1,faceQuadratureScheme%NUMBER_OF_GAUSS
            gaussWeight=faceQuadratureScheme%GAUSS_WEIGHTS(gaussIdx)
            !Get interpolated Darcy pressure
            CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussIdx, &
              & dependentInterpolatedPoint,err,error,*999)
            pressureGauss=dependentInterpolatedPoint%values(4,1) !(component,derivative)

            !Use the geometric field to find the face normal and the Jacobian for the face integral
            geometricInterpolationParameters=>equations%interpolation%geometricInterpParameters( &
              & FIELD_U_VARIABLE_TYPE)%ptr
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,elementNumber, &
              & geometricInterpolationParameters,err,error,*999)
            geometricInterpolatedPoint=>equations%interpolation%geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr
            CALL FIELD_INTERPOLATE_LOCAL_FACE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,faceIdx,gaussIdx, &
              & geometricInterpolatedPoint,err,error,*999)
            !Calculate the metric tensors and Jacobian
            pointMetrics=>equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(COORDINATE_JACOBIAN_VOLUME_TYPE,pointMetrics,err,error,*999)

            DO componentIdx=1,dependentVariable%NUMBER_OF_COMPONENTS-1
              normalProjection=DOT_PRODUCT(pointMetrics%GU(normalComponentIdx,:),pointMetrics%DX_DXI(componentIdx,:))
              IF(face%XI_DIRECTION<0) THEN
                normalProjection=-normalProjection
              END IF
              IF(ABS(normalProjection)<ZERO_TOLERANCE) CYCLE
              !Work out the first index of the rhs vector for this element - 1
              elementBaseDofIdx=dependentBasis%NUMBER_OF_ELEMENT_PARAMETERS*(componentIdx-1)
              DO faceNodeIdx=1,faceBasis%NUMBER_OF_NODES
                elementNodeIdx=dependentBasis%NODE_NUMBERS_IN_LOCAL_FACE(faceNodeIdx,faceIdx)
                DO faceNodeDerivativeIdx=1,faceBasis%NUMBER_OF_DERIVATIVES(faceNodeIdx)
                  nodeDerivativeIdx=dependentBasis%DERIVATIVE_NUMBERS_IN_LOCAL_FACE(faceNodeDerivativeIdx,faceNodeIdx,faceIdx)
                  parameterIdx=dependentBasis%ELEMENT_PARAMETER_INDEX(nodeDerivativeIdx,elementNodeIdx)
                  faceParameterIdx=faceBasis%ELEMENT_PARAMETER_INDEX(faceNodeDerivativeIdx,faceNodeIdx)
                  elementDofIdx=elementBaseDofIdx+parameterIdx
                  rhsVector%elementVector%vector(elementDofIdx) = rhsVector%elementVector%vector(elementDofIdx) - &
                    & gaussWeight*pressureGauss*normalProjection* &
                    & faceQuadratureScheme%GAUSS_BASIS_FNS(faceParameterIdx,NO_PART_DERIV,gaussIdx)* &
                    & pointMetrics%JACOBIAN
                END DO !nodeDerivativeIdx
              END DO !faceNodeIdx
            END DO !componentIdx
          END DO !gaussIdx
        END DO !faceIdx
      END IF !decomposition%calculate_faces

    CASE DEFAULT
      ! Do nothing for other equation set subtypes
    END SELECT

    EXITS("Darcy_FiniteElementFaceIntegrate")
    RETURN
999 ERRORSEXITS("Darcy_FiniteElementFaceIntegrate",err,error)
    RETURN 1
  END SUBROUTINE

  !
  !================================================================================================================================
  !

  !>Sets the equation specification for a Darcy equation type of a fluid mechanics equations set class.
  SUBROUTINE Darcy_EquationsSetSpecificationSet(equationsSet,specification,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to set the specification for
    INTEGER(INTG), INTENT(IN) :: specification(:) !<The equations set specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: subtype

    ENTERS("Darcy_EquationsSetSpecificationSet",err,error,*999)

    IF(ASSOCIATED(equationsSet)) THEN
      IF(SIZE(specification,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a Darcy type equations set.", &
          & err,error,*999)
      END IF
      subtype=specification(3)
      SELECT CASE(subtype)
      CASE(EQUATIONS_SET_STANDARD_DARCY_SUBTYPE, &
          & EQUATIONS_SET_QUASISTATIC_DARCY_SUBTYPE, &
          & EQUATIONS_SET_ALE_DARCY_SUBTYPE, &
          & EQUATIONS_SET_TRANSIENT_DARCY_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE, &
          & EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE, &
          & EQUATIONS_SET_TRANSIENT_ALE_DARCY_SUBTYPE, &
          & EQUATIONS_SET_MULTI_COMPARTMENT_DARCY_SUBTYPE, &
          & EQUATIONS_SET_ELASTICITY_MULTI_COMPARTMENT_DARCY_INRIA_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
        !ok
      CASE DEFAULT
        localError="The third equations set specification of "//TRIM(NumberToVstring(specification(3),"*",err,error))// &
          & " is not valid for a Darcy type of a fluid mechanics equations set."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      !Set full specification
      IF(ALLOCATED(equationsSet%specification)) THEN
        CALL FlagError("Equations set specification is already allocated.",err,error,*999)
      ELSE
        ALLOCATE(equationsSet%specification(3),stat=err)
        IF(err/=0) CALL FlagError("Could not allocate equations set specification.",err,error,*999)
      END IF
      equationsSet%specification(1:3)=[EQUATIONS_SET_FLUID_MECHANICS_CLASS,EQUATIONS_SET_DARCY_EQUATION_TYPE,subtype]
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    END IF

    EXITS("Darcy_EquationsSetSpecificationSet")
    RETURN
999 ERRORS("Darcy_EquationsSetSpecificationSet",err,error)
    EXITS("Darcy_EquationsSetSpecificationSet")
    RETURN 1
    
  END SUBROUTINE Darcy_EquationsSetSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets the problem specification for a Darcy problem.
  SUBROUTINE Darcy_ProblemSpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: problem !<A pointer to the problem to set the problem specification for
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The problem specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: problemSubtype

    ENTERS("Darcy_ProblemSpecificationSet",err,error,*998)

    IF(ASSOCIATED(problem)) THEN
      IF(SIZE(problemSpecification,1)==3) THEN
        problemSubtype=problemSpecification(3)
        SELECT CASE(problemSubtype)
        CASE(PROBLEM_STANDARD_DARCY_SUBTYPE, &
            & PROBLEM_QUASISTATIC_DARCY_SUBTYPE, &
            & PROBLEM_ALE_DARCY_SUBTYPE, &
            & PROBLEM_TRANSIENT_DARCY_SUBTYPE, &
            & PROBLEM_PGM_DARCY_SUBTYPE, &
            & PROBLEM_PGM_TRANSIENT_DARCY_SUBTYPE)
          !All ok
        CASE DEFAULT
          localError="The third problem subtype of "//TRIM(NumberToVstring(problemSubtype,"*",err,error))// &
            & " is not valid for a Darcy type of a fluid mechanics problem."
          CALL FlagError(localError,err,error,*998)
        END SELECT
        IF(ALLOCATED(problem%specification)) THEN
          CALL FlagError("Problem specification is already allocated.",err,error,*998)
        ELSE
          ALLOCATE(problem%specification(3),stat=err)
          IF(err/=0) CALL FlagError("Could not allocate problem specification.",err,error,*999)
        END IF
        problem%specification(1:3)=[PROBLEM_FLUID_MECHANICS_CLASS,PROBLEM_DARCY_EQUATION_TYPE,problemSubtype]
      ELSE
        CALL FlagError("Darcy problem specification must have three entries.",err,error,*998)
      END IF
    ELSE
      CALL FlagError("Problem is not associated.",err,error,*998)
    END IF

    EXITS("Darcy_ProblemSpecificationSet")
    RETURN
999 IF(ALLOCATED(problem%specification)) DEALLOCATE(problem%specification)
998 ERRORS("Darcy_ProblemSpecificationSet",err,error)
    EXITS("Darcy_ProblemSpecificationSet")
    RETURN 1
    
  END SUBROUTINE Darcy_ProblemSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets up the Darcy equations problem.
!   SUBROUTINE DARCY_EQUATION_PROBLEM_STANDARD_SETUP(PROBLEM,PROBLEM_SETUP,err,error,*)
  SUBROUTINE DARCY_EQUATION_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to setup
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP,CONTROL_LOOP_ROOT
    TYPE(SOLVER_TYPE), POINTER :: SOLVER, SOLVER_MAT_PROPERTIES
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS, SOLVER_EQUATIONS_MAT_PROPERTIES
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: localError

!     ENTERS("DARCY_EQUATION_PROBLEM_STANDARD_SETUP",err,error,*999)
    ENTERS("DARCY_EQUATION_PROBLEM_SETUP",err,error,*999)

    NULLIFY(CONTROL_LOOP)
    NULLIFY(SOLVERS)
    NULLIFY(SOLVER)
    NULLIFY(SOLVER_MAT_PROPERTIES)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(SOLVER_EQUATIONS_MAT_PROPERTIES)
    IF(ASSOCIATED(PROBLEM)) THEN
      IF(.NOT.ALLOCATED(PROBLEM%SPECIFICATION)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(PROBLEM%SPECIFICATION,1)<3) THEN
        CALL FlagError("Problem specification must have three entries for a Darcy equation problem.",err,error,*999)
      END IF
      SELECT CASE(PROBLEM%SPECIFICATION(3))

      !-----------------------------------------------------------------
      !   s t a n d a r d   D a r c y
      !-----------------------------------------------------------------
      CASE(PROBLEM_STANDARD_DARCY_SUBTYPE)
        SELECT CASE(PROBLEM_SETUP%SETUP_TYPE)
        CASE(PROBLEM_SETUP_INITIAL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Do nothing????
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Do nothing???
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a standard Darcy equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_CONTROL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Set up a simple control loop
            CALL CONTROL_LOOP_CREATE_START(PROBLEM,CONTROL_LOOP,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Finish the control loops
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            CALL CONTROL_LOOP_CREATE_FINISH(CONTROL_LOOP,err,error,*999)
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a standard Darcy equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVERS_TYPE)
          !Get the control loop
          CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
          CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Start the solvers creation
            CALL SOLVERS_CREATE_START(CONTROL_LOOP,SOLVERS,err,error,*999)
            CALL SOLVERS_NUMBER_SET(SOLVERS,1,err,error,*999)
            !Set the solver to be a linear solver
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
            CALL SOLVER_TYPE_SET(SOLVER,SOLVER_LINEAR_TYPE,err,error,*999)
            !Set solver defaults
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_PETSC_LIBRARY,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the solvers
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            !Finish the solvers creation
            CALL SOLVERS_CREATE_FINISH(SOLVERS,err,error,*999)
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
                & " is invalid for a standard Darcy equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            !Get the solver
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
            !Create the solver equations
            CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,err,error,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_LINEAR,err,error,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_STATIC,err,error,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            !Get the solver equations
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,err,error,*999)
            !Finish the solver equations creation
            CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,err,error,*999)
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a standard Darcy equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
            & " is invalid for a standard Darcy equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT

      !-----------------------------------------------------------------
      !   q u a s i s t a t i c   D a r c y
      !-----------------------------------------------------------------
      CASE(PROBLEM_QUASISTATIC_DARCY_SUBTYPE)
        SELECT CASE(PROBLEM_SETUP%SETUP_TYPE)
        CASE(PROBLEM_SETUP_INITIAL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Do nothing????
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Do nothing???
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a quasistatic Darcy equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_CONTROL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Set up a time control loop
            CALL CONTROL_LOOP_CREATE_START(PROBLEM,CONTROL_LOOP,err,error,*999)
            CALL CONTROL_LOOP_TYPE_SET(CONTROL_LOOP,PROBLEM_CONTROL_TIME_LOOP_TYPE,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Finish the control loops
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            CALL CONTROL_LOOP_CREATE_FINISH(CONTROL_LOOP,err,error,*999)
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a quasistatic Darcy equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVERS_TYPE)
          !Get the control loop
          CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
          CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Start the solvers creation
            CALL SOLVERS_CREATE_START(CONTROL_LOOP,SOLVERS,err,error,*999)
            CALL SOLVERS_NUMBER_SET(SOLVERS,1,err,error,*999)
            !Set the solver to be a linear solver
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
            CALL SOLVER_TYPE_SET(SOLVER,SOLVER_LINEAR_TYPE,err,error,*999)
            !Set solver defaults
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_PETSC_LIBRARY,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the solvers
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            !Finish the solvers creation
            CALL SOLVERS_CREATE_FINISH(SOLVERS,err,error,*999)
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
                & " is invalid for a quasistatic Darcy equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            !Get the solver
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
            !Create the solver equations
            CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,err,error,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_LINEAR,err,error,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_QUASISTATIC,err,error,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            !Get the solver equations
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,err,error,*999)
            !Finish the solver equations creation
            CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,err,error,*999)
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a quasistatic Darcy equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
            & " is invalid for a quasistatic Darcy equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT

      !-----------------------------------------------------------------
      !   A L E / P G M   D a r c y
      !-----------------------------------------------------------------
      CASE(PROBLEM_ALE_DARCY_SUBTYPE,PROBLEM_PGM_DARCY_SUBTYPE)
        SELECT CASE(PROBLEM_SETUP%SETUP_TYPE)
        CASE(PROBLEM_SETUP_INITIAL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Do nothing????
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Do nothing???
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for an ALE Darcy equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_CONTROL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Set up a time control loop
            CALL CONTROL_LOOP_CREATE_START(PROBLEM,CONTROL_LOOP,err,error,*999)
            CALL CONTROL_LOOP_TYPE_SET(CONTROL_LOOP,PROBLEM_CONTROL_TIME_LOOP_TYPE,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Finish the control loops
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            CALL CONTROL_LOOP_CREATE_FINISH(CONTROL_LOOP,err,error,*999)
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for an ALE Darcy equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVERS_TYPE)
          !Get the control loop
          CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
          CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Start the solvers creation
            CALL SOLVERS_CREATE_START(CONTROL_LOOP,SOLVERS,err,error,*999)
            CALL SOLVERS_NUMBER_SET(SOLVERS,2,err,error,*999)
            !
            !Set the first solver to be a linear solver for the material update
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER_MAT_PROPERTIES,err,error,*999)
            CALL SOLVER_TYPE_SET(SOLVER_MAT_PROPERTIES,SOLVER_LINEAR_TYPE,err,error,*999)
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER_MAT_PROPERTIES,SOLVER_PETSC_LIBRARY,err,error,*999)
            !
            !Set the second solver to be a linear solver for the ALE Darcy
            CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
            CALL SOLVER_TYPE_SET(SOLVER,SOLVER_LINEAR_TYPE,err,error,*999)
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_PETSC_LIBRARY,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the solvers
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            !Finish the solvers creation
            CALL SOLVERS_CREATE_FINISH(SOLVERS,err,error,*999)
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
                & " is invalid for an ALE Darcy equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Get the control loop and solvers
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            !Get the material-properties solver and create the material-properties solver equations
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER_MAT_PROPERTIES,err,error,*999)
            CALL SOLVER_EQUATIONS_CREATE_START(SOLVER_MAT_PROPERTIES,SOLVER_EQUATIONS_MAT_PROPERTIES,err,error,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS_MAT_PROPERTIES,SOLVER_EQUATIONS_LINEAR,err,error,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS_MAT_PROPERTIES,SOLVER_EQUATIONS_QUASISTATIC, &
              & err,error,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS_MAT_PROPERTIES,SOLVER_SPARSE_MATRICES,err,error,*999)
            !Get the Darcy-ALE solver and create the Darcy-ALE solver equations
            CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
            CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,err,error,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_LINEAR,err,error,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_QUASISTATIC,err,error,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            !Finish the creation of the material-properties solver equations
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER_MAT_PROPERTIES,err,error,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER_MAT_PROPERTIES,SOLVER_EQUATIONS_MAT_PROPERTIES,err,error,*999)
            CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS_MAT_PROPERTIES,err,error,*999)
            !Finish the creation of the Darcy-ALE solver equations
            CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,err,error,*999)
            CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,err,error,*999)
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for an ALE Darcy equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
            & " is invalid for an ALE Darcy equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT

      !-----------------------------------------------------------------
      !   D Y N A M I C   A L E / P G M   D a r c y
      !-----------------------------------------------------------------
      CASE(PROBLEM_PGM_TRANSIENT_DARCY_SUBTYPE)
        SELECT CASE(PROBLEM_SETUP%SETUP_TYPE)
        CASE(PROBLEM_SETUP_INITIAL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Do nothing????
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Do nothing???
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for an ALE Darcy equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_CONTROL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Set up a time control loop
            CALL CONTROL_LOOP_CREATE_START(PROBLEM,CONTROL_LOOP,err,error,*999)
            CALL CONTROL_LOOP_TYPE_SET(CONTROL_LOOP,PROBLEM_CONTROL_TIME_LOOP_TYPE,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Finish the control loops
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            CALL CONTROL_LOOP_CREATE_FINISH(CONTROL_LOOP,err,error,*999)
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for an ALE Darcy equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVERS_TYPE)
          !Get the control loop
          CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
          CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Start the solvers creation
            CALL SOLVERS_CREATE_START(CONTROL_LOOP,SOLVERS,err,error,*999)
            CALL SOLVERS_NUMBER_SET(SOLVERS,2,err,error,*999)
            !
            !Set the first solver to be a linear solver for the material update
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER_MAT_PROPERTIES,err,error,*999)
            CALL SOLVER_TYPE_SET(SOLVER_MAT_PROPERTIES,SOLVER_LINEAR_TYPE,err,error,*999)
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER_MAT_PROPERTIES,SOLVER_PETSC_LIBRARY,err,error,*999)
            !
            !Set the second solver to be a first order dynamic solver for the ALE Darcy
            CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
            CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DYNAMIC_TYPE,err,error,*999)
            CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
            !Set solver defaults
            CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
            CALL SOLVER_DYNAMIC_SCHEME_SET(SOLVER,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,err,error,*999)
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,err,error,*999)
!             CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_PETSC_LIBRARY,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the solvers
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            !Finish the solvers creation
            CALL SOLVERS_CREATE_FINISH(SOLVERS,err,error,*999)
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
                & " is invalid for an ALE Darcy equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Get the control loop and solvers
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            !Get the material-properties solver and create the material-properties solver equations
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER_MAT_PROPERTIES,err,error,*999)
            CALL SOLVER_EQUATIONS_CREATE_START(SOLVER_MAT_PROPERTIES,SOLVER_EQUATIONS_MAT_PROPERTIES,err,error,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS_MAT_PROPERTIES,SOLVER_EQUATIONS_LINEAR,err,error,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS_MAT_PROPERTIES,SOLVER_EQUATIONS_QUASISTATIC, &
              & err,error,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS_MAT_PROPERTIES,SOLVER_SPARSE_MATRICES,err,error,*999)
            !Get the Darcy-ALE solver and create the Darcy-ALE solver equations
            CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
            CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,err,error,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_LINEAR,err,error,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            !Finish the creation of the material-properties solver equations
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER_MAT_PROPERTIES,err,error,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER_MAT_PROPERTIES,SOLVER_EQUATIONS_MAT_PROPERTIES,err,error,*999)
            CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS_MAT_PROPERTIES,err,error,*999)
            !Finish the creation of the Darcy-ALE solver equations
            CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,err,error,*999)
            CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,err,error,*999)
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for an ALE Darcy equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
            & " is invalid for an ALE Darcy equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT

      !-----------------------------------------------------------------
      !   d y n a m i c   D a r c y
      !-----------------------------------------------------------------
      CASE(PROBLEM_TRANSIENT_DARCY_SUBTYPE)
        SELECT CASE(PROBLEM_SETUP%SETUP_TYPE)
          CASE(PROBLEM_SETUP_INITIAL_TYPE)
            SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
              CASE(PROBLEM_SETUP_START_ACTION)
                !Do nothing????
              CASE(PROBLEM_SETUP_FINISH_ACTION)
                !Do nothing????
              CASE DEFAULT
                localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
                  & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
                  & " is invalid for a transient Darcy fluid."
                CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE(PROBLEM_SETUP_CONTROL_TYPE)
            SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
              CASE(PROBLEM_SETUP_START_ACTION)
                !Set up a time control loop
                CALL CONTROL_LOOP_CREATE_START(PROBLEM,CONTROL_LOOP,err,error,*999)
                CALL CONTROL_LOOP_TYPE_SET(CONTROL_LOOP,PROBLEM_CONTROL_TIME_LOOP_TYPE,err,error,*999)
              CASE(PROBLEM_SETUP_FINISH_ACTION)
                !Finish the control loops
                CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
                CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
                CALL CONTROL_LOOP_CREATE_FINISH(CONTROL_LOOP,err,error,*999)
              CASE DEFAULT
                localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
                  & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
                  & " is invalid for a transient Darcy fluid."
                CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE(PROBLEM_SETUP_SOLVERS_TYPE)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
              CASE(PROBLEM_SETUP_START_ACTION)
                !Start the solvers creation
                CALL SOLVERS_CREATE_START(CONTROL_LOOP,SOLVERS,err,error,*999)
                CALL SOLVERS_NUMBER_SET(SOLVERS,1,err,error,*999)
                !Set the solver to be a first order dynamic solver
                CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
                CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DYNAMIC_TYPE,err,error,*999)
                CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
                !Set solver defaults
                CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
                CALL SOLVER_DYNAMIC_SCHEME_SET(SOLVER,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,err,error,*999)
                CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,err,error,*999)
              CASE(PROBLEM_SETUP_FINISH_ACTION)
                !Get the solvers
                CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
                !Finish the solvers creation
                CALL SOLVERS_CREATE_FINISH(SOLVERS,err,error,*999)
              CASE DEFAULT
                localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
                 & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
                 & " is invalid for a transient Darcy fluid."
                CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
            SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
              CASE(PROBLEM_SETUP_START_ACTION)
                !Get the control loop
                CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
                CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
                !Get the solver
                CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
                CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
                !Create the solver equations
                CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,err,error,*999)
                CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_LINEAR,err,error,*999)
                CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,&
                & err,error,*999)
                CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,err,error,*999)
              CASE(PROBLEM_SETUP_FINISH_ACTION)
                !Get the control loop
                CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
                CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
                !Get the solver equations
                CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
                CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
                CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,err,error,*999)
                !Finish the solver equations creation
                CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,err,error,*999)
              CASE DEFAULT
                localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
                  & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
                  & " is invalid for a transient Darcy fluid."
                CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            localError="The setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a transient Darcy fluid."
            CALL FlagError(localError,err,error,*999)
        END SELECT

      !-----------------------------------------------------------------
      !   c a s e   d e f a u l t
      !-----------------------------------------------------------------
      CASE DEFAULT
        localError="The problem subtype of "//TRIM(NumberToVString(PROBLEM%SPECIFICATION(3),"*",err,error))// &
          & " does not equal a standard, quasistatic or ALE Darcy equation subtype."
        CALL FlagError(localError,err,error,*999)

      END SELECT
    ELSE
      CALL FlagError("Problem is not associated.",err,error,*999)
    ENDIF

!     EXITS("DARCY_EQUATION_PROBLEM_STANDARD_SETUP")
    EXITS("DARCY_EQUATION_PROBLEM_SETUP")
    RETURN
! 999 ERRORSEXITS("DARCY_EQUATION_PROBLEM_STANDARD_SETUP",err,error)
999 ERRORSEXITS("DARCY_EQUATION_PROBLEM_SETUP",err,error)
    RETURN 1
!   END SUBROUTINE DARCY_EQUATION_PROBLEM_STANDARD_SETUP
  END SUBROUTINE DARCY_EQUATION_PROBLEM_SETUP

  !
  !================================================================================================================================
  !

  !>Sets up the Darcy problem pre-solve.
  SUBROUTINE DARCY_EQUATION_PRE_SOLVE(CONTROL_LOOP,SOLVER,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS  !<A pointer to the solver equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES
    TYPE(SOLVER_MATRIX_TYPE), POINTER :: SOLVER_MATRIX
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equations_SET !<A pointer to the equations set
    TYPE(EquationsType), POINTER :: equations

    TYPE(SOLVER_TYPE), POINTER :: SOLVER_ALE_DARCY  !<A pointer to the solvers
    TYPE(VARYING_STRING) :: localError

    INTEGER(INTG) :: solver_matrix_idx


    ENTERS("DARCY_EQUATION_PRE_SOLVE",err,error,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a Darcy equation problem.",err,error,*999)
          END IF
          !--- Set 'SOLVER_NUMBER' depending on CONTROL_LOOP%PROBLEM%SPECIFICATION(3)
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_STANDARD_DARCY_SUBTYPE,PROBLEM_QUASISTATIC_DARCY_SUBTYPE,PROBLEM_TRANSIENT_DARCY_SUBTYPE)
              SOLVER_NUMBER_DARCY=1
            CASE(PROBLEM_ALE_DARCY_SUBTYPE,PROBLEM_PGM_DARCY_SUBTYPE,PROBLEM_PGM_TRANSIENT_DARCY_SUBTYPE)
              SOLVER_NUMBER_MAT_PROPERTIES=1
              SOLVER_NUMBER_DARCY=2
            CASE(PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE)
              SOLVER_NUMBER_SOLID=1
              SOLVER_NUMBER_DARCY=1
            CASE(PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE,PROBLEM_PGM_ELASTICITY_DARCY_SUBTYPE, &
              & PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)
              SOLVER_NUMBER_SOLID=1
              SOLVER_NUMBER_MAT_PROPERTIES=1
              SOLVER_NUMBER_DARCY=2
          END SELECT

          !--- Set explicitly 'SOLVER_MATRIX%UPDATE_MATRIX=.TRUE.'
          SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
          IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
            SOLVER_MAPPING=>SOLVER_equations%SOLVER_MAPPING
            IF(ASSOCIATED(SOLVER_MAPPING)) THEN
              SOLVER_MATRICES=>SOLVER_equations%SOLVER_MATRICES
              IF(ASSOCIATED(SOLVER_MATRICES)) THEN
                DO solver_matrix_idx=1,SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES
                  SOLVER_MATRIX=>SOLVER_MATRICES%matrices(solver_matrix_idx)%ptr
                  IF(ASSOCIATED(SOLVER_MATRIX)) THEN
                    SOLVER_MATRIX%UPDATE_MATRIX=.TRUE.
                  ELSE
                    CALL FlagError("Solver Matrix is not associated.",err,error,*999)
                  ENDIF
                ENDDO
              ELSE
                CALL FlagError("Solver Matrices is not associated.",err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Solver mapping is not associated.",err,error,*999)
            ENDIF
          ELSE
            CALL FlagError("Solver equations is not associated.",err,error,*999)
          ENDIF

          !--- pre_solve calls for various actions
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_STANDARD_DARCY_SUBTYPE)
              ! do nothing
            CASE(PROBLEM_QUASISTATIC_DARCY_SUBTYPE)
              ! do nothing
            CASE(PROBLEM_TRANSIENT_DARCY_SUBTYPE)
              CALL Darcy_PreSolveUpdateBoundaryConditions(CONTROL_LOOP,SOLVER,err,error,*999)
            CASE(PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE)

              IF((CONTROL_LOOP%LOOP_TYPE==PROBLEM_CONTROL_SIMPLE_TYPE.OR.CONTROL_LOOP%LOOP_TYPE==PROBLEM_CONTROL_TIME_LOOP_TYPE) &
                  & .AND.SOLVER%GLOBAL_NUMBER==SOLVER_NUMBER_DARCY) THEN
                !--- flags to ensure once-per-time-step output in conjunction with diagnostics
                idebug1 = .TRUE.
                idebug2 = .TRUE.
                idebug3 = .TRUE.

                NULLIFY(SOLVER_ALE_DARCY)
                CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,SOLVER_NUMBER_DARCY,SOLVER_ALE_DARCY,err,error,*999)
                EQUATIONS=>SOLVER_ALE_DARCY%SOLVER_equations%SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(1)%EQUATIONS
                IF(ASSOCIATED(EQUATIONS)) THEN
                 EQUATIONS_SET=>equations%equationsSet
                 IF(ASSOCIATED(EQUATIONS_SET)) THEN
                  IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
                   IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_INCOMP_ELAST_DARCY_ANALYTIC_DARCY)THEN
                   !call only analytic update and DO NOT call the other standard pre-solve routines as the mesh does not require deformation
                     CALL Darcy_PreSolveUpdateAnalyticValues(CONTROL_LOOP,SOLVER_ALE_DARCY,err,error,*999)
                   ENDIF
                  ELSE
                   !default
                     !--- 1.1 Transfer solid displacement to Darcy geometric field
                     CALL Darcy_PreSolveGetSolidDisplacement(CONTROL_LOOP,SOLVER_ALE_DARCY,err,error,*999)

                     !--- 1.2 Update the mesh (and calculate boundary velocities) PRIOR to solving for new material properties
                     CALL DARCY_EQUATION_PRE_SOLVE_ALE_UPDATE_MESH(CONTROL_LOOP,SOLVER_ALE_DARCY,err,error,*999)

  ! ! !                 !i n   p r i n c i p l e   c u r r e n t l y   d o   n o t   n e e d   t o   u p d a t e   B C s
  ! ! !              !unless:
  ! ! !              !--- 1.3 Apply both normal and moving mesh boundary conditions, OR:
  ! ! !              !--- 1.3 (Iteratively) Render the boundary impermeable (ellipsoid, general curvilinear mesh)
  ! ! !              CALL Darcy_PreSolveUpdateBoundaryConditions(CONTROL_LOOP,SOLVER_ALE_DARCY,err,error,*999)
                  ENDIF
                 ENDIF
                ENDIF
              END IF
            CASE(PROBLEM_ALE_DARCY_SUBTYPE,PROBLEM_PGM_DARCY_SUBTYPE,PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE, &
              & PROBLEM_PGM_ELASTICITY_DARCY_SUBTYPE,PROBLEM_PGM_TRANSIENT_DARCY_SUBTYPE, &
              & PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)

              IF((CONTROL_LOOP%LOOP_TYPE==PROBLEM_CONTROL_SIMPLE_TYPE.OR.CONTROL_LOOP%LOOP_TYPE==PROBLEM_CONTROL_TIME_LOOP_TYPE) &
                  & .AND.SOLVER%GLOBAL_NUMBER==SOLVER_NUMBER_MAT_PROPERTIES) THEN
                !--- flags to ensure once-per-time-step output in conjunction with diagnostics
                idebug1 = .TRUE.
                idebug2 = .TRUE.
                idebug3 = .TRUE.

                NULLIFY(SOLVER_ALE_DARCY)
                CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,SOLVER_NUMBER_DARCY,SOLVER_ALE_DARCY,err,error,*999)
                EQUATIONS=>SOLVER_ALE_DARCY%SOLVER_equations%SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(1)%EQUATIONS
                IF(ASSOCIATED(EQUATIONS)) THEN
                 EQUATIONS_SET=>equations%equationsSet
                 IF(ASSOCIATED(EQUATIONS_SET)) THEN
                  IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
                   IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_INCOMP_ELAST_DARCY_ANALYTIC_DARCY)THEN
                   !call only analytic update and DO NOT call the other standard pre-solve routines as the mesh does not require deformation
                     CALL Darcy_PreSolveUpdateAnalyticValues(CONTROL_LOOP,SOLVER_ALE_DARCY,err,error,*999)
                   ENDIF
                  ELSE
                   !default
                     !--- 1.1 Transfer solid displacement to Darcy geometric field
                     CALL Darcy_PreSolveGetSolidDisplacement(CONTROL_LOOP,SOLVER_ALE_DARCY,err,error,*999)

                     !--- 1.2 Update the mesh (and calculate boundary velocities) PRIOR to solving for new material properties
                     CALL DARCY_EQUATION_PRE_SOLVE_ALE_UPDATE_MESH(CONTROL_LOOP,SOLVER_ALE_DARCY,err,error,*999)

  ! ! !                 !i n   p r i n c i p l e   c u r r e n t l y   d o   n o t   n e e d   t o   u p d a t e   B C s
  ! ! !                 !--- 1.3 Apply both normal and moving mesh boundary conditions
  ! ! !                 CALL Darcy_PreSolveUpdateBoundaryConditions(CONTROL_LOOP,SOLVER_ALE_DARCY,err,error,*999)
                  ENDIF
                 ENDIF
                ENDIF
              ELSE IF((CONTROL_LOOP%LOOP_TYPE==PROBLEM_CONTROL_SIMPLE_TYPE.OR. &
                    & CONTROL_LOOP%LOOP_TYPE==PROBLEM_CONTROL_TIME_LOOP_TYPE).AND.SOLVER%GLOBAL_NUMBER==SOLVER_NUMBER_DARCY) THEN
! ! !                 !n o t   f o r   n o w   ! ! !
! ! !                 !--- 2.1 Update the material field
! ! !                 CALL Darcy_PreSolveUpdateMatrixProperties(CONTROL_LOOP,SOLVER,err,error,*999)
              END IF
            CASE DEFAULT
              localError="Problem subtype "//TRIM(NumberToVString(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",err,error))// &
                & " is not valid for a Darcy fluid type of a fluid mechanics problem class."
              CALL FlagError(localError,err,error,*999)
          END SELECT
        ELSE
          CALL FlagError("Problem is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Solver is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",err,error,*999)
    ENDIF

    EXITS("DARCY_EQUATION_PRE_SOLVE")
    RETURN
999 ERRORSEXITS("DARCY_EQUATION_PRE_SOLVE",err,error)
    RETURN 1
  END SUBROUTINE DARCY_EQUATION_PRE_SOLVE

  !
  !================================================================================================================================
  !

  SUBROUTINE DARCY_CONTROL_TIME_LOOP_PRE_LOOP(CONTROL_LOOP,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the time control loop for the Darcy  problem
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP_DARCY
    TYPE(SOLVER_TYPE), POINTER :: SOLVER_DARCY

    ENTERS("DARCY_CONTROL_TIME_LOOP_PRE_LOOP",err,error,*999)

    !Get the solver for the Darcy problem
    NULLIFY(SOLVER_DARCY)
    NULLIFY(CONTROL_LOOP_DARCY)
    IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
      CALL FlagError("Problem specification is not allocated.",err,error,*999)
    ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
      CALL FlagError("Problem specification must have three entries for a Darcy equation problem.",err,error,*999)
    END IF
    SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
      CASE(PROBLEM_STANDARD_DARCY_SUBTYPE,PROBLEM_QUASISTATIC_DARCY_SUBTYPE,PROBLEM_TRANSIENT_DARCY_SUBTYPE)
        !SOLVER_NUMBER_DARCY has to be set here so that store_reference_data and store_previous_data have access to it
        SOLVER_NUMBER_DARCY=1
        CALL CONTROL_LOOP_GET(CONTROL_LOOP,[CONTROL_LOOP_NODE],CONTROL_LOOP_DARCY,err,error,*999)
        CALL SOLVERS_SOLVER_GET(CONTROL_LOOP_DARCY%SOLVERS,SOLVER_NUMBER_DARCY,SOLVER_DARCY,err,error,*999)
      CASE(PROBLEM_ALE_DARCY_SUBTYPE,PROBLEM_PGM_DARCY_SUBTYPE,PROBLEM_PGM_TRANSIENT_DARCY_SUBTYPE)
        SOLVER_NUMBER_MAT_PROPERTIES=1
        SOLVER_NUMBER_DARCY=2
        CALL CONTROL_LOOP_GET(CONTROL_LOOP,[CONTROL_LOOP_NODE],CONTROL_LOOP_DARCY,err,error,*999)
        CALL SOLVERS_SOLVER_GET(CONTROL_LOOP_DARCY%SOLVERS,SOLVER_NUMBER_DARCY,SOLVER_DARCY,err,error,*999)
      CASE(PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE,PROBLEM_PGM_ELASTICITY_DARCY_SUBTYPE)
        SOLVER_NUMBER_SOLID=1
        SOLVER_NUMBER_MAT_PROPERTIES=1
        SOLVER_NUMBER_DARCY=2
        CALL CONTROL_LOOP_GET(CONTROL_LOOP,[2,CONTROL_LOOP_NODE],CONTROL_LOOP_DARCY,err,error,*999)
        CALL SOLVERS_SOLVER_GET(CONTROL_LOOP_DARCY%SOLVERS,SOLVER_NUMBER_DARCY,SOLVER_DARCY,err,error,*999)
      CASE(PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE)
        SOLVER_NUMBER_SOLID=1
        SOLVER_NUMBER_DARCY=1
        CALL CONTROL_LOOP_GET(CONTROL_LOOP,[1,2,CONTROL_LOOP_NODE],CONTROL_LOOP_DARCY,err,error,*999)
        CALL SOLVERS_SOLVER_GET(CONTROL_LOOP_DARCY%SOLVERS,SOLVER_NUMBER_DARCY,SOLVER_DARCY,err,error,*999)
      CASE(PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)
        SOLVER_NUMBER_SOLID=1
        SOLVER_NUMBER_MAT_PROPERTIES=1
        SOLVER_NUMBER_DARCY=2
        CALL CONTROL_LOOP_GET(CONTROL_LOOP,[1,2,CONTROL_LOOP_NODE],CONTROL_LOOP_DARCY,err,error,*999)
        CALL SOLVERS_SOLVER_GET(CONTROL_LOOP_DARCY%SOLVERS,SOLVER_NUMBER_DARCY,SOLVER_DARCY,err,error,*999)
    END SELECT

    !If this is the first time step then store reference data
    IF(CONTROL_LOOP%TIME_LOOP%ITERATION_NUMBER==0) THEN
      IF(CONTROL_LOOP%outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
        CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,'== Storing reference data',err,error,*999)
      ENDIF
      CALL Darcy_PreSolveStoreReferenceData(CONTROL_LOOP,SOLVER_DARCY,err,error,*999)
    ENDIF

    !Store data of previous time step (mesh position); executed once per time step before subiteration
    IF(CONTROL_LOOP%outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
      CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,'== Storing previous data',err,error,*999)
    ENDIF
    CALL Darcy_PreSolveStorePreviousData(CONTROL_LOOP,SOLVER_DARCY,err,error,*999)

    EXITS("DARCY_CONTROL_TIME_LOOP_PRE_LOOP")
    RETURN
999 ERRORSEXITS("DARCY_CONTROL_TIME_LOOP_PRE_LOOP",err,error)
    RETURN 1
  END SUBROUTINE DARCY_CONTROL_TIME_LOOP_PRE_LOOP

  !
  !================================================================================================================================
  !

  !>Store some reference data for ALE Darcy problem
  SUBROUTINE Darcy_PreSolveStoreReferenceData(CONTROL_LOOP,SOLVER,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solvers
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(FIELD_TYPE), POINTER :: dependentField, geometricField
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS  !<A pointer to the solver equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equations_SET !<A pointer to the equations set
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: localError

    REAL(DP) :: ALPHA
    REAL(DP), POINTER :: INITIAL_VALUES(:)

    INTEGER(INTG) :: FIELD_VAR_TYPE
    INTEGER(INTG) :: NDOFS_TO_PRINT,equations_set_idx


    ENTERS("Darcy_PreSolveStoreReferenceData",err,error,*999)

    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(SOLVER_MAPPING)
    NULLIFY(EQUATIONS_SET)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(SOLVER%GLOBAL_NUMBER==SOLVER_NUMBER_DARCY) THEN
          IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
            IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
              CALL FlagError("Problem specification is not allocated.",err,error,*999)
            ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
              CALL FlagError("Problem specification must have three entries for a Darcy equation problem.",err,error,*999)
            END IF
            SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
              CASE(PROBLEM_STANDARD_DARCY_SUBTYPE)
                ! do nothing
              CASE(PROBLEM_QUASISTATIC_DARCY_SUBTYPE)
                ! do nothing
              CASE(PROBLEM_TRANSIENT_DARCY_SUBTYPE)
                ! do nothing
              CASE(PROBLEM_ALE_DARCY_SUBTYPE,PROBLEM_PGM_DARCY_SUBTYPE,PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE, &
                & PROBLEM_PGM_ELASTICITY_DARCY_SUBTYPE,PROBLEM_PGM_TRANSIENT_DARCY_SUBTYPE, &
                & PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE,PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)
                SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                  SOLVER_MAPPING=>SOLVER_equations%SOLVER_MAPPING
                  IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                   DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                    EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%ptr
                    IF(ASSOCIATED(EQUATIONS_SET)) THEN
                      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
                        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
                      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
                        CALL FlagError("Equations set specification must have three entries for a Darcy type equations set.", &
                          & err,error,*999)
                      END IF
                      SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
                        CASE(EQUATIONS_SET_STANDARD_DARCY_SUBTYPE)
                          ! do nothing
                        CASE(EQUATIONS_SET_QUASISTATIC_DARCY_SUBTYPE)
                          ! do nothing
                        CASE(EQUATIONS_SET_TRANSIENT_DARCY_SUBTYPE)
                          ! do nothing
                        CASE(EQUATIONS_SET_ALE_DARCY_SUBTYPE,EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE, &
                          & EQUATIONS_SET_TRANSIENT_ALE_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
                          & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE, &
                          & EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
                          dependentField=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                          geometricField=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                          IF(ASSOCIATED(dependentField).AND.ASSOCIATED(geometricField)) THEN
                            !--- Store the initial (= reference) GEOMETRY field values
                            ALPHA = 1.0_DP
                            CALL FIELD_PARAMETER_SETS_COPY(geometricField,FIELD_U_VARIABLE_TYPE, &
                              & FIELD_VALUES_SET_TYPE,FIELD_INITIAL_VALUES_SET_TYPE,ALPHA,err,error,*999)
                            NULLIFY(vectorEquations)
                            CALL Equations_VectorEquationsGet(EQUATIONS_SET%equations,vectorEquations,err,error,*999)
                            vectorMapping=>vectorEquations%vectorMapping
                            IF(ASSOCIATED(vectorMapping)) THEN

                              SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
                              CASE(EQUATIONS_SET_ALE_DARCY_SUBTYPE,EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE)
                              FIELD_VARIABLE=>vectorMapping%linearMapping%equationsMatrixToVarMaps(1)%VARIABLE
                              ! '1' associated with linear matrix
                              CASE(EQUATIONS_SET_TRANSIENT_ALE_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
                                  & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE, &
                                  & EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
                                FIELD_VARIABLE=>vectorMapping%dynamicMapping%equationsMatrixToVarMaps(1)%VARIABLE
                              END SELECT

                              IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                                FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
                                !--- Store the initial DEPENDENT field values
                                ALPHA = 1.0_DP
                                CALL FIELD_PARAMETER_SETS_COPY(dependentField,FIELD_VAR_TYPE, &
                                  & FIELD_VALUES_SET_TYPE,FIELD_INITIAL_VALUES_SET_TYPE,ALPHA,err,error,*999)

                                IF(DIAGNOSTICS1) THEN
                                  NULLIFY(INITIAL_VALUES)
                                  CALL FIELD_PARAMETER_SET_DATA_GET(dependentField,FIELD_VAR_TYPE, &
                                    & FIELD_INITIAL_VALUES_SET_TYPE,INITIAL_VALUES,err,error,*999)
                                  NDOFS_TO_PRINT = SIZE(INITIAL_VALUES,1)
                                  CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NDOFS_TO_PRINT,NDOFS_TO_PRINT,NDOFS_TO_PRINT,&
                                    & INITIAL_VALUES, &
                                    & '(" dependentField,FIELD_U_VARIABLE_TYPE,FIELD_INITIAL_VALUES_SET_TYPE = ",4(X,E13.6))', &
                                    & '4(4(X,E13.6))',err,error,*999)
                                  CALL FIELD_PARAMETER_SET_DATA_RESTORE(dependentField,FIELD_VAR_TYPE, &
                                    & FIELD_INITIAL_VALUES_SET_TYPE,INITIAL_VALUES,err,error,*999)
                                ENDIF
                              ELSE
                                CALL FlagError("FIELD_VAR_TYPE is not associated.",err,error,*999)
                              ENDIF
                            ELSE
                              CALL FlagError("vectorMapping is not associated.",err,error,*999)
                            ENDIF
                          ELSE
                            CALL FlagError("Dependent field and / or geometric field is / are not associated.",err,error,*999)
                          ENDIF
                        CASE DEFAULT
                          localError="Equations set subtype " &
                            & //TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
                            & " is not valid for a Darcy equation fluid type of a fluid mechanics problem class."
                          CALL FlagError(localError,err,error,*999)
                      END SELECT
                    ELSE
                      CALL FlagError("Equations set is not associated.",err,error,*999)
                    ENDIF
                   ENDDO
                  ELSE
                    CALL FlagError("Solver mapping is not associated.",err,error,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Solver equations is not associated.",err,error,*999)
                ENDIF
              CASE DEFAULT
                localError="Problem subtype "//TRIM(NumberToVString(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",err,error))// &
                  & " is not valid for a Darcy equation fluid type of a fluid mechanics problem class."
                CALL FlagError(localError,err,error,*999)
            END SELECT
          ELSE
            CALL FlagError("Problem is not associated.",err,error,*999)
          ENDIF
        ELSE
          ! do nothing
        ENDIF
      ELSE
        CALL FlagError("Solver is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",err,error,*999)
    ENDIF


    EXITS("Darcy_PreSolveStoreReferenceData")
    RETURN
999 ERRORSEXITS("Darcy_PreSolveStoreReferenceData",err,error)
    RETURN 1
    
  END SUBROUTINE Darcy_PreSolveStoreReferenceData

  !
  !================================================================================================================================
  !

  !>Store data of previous time step (mesh position) for ALE Darcy problem
  SUBROUTINE Darcy_PreSolveStorePreviousData(CONTROL_LOOP,SOLVER,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solvers
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(FIELD_TYPE), POINTER :: geometricField
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS  !<A pointer to the solver equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equations_SET !<A pointer to the equations set
    TYPE(VARYING_STRING) :: localError

    REAL(DP) :: ALPHA


    ENTERS("Darcy_PreSolveStorePreviousData",err,error,*999)

    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(SOLVER_MAPPING)
    NULLIFY(EQUATIONS_SET)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(SOLVER%GLOBAL_NUMBER==SOLVER_NUMBER_DARCY) THEN
          IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
            IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
              CALL FlagError("Problem specification is not allocated.",err,error,*999)
            ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
              CALL FlagError("Problem specification must have three entries for a Darcy equation problem.",err,error,*999)
            END IF
            SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
              CASE(PROBLEM_STANDARD_DARCY_SUBTYPE)
                ! do nothing
              CASE(PROBLEM_QUASISTATIC_DARCY_SUBTYPE)
                ! do nothing
              CASE(PROBLEM_TRANSIENT_DARCY_SUBTYPE)
                ! do nothing
              CASE(PROBLEM_ALE_DARCY_SUBTYPE,PROBLEM_PGM_DARCY_SUBTYPE,PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE, &
                & PROBLEM_PGM_ELASTICITY_DARCY_SUBTYPE,PROBLEM_PGM_TRANSIENT_DARCY_SUBTYPE, &
                & PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE,PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)
                SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                  SOLVER_MAPPING=>SOLVER_equations%SOLVER_MAPPING
                  IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                    EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(1)%ptr
                    IF(ASSOCIATED(EQUATIONS_SET)) THEN
                      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
                        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
                      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
                        CALL FlagError("Equations set specification must have three entries for a Darcy type equations set.", &
                          & err,error,*999)
                      END IF
                      SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
                        CASE(EQUATIONS_SET_STANDARD_DARCY_SUBTYPE)
                          ! do nothing
                        CASE(EQUATIONS_SET_QUASISTATIC_DARCY_SUBTYPE)
                          ! do nothing
                        CASE(EQUATIONS_SET_TRANSIENT_DARCY_SUBTYPE)
                          ! do nothing
                        CASE(EQUATIONS_SET_ALE_DARCY_SUBTYPE,EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE, &
                          & EQUATIONS_SET_TRANSIENT_ALE_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
                          & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE, &
                          & EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
                          geometricField=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                          IF(ASSOCIATED(geometricField)) THEN
                            !--- Store the GEOMETRY field values of the previous time step
                            ALPHA = 1.0_DP
                            CALL FIELD_PARAMETER_SETS_COPY(geometricField,FIELD_U_VARIABLE_TYPE, &
                              & FIELD_VALUES_SET_TYPE,FIELD_PREVIOUS_VALUES_SET_TYPE,ALPHA,err,error,*999)
                          ELSE
                            CALL FlagError("Dependent field and / or geometric field is / are not associated.",err,error,*999)
                          ENDIF
                        CASE DEFAULT
                          localError="Equations set subtype " &
                            & //TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
                            & " is not valid for a Darcy equation fluid type of a fluid mechanics problem class."
                          CALL FlagError(localError,err,error,*999)
                      END SELECT
                    ELSE
                      CALL FlagError("Equations set is not associated.",err,error,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Solver mapping is not associated.",err,error,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Solver equations is not associated.",err,error,*999)
                ENDIF
              CASE DEFAULT
                localError="Problem subtype "//TRIM(NumberToVString(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",err,error))// &
                  & " is not valid for a Darcy equation fluid type of a fluid mechanics problem class."
                CALL FlagError(localError,err,error,*999)
            END SELECT
          ELSE
            CALL FlagError("Problem is not associated.",err,error,*999)
          ENDIF
        ELSE
          ! do nothing
        ENDIF
      ELSE
        CALL FlagError("Solver is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",err,error,*999)
    ENDIF

    EXITS("Darcy_PreSolveStorePreviousData")
    RETURN
999 ERRORSEXITS("Darcy_PreSolveStorePreviousData",err,error)
    RETURN 1
    
  END SUBROUTINE Darcy_PreSolveStorePreviousData

  !
  !================================================================================================================================
  !

  !>Update mesh position and velocity for ALE Darcy problem
  SUBROUTINE DARCY_EQUATION_PRE_SOLVE_ALE_UPDATE_MESH(CONTROL_LOOP,SOLVER,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solvers
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(FIELD_TYPE), POINTER :: geometricField
    TYPE(SOLVER_TYPE), POINTER :: SOLVER_ALE_DARCY !<A pointer to the solvers
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS  !<A pointer to the solver equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equations_SET !<A pointer to the equations set
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_TIME_LOOP !<A pointer to the control time loop
    TYPE(VARYING_STRING) :: localError

    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT,ALPHA
    REAL(DP), POINTER :: MESH_DISPLACEMENT_VALUES(:)

    INTEGER(INTG) :: dof_number,NUMBER_OF_DOFS,NDOFS_TO_PRINT,loop_idx
    INTEGER(INTG) :: PROBLEM_SUBTYPE

    ENTERS("DARCY_EQUATION_PRE_SOLVE_ALE_UPDATE_MESH",err,error,*999)

    NULLIFY(SOLVER_ALE_DARCY)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(SOLVER_MAPPING)
    NULLIFY(EQUATIONS_SET)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      CONTROL_TIME_LOOP=>CONTROL_LOOP
      DO loop_idx=1,CONTROL_LOOP%CONTROL_LOOP_LEVEL
        IF(CONTROL_TIME_LOOP%LOOP_TYPE==PROBLEM_CONTROL_TIME_LOOP_TYPE) THEN
          CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_TIME_LOOP,CURRENT_TIME,TIME_INCREMENT,err,error,*999)
          EXIT
        ENDIF
        IF (ASSOCIATED(CONTROL_LOOP%PARENT_LOOP)) THEN
          CONTROL_TIME_LOOP=>CONTROL_TIME_LOOP%PARENT_LOOP
        ELSE
          CALL FlagError("Could not find a time control loop.",err,error,*999)
        ENDIF
      ENDDO

      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a Darcy equation problem.",err,error,*999)
          END IF
          PROBLEM_SUBTYPE=CONTROL_LOOP%PROBLEM%SPECIFICATION(3)
          SELECT CASE(PROBLEM_SUBTYPE)
            CASE(PROBLEM_STANDARD_DARCY_SUBTYPE)
              ! do nothing
            CASE(PROBLEM_QUASISTATIC_DARCY_SUBTYPE)
              ! do nothing
            CASE(PROBLEM_TRANSIENT_DARCY_SUBTYPE)
              ! do nothing
            CASE(PROBLEM_ALE_DARCY_SUBTYPE,PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE,PROBLEM_PGM_DARCY_SUBTYPE, &
              & PROBLEM_PGM_ELASTICITY_DARCY_SUBTYPE,PROBLEM_PGM_TRANSIENT_DARCY_SUBTYPE, &
              & PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE,PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)
              IF(SOLVER%GLOBAL_NUMBER==SOLVER_NUMBER_DARCY) THEN
                SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                  SOLVER_MAPPING=>SOLVER_equations%SOLVER_MAPPING
                  IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                    EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(1)%ptr
                    IF(ASSOCIATED(EQUATIONS_SET)) THEN
                      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
                        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
                      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
                        CALL FlagError("Equations set specification must have three entries for a Darcy type equations set.", &
                          & err,error,*999)
                      END IF
                      SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
                        CASE(EQUATIONS_SET_STANDARD_DARCY_SUBTYPE)
                          ! do nothing
                        CASE(EQUATIONS_SET_QUASISTATIC_DARCY_SUBTYPE)
                          ! do nothing
                        CASE(EQUATIONS_SET_ALE_DARCY_SUBTYPE,EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE, &
                          & EQUATIONS_SET_TRANSIENT_ALE_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
                          & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE, &
                          & EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
                          IF(SOLVER%outputType>=SOLVER_PROGRESS_OUTPUT) THEN
                            CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Darcy update mesh ... ",err,error,*999)
                          ENDIF
                          geometricField=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                          IF(ASSOCIATED(geometricField)) THEN
                            !--- First, get pointer to mesh displacement values
                            NULLIFY(MESH_DISPLACEMENT_VALUES)
                            CALL FIELD_PARAMETER_SET_DATA_GET(geometricField,FIELD_U_VARIABLE_TYPE, &
                              & FIELD_MESH_DISPLACEMENT_SET_TYPE,MESH_DISPLACEMENT_VALUES,err,error,*999)
                            IF(DIAGNOSTICS1) THEN
                              NDOFS_TO_PRINT = SIZE(MESH_DISPLACEMENT_VALUES,1)
                              CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NDOFS_TO_PRINT,NDOFS_TO_PRINT,NDOFS_TO_PRINT,&
                                & MESH_DISPLACEMENT_VALUES,'(" MESH_DISPLACEMENT_VALUES = ",3(X,E13.6))','3(3(X,E13.6))', &
                                & err,error,*999)
                            ENDIF

                            NUMBER_OF_DOFS = geometricField%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%ptr%NUMBER_OF_DOFS

                            IF(PROBLEM_SUBTYPE==PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE &
                                & .OR. PROBLEM_SUBTYPE==PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE &
                                & .OR. PROBLEM_SUBTYPE==PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE) THEN
                              !--- Don't update geometric field here, this is done in
                              !    darcy_equation_pre_solve_get_solid_displacement for these problems, but
                              !    needs to be made consistent between the different problem types
                            ELSE
                              !--- Second, update geometric field
                              DO dof_number=1,NUMBER_OF_DOFS
                                CALL FIELD_PARAMETER_SET_ADD_LOCAL_DOF(geometricField, &
                                  & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,dof_number, &
                                  & MESH_DISPLACEMENT_VALUES(dof_number), &
                                  & err,error,*999)
                              END DO
                              CALL FIELD_PARAMETER_SET_UPDATE_START(geometricField, &
                                & FIELD_U_VARIABLE_TYPE, FIELD_VALUES_SET_TYPE,err,error,*999)
                              CALL FIELD_PARAMETER_SET_UPDATE_FINISH(geometricField, &
                                & FIELD_U_VARIABLE_TYPE, FIELD_VALUES_SET_TYPE,err,error,*999)
                            ENDIF

                            !--- Third, use displacement values to calculate velocity values
                            ALPHA=1.0_DP/TIME_INCREMENT
                            CALL FIELD_PARAMETER_SETS_COPY(geometricField,FIELD_U_VARIABLE_TYPE, &
                              & FIELD_MESH_DISPLACEMENT_SET_TYPE,FIELD_MESH_VELOCITY_SET_TYPE,ALPHA,err,error,*999)
                            CALL FIELD_PARAMETER_SET_DATA_RESTORE(geometricField,FIELD_U_VARIABLE_TYPE, &
                              & FIELD_MESH_DISPLACEMENT_SET_TYPE,MESH_DISPLACEMENT_VALUES,err,error,*999)
                          ELSE
                            CALL FlagError("Geometric field is not associated.",err,error,*999)
                          ENDIF
                        CASE DEFAULT
                          localError="Equations set subtype " &
                            & //TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
                            & " is not valid for a Darcy equation fluid type of a fluid mechanics problem class."
                          CALL FlagError(localError,err,error,*999)
                      END SELECT
                    ELSE
                      CALL FlagError("Equations set is not associated.",err,error,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Solver mapping is not associated.",err,error,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Solver equations is not associated.",err,error,*999)
                ENDIF
              ELSE
                ! do nothing
              ENDIF
            CASE DEFAULT
              localError="Problem subtype "//TRIM(NumberToVString(PROBLEM_SUBTYPE,"*",err,error))// &
                & " is not valid for a Darcy equation fluid type of a fluid mechanics problem class."
              CALL FlagError(localError,err,error,*999)
          END SELECT
        ELSE
          CALL FlagError("Problem is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Solver is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",err,error,*999)
    ENDIF

    EXITS("DARCY_EQUATION_PRE_SOLVE_ALE_UPDATE_MESH")
    RETURN
999 ERRORSEXITS("DARCY_EQUATION_PRE_SOLVE_ALE_UPDATE_MESH",err,error)
    RETURN 1
  END SUBROUTINE DARCY_EQUATION_PRE_SOLVE_ALE_UPDATE_MESH

  !
  !================================================================================================================================
  !

  !>Update boundary conditions for Darcy equation pre solve
  SUBROUTINE Darcy_PreSolveUpdateBoundaryConditions(CONTROL_LOOP,SOLVER,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS  !<A pointer to the solver equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equations_SET !<A pointer to the equations set
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: localError
    TYPE(FIELD_TYPE), POINTER :: dependentField, geometricField
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_TIME_LOOP

    REAL(DP), POINTER :: MESH_VELOCITY_VALUES(:)
    REAL(DP), POINTER :: INITIAL_VALUES(:)
    REAL(DP), POINTER :: DUMMY_VALUES1(:)
    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT
    REAL(DP) :: PRESSURE

    INTEGER(INTG) :: FIELD_VAR_TYPE
    INTEGER(INTG) :: BOUNDARY_CONDITION_CHECK_VARIABLE
    INTEGER(INTG) :: dof_number,NUMBER_OF_DOFS,loop_idx
    INTEGER(INTG) :: NDOFS_TO_PRINT

    ENTERS("Darcy_PreSolveUpdateBoundaryConditions",err,error,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      CONTROL_TIME_LOOP=>CONTROL_LOOP
      DO loop_idx=1,CONTROL_LOOP%CONTROL_LOOP_LEVEL
        IF(CONTROL_TIME_LOOP%LOOP_TYPE==PROBLEM_CONTROL_TIME_LOOP_TYPE) THEN
          CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_TIME_LOOP,CURRENT_TIME,TIME_INCREMENT,err,error,*999)
          EXIT
        ENDIF
        IF (ASSOCIATED(CONTROL_LOOP%PARENT_LOOP)) THEN
          CONTROL_TIME_LOOP=>CONTROL_TIME_LOOP%PARENT_LOOP
        ELSE
          CALL FlagError("Could not find a time control loop.",err,error,*999)
        ENDIF
      ENDDO

      IF(ASSOCIATED(SOLVER)) THEN
        IF(SOLVER%GLOBAL_NUMBER==SOLVER_NUMBER_DARCY) THEN
          IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
            IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
              CALL FlagError("Problem specification is not allocated.",err,error,*999)
            ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
              CALL FlagError("Problem specification must have three entries for a Darcy equation problem.",err,error,*999)
            END IF
            SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
              CASE(PROBLEM_STANDARD_DARCY_SUBTYPE)
                ! do nothing
              CASE(PROBLEM_QUASISTATIC_DARCY_SUBTYPE)
                ! do nothing
              CASE(PROBLEM_TRANSIENT_DARCY_SUBTYPE)
                ! do nothing
              CASE(PROBLEM_ALE_DARCY_SUBTYPE,PROBLEM_PGM_DARCY_SUBTYPE,PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE, &
                & PROBLEM_PGM_ELASTICITY_DARCY_SUBTYPE,PROBLEM_PGM_TRANSIENT_DARCY_SUBTYPE, &
                & PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE,PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)
                SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                  SOLVER_MAPPING=>SOLVER_equations%SOLVER_MAPPING
                  IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                    EQUATIONS=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(1)%EQUATIONS
                    IF(ASSOCIATED(EQUATIONS)) THEN
                      EQUATIONS_SET=>equations%equationsSet
                      IF(ASSOCIATED(EQUATIONS_SET)) THEN
                        IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
                          CALL FlagError("Equations set specification is not allocated.",err,error,*999)
                        ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
                          CALL FlagError("Equations set specification must have three entries for a Darcy type equations set.", &
                            & err,error,*999)
                        END IF
                        SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
                          CASE(EQUATIONS_SET_STANDARD_DARCY_SUBTYPE)
                            ! do nothing
                          CASE(EQUATIONS_SET_QUASISTATIC_DARCY_SUBTYPE)
                            ! do nothing
                          CASE(EQUATIONS_SET_ALE_DARCY_SUBTYPE,EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE, &
                            & EQUATIONS_SET_TRANSIENT_ALE_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
                            & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE, &
                            & EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
                            IF(SOLVER%outputType>=SOLVER_PROGRESS_OUTPUT) THEN
                              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Darcy update boundary conditions ... ",err,error,*999)
                            ENDIF
                            dependentField=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                            geometricField=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                            IF(ASSOCIATED(dependentField).AND.ASSOCIATED(geometricField)) THEN
                              BOUNDARY_CONDITIONS=>SOLVER_equations%BOUNDARY_CONDITIONS
                              IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
                                vectorMapping=>EQUATIONS_SET%equations%vectorEquations%vectorMapping
                                IF(ASSOCIATED(vectorMapping)) THEN

                                  SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
                                  CASE(EQUATIONS_SET_ALE_DARCY_SUBTYPE, &
                                    & EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE)
                                  FIELD_VARIABLE=>vectorMapping%linearMapping%equationsMatrixToVarMaps(1)%VARIABLE
                                  ! '1' associated with linear matrix
                                  CASE(EQUATIONS_SET_TRANSIENT_ALE_DARCY_SUBTYPE, &
                                      & EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
                                      & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE, &
                                      & EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
                                    FIELD_VARIABLE=>vectorMapping%dynamicMapping%equationsMatrixToVarMaps(1)%VARIABLE
                                  END SELECT

                                  IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                                    FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE

                                    CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,FIELD_VARIABLE, &
                                      & BOUNDARY_CONDITIONS_VARIABLE,err,error,*999)
                                    IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
                                      NULLIFY(MESH_VELOCITY_VALUES)
                                      CALL FIELD_PARAMETER_SET_DATA_GET(geometricField,FIELD_U_VARIABLE_TYPE, &
                                        & FIELD_MESH_VELOCITY_SET_TYPE,MESH_VELOCITY_VALUES,err,error,*999)
                                      NULLIFY(INITIAL_VALUES)
                                      CALL FIELD_PARAMETER_SET_DATA_GET(dependentField,FIELD_VAR_TYPE, &
                                        & FIELD_INITIAL_VALUES_SET_TYPE,INITIAL_VALUES,err,error,*999)
                                      IF(DIAGNOSTICS1) THEN
                                        NULLIFY( DUMMY_VALUES1 )
                                        CALL FIELD_PARAMETER_SET_DATA_GET(dependentField,FIELD_VAR_TYPE, &
                                          & FIELD_VALUES_SET_TYPE,DUMMY_VALUES1,err,error,*999)
                                        NDOFS_TO_PRINT = SIZE(DUMMY_VALUES1,1)
                                        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NDOFS_TO_PRINT,NDOFS_TO_PRINT, &
                                          & NDOFS_TO_PRINT,DUMMY_VALUES1, &
                                          & '(" dependentField,FIELD_VAR_TYPE,FIELD_VALUES_SET_TYPE (before) = ",4(X,E13.6))', &
                                          & '4(4(X,E13.6))',err,error,*999)
                                      ENDIF
                                      NUMBER_OF_DOFS = dependentField%VARIABLE_TYPE_MAP(FIELD_VAR_TYPE)%ptr%NUMBER_OF_DOFS
                                      DO dof_number=1,NUMBER_OF_DOFS
                                        BOUNDARY_CONDITION_CHECK_VARIABLE=BOUNDARY_CONDITIONS_VARIABLE% &
                                          & CONDITION_TYPES(dof_number)
                                        IF(BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_MOVED_WALL) THEN
                                          !--- Reset boundary condition to the initial normal-velocity boundary condition
                                          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(dependentField, &
                                            & FIELD_VAR_TYPE,FIELD_VALUES_SET_TYPE,dof_number, &
                                            & INITIAL_VALUES(dof_number),err,error,*999)
                                          !--- Add the velocity of the moving boundary on top of the initial boundary condition
                                          !! === If we solve in terms of Darcy flow vector, then do not add mesh velocity === !!
                                          !! === The BC is kept to the initial BC, for instance: null-flux                === !!
!                                          CALL FIELD_PARAMETER_SET_ADD_LOCAL_DOF(dependentField, &
!                                            & FIELD_VAR_TYPE,FIELD_VALUES_SET_TYPE,dof_number, &
!                                            & MESH_VELOCITY_VALUES(dof_number),err,error,*999)
!                                            ! dependent field      ( V_u, V_v, V_w, P_p )
!                                            ! MESH_VELOCITY_VALUES ( V_u, V_v, V_w )

                                        ELSE IF( BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_FIXED .AND. &
                                          & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE) THEN
                                          !\ToDo: Check component number; this way we can also apply it to velocity
                                          !--- Set the time-dependent pressure BC
                                          PRESSURE = INITIAL_VALUES(dof_number) * (1.0_DP - exp(- CURRENT_TIME**2.0_DP / 0.25_DP))

                                          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(dependentField, &
                                            & FIELD_VAR_TYPE,FIELD_VALUES_SET_TYPE,dof_number, &
                                            & PRESSURE,err,error,*999)
                                        ELSE
                                          ! do nothing
                                        END IF
                                      END DO
                                      CALL FIELD_PARAMETER_SET_UPDATE_START(dependentField, &
                                        & FIELD_VAR_TYPE, FIELD_VALUES_SET_TYPE,err,error,*999)
                                      CALL FIELD_PARAMETER_SET_UPDATE_FINISH(dependentField, &
                                        & FIELD_VAR_TYPE, FIELD_VALUES_SET_TYPE,err,error,*999)
                                      IF(DIAGNOSTICS1) THEN
                                        NDOFS_TO_PRINT = SIZE(MESH_VELOCITY_VALUES,1)
                                        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NDOFS_TO_PRINT,NDOFS_TO_PRINT, &
                                          & NDOFS_TO_PRINT,MESH_VELOCITY_VALUES, &
                                          & '(" MESH_VELOCITY_VALUES = ",4(X,E13.6))','4(4(X,E13.6))',err,error,*999)
                                          CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE," ",err,error,*999)
                                        !
                                        NDOFS_TO_PRINT = SIZE(INITIAL_VALUES,1)
                                        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NDOFS_TO_PRINT,NDOFS_TO_PRINT, &
                                          & NDOFS_TO_PRINT,INITIAL_VALUES, &
                                          & '(" INITIAL_VALUES = ",4(X,E13.6))', &
                                          & '4(4(X,E13.6))',err,error,*999)
                                        !
                                        NULLIFY( DUMMY_VALUES1 )
                                        CALL FIELD_PARAMETER_SET_DATA_GET(dependentField,FIELD_VAR_TYPE, &
                                          & FIELD_VALUES_SET_TYPE,DUMMY_VALUES1,err,error,*999)
                                        NDOFS_TO_PRINT = SIZE(DUMMY_VALUES1,1)
                                        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NDOFS_TO_PRINT,NDOFS_TO_PRINT, &
                                          & NDOFS_TO_PRINT,DUMMY_VALUES1, &
                                          & '(" dependentField,FIELD_VAR_TYPE,FIELD_VALUES_SET_TYPE (after) = ",4(X,E13.6))', &
                                          & '4(4(X,E13.6))',err,error,*999)
                                        CALL FIELD_PARAMETER_SET_DATA_RESTORE(dependentField,FIELD_VAR_TYPE, &
                                          & FIELD_VALUES_SET_TYPE,DUMMY_VALUES1,err,error,*999)
                                      ENDIF
                                      CALL FIELD_PARAMETER_SET_DATA_RESTORE(geometricField,FIELD_U_VARIABLE_TYPE, &
                                        & FIELD_MESH_VELOCITY_SET_TYPE,MESH_VELOCITY_VALUES,err,error,*999)
                                      CALL FIELD_PARAMETER_SET_DATA_RESTORE(dependentField,FIELD_VAR_TYPE, &
                                        & FIELD_INITIAL_VALUES_SET_TYPE,INITIAL_VALUES,err,error,*999)
                                    ELSE
                                      CALL FlagError("Boundary condition variable is not associated.",err,error,*999)
                                    END IF

                                    CALL FIELD_PARAMETER_SET_UPDATE_START(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_VAR_TYPE, &
                                      & FIELD_VALUES_SET_TYPE,err,error,*999)
                                    CALL FIELD_PARAMETER_SET_UPDATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_VAR_TYPE, &
                                      & FIELD_VALUES_SET_TYPE,err,error,*999)

                                  ELSE
                                    CALL FlagError("FIELD_VAR_TYPE is not associated.",err,error,*999)
                                  ENDIF
                                ELSE
                                  CALL FlagError("vectorMapping is not associated.",err,error,*999)
                                ENDIF
                              ELSE
                                CALL FlagError("Boundary conditions are not associated.",err,error,*999)
                              END IF
                            ELSE
                              CALL FlagError("Dependent field and/or geometric field is/are not associated.",err,error,*999)
                            END IF
                          CASE DEFAULT
                            localError="Equations set subtype " &
                              & //TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
                              & " is not valid for a Darcy equation fluid type of a fluid mechanics problem class."
                            CALL FlagError(localError,err,error,*999)
                        END SELECT
                      ELSE
                        CALL FlagError("Equations set is not associated.",err,error,*999)
                      END IF
                    ELSE
                      CALL FlagError("Equations are not associated.",err,error,*999)
                    END IF
                  ELSE
                    CALL FlagError("Solver mapping is not associated.",err,error,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Solver equations are not associated.",err,error,*999)
                END IF
              CASE DEFAULT
                localError="Problem subtype "//TRIM(NumberToVString(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",err,error))// &
                  & " is not valid for a Darcy equation fluid type of a fluid mechanics problem class."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          ELSE
            CALL FlagError("Problem is not associated.",err,error,*999)
          ENDIF
        ELSE
          ! do nothing
        ENDIF
      ELSE
        CALL FlagError("Solver is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",err,error,*999)
    ENDIF

    EXITS("Darcy_PreSolveUpdateBoundaryConditions")
    RETURN
999 ERRORS("Darcy_PreSolveUpdateBoundaryConditions",err,error)
    EXITS("Darcy_PreSolveUpdateBoundaryConditions")
    RETURN 1
    
  END SUBROUTINE Darcy_PreSolveUpdateBoundaryConditions

  !
  !================================================================================================================================
  !

  !>Update materials field for ALE Darcy problem
  SUBROUTINE Darcy_PreSolveUpdateMatrixProperties(CONTROL_LOOP,SOLVER,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solvers
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER_MAT_PROPERTIES, SOLVER_ALE_DARCY  !<A pointer to the solvers
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD_MAT_PROPERTIES, MATERIALS_FIELD_ALE_DARCY
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS_MAT_PROPERTIES, SOLVER_EQUATIONS_ALE_DARCY  !<A pointer to the solver equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING_MAT_PROPERTIES, SOLVER_MAPPING_ALE_DARCY !<A pointer to the solver mapping
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equations_SET_MAT_PROPERTIES, EQUATIONS_SET_ALE_DARCY !<A pointer to the equations set
    TYPE(VARYING_STRING) :: localError

    REAL(DP), POINTER :: DUMMY_VALUES2(:)

    INTEGER(INTG) :: NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_MAT_PROPERTIES,NUMBER_OF_COMPONENTS_MATERIALS_FIELD_ALE_DARCY
    INTEGER(INTG) :: NDOFS_TO_PRINT, I


    ENTERS("Darcy_PreSolveUpdateMatrixProperties",err,error,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN

      NULLIFY(SOLVER_MAT_PROPERTIES)
      NULLIFY(SOLVER_ALE_DARCY)

      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a Darcy equation problem.",err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_STANDARD_DARCY_SUBTYPE)
              ! do nothing
            CASE(PROBLEM_QUASISTATIC_DARCY_SUBTYPE)
              ! do nothing
            CASE(PROBLEM_TRANSIENT_DARCY_SUBTYPE)
              ! do nothing
            CASE(PROBLEM_ALE_DARCY_SUBTYPE,PROBLEM_PGM_DARCY_SUBTYPE,PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE, &
              & PROBLEM_PGM_ELASTICITY_DARCY_SUBTYPE,PROBLEM_PGM_TRANSIENT_DARCY_SUBTYPE, &
              & PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)
              IF((CONTROL_LOOP%LOOP_TYPE==PROBLEM_CONTROL_SIMPLE_TYPE.OR.CONTROL_LOOP%LOOP_TYPE==PROBLEM_CONTROL_TIME_LOOP_TYPE) &
                  & .AND.SOLVER%GLOBAL_NUMBER==SOLVER_NUMBER_DARCY) THEN
                !--- Get the dependent field of the Material-Properties Galerkin-Projection equations
                IF(SOLVER%outputType>=SOLVER_PROGRESS_OUTPUT) THEN
                  CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Darcy update materials ... ",err,error,*999)
                ENDIF
                CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,SOLVER_NUMBER_MAT_PROPERTIES,SOLVER_MAT_PROPERTIES,err,error,*999)
                SOLVER_EQUATIONS_MAT_PROPERTIES=>SOLVER_MAT_PROPERTIES%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS_MAT_PROPERTIES)) THEN
                  SOLVER_MAPPING_MAT_PROPERTIES=>SOLVER_EQUATIONS_MAT_PROPERTIES%SOLVER_MAPPING
                  IF(ASSOCIATED(SOLVER_MAPPING_MAT_PROPERTIES)) THEN
                    EQUATIONS_SET_MAT_PROPERTIES=>SOLVER_MAPPING_MAT_PROPERTIES%EQUATIONS_SETS(1)%ptr
                    IF(ASSOCIATED(EQUATIONS_SET_MAT_PROPERTIES)) THEN
                      DEPENDENT_FIELD_MAT_PROPERTIES=>EQUATIONS_SET_MAT_PROPERTIES%DEPENDENT%DEPENDENT_FIELD
                      IF(ASSOCIATED(DEPENDENT_FIELD_MAT_PROPERTIES)) THEN
                        CALL FIELD_NUMBER_OF_COMPONENTS_GET(DEPENDENT_FIELD_MAT_PROPERTIES, &
                          & FIELD_U_VARIABLE_TYPE,NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_MAT_PROPERTIES,err,error,*999)
                      ELSE
                        CALL FlagError("DEPENDENT_FIELD_MAT_PROPERTIES is not associated.",err,error,*999)
                      END IF
                    ELSE
                      CALL FlagError("Galerkin Projection equations set is not associated.",err,error,*999)
                    END IF
                  ELSE
                    CALL FlagError("Galerkin Projection solver mapping is not associated.",err,error,*999)
                  END IF
                ELSE
                  CALL FlagError("Galerkin Projection solver equations are not associated.",err,error,*999)
                END IF

                !--- Get the materials field for the ALE Darcy equations
                CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,SOLVER_NUMBER_DARCY,SOLVER_ALE_DARCY,err,error,*999)
                SOLVER_EQUATIONS_ALE_DARCY=>SOLVER_ALE_DARCY%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS_ALE_DARCY)) THEN
                  SOLVER_MAPPING_ALE_DARCY=>SOLVER_EQUATIONS_ALE_DARCY%SOLVER_MAPPING
                  IF(ASSOCIATED(SOLVER_MAPPING_ALE_DARCY)) THEN
                    EQUATIONS_SET_ALE_DARCY=>SOLVER_MAPPING_ALE_DARCY%EQUATIONS_SETS(1)%ptr
                    IF(ASSOCIATED(EQUATIONS_SET_ALE_DARCY)) THEN
                      MATERIALS_FIELD_ALE_DARCY=>EQUATIONS_SET_ALE_DARCY%MATERIALS%MATERIALS_FIELD
                      IF(ASSOCIATED(MATERIALS_FIELD_ALE_DARCY)) THEN
                        CALL FIELD_NUMBER_OF_COMPONENTS_GET(MATERIALS_FIELD_ALE_DARCY, &
                          & FIELD_U_VARIABLE_TYPE,NUMBER_OF_COMPONENTS_MATERIALS_FIELD_ALE_DARCY,err,error,*999)
                      ELSE
                        CALL FlagError("MATERIALS_FIELD_ALE_DARCY is not associated.",err,error,*999)
                      END IF
                    ELSE
                      CALL FlagError("ALE Darcy equations set is not associated.",err,error,*999)
                    END IF
                  ELSE
                    CALL FlagError("ALE Darcy solver mapping is not associated.",err,error,*999)
                  END IF
                ELSE
                  CALL FlagError("ALE Darcy solver equations are not associated.",err,error,*999)
                END IF

                !--- Copy the result from Galerkin-Projection's dependent field to ALE Darcy's material field
                IF(NUMBER_OF_COMPONENTS_MATERIALS_FIELD_ALE_DARCY==NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_MAT_PROPERTIES) THEN
                  DO I=1,NUMBER_OF_COMPONENTS_MATERIALS_FIELD_ALE_DARCY
                    CALL Field_ParametersToFieldParametersCopy(DEPENDENT_FIELD_MAT_PROPERTIES, &
                      & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,I,MATERIALS_FIELD_ALE_DARCY, &
                      & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,I,err,error,*999)
                  END DO
                ELSE
!                   CALL FlagError("Dimension of Galerkin Projection and ALE Darcy equations set is not consistent",err,error,*999)
                  localError="Number of components of Galerkin-Projection dependent field "// &
                    & "is not consistent with ALE-Darcy-equation material field."
                  CALL FlagError(localError,err,error,*999)
                END IF

                IF(DIAGNOSTICS3) THEN
                  NULLIFY( DUMMY_VALUES2 )
                  CALL FIELD_PARAMETER_SET_DATA_GET(DEPENDENT_FIELD_MAT_PROPERTIES,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,DUMMY_VALUES2,err,error,*999)
                  NDOFS_TO_PRINT = SIZE(DUMMY_VALUES2,1)
                  CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NDOFS_TO_PRINT,NDOFS_TO_PRINT,NDOFS_TO_PRINT,DUMMY_VALUES2, &
                    & '(" DEPENDENT_FIELD_MAT_PROPERTIES,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE = ",4(X,E13.6))',&
                    & '4(4(X,E13.6))',err,error,*999)
                  CALL FIELD_PARAMETER_SET_DATA_RESTORE(DEPENDENT_FIELD_MAT_PROPERTIES,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,DUMMY_VALUES2,err,error,*999)
                ENDIF

              ELSE
                ! do nothing
              END IF
            CASE DEFAULT
              localError="The third problem specification of "// &
                & TRIM(NumberToVString(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",err,error))// &
                & " is not valid for Darcy_PreSolveUpdateMatrixProperties."
            CALL FLAG_ERROR(localError,err,error,*999)
          END SELECT
        ELSE
          CALL FlagError("Problem is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Solver is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",err,error,*999)
    ENDIF

    EXITS("Darcy_PreSolveUpdateMatrixProperties")
    RETURN
999 ERRORS("Darcy_PreSolveUpdateMatrixProperties",err,error)
    EXITS("Darcy_PreSolveUpdateMatrixProperties")
    RETURN 1
    
  END SUBROUTINE Darcy_PreSolveUpdateMatrixProperties

  !
  !================================================================================================================================
  !

  !>Sets up the Darcy problem post solve.
  SUBROUTINE DARCY_EQUATION_POST_SOLVE(CONTROL_LOOP,SOLVER,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER!<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("DARCY_EQUATION_POST_SOLVE",err,error,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a Darcy equation problem.",err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_STANDARD_DARCY_SUBTYPE)
              CALL DARCY_EQUATION_POST_SOLVE_OUTPUT_DATA(CONTROL_LOOP,SOLVER,err,error,*999)
            CASE(PROBLEM_QUASISTATIC_DARCY_SUBTYPE)
              CALL DARCY_EQUATION_POST_SOLVE_OUTPUT_DATA(CONTROL_LOOP,SOLVER,err,error,*999)
            CASE(PROBLEM_ALE_DARCY_SUBTYPE,PROBLEM_PGM_ELASTICITY_DARCY_SUBTYPE)
              IF(SOLVER%GLOBAL_NUMBER==SOLVER_NUMBER_DARCY) THEN
                CALL DARCY_EQUATION_POST_SOLVE_OUTPUT_DATA(CONTROL_LOOP,SOLVER,err,error,*999)
              END IF
            CASE(PROBLEM_PGM_DARCY_SUBTYPE,PROBLEM_PGM_TRANSIENT_DARCY_SUBTYPE,PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE, &
              & PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE,PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)
              IF(SOLVER%GLOBAL_NUMBER==SOLVER_NUMBER_DARCY) THEN
                CALL DARCY_EQUATION_POST_SOLVE_OUTPUT_DATA(CONTROL_LOOP,SOLVER,err,error,*999)

              ! The following command only when setting the Darcy mass increase explicitly to test finite elasticity !!!
! ! !               CALL DARCY_EQUATION_POST_SOLVE_SET_MASS_INCREASE(CONTROL_LOOP,SOLVER,err,error,*999)

              END IF
            CASE(PROBLEM_TRANSIENT_DARCY_SUBTYPE)
              CALL DARCY_EQUATION_POST_SOLVE_OUTPUT_DATA(CONTROL_LOOP,SOLVER,err,error,*999)
            CASE DEFAULT
              localError="Problem subtype "//TRIM(NumberToVString(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",err,error))// &
                & " is not valid for a Darcy fluid type of a fluid mechanics problem class."
              CALL FlagError(localError,err,error,*999)
          END SELECT
        ELSE
          CALL FlagError("Problem is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Solver is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",err,error,*999)
    ENDIF

    EXITS("DARCY_EQUATION_POST_SOLVE")
    RETURN
999 ERRORSEXITS("DARCY_EQUATION_POST_SOLVE",err,error)
    RETURN 1
  END SUBROUTINE DARCY_EQUATION_POST_SOLVE

  !
  !================================================================================================================================
  !

  !>Sets up the Darcy problem post solve output data.
  SUBROUTINE DARCY_EQUATION_POST_SOLVE_OUTPUT_DATA(CONTROL_LOOP,SOLVER,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS  !<A pointer to the solver equations
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equations_SET !<A pointer to the equations set
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_TIME_LOOP !<A pointer to the control time loop.
    TYPE(VARYING_STRING) :: localError
    TYPE(VARYING_STRING) :: METHOD !,FILE
    CHARACTER(14) :: FILE
    CHARACTER(14) :: OUTPUT_FILE
    LOGICAL :: EXPORT_FIELD
    INTEGER(INTG) :: CURRENT_LOOP_ITERATION,SUBITERATION_NUMBER
    INTEGER(INTG) :: OUTPUT_ITERATION_NUMBER
    INTEGER(INTG) :: equations_set_idx,loop_idx

    ENTERS("DARCY_EQUATION_POST_SOLVE_OUTPUT_DATA",err,error,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a Darcy equation problem.",err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_STANDARD_DARCY_SUBTYPE)
              SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                  SOLVER_MAPPING=>SOLVER_equations%SOLVER_MAPPING
                  IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                    !Make sure the equations sets are up to date
                    DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                      EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%ptr
                      METHOD="FORTRAN"
                      EXPORT_FIELD=.TRUE.
                      IF(EXPORT_FIELD) THEN
                        IF(SOLVER%outputType>=SOLVER_PROGRESS_OUTPUT) THEN
                          CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Darcy export fields ... ",err,error,*999)
                          CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"STATICSOLUTION",err,error,*999)
                        ENDIF
                        CALL FLUID_MECHANICS_IO_WRITE_CMGUI(EQUATIONS_SET%REGION,EQUATIONS_SET%GLOBAL_NUMBER,"STATICSOLUTION", &
                          & err,error,*999)
                      ENDIF
                    ENDDO
                  ENDIF
                ENDIF
            CASE(PROBLEM_QUASISTATIC_DARCY_SUBTYPE, PROBLEM_ALE_DARCY_SUBTYPE, PROBLEM_PGM_DARCY_SUBTYPE, &
              & PROBLEM_TRANSIENT_DARCY_SUBTYPE, PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE, &
              & PROBLEM_PGM_ELASTICITY_DARCY_SUBTYPE,PROBLEM_PGM_TRANSIENT_DARCY_SUBTYPE, &
              & PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE,PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)
              SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
              IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                SOLVER_MAPPING=>SOLVER_equations%SOLVER_MAPPING
                IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                  !Make sure the equations sets are up to date
                  DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                  EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%ptr
                   IF(EQUATIONS_SET%SPECIFICATION(2)==EQUATIONS_SET_DARCY_EQUATION_TYPE)THEN
                    IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)THEN
                      CONTROL_TIME_LOOP=>CONTROL_LOOP
                      DO loop_idx=1,CONTROL_LOOP%CONTROL_LOOP_LEVEL
                        IF(CONTROL_TIME_LOOP%LOOP_TYPE==PROBLEM_CONTROL_TIME_LOOP_TYPE) THEN
                          CURRENT_LOOP_ITERATION=CONTROL_TIME_LOOP%TIME_LOOP%ITERATION_NUMBER
                          OUTPUT_ITERATION_NUMBER=CONTROL_TIME_LOOP%TIME_LOOP%OUTPUT_NUMBER
                          EXIT
                        ENDIF
                        IF (ASSOCIATED(CONTROL_LOOP%PARENT_LOOP)) THEN
                          CONTROL_TIME_LOOP=>CONTROL_TIME_LOOP%PARENT_LOOP
                        ELSE
                          CURRENT_LOOP_ITERATION=0
                          OUTPUT_ITERATION_NUMBER=0
                        ENDIF
                      ENDDO
                      IF(CONTROL_LOOP%PARENT_LOOP%LOOP_TYPE==PROBLEM_CONTROL_WHILE_LOOP_TYPE) THEN
                        SUBITERATION_NUMBER=CONTROL_LOOP%PARENT_LOOP%WHILE_LOOP%ITERATION_NUMBER
                      ENDIF

                       IF(OUTPUT_ITERATION_NUMBER/=0) THEN
                        IF(CONTROL_TIME_LOOP%TIME_LOOP%CURRENT_TIME<=CONTROL_TIME_LOOP%TIME_LOOP%STOP_TIME) THEN
                          IF(CURRENT_LOOP_ITERATION<10) THEN
                            WRITE(OUTPUT_FILE, '("T_STEP_000",I0,"_C",I0)') CURRENT_LOOP_ITERATION, equations_set_idx
                          ELSE IF(CURRENT_LOOP_ITERATION<100) THEN
                            WRITE(OUTPUT_FILE,'("T_STEP_00",I0,"_C",I0)') CURRENT_LOOP_ITERATION, equations_set_idx
                          ELSE IF(CURRENT_LOOP_ITERATION<1000) THEN
                            WRITE(OUTPUT_FILE,'("T_STEP_0",I0,"_C",I0)') CURRENT_LOOP_ITERATION, equations_set_idx
                          ELSE IF(CURRENT_LOOP_ITERATION<10000) THEN
                            WRITE(OUTPUT_FILE,'("T_STEP_",I0,"_C",I0)') CURRENT_LOOP_ITERATION, equations_set_idx
                          END IF
                          FILE=OUTPUT_FILE
                          METHOD="FORTRAN"
                          EXPORT_FIELD=.TRUE.
                          IF(EXPORT_FIELD) THEN
                            IF(MOD(CURRENT_LOOP_ITERATION,OUTPUT_ITERATION_NUMBER)==0)  THEN
                              IF(SOLVER%outputType>=SOLVER_PROGRESS_OUTPUT) THEN
                                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Darcy export fields ...",err,error,*999)
                                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,OUTPUT_FILE,err,error,*999)
                              ENDIF
                              CALL FLUID_MECHANICS_IO_WRITE_CMGUI(EQUATIONS_SET%REGION,EQUATIONS_SET%GLOBAL_NUMBER,FILE, &
                                & err,error,*999)
                              IF(SOLVER%outputType>=SOLVER_PROGRESS_OUTPUT) THEN
                                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Darcy all fields exported ...",err,error,*999)
                              ENDIF
                            ENDIF
                          ENDIF
                        ENDIF
                       ENDIF


                      !Subiteration intermediate solutions / iterates output:
!                        IF(CONTROL_LOOP%PARENT_LOOP%LOOP_TYPE==PROBLEM_CONTROL_WHILE_LOOP_TYPE) THEN  !subiteration exists
!                         IF(CURRENT_LOOP_ITERATION<10) THEN
!                           IF(SUBITERATION_NUMBER<10) THEN
!                             WRITE(OUTPUT_FILE,'("T_00",I0,"_SB_0",I0,"_C",I0)') CURRENT_LOOP_ITERATION,SUBITERATION_NUMBER, &
!                               & equations_set_idx
!                           ELSE IF(SUBITERATION_NUMBER<100) THEN
!                             WRITE(OUTPUT_FILE,'("T_00",I0,"_SB_",I0,"_C",I0)') CURRENT_LOOP_ITERATION,SUBITERATION_NUMBER, &
!                               & equations_set_idx
!                           END IF
!                           FILE=OUTPUT_FILE
!                           METHOD="FORTRAN"
!                           EXPORT_FIELD=.TRUE.
!                           IF(EXPORT_FIELD) THEN
!                             CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Darcy export subiterates ...",err,error,*999)
!                             CALL FLUID_MECHANICS_IO_WRITE_CMGUI(EQUATIONS_SET%REGION,EQUATIONS_SET%GLOBAL_NUMBER,FILE, &
!                               & err,error,*999)
!                             CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,OUTPUT_FILE,err,error,*999)
!                           ENDIF
!                         ENDIF
!                        ENDIF

                    ELSE !for single compartment (i.e. standary Darcy flow) equations sets
                      !Find the time loop
                      CONTROL_TIME_LOOP=>CONTROL_LOOP
                      DO loop_idx=1,CONTROL_LOOP%CONTROL_LOOP_LEVEL
                        IF(CONTROL_TIME_LOOP%LOOP_TYPE==PROBLEM_CONTROL_TIME_LOOP_TYPE) THEN
                          CURRENT_LOOP_ITERATION=CONTROL_TIME_LOOP%TIME_LOOP%ITERATION_NUMBER
                          OUTPUT_ITERATION_NUMBER=CONTROL_TIME_LOOP%TIME_LOOP%OUTPUT_NUMBER
                          EXIT
                        ENDIF
                        IF (ASSOCIATED(CONTROL_LOOP%PARENT_LOOP)) THEN
                          CONTROL_TIME_LOOP=>CONTROL_TIME_LOOP%PARENT_LOOP
                        ELSE
                          CURRENT_LOOP_ITERATION=0
                          OUTPUT_ITERATION_NUMBER=0
                        ENDIF
                      ENDDO
                      !If coupled with finite elasticity and using subiterations, get the while loop iteration number
                      IF (ASSOCIATED(CONTROL_LOOP%PARENT_LOOP)) THEN
                        IF(CONTROL_LOOP%PARENT_LOOP%LOOP_TYPE==PROBLEM_CONTROL_WHILE_LOOP_TYPE) THEN
                          SUBITERATION_NUMBER=CONTROL_LOOP%PARENT_LOOP%WHILE_LOOP%ITERATION_NUMBER
                        ELSE
                          SUBITERATION_NUMBER=0
                        ENDIF
                      ENDIF

                       IF(OUTPUT_ITERATION_NUMBER/=0) THEN
                        IF(CONTROL_TIME_LOOP%TIME_LOOP%CURRENT_TIME<=CONTROL_TIME_LOOP%TIME_LOOP%STOP_TIME) THEN
                          IF(CURRENT_LOOP_ITERATION<10) THEN
                            WRITE(OUTPUT_FILE,'("TIME_STEP_000",I0)') CURRENT_LOOP_ITERATION
                          ELSE IF(CURRENT_LOOP_ITERATION<100) THEN
                            WRITE(OUTPUT_FILE,'("TIME_STEP_00",I0)') CURRENT_LOOP_ITERATION
                          ELSE IF(CURRENT_LOOP_ITERATION<1000) THEN
                            WRITE(OUTPUT_FILE,'("TIME_STEP_0",I0)') CURRENT_LOOP_ITERATION
                          ELSE IF(CURRENT_LOOP_ITERATION<10000) THEN
                            WRITE(OUTPUT_FILE,'("TIME_STEP_",I0)') CURRENT_LOOP_ITERATION
                          END IF
                          FILE=OUTPUT_FILE
                          METHOD="FORTRAN"
                          EXPORT_FIELD=.TRUE.
                          IF(EXPORT_FIELD) THEN
                            IF(MOD(CURRENT_LOOP_ITERATION,OUTPUT_ITERATION_NUMBER)==0)  THEN
                              IF(SOLVER%outputType>=SOLVER_PROGRESS_OUTPUT) THEN
                                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Darcy export fields ...",err,error,*999)
                                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,OUTPUT_FILE,err,error,*999)
                              ENDIF
                              CALL FLUID_MECHANICS_IO_WRITE_CMGUI(EQUATIONS_SET%REGION,EQUATIONS_SET%GLOBAL_NUMBER,FILE, &
                                & err,error,*999)
                              IF(SOLVER%outputType>=SOLVER_PROGRESS_OUTPUT) THEN
                                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Darcy all fields exported ...",err,error,*999)
                              ENDIF
                            ENDIF
                          ENDIF
                        ENDIF
                       ENDIF


!                       !Subiteration intermediate solutions / iterates output:
!                        IF(CONTROL_LOOP%PARENT_LOOP%LOOP_TYPE==PROBLEM_CONTROL_WHILE_LOOP_TYPE) THEN  !subiteration exists
!                         IF(CURRENT_LOOP_ITERATION<10) THEN
!                           IF(SUBITERATION_NUMBER<10) THEN
!                             WRITE(OUTPUT_FILE,'("T_00",I0,"_SUB_000",I0)') CURRENT_LOOP_ITERATION,SUBITERATION_NUMBER
!                           ELSE IF(SUBITERATION_NUMBER<100) THEN
!                             WRITE(OUTPUT_FILE,'("T_00",I0,"_SUB_00",I0)') CURRENT_LOOP_ITERATION,SUBITERATION_NUMBER
!                           END IF
!                           FILE=OUTPUT_FILE
!                           METHOD="FORTRAN"
!                           EXPORT_FIELD=.TRUE.
!                           IF(EXPORT_FIELD) THEN
!                             CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Darcy export subiterates ...",err,error,*999)
!                             CALL FLUID_MECHANICS_IO_WRITE_CMGUI(EQUATIONS_SET%REGION,EQUATIONS_SET%GLOBAL_NUMBER,FILE, &
!                               & err,error,*999)
!                             CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,OUTPUT_FILE,err,error,*999)
!                           ENDIF
!                         ENDIF
!                        ENDIF

                    ENDIF
                   ENDIF
                  ENDDO
                ENDIF
              ENDIF
            CASE DEFAULT
              localError="Problem subtype "//TRIM(NumberToVString(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",err,error))// &
                & " is not valid for a Darcy fluid type of a fluid mechanics problem class."
              CALL FlagError(localError,err,error,*999)
          END SELECT
        ELSE
          CALL FlagError("Problem is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Solver is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",err,error,*999)
    ENDIF

    EXITS("DARCY_EQUATION_POST_SOLVE_OUTPUT_DATA")
    RETURN
999 ERRORSEXITS("DARCY_EQUATION_POST_SOLVE_OUTPUT_DATA",err,error)
    RETURN 1
  END SUBROUTINE DARCY_EQUATION_POST_SOLVE_OUTPUT_DATA

  !
  !================================================================================================================================
  !

  !>Calculates the analytic solution and sets the boundary conditions for an analytic problem.
  SUBROUTINE Darcy_BoundaryConditionsAnalyticCalculate(EQUATIONS_SET,BOUNDARY_CONDITIONS,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equations_SET
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,deriv_idx,dim_idx,local_ny,node_idx,NUMBER_OF_DIMENSIONS,variable_idx,variable_type,I,J,K
    INTEGER(INTG) :: number_of_nodes_xic(3),element_idx,en_idx,BOUND_COUNT
    REAL(DP) :: VALUE,X(3),ARG(3),L,XI_COORDINATES(3),FACT,PERM_OVER_VIS_PARAM
    REAL(DP) :: BOUNDARY_TOLERANCE, BOUNDARY_X(3,2), T_COORDINATES(20,3)
    REAL(DP), POINTER :: GEOMETRIC_PARAMETERS(:)
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,geometricField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE,GEOMETRIC_VARIABLE
    TYPE(FIELD_INTERPOLATED_POINT_PTR_TYPE), POINTER :: INTERPOLATED_POINT (:)
    TYPE(FIELD_INTERPOLATION_PARAMETERS_PTR_TYPE), POINTER :: INTERPOLATION_PARAMETERS(:)
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: GLOBAL_DERIV_INDEX,ANALYTIC_FUNCTION_TYPE
    REAL(DP) :: CURRENT_TIME
    !Temp variables
    INTEGER(INTG) :: number_of_element_nodes,temp_local_ny,temp_node_number,velocity_DOF_check,temp_local_node_number

    ENTERS("Darcy_BoundaryConditionsAnalyticCalculate",err,error,*999)

    BOUND_COUNT=0

    PERM_OVER_VIS_PARAM = 1.0_DP  !temporarily hard-coded: Should rather be determined by interpolating materials field

    L=10.0_DP
    XI_COORDINATES(3)=0.0_DP
    BOUNDARY_TOLERANCE=0.000000001_DP
    BOUNDARY_X=0.0_DP
    T_COORDINATES=0.0_DP

    number_of_element_nodes=0
    temp_local_node_number=0
    temp_local_ny=0
    temp_node_number=0
    velocity_DOF_check=0


    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
      SELECT CASE(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE)
      CASE(EQUATIONS_SET_INCOMP_ELAST_DARCY_ANALYTIC_DARCY)
        DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
        IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
          geometricField=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
          IF(ASSOCIATED(geometricField)) THEN
            CALL FIELD_NUMBER_OF_COMPONENTS_GET(geometricField,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS,err,error,*999)
            NULLIFY(GEOMETRIC_VARIABLE)
            NULLIFY(GEOMETRIC_PARAMETERS)
            CALL Field_VariableGet(geometricField,FIELD_U_VARIABLE_TYPE,GEOMETRIC_VARIABLE,err,error,*999)
            CALL FIELD_PARAMETER_SET_DATA_GET(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,GEOMETRIC_PARAMETERS, &
              & err,error,*999)
            IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
              CURRENT_TIME=EQUATIONS_SET%ANALYTIC%ANALYTIC_USER_PARAMS(1)
              DO variable_idx=3,DEPENDENT_FIELD%NUMBER_OF_VARIABLES
                variable_type=DEPENDENT_FIELD%VARIABLES(variable_idx)%VARIABLE_TYPE
                FIELD_VARIABLE=>DEPENDENT_FIELD%VARIABLE_TYPE_MAP(variable_type)%ptr
                IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                  CALL FIELD_PARAMETER_SET_CREATE(DEPENDENT_FIELD,variable_type,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
                  DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                    IF(FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION) THEN
                      DOMAIN=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN
                      IF(ASSOCIATED(DOMAIN)) THEN
                        IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
                          DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
                          IF(ASSOCIATED(DOMAIN_NODES)) THEN
                            !Loop over the local nodes excluding the ghosts.
                            DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
                              !!TODO \todo We should interpolate the geometric field here and the node position.
                              DO dim_idx=1,NUMBER_OF_DIMENSIONS
                                !Default to version 1 of each node derivative
                                local_ny=GEOMETRIC_VARIABLE%COMPONENTS(dim_idx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                                  & NODES(node_idx)%DERIVATIVES(1)%VERSIONS(1)
                                X(dim_idx)=GEOMETRIC_PARAMETERS(local_ny)
                              ENDDO !dim_idx
                              !Loop over the derivatives
                              DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
                                ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE
                                GLOBAL_DERIV_INDEX=DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX
!                                 CALL DIFFUSION_EQUATION_ANALYTIC_FUNCTIONS(VALUE,X,CURRENT_TIME,variable_type, &
!                                   & GLOBAL_DERIV_INDEX,ANALYTIC_FUNCTION_TYPE,err,error,*999)
!!!!!!!!!!!!NEED TO SET APPROPRIATE VALUE DEPENDING ON WHETHER IT IS A VELOCITY COMPONENT OR THE MASS INCREASE COMPONENT
                                VALUE=0.0_DP
                                !Default to version 1 of each node derivative
                                local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                  & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(deriv_idx)%VERSIONS(1)
                                CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD,variable_type, &
                                  & FIELD_ANALYTIC_VALUES_SET_TYPE,local_ny,VALUE,err,error,*999)
                                IF(variable_type==FIELD_V_VARIABLE_TYPE) THEN
                                  IF(DOMAIN_NODES%NODES(node_idx)%BOUNDARY_NODE) THEN
                                    !If we are a boundary node then set the analytic value on the boundary
                                    CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(BOUNDARY_CONDITIONS,DEPENDENT_FIELD,variable_type, &
                                      & local_ny,BOUNDARY_CONDITION_FIXED,VALUE,err,error,*999)
                                  ELSE
                                    CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD,variable_type, &
                                    & FIELD_VALUES_SET_TYPE,local_ny,VALUE,err,error,*999)
                                  ENDIF
                                ENDIF
                              ENDDO !deriv_idx
                            ENDDO !node_idx
                          ELSE
                            CALL FlagError("Domain topology nodes is not associated.",err,error,*999)
                          ENDIF
                        ELSE
                          CALL FlagError("Domain topology is not associated.",err,error,*999)
                        ENDIF
                      ELSE
                        CALL FlagError("Domain is not associated.",err,error,*999)
                      ENDIF
                    ELSE
                      CALL FlagError("Only node based interpolation is implemented.",err,error,*999)
                    ENDIF
                  ENDDO !component_idx
                  CALL FIELD_PARAMETER_SET_UPDATE_START(DEPENDENT_FIELD,variable_type,FIELD_ANALYTIC_VALUES_SET_TYPE, &
                    & err,error,*999)
                  CALL FIELD_PARAMETER_SET_UPDATE_FINISH(DEPENDENT_FIELD,variable_type,FIELD_ANALYTIC_VALUES_SET_TYPE, &
                    & err,error,*999)
                ELSE
                  CALL FlagError("Field variable is not associated.",err,error,*999)
                ENDIF
              ENDDO !variable_idx
              CALL FIELD_PARAMETER_SET_DATA_RESTORE(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
               & GEOMETRIC_PARAMETERS,err,error,*999)
            ELSE
              CALL FlagError("Equations set boundary conditions is not associated.",err,error,*999)
            ENDIF
          ELSE
            CALL FlagError("Equations set geometric field is not associated.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Equations set dependent field is not associated.",err,error,*999)
        ENDIF
      CASE DEFAULT
        DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
        IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
          geometricField=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
          IF(ASSOCIATED(geometricField)) THEN
            NULLIFY(INTERPOLATION_PARAMETERS)
            NULLIFY(INTERPOLATED_POINT)
            CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(geometricField,INTERPOLATION_PARAMETERS,err,error,*999)
            CALL FIELD_INTERPOLATED_POINTS_INITIALISE(INTERPOLATION_PARAMETERS,INTERPOLATED_POINT,err,error,*999)

            CALL FIELD_NUMBER_OF_COMPONENTS_GET(geometricField,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS,err,error,*999)

            IF(NUMBER_OF_DIMENSIONS==2) THEN
              BOUNDARY_X(1,1)=0.0_DP
              BOUNDARY_X(1,2)=10.0_DP
              BOUNDARY_X(2,1)=0.0_DP
              BOUNDARY_X(2,2)=10.0_DP
            ELSE IF(NUMBER_OF_DIMENSIONS==3) THEN
              BOUNDARY_X(1,1)=-5.0_DP
              BOUNDARY_X(1,2)=5.0_DP
              BOUNDARY_X(2,1)=-5.0_DP
              BOUNDARY_X(2,2)=5.0_DP
              BOUNDARY_X(3,1)=-5.0_DP
              BOUNDARY_X(3,2)=5.0_DP
            ENDIF

            NULLIFY(GEOMETRIC_VARIABLE)
            CALL Field_VariableGet(geometricField,FIELD_U_VARIABLE_TYPE,GEOMETRIC_VARIABLE,err,error,*999)
            CALL FIELD_PARAMETER_SET_DATA_GET(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,GEOMETRIC_PARAMETERS, &
              & err,error,*999)
            IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
              DO variable_idx=1,DEPENDENT_FIELD%NUMBER_OF_VARIABLES
                variable_type=DEPENDENT_FIELD%VARIABLES(variable_idx)%VARIABLE_TYPE
                FIELD_VARIABLE=>DEPENDENT_FIELD%VARIABLE_TYPE_MAP(variable_type)%ptr
                IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                  CALL FIELD_PARAMETER_SET_CREATE(DEPENDENT_FIELD,variable_type,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
                  DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                    BOUND_COUNT=0
                    IF(FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION) THEN
                      DOMAIN=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN
                      IF(ASSOCIATED(DOMAIN)) THEN
                        IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
                          DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
                          IF(ASSOCIATED(DOMAIN_NODES)) THEN
                            !Loop over the local nodes excluding the ghosts.
                            DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES

                              element_idx=DOMAIN%topology%nodes%nodes(node_idx)%surrounding_elements(1)
                              CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,element_idx, &
                                & INTERPOLATION_PARAMETERS(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)

!                             DO I=1,DOMAIN%topology%elements%maximum_number_of_element_parameters
!                               IF(DOMAIN%topology%elements%elements(element_idx)%element_nodes(I)=node_idx THEN

                              en_idx=0
                              XI_COORDINATES=0.0_DP
                              number_of_nodes_xic(1)=DOMAIN%topology%elements%elements(element_idx)%basis%number_of_nodes_xic(1)
                              number_of_nodes_xic(2)=DOMAIN%topology%elements%elements(element_idx)%basis%number_of_nodes_xic(2)
                              IF(NUMBER_OF_DIMENSIONS==3) THEN
                                number_of_nodes_xic(3)=DOMAIN%topology%elements%elements(element_idx)%basis%number_of_nodes_xic(3)
                              ELSE
                                number_of_nodes_xic(3)=1
                              ENDIF

                              IF(DOMAIN%topology%elements%maximum_number_of_element_parameters==4.AND.NUMBER_OF_DIMENSIONS==2 .OR. &
                                & DOMAIN%topology%elements%maximum_number_of_element_parameters==9.OR. &
                                & DOMAIN%topology%elements%maximum_number_of_element_parameters==16.OR. &
                                & DOMAIN%topology%elements%maximum_number_of_element_parameters==8.OR. &
                                & DOMAIN%topology%elements%maximum_number_of_element_parameters==27.OR. &
                                & DOMAIN%topology%elements%maximum_number_of_element_parameters==64) THEN

                                DO K=1,number_of_nodes_xic(3)
                                  DO J=1,number_of_nodes_xic(2)
                                    DO I=1,number_of_nodes_xic(1)
                                      en_idx=en_idx+1
                                      IF(DOMAIN%topology%elements%elements(element_idx)%element_nodes(en_idx)==node_idx) EXIT
                                      XI_COORDINATES(1)=XI_COORDINATES(1)+(1.0_DP/(number_of_nodes_xic(1)-1))
                                    ENDDO
                                      IF(DOMAIN%topology%elements%elements(element_idx)%element_nodes(en_idx)==node_idx) EXIT
                                      XI_COORDINATES(1)=0.0_DP
                                      XI_COORDINATES(2)=XI_COORDINATES(2)+(1.0_DP/(number_of_nodes_xic(2)-1))
                                  ENDDO
                                  IF(DOMAIN%topology%elements%elements(element_idx)%element_nodes(en_idx)==node_idx) EXIT
                                  XI_COORDINATES(1)=0.0_DP
                                  XI_COORDINATES(2)=0.0_DP
                                  IF(number_of_nodes_xic(3)/=1) THEN
                                    XI_COORDINATES(3)=XI_COORDINATES(3)+(1.0_DP/(number_of_nodes_xic(3)-1))
                                  ENDIF
                                ENDDO
                                CALL FIELD_INTERPOLATE_XI(NO_PART_DERIV,XI_COORDINATES, &
                                  & INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                              ELSE
                                IF(DOMAIN%topology%elements%maximum_number_of_element_parameters==3) THEN
                                  T_COORDINATES(1,1:2)=[0.0_DP,1.0_DP]
                                  T_COORDINATES(2,1:2)=[1.0_DP,0.0_DP]
                                  T_COORDINATES(3,1:2)=[1.0_DP,1.0_DP]
                                ELSE IF(DOMAIN%topology%elements%maximum_number_of_element_parameters==6) THEN
                                  T_COORDINATES(1,1:2)=[0.0_DP,1.0_DP]
                                  T_COORDINATES(2,1:2)=[1.0_DP,0.0_DP]
                                  T_COORDINATES(3,1:2)=[1.0_DP,1.0_DP]
                                  T_COORDINATES(4,1:2)=[0.5_DP,0.5_DP]
                                  T_COORDINATES(5,1:2)=[1.0_DP,0.5_DP]
                                  T_COORDINATES(6,1:2)=[0.5_DP,1.0_DP]
                                ELSE IF(DOMAIN%topology%elements%maximum_number_of_element_parameters==10.AND. &
                                  & NUMBER_OF_DIMENSIONS==2) THEN
                                  T_COORDINATES(1,1:2)=[0.0_DP,1.0_DP]
                                  T_COORDINATES(2,1:2)=[1.0_DP,0.0_DP]
                                  T_COORDINATES(3,1:2)=[1.0_DP,1.0_DP]
                                  T_COORDINATES(4,1:2)=[1.0_DP/3.0_DP,2.0_DP/3.0_DP]
                                  T_COORDINATES(5,1:2)=[2.0_DP/3.0_DP,1.0_DP/3.0_DP]
                                  T_COORDINATES(6,1:2)=[1.0_DP,1.0_DP/3.0_DP]
                                  T_COORDINATES(7,1:2)=[1.0_DP,2.0_DP/3.0_DP]
                                  T_COORDINATES(8,1:2)=[2.0_DP/3.0_DP,1.0_DP]
                                  T_COORDINATES(9,1:2)=[1.0_DP/3.0_DP,1.0_DP]
                                  T_COORDINATES(10,1:2)=[2.0_DP/3.0_DP,2.0_DP/3.0_DP]
                                ELSE IF(DOMAIN%topology%elements%maximum_number_of_element_parameters==4) THEN
                                  T_COORDINATES(1,1:3)=[0.0_DP,1.0_DP,1.0_DP]
                                  T_COORDINATES(2,1:3)=[1.0_DP,0.0_DP,1.0_DP]
                                  T_COORDINATES(3,1:3)=[1.0_DP,1.0_DP,0.0_DP]
                                  T_COORDINATES(4,1:3)=[1.0_DP,1.0_DP,1.0_DP]
                                ELSE IF(DOMAIN%topology%elements%maximum_number_of_element_parameters==10.AND. &
                                  & NUMBER_OF_DIMENSIONS==3) THEN
                                  T_COORDINATES(1,1:3)=[0.0_DP,1.0_DP,1.0_DP]
                                  T_COORDINATES(2,1:3)=[1.0_DP,0.0_DP,1.0_DP]
                                  T_COORDINATES(3,1:3)=[1.0_DP,1.0_DP,0.0_DP]
                                  T_COORDINATES(4,1:3)=[1.0_DP,1.0_DP,1.0_DP]
                                  T_COORDINATES(5,1:3)=[0.5_DP,0.5_DP,1.0_DP]
                                  T_COORDINATES(6,1:3)=[0.5_DP,1.0_DP,0.5_DP]
                                  T_COORDINATES(7,1:3)=[0.5_DP,1.0_DP,1.0_DP]
                                  T_COORDINATES(8,1:3)=[1.0_DP,0.5_DP,0.5_DP]
                                  T_COORDINATES(9,1:3)=[1.0_DP,1.0_DP,0.5_DP]
                                  T_COORDINATES(10,1:3)=[1.0_DP,0.5_DP,1.0_DP]
                                ELSE IF(DOMAIN%topology%elements%maximum_number_of_element_parameters==20) THEN
                                  T_COORDINATES(1,1:3)=[0.0_DP,1.0_DP,1.0_DP]
                                  T_COORDINATES(2,1:3)=[1.0_DP,0.0_DP,1.0_DP]
                                  T_COORDINATES(3,1:3)=[1.0_DP,1.0_DP,0.0_DP]
                                  T_COORDINATES(4,1:3)=[1.0_DP,1.0_DP,1.0_DP]
                                  T_COORDINATES(5,1:3)=[1.0_DP/3.0_DP,2.0_DP/3.0_DP,1.0_DP]
                                  T_COORDINATES(6,1:3)=[2.0_DP/3.0_DP,1.0_DP/3.0_DP,1.0_DP]
                                  T_COORDINATES(7,1:3)=[1.0_DP/3.0_DP,1.0_DP,2.0_DP/3.0_DP]
                                  T_COORDINATES(8,1:3)=[2.0_DP/3.0_DP,1.0_DP,1.0_DP/3.0_DP]
                                  T_COORDINATES(9,1:3)=[1.0_DP/3.0_DP,1.0_DP,1.0_DP]
                                  T_COORDINATES(10,1:3)=[2.0_DP/3.0_DP,1.0_DP,1.0_DP]
                                  T_COORDINATES(11,1:3)=[1.0_DP,1.0_DP/3.0_DP,2.0_DP/3.0_DP]
                                  T_COORDINATES(12,1:3)=[1.0_DP,2.0_DP/3.0_DP,1.0_DP/3.0_DP]
                                  T_COORDINATES(13,1:3)=[1.0_DP,1.0_DP,1.0_DP/3.0_DP]
                                  T_COORDINATES(14,1:3)=[1.0_DP,1.0_DP,2.0_DP/3.0_DP]
                                  T_COORDINATES(15,1:3)=[1.0_DP,1.0_DP/3.0_DP,1.0_DP]
                                  T_COORDINATES(16,1:3)=[1.0_DP,2.0_DP/3.0_DP,1.0_DP]
                                  T_COORDINATES(17,1:3)=[2.0_DP/3.0_DP,2.0_DP/3.0_DP,2.0_DP/3.0_DP]
                                  T_COORDINATES(18,1:3)=[2.0_DP/3.0_DP,2.0_DP/3.0_DP,1.0_DP]
                                  T_COORDINATES(19,1:3)=[2.0_DP/3.0_DP,1.0_DP,2.0_DP/3.0_DP]
                                  T_COORDINATES(20,1:3)=[1.0_DP,2.0_DP/3.0_DP,2.0_DP/3.0_DP]
                                ENDIF

                                DO K=1,DOMAIN%topology%elements%maximum_number_of_element_parameters
                                  IF(DOMAIN%topology%elements%elements(element_idx)%element_nodes(K)==node_idx) EXIT
                                ENDDO

                                IF(NUMBER_OF_DIMENSIONS==2) THEN
                                  CALL FIELD_INTERPOLATE_XI(NO_PART_DERIV,T_COORDINATES(K,1:2), &
                                    & INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                                ELSE IF(NUMBER_OF_DIMENSIONS==3) THEN
                                  CALL FIELD_INTERPOLATE_XI(NO_PART_DERIV,T_COORDINATES(K,1:3), &
                                    & INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                                ENDIF
                              ENDIF

                              X=0.0_DP
                              DO dim_idx=1,NUMBER_OF_DIMENSIONS
                                X(dim_idx)=INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(dim_idx,1)
                              ENDDO !dim_idx

                              !Loop over the derivatives
                              DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
                                SELECT CASE(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE)
                                CASE(EQUATIONS_SET_DARCY_EQUATION_TWO_DIM_1)
                                  IF(NUMBER_OF_DIMENSIONS==2.AND.FIELD_VARIABLE%NUMBER_OF_COMPONENTS==3) THEN
!POLYNOM
                                    SELECT CASE(variable_type)
                                      CASE(FIELD_U_VARIABLE_TYPE)
                                        SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                          CASE(NO_GLOBAL_DERIV)
                                            FACT = PERM_OVER_VIS_PARAM
                                            IF(component_idx==1) THEN
                                              !calculate u
                                              VALUE = - FACT * ( 2.0_DP*X(1) + 2.0_DP*X(2) )
                                            ELSE IF(component_idx==2) THEN
                                              !calculate v
                                              VALUE = - FACT * ( 2.0_DP*X(1) - 2.0_DP*X(2) )
                                            ELSE IF(component_idx==3) THEN
                                              !calculate p
                                              VALUE = X(1)**2.0_DP + 2.0_DP*X(1)*X(2) - X(2)**2.0_DP
                                            ELSE
                                              CALL FlagError("Not implemented.",err,error,*999)
                                            ENDIF
                                          CASE(GLOBAL_DERIV_S1)
                                            CALL FlagError("Not implemented.",err,error,*999)
                                          CASE(GLOBAL_DERIV_S2)
                                            CALL FlagError("Not implemented.",err,error,*999)
                                          CASE(GLOBAL_DERIV_S1_S2)
                                            CALL FlagError("Not implemented.",err,error,*999)
                                          CASE DEFAULT
                                            localError="The global derivative index of "//TRIM(NumberToVString( &
                                              & DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                              & err,error))//" is invalid."
                                            CALL FlagError(localError,err,error,*999)
                                        END SELECT
                                      CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                                        SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                          CASE(NO_GLOBAL_DERIV)
                                              VALUE= 0.0_DP
                                          CASE(GLOBAL_DERIV_S1)
                                            CALL FlagError("Not implemented.",err,error,*999)
                                          CASE(GLOBAL_DERIV_S2)
                                            CALL FlagError("Not implemented.",err,error,*999)
                                          CASE(GLOBAL_DERIV_S1_S2)
                                            CALL FlagError("Not implemented.",err,error,*999)
                                          CASE DEFAULT
                                            localError="The global derivative index of "//TRIM(NumberToVString( &
                                              & DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                              & err,error))//" is invalid."
                                            CALL FlagError(localError,err,error,*999)
                                        END SELECT
                                      CASE DEFAULT
                                        localError="The variable type of "//TRIM(NumberToVString(variable_type,"*",err,error))//&
                                          & " is invalid."
                                        CALL FlagError(localError,err,error,*999)
                                    END SELECT
                                  ELSE
                                    localError="The number of components does not correspond to the number of dimensions."
                                    CALL FlagError(localError,err,error,*999)
                                  ENDIF


                                CASE(EQUATIONS_SET_DARCY_EQUATION_TWO_DIM_2)
                                  IF(NUMBER_OF_DIMENSIONS==2.AND.FIELD_VARIABLE%NUMBER_OF_COMPONENTS==3) THEN
!EXPONENTIAL
                                    SELECT CASE(variable_type)
                                      CASE(FIELD_U_VARIABLE_TYPE)
                                        SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                          CASE(NO_GLOBAL_DERIV)
                                            FACT   = PERM_OVER_VIS_PARAM / L
                                            ARG(1) = X(1) / L
                                            ARG(2) = X(2) / L
                                            IF(component_idx==1) THEN
                                              !calculate u
                                              VALUE = - FACT * EXP( ARG(1) ) * EXP( ARG(2) )
                                            ELSE IF(component_idx==2) THEN
                                              !calculate v
                                              VALUE = - FACT * EXP( ARG(1) ) * EXP( ARG(2) )
                                            ELSE IF(component_idx==3) THEN
                                              !calculate p
                                              VALUE =          EXP( ARG(1) ) * EXP( ARG(2) )
                                            ELSE
                                              CALL FlagError("Not implemented.",err,error,*999)
                                            ENDIF
                                          CASE(GLOBAL_DERIV_S1)
                                            CALL FlagError("Not implemented.",err,error,*999)
                                          CASE(GLOBAL_DERIV_S2)
                                            CALL FlagError("Not implemented.",err,error,*999)
                                          CASE(GLOBAL_DERIV_S1_S2)
                                            CALL FlagError("Not implemented.",err,error,*999)
                                          CASE DEFAULT
                                            localError="The global derivative index of "//TRIM(NumberToVString( &
                                              & DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                              & err,error))//" is invalid."
                                            CALL FlagError(localError,err,error,*999)
                                        END SELECT
                                      CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                                        SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                          CASE(NO_GLOBAL_DERIV)
                                            IF(component_idx==1) THEN
                                              !calculate u
                                              VALUE= 0.0_DP
                                            ELSE IF(component_idx==2) THEN
                                              !calculate v
                                              VALUE= 0.0_DP
                                            ELSE IF(component_idx==3) THEN
                                              !calculate p
                                              VALUE= 0.0_DP
                                            ELSE
                                              CALL FlagError("Not implemented.",err,error,*999)
                                            ENDIF
                                          CASE(GLOBAL_DERIV_S1)
                                            CALL FlagError("Not implemented.",err,error,*999)
                                          CASE(GLOBAL_DERIV_S2)
                                            CALL FlagError("Not implemented.",err,error,*999)
                                          CASE(GLOBAL_DERIV_S1_S2)
                                            CALL FlagError("Not implemented.",err,error,*999)

                                          CASE DEFAULT
                                            localError="The global derivative index of "//TRIM(NumberToVString( &
                                              & DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                              & err,error))//" is invalid."
                                            CALL FlagError(localError,err,error,*999)
                                        END SELECT


                                      CASE DEFAULT
                                        localError="The variable type of "//TRIM(NumberToVString(variable_type,"*",err,error))//&
                                          & " is invalid."
                                        CALL FlagError(localError,err,error,*999)
                                    END SELECT
                                  ELSE
                                    localError="The number of components does not correspond to the number of dimensions."
                                    CALL FlagError(localError,err,error,*999)
                                  ENDIF


                                CASE(EQUATIONS_SET_DARCY_EQUATION_TWO_DIM_3)
                                  IF(NUMBER_OF_DIMENSIONS==2.AND.FIELD_VARIABLE%NUMBER_OF_COMPONENTS==3) THEN
!SINUS/COSINUS
                                    SELECT CASE(variable_type)
                                      CASE(FIELD_U_VARIABLE_TYPE)
                                        SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                          CASE(NO_GLOBAL_DERIV)
                                            FACT = 2.0_DP * PI * PERM_OVER_VIS_PARAM / L
                                            ARG(1) = 2.0_DP * PI * X(1) / L
                                            ARG(2) = 2.0_DP * PI * X(2) / L
                                            IF(component_idx==1) THEN
                                              !calculate u
                                              VALUE = - FACT * COS( ARG(1) ) * SIN( ARG(2) )
                                            ELSE IF(component_idx==2) THEN
                                              !calculate v
                                              VALUE = - FACT * SIN( ARG(1) ) * COS( ARG(2) )
                                            ELSE IF(component_idx==3) THEN
                                              !calculate p
                                              VALUE =          SIN( ARG(1) ) * SIN( ARG(2) )
                                            ELSE
                                              CALL FlagError("Not implemented.",err,error,*999)
                                            ENDIF
                                          CASE(GLOBAL_DERIV_S1)
                                            CALL FlagError("Not implemented.",err,error,*999)
                                          CASE(GLOBAL_DERIV_S2)
                                            CALL FlagError("Not implemented.",err,error,*999)
                                          CASE(GLOBAL_DERIV_S1_S2)
                                            CALL FlagError("Not implemented.",err,error,*999)
                                          CASE DEFAULT
                                            localError="The global derivative index of "//TRIM(NumberToVString( &
                                              & DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                              & err,error))//" is invalid."
                                            CALL FlagError(localError,err,error,*999)
                                        END SELECT
                                      CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                                        SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                          CASE(NO_GLOBAL_DERIV)
                                            IF(component_idx==1) THEN
                                              !calculate u
                                              VALUE=0.0_DP
                                            ELSE IF(component_idx==2) THEN
                                              !calculate v
                                              VALUE=0.0_DP
                                            ELSE IF(component_idx==3) THEN
                                              !calculate p
                                              VALUE=0.0_DP
                                            ELSE
                                              CALL FlagError("Not implemented.",err,error,*999)
                                            ENDIF
                                          CASE(GLOBAL_DERIV_S1)
                                            CALL FlagError("Not implemented.",err,error,*999)
                                          CASE(GLOBAL_DERIV_S2)
                                            CALL FlagError("Not implemented.",err,error,*999)
                                          CASE(GLOBAL_DERIV_S1_S2)
                                            CALL FlagError("Not implemented.",err,error,*999)
                                          CASE DEFAULT
                                            localError="The global derivative index of "//TRIM(NumberToVString( &
                                              & DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                              & err,error))//" is invalid."
                                            CALL FlagError(localError,err,error,*999)
                                        END SELECT
                                      CASE DEFAULT
                                        localError="The variable type of "//TRIM(NumberToVString(variable_type,"*",err,error))//&
                                          & " is invalid."
                                        CALL FlagError(localError,err,error,*999)
                                    END SELECT
                                  ELSE
                                    localError="The number of components does not correspond to the number of dimensions."
                                    CALL FlagError(localError,err,error,*999)
                                  ENDIF

                                CASE(EQUATIONS_SET_DARCY_EQUATION_THREE_DIM_1)
                                  IF(NUMBER_OF_DIMENSIONS==3.AND.FIELD_VARIABLE%NUMBER_OF_COMPONENTS==4) THEN
!POLYNOM
                                    SELECT CASE(variable_type)
                                      CASE(FIELD_U_VARIABLE_TYPE)
                                        SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                          CASE(NO_GLOBAL_DERIV)
                                            FACT = PERM_OVER_VIS_PARAM
                                            IF(component_idx==1) THEN
                                              !calculate u
                                              VALUE = - FACT * ( 2.0_DP*X(1) + 2.0_DP*X(2) + X(3) )
                                            ELSE IF(component_idx==2) THEN
                                              !calculate v
                                              VALUE = - FACT * ( 2.0_DP*X(1) - 2.0_DP*X(2) + X(3) )
                                            ELSE IF(component_idx==3) THEN
                                              !calculate w
                                              VALUE = - FACT * ( 3.0_DP + X(1) + X(2) )
                                            ELSE IF(component_idx==4) THEN
                                              !calculate p
                                              VALUE = X(1)**2.0_DP + 2.0_DP*X(1)*X(2) - X(2)**2.0_DP + &
                                                       & 3.0_DP*X(3) + X(3)*X(1) + X(3)*X(2)
                                            ELSE
                                              CALL FlagError("Not implemented.",err,error,*999)
                                            ENDIF
                                          CASE(GLOBAL_DERIV_S1)
                                            CALL FlagError("Not implemented.",err,error,*999)
                                          CASE(GLOBAL_DERIV_S2)
                                            CALL FlagError("Not implemented.",err,error,*999)
                                          CASE(GLOBAL_DERIV_S1_S2)
                                            CALL FlagError("Not implemented.",err,error,*999)
                                          CASE DEFAULT
                                            localError="The global derivative index of "//TRIM(NumberToVString( &
                                              & DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                              & err,error))//" is invalid."
                                            CALL FlagError(localError,err,error,*999)
                                        END SELECT
                                      CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                                        SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                          CASE(NO_GLOBAL_DERIV)
                                            VALUE=0.0_DP
                                          CASE(GLOBAL_DERIV_S1)
                                            CALL FlagError("Not implemented.",err,error,*999)
                                          CASE(GLOBAL_DERIV_S2)
                                            CALL FlagError("Not implemented.",err,error,*999)
                                          CASE(GLOBAL_DERIV_S1_S2)
                                            CALL FlagError("Not implemented.",err,error,*999)
                                          CASE DEFAULT
                                            localError="The global derivative index of "//TRIM(NumberToVString( &
                                              & DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                              & err,error))//" is invalid."
                                            CALL FlagError(localError,err,error,*999)
                                        END SELECT
                                      CASE DEFAULT
                                        localError="The variable type of "//TRIM(NumberToVString(variable_type,"*",err,error))//&
                                          & " is invalid."
                                        CALL FlagError(localError,err,error,*999)
                                    END SELECT
                                  ELSE
                                    localError="The number of components does not correspond to the number of dimensions."
                                    CALL FlagError(localError,err,error,*999)
                                  ENDIF


                                CASE(EQUATIONS_SET_DARCY_EQUATION_THREE_DIM_2)
                                  IF(NUMBER_OF_DIMENSIONS==3.AND.FIELD_VARIABLE%NUMBER_OF_COMPONENTS==4) THEN
!EXPONENTIAL
                                    SELECT CASE(variable_type)
                                      CASE(FIELD_U_VARIABLE_TYPE)
                                        SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                          CASE(NO_GLOBAL_DERIV)
                                            FACT   = PERM_OVER_VIS_PARAM / L
                                            ARG(1) = X(1) / L
                                            ARG(2) = X(2) / L
                                            ARG(3) = X(3) / L
                                            IF(component_idx==1) THEN
                                              !calculate u
                                              VALUE = - FACT * EXP( ARG(1) ) * EXP( ARG(2) ) * EXP( ARG(3) )
                                            ELSE IF(component_idx==2) THEN
                                              !calculate v
                                              VALUE = - FACT * EXP( ARG(1) ) * EXP( ARG(2) ) * EXP( ARG(3) )
                                            ELSE IF(component_idx==3) THEN
                                              !calculate w
                                              VALUE = - FACT * EXP( ARG(1) ) * EXP( ARG(2) ) * EXP( ARG(3) )
                                            ELSE IF(component_idx==4) THEN
                                              !calculate p
                                              VALUE =          EXP( ARG(1) ) * EXP( ARG(2) ) * EXP( ARG(3) )
                                            ELSE
                                              CALL FlagError("Not implemented.",err,error,*999)
                                            ENDIF
                                          CASE(GLOBAL_DERIV_S1)
                                            CALL FlagError("Not implemented.",err,error,*999)
                                          CASE(GLOBAL_DERIV_S2)
                                            CALL FlagError("Not implemented.",err,error,*999)
                                          CASE(GLOBAL_DERIV_S1_S2)
                                            CALL FlagError("Not implemented.",err,error,*999)
                                          CASE DEFAULT
                                            localError="The global derivative index of "//TRIM(NumberToVString( &
                                              & DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                              & err,error))//" is invalid."
                                            CALL FlagError(localError,err,error,*999)
                                        END SELECT
                                      CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                                        SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                          CASE(NO_GLOBAL_DERIV)
                                            IF(component_idx==1) THEN
                                              !calculate u
                                              VALUE=0.0_DP
                                            ELSE IF(component_idx==2) THEN
                                              !calculate v
                                              VALUE=0.0_DP
                                            ELSE IF(component_idx==3) THEN
                                              !calculate w
                                              VALUE=0.0_DP
                                            ELSE IF(component_idx==4) THEN
                                              !calculate p
                                              VALUE=0.0_DP
                                            ELSE
                                              CALL FlagError("Not implemented.",err,error,*999)
                                            ENDIF

                                          CASE(GLOBAL_DERIV_S1)
                                            CALL FlagError("Not implemented.",err,error,*999)
                                          CASE(GLOBAL_DERIV_S2)
                                            CALL FlagError("Not implemented.",err,error,*999)
                                          CASE(GLOBAL_DERIV_S1_S2)
                                            CALL FlagError("Not implemented.",err,error,*999)
                                          CASE DEFAULT
                                            localError="The global derivative index of "//TRIM(NumberToVString( &
                                              & DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                              & err,error))//" is invalid."
                                            CALL FlagError(localError,err,error,*999)
                                        END SELECT
                                      CASE DEFAULT
                                        localError="The variable type of "//TRIM(NumberToVString(variable_type,"*",err,error))//&
                                          & " is invalid."
                                        CALL FlagError(localError,err,error,*999)
                                    END SELECT
                                  ELSE
                                    localError="The number of components does not correspond to the number of dimensions."
                                    CALL FlagError(localError,err,error,*999)
                                  ENDIF

                                CASE(EQUATIONS_SET_DARCY_EQUATION_THREE_DIM_3)
                                  IF(NUMBER_OF_DIMENSIONS==3.AND.FIELD_VARIABLE%NUMBER_OF_COMPONENTS==4) THEN
  !SINE/COSINE
                                    SELECT CASE(variable_type)
                                      CASE(FIELD_U_VARIABLE_TYPE)
                                        SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                          CASE(NO_GLOBAL_DERIV)
                                            FACT = 2.0_DP * PI * PERM_OVER_VIS_PARAM / L
                                            ARG(1) = 2.0_DP * PI * X(1) / L
                                            ARG(2) = 2.0_DP * PI * X(2) / L
                                            ARG(3) = 2.0_DP * PI * X(3) / L
                                            IF(component_idx==1) THEN
                                              !calculate u
                                              VALUE = - FACT * COS( ARG(1) ) * SIN( ARG(2) )  * SIN( ARG(3) )
                                            ELSE IF(component_idx==2) THEN
                                              !calculate v
                                              VALUE = - FACT * SIN( ARG(1) ) * COS( ARG(2) )  * SIN( ARG(3) )
                                            ELSE IF(component_idx==3) THEN
                                              !calculate w
                                              VALUE = - FACT * SIN( ARG(1) ) * SIN( ARG(2) )  * COS( ARG(3) )
                                            ELSE IF(component_idx==4) THEN
                                              !calculate p
                                              VALUE =          SIN( ARG(1) ) * SIN( ARG(2) )  * SIN( ARG(3) )
                                            ELSE
                                              CALL FlagError("Not implemented.",err,error,*999)
                                            ENDIF
                                          CASE(GLOBAL_DERIV_S1)
                                            CALL FlagError("Not implemented.",err,error,*999)
                                          CASE(GLOBAL_DERIV_S2)
                                            CALL FlagError("Not implemented.",err,error,*999)
                                          CASE(GLOBAL_DERIV_S1_S2)
                                            CALL FlagError("Not implemented.",err,error,*999)
                                          CASE DEFAULT
                                            localError="The global derivative index of "//TRIM(NumberToVString( &
                                              & DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                              & err,error))//" is invalid."
                                            CALL FlagError(localError,err,error,*999)
                                        END SELECT
                                      CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                                        SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                          CASE(NO_GLOBAL_DERIV)
                                            IF(component_idx==1) THEN
                                              !calculate u
                                              VALUE=0.0_DP
                                            ELSE IF(component_idx==2) THEN
                                              !calculate v
                                              VALUE=0.0_DP
                                            ELSE IF(component_idx==3) THEN
                                              !calculate w
                                              VALUE=0.0_DP
                                            ELSE IF(component_idx==4) THEN
                                              !calculate p
                                              VALUE=0.0_DP
                                            ELSE
                                              CALL FlagError("Not implemented.",err,error,*999)
                                            ENDIF
                                          CASE(GLOBAL_DERIV_S1)
                                            CALL FlagError("Not implemented.",err,error,*999)
                                          CASE(GLOBAL_DERIV_S2)
                                            CALL FlagError("Not implemented.",err,error,*999)
                                          CASE(GLOBAL_DERIV_S1_S2)
                                            CALL FlagError("Not implemented.",err,error,*999)
                                          CASE DEFAULT
                                            localError="The global derivative index of "//TRIM(NumberToVString( &
                                              & DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                              & err,error))//" is invalid."
                                            CALL FlagError(localError,err,error,*999)
                                        END SELECT
                                      CASE DEFAULT
                                        localError="The variable type of "//TRIM(NumberToVString(variable_type,"*",err,error))//&
                                          & " is invalid."
                                        CALL FlagError(localError,err,error,*999)
                                    END SELECT
                                  ELSE
                                    localError="The number of components does not correspond to the number of dimensions."
                                    CALL FlagError(localError,err,error,*999)
                                  ENDIF
                                CASE DEFAULT
                                  localError="The analytic function type of "// &
                                    & TRIM(NumberToVString(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE,"*",err,error))// &
                                    & " is invalid."
                                  CALL FlagError(localError,err,error,*999)
                                END SELECT
                                !Default to version 1 of each node derivative
                                local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                  & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(deriv_idx)%VERSIONS(1)
                                CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD,variable_type, &
                                  & FIELD_ANALYTIC_VALUES_SET_TYPE,local_ny,VALUE,err,error,*999)
                                IF(variable_type==FIELD_U_VARIABLE_TYPE) THEN


! ! !                                 IF(DOMAIN_NODES%NODES(node_idx)%BOUNDARY_NODE) THEN
! ! !                                   !If we are a boundary node then set the analytic value on the boundary
! ! !                                   IF(component_idx<=NUMBER_OF_DIMENSIONS) THEN
! ! !                                     CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(BOUNDARY_CONDITIONS,variable_type,local_ny, &
! ! !                                       & BOUNDARY_CONDITION_FIXED,VALUE,err,error,*999)
! ! !                                   BOUND_COUNT=BOUND_COUNT+1
! ! !                                   ENDIF
! ! !                                 ELSE
! ! !                                   IF(component_idx<=NUMBER_OF_DIMENSIONS) THEN
! ! !                                     CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD,variable_type, &
! ! !                                       & FIELD_VALUES_SET_TYPE,local_ny,VALUE,err,error,*999)
! ! !                                   ENDIF
! ! !                                 ENDIF




                                  !If we are a boundary node then set the analytic value on the boundary
                                  IF(NUMBER_OF_DIMENSIONS==2) THEN
                                    IF(ABS(X(1)-BOUNDARY_X(1,1))<BOUNDARY_TOLERANCE.OR. &
                                      & ABS(X(1)-BOUNDARY_X(1,2))<BOUNDARY_TOLERANCE.OR. &
                                      & ABS(X(2)-BOUNDARY_X(2,1))<BOUNDARY_TOLERANCE.OR. &
                                      & ABS(X(2)-BOUNDARY_X(2,2))<BOUNDARY_TOLERANCE) THEN
                                      IF(component_idx<=NUMBER_OF_DIMENSIONS) THEN
                                        CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(BOUNDARY_CONDITIONS,DEPENDENT_FIELD,variable_type, &
                                          & local_ny,BOUNDARY_CONDITION_FIXED,VALUE,err,error,*999)
                                      BOUND_COUNT=BOUND_COUNT+1
!Apply boundary conditions check for pressure nodes
                                      ELSE IF(component_idx>NUMBER_OF_DIMENSIONS) THEN
                                        IF(DOMAIN%topology%elements%maximum_number_of_element_parameters==3.OR. &
                                          & DOMAIN%topology%elements%maximum_number_of_element_parameters==6.OR. &
                                          & DOMAIN%topology%elements%maximum_number_of_element_parameters==10) THEN

                                        IF(ABS(X(1)-BOUNDARY_X(1,1))<BOUNDARY_TOLERANCE.AND. &
                                          & ABS(X(2)-BOUNDARY_X(2,1))<BOUNDARY_TOLERANCE.OR. &
                                          & ABS(X(1)-BOUNDARY_X(1,1))<BOUNDARY_TOLERANCE.AND.&
                                          & ABS(X(2)-BOUNDARY_X(2,2))<BOUNDARY_TOLERANCE.OR. &
                                          & ABS(X(1)-BOUNDARY_X(1,2))<BOUNDARY_TOLERANCE.AND.&
                                          & ABS(X(2)-BOUNDARY_X(2,1))<BOUNDARY_TOLERANCE.OR. &
                                          & ABS(X(1)-BOUNDARY_X(1,2))<BOUNDARY_TOLERANCE.AND.&
                                          & ABS(X(2)-BOUNDARY_X(2,2))<BOUNDARY_TOLERANCE) &
                                          & THEN
                                            CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(BOUNDARY_CONDITIONS,DEPENDENT_FIELD, &
                                              & variable_type,local_ny,BOUNDARY_CONDITION_FIXED,VALUE,err,error,*999)
                                            BOUND_COUNT=BOUND_COUNT+1
                                        ENDIF
                                        ENDIF
                                      ENDIF
                                    ELSE
                                      IF(component_idx<=NUMBER_OF_DIMENSIONS) THEN
                                        CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD,variable_type, &
                                          & FIELD_VALUES_SET_TYPE,local_ny,VALUE,err,error,*999)
                                      ENDIF
                                    ENDIF
                                  ELSE IF(NUMBER_OF_DIMENSIONS==3) THEN
                                    IF(ABS(X(1)-BOUNDARY_X(1,1))<BOUNDARY_TOLERANCE.OR. &
                                      & ABS(X(1)-BOUNDARY_X(1,2))<BOUNDARY_TOLERANCE.OR. &
                                      & ABS(X(2)-BOUNDARY_X(2,1))<BOUNDARY_TOLERANCE.OR. &
                                      & ABS(X(2)-BOUNDARY_X(2,2))<BOUNDARY_TOLERANCE.OR. &
                                      & ABS(X(3)-BOUNDARY_X(3,1))<BOUNDARY_TOLERANCE.OR. &
                                      & ABS(X(3)-BOUNDARY_X(3,2))<BOUNDARY_TOLERANCE) THEN
                                      IF(component_idx<=NUMBER_OF_DIMENSIONS) THEN
                                        CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(BOUNDARY_CONDITIONS,DEPENDENT_FIELD,variable_type, &
                                          & local_ny,BOUNDARY_CONDITION_FIXED,VALUE,err,error,*999)
                                      BOUND_COUNT=BOUND_COUNT+1
!Apply boundary conditions check for pressure nodes
                                      ELSE IF(component_idx>NUMBER_OF_DIMENSIONS) THEN
                                        IF(DOMAIN%topology%elements%maximum_number_of_element_parameters==4.OR. &
                                          & DOMAIN%topology%elements%maximum_number_of_element_parameters==10.OR. &
                                          & DOMAIN%topology%elements%maximum_number_of_element_parameters==20) THEN
                                        IF(ABS(X(1)-BOUNDARY_X(1,1))<BOUNDARY_TOLERANCE.AND. &
                                         & ABS(X(2)-BOUNDARY_X(2,1))<BOUNDARY_TOLERANCE.AND. &
                                         & ABS(X(3)-BOUNDARY_X(3,1))<BOUNDARY_TOLERANCE.OR. &
                                         & ABS(X(1)-BOUNDARY_X(1,1))<BOUNDARY_TOLERANCE.AND. &
                                         & ABS(X(2)-BOUNDARY_X(2,1))<BOUNDARY_TOLERANCE.AND. &
                                         & ABS(X(3)-BOUNDARY_X(3,2))<BOUNDARY_TOLERANCE.OR. &
                                         & ABS(X(1)-BOUNDARY_X(1,1))<BOUNDARY_TOLERANCE.AND. &
                                         & ABS(X(2)-BOUNDARY_X(2,2))<BOUNDARY_TOLERANCE.AND. &
                                         & ABS(X(3)-BOUNDARY_X(3,1))<BOUNDARY_TOLERANCE.OR. &
                                         & ABS(X(1)-BOUNDARY_X(1,1))<BOUNDARY_TOLERANCE.AND. &
                                         & ABS(X(2)-BOUNDARY_X(2,2))<BOUNDARY_TOLERANCE.AND. &
                                         & ABS(X(3)-BOUNDARY_X(3,2))<BOUNDARY_TOLERANCE.OR. &
                                         & ABS(X(1)-BOUNDARY_X(1,2))<BOUNDARY_TOLERANCE.AND. &
                                         & ABS(X(2)-BOUNDARY_X(2,1))<BOUNDARY_TOLERANCE.AND. &
                                         & ABS(X(3)-BOUNDARY_X(3,1))<BOUNDARY_TOLERANCE.OR. &
                                         & ABS(X(1)-BOUNDARY_X(1,2))<BOUNDARY_TOLERANCE.AND. &
                                         & ABS(X(2)-BOUNDARY_X(2,1))<BOUNDARY_TOLERANCE.AND. &
                                         & ABS(X(3)-BOUNDARY_X(3,2))<BOUNDARY_TOLERANCE.OR. &
                                         & ABS(X(1)-BOUNDARY_X(1,2))<BOUNDARY_TOLERANCE.AND. &
                                         & ABS(X(2)-BOUNDARY_X(2,2))<BOUNDARY_TOLERANCE.AND. &
                                         & ABS(X(3)-BOUNDARY_X(3,1))<BOUNDARY_TOLERANCE.OR. &
                                         & ABS(X(1)-BOUNDARY_X(1,2))<BOUNDARY_TOLERANCE.AND. &
                                         & ABS(X(2)-BOUNDARY_X(2,2))<BOUNDARY_TOLERANCE.AND. &
                                         & ABS(X(3)-BOUNDARY_X(3,2))<BOUNDARY_TOLERANCE) THEN
                                           CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(BOUNDARY_CONDITIONS,DEPENDENT_FIELD, &
                                             & variable_type,local_ny,BOUNDARY_CONDITION_FIXED,VALUE,err,error,*999)
                                           BOUND_COUNT=BOUND_COUNT+1
                                        ENDIF
                                        ENDIF
                                      ENDIF
                                    ELSE
                                      IF(component_idx<=NUMBER_OF_DIMENSIONS) THEN
                                        CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD,variable_type, &
                                          & FIELD_VALUES_SET_TYPE,local_ny,VALUE,err,error,*999)
                                      ENDIF
                                    ENDIF
                                  ENDIF
                                ENDIF
                              ENDDO !deriv_idx
                            ENDDO !node_idx
                          ELSE
                            CALL FlagError("Domain topology nodes is not associated.",err,error,*999)
                          ENDIF
                        ELSE
                          CALL FlagError("Domain topology is not associated.",err,error,*999)
                        ENDIF
                      ELSE
                        CALL FlagError("Domain is not associated.",err,error,*999)
                      ENDIF
                    ELSE
                      CALL FlagError("Only node based interpolation is implemented.",err,error,*999)
                    ENDIF
                  ENDDO !component_idx
                  CALL FIELD_PARAMETER_SET_UPDATE_START(DEPENDENT_FIELD,variable_type,FIELD_ANALYTIC_VALUES_SET_TYPE, &
                    & err,error,*999)
                  CALL FIELD_PARAMETER_SET_UPDATE_FINISH(DEPENDENT_FIELD,variable_type,FIELD_ANALYTIC_VALUES_SET_TYPE, &
                    & err,error,*999)
                  CALL FIELD_PARAMETER_SET_UPDATE_START(DEPENDENT_FIELD,variable_type,FIELD_VALUES_SET_TYPE, &
                    & err,error,*999)
                  CALL FIELD_PARAMETER_SET_UPDATE_FINISH(DEPENDENT_FIELD,variable_type,FIELD_VALUES_SET_TYPE, &
                    & err,error,*999)
                ELSE
                  CALL FlagError("Field variable is not associated.",err,error,*999)
                ENDIF
              ENDDO !variable_idx
              CALL FIELD_PARAMETER_SET_DATA_RESTORE(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & GEOMETRIC_PARAMETERS,err,error,*999)
            ELSE
              CALL FlagError("Equations set boundary conditions is not associated.",err,error,*999)
            ENDIF
          ELSE
            CALL FlagError("Equations set geometric field is not associated.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Equations set dependent field is not associated.",err,error,*999)
        ENDIF
       END SELECT
      ELSE
        CALL FlagError("Equations set analytic is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF

    EXITS("Darcy_BoundaryConditionsAnalyticCalculate")
    RETURN
999 ERRORSEXITS("Darcy_BoundaryConditionsAnalyticCalculate",err,error)
    RETURN 1
  END SUBROUTINE Darcy_BoundaryConditionsAnalyticCalculate

  !
  !================================================================================================================================
  !

  !>Update geometric field for ALE Darcy problem
  SUBROUTINE Darcy_PreSolveGetSolidDisplacement(CONTROL_LOOP,SOLVER,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solvers
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER_FINITE_ELASTICITY, SOLVER_DARCY  !<A pointer to the solvers
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD_FINITE_ELASTICITY, GEOMETRIC_FIELD_DARCY
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS_FINITE_ELASTICITY, SOLVER_EQUATIONS_DARCY  !<A pointer to the solver equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING_FINITE_ELASTICITY, SOLVER_MAPPING_DARCY !<A pointer to the solver mapping
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equations_SET_FINITE_ELASTICITY, EQUATIONS_SET_DARCY !<A pointer to the equations set
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_TIME_LOOP !<A pointer to the control time loop
    TYPE(CONTROL_LOOP_TYPE), POINTER :: ROOT_CONTROL_LOOP, CONTROL_LOOP_SOLID
    TYPE(VARYING_STRING) :: localError

    REAL(DP), POINTER :: MESH_DISPLACEMENT_VALUES(:),SOLUTION_VALUES_SOLID(:)
    REAL(DP), POINTER :: DUMMY_VALUES2(:)
    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT,ALPHA

!     INTEGER(INTG) :: NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_FINITE_ELASTICITY,NUMBER_OF_COMPONENTS_GEOMETRIC_FIELD_DARCY
    INTEGER(INTG) :: NUMBER_OF_DIMENSIONS,NUMBER_OF_DOFS,NDOFS_TO_PRINT,dof_number,loop_idx
    INTEGER(INTG) :: INPUT_TYPE,INPUT_OPTION


    ENTERS("Darcy_PreSolveGetSolidDisplacement",err,error,*999)

!--- \todo : Do we need for each case a FIELD_PARAMETER_SET_UPDATE_START / FINISH on FIELD_MESH_DISPLACEMENT_SET_TYPE ?

    NULLIFY(SOLVER_FINITE_ELASTICITY)
    NULLIFY(SOLVER_DARCY)
    NULLIFY(MESH_DISPLACEMENT_VALUES)
    NULLIFY(SOLUTION_VALUES_SOLID)
    NULLIFY(ROOT_CONTROL_LOOP)
    NULLIFY(CONTROL_LOOP_SOLID)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      CONTROL_TIME_LOOP=>CONTROL_LOOP
      DO loop_idx=1,CONTROL_LOOP%CONTROL_LOOP_LEVEL
        IF(CONTROL_TIME_LOOP%LOOP_TYPE==PROBLEM_CONTROL_TIME_LOOP_TYPE) THEN
          CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_TIME_LOOP,CURRENT_TIME,TIME_INCREMENT,err,error,*999)
          EXIT
        ENDIF
        IF (ASSOCIATED(CONTROL_LOOP%PARENT_LOOP)) THEN
          CONTROL_TIME_LOOP=>CONTROL_TIME_LOOP%PARENT_LOOP
        ELSE
          CALL FlagError("Could not find a time control loop.",err,error,*999)
        ENDIF
      ENDDO

      IF(DIAGNOSTICS1) THEN
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE, &
          & "*******************************************************************************************************", &
          & err,error,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"CURRENT_TIME   = ",CURRENT_TIME,err,error,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"TIME_INCREMENT = ",TIME_INCREMENT,err,error,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE, &
          & "*******************************************************************************************************", &
          & err,error,*999)
      ENDIF

      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          ROOT_CONTROL_LOOP=>CONTROL_LOOP%PROBLEM%CONTROL_LOOP
          IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a Darcy equation problem.",err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_STANDARD_DARCY_SUBTYPE)
              ! do nothing
            CASE(PROBLEM_QUASISTATIC_DARCY_SUBTYPE)
              ! do nothing
            CASE(PROBLEM_TRANSIENT_DARCY_SUBTYPE)
              ! do nothing
            CASE(PROBLEM_ALE_DARCY_SUBTYPE)
              !--- Motion: specified
              IF(SOLVER%GLOBAL_NUMBER==SOLVER_NUMBER_DARCY) THEN
                SOLVER_EQUATIONS_DARCY=>SOLVER%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS_DARCY)) THEN
                  SOLVER_MAPPING_DARCY=>SOLVER_EQUATIONS_DARCY%SOLVER_MAPPING
                  IF(ASSOCIATED(SOLVER_MAPPING_DARCY)) THEN
                    EQUATIONS_SET_DARCY=>SOLVER_MAPPING_DARCY%EQUATIONS_SETS(1)%ptr
                    IF(ASSOCIATED(EQUATIONS_SET_DARCY)) THEN
                      IF(.NOT.ALLOCATED(equations_set_darcy%specification)) THEN
                        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
                      ELSE IF(SIZE(equations_set_darcy%specification,1)/=3) THEN
                        CALL FlagError("Equations set specification must have three entries for a Darcy type equations set.", &
                          & err,error,*999)
                      END IF
                      SELECT CASE(EQUATIONS_SET_DARCY%SPECIFICATION(3))
                        CASE(EQUATIONS_SET_STANDARD_DARCY_SUBTYPE)
                          ! do nothing
                        CASE(EQUATIONS_SET_QUASISTATIC_DARCY_SUBTYPE)
                          ! do nothing
                        CASE(EQUATIONS_SET_ALE_DARCY_SUBTYPE,EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE, &
                          & EQUATIONS_SET_TRANSIENT_ALE_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
                          & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE)
                          CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Darcy motion specified ... ",err,error,*999)
                          GEOMETRIC_FIELD_DARCY=>EQUATIONS_SET_DARCY%GEOMETRY%GEOMETRIC_FIELD
                          IF(ASSOCIATED(GEOMETRIC_FIELD_DARCY)) THEN
                            ALPHA = 0.085_DP * sin( 2.0_DP * PI * CURRENT_TIME / 4.0_DP )

                            CALL FIELD_PARAMETER_SETS_COPY(GEOMETRIC_FIELD_DARCY,FIELD_U_VARIABLE_TYPE, &
                              & FIELD_INITIAL_VALUES_SET_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE,ALPHA,err,error,*999)
                          ELSE
                            CALL FlagError("Geometric field is not associated.",err,error,*999)
                          ENDIF
                        CASE DEFAULT
                          localError="Equations set subtype " &
                            & //TRIM(NumberToVString(EQUATIONS_SET_DARCY%SPECIFICATION(3),"*",err,error))// &
                            & " is not valid for a Darcy equation fluid type of a fluid mechanics problem class."
                          CALL FlagError(localError,err,error,*999)
                      END SELECT
                    ELSE
                      CALL FlagError("Equations set is not associated.",err,error,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Solver mapping is not associated.",err,error,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Solver equations is not associated.",err,error,*999)
                ENDIF
              ELSE
                ! do nothing
              ENDIF
            CASE(PROBLEM_PGM_DARCY_SUBTYPE,PROBLEM_PGM_TRANSIENT_DARCY_SUBTYPE,PROBLEM_PGM_ELASTICITY_DARCY_SUBTYPE)
              !--- Motion: read in from a file
              IF(SOLVER%GLOBAL_NUMBER==SOLVER_NUMBER_DARCY) THEN
                CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,SOLVER_NUMBER_DARCY,SOLVER_DARCY,err,error,*999)
                SOLVER_EQUATIONS_DARCY=>SOLVER_DARCY%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS_DARCY)) THEN
                  SOLVER_MAPPING_DARCY=>SOLVER_EQUATIONS_DARCY%SOLVER_MAPPING
                  IF(ASSOCIATED(SOLVER_MAPPING_DARCY)) THEN
                    EQUATIONS_SET_DARCY=>SOLVER_MAPPING_DARCY%EQUATIONS_SETS(1)%ptr
                    IF(ASSOCIATED(EQUATIONS_SET_DARCY)) THEN
                      GEOMETRIC_FIELD_DARCY=>EQUATIONS_SET_DARCY%GEOMETRY%GEOMETRIC_FIELD
                    ELSE
                      CALL FlagError("Darcy equations set is not associated.",err,error,*999)
                    END IF
                    IF(SOLVER%outputType>=SOLVER_PROGRESS_OUTPUT) THEN
                      CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Darcy motion read from a file ... ",err,error,*999)
                    ENDIF

                    CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET_DARCY%GEOMETRY%GEOMETRIC_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS,err,error,*999)

                    !Copy input to Darcy' geometric field
                    INPUT_TYPE=42
                    INPUT_OPTION=2
                    CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET_DARCY%GEOMETRY%GEOMETRIC_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE,MESH_DISPLACEMENT_VALUES,err,error,*999)
                    CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,MESH_DISPLACEMENT_VALUES, &
                      & NUMBER_OF_DIMENSIONS,INPUT_TYPE,INPUT_OPTION,CONTROL_LOOP%TIME_LOOP%ITERATION_NUMBER,1.0_DP, &
                      & err,error,*999)
                    CALL FIELD_PARAMETER_SET_UPDATE_START(EQUATIONS_SET_DARCY%GEOMETRY%GEOMETRIC_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE,err,error,*999)
                    CALL FIELD_PARAMETER_SET_UPDATE_FINISH(EQUATIONS_SET_DARCY%GEOMETRY%GEOMETRIC_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE,err,error,*999)
                  ELSE
                    CALL FlagError("Darcy solver mapping is not associated.",err,error,*999)
                  END IF
                ELSE
                  CALL FlagError("Darcy solver equations are not associated.",err,error,*999)
                END IF

               IF(DIAGNOSTICS1) THEN
                 NDOFS_TO_PRINT = SIZE(MESH_DISPLACEMENT_VALUES,1)
                 CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NDOFS_TO_PRINT,NDOFS_TO_PRINT,NDOFS_TO_PRINT,&
                   & MESH_DISPLACEMENT_VALUES,'(" MESH_DISPLACEMENT_VALUES = ",4(X,E13.6))','4(4(X,E13.6))', &
                   & err,error,*999)
               ENDIF
              ELSE
                ! in case of a solver number different from 2: do nothing
              ENDIF
            CASE(PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE,PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE, &
              & PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)
              !--- Motion: defined by fluid-solid interaction (thus read from solid's dependent field)
              IF(SOLVER%GLOBAL_NUMBER==SOLVER_NUMBER_DARCY) THEN  !It is called with 'SOLVER%GLOBAL_NUMBER=SOLVER_NUMBER_DARCY', otherwise it doesn't work
                !--- Get the dependent field of the finite elasticity equations
                IF(SOLVER%outputType>=SOLVER_PROGRESS_OUTPUT) THEN
                  CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Darcy motion read from solid's dependent field ... ",err,error,*999)
                ENDIF
                SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
                CASE(PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE)
                  CALL CONTROL_LOOP_GET(ROOT_CONTROL_LOOP,[1,CONTROL_LOOP_NODE],CONTROL_LOOP_SOLID,err,error,*999)
                CASE(PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE,PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)
                  CALL CONTROL_LOOP_GET(ROOT_CONTROL_LOOP,[1,1,CONTROL_LOOP_NODE],CONTROL_LOOP_SOLID,err,error,*999)
                END SELECT
                CALL SOLVERS_SOLVER_GET(CONTROL_LOOP_SOLID%SOLVERS,SOLVER_NUMBER_SOLID, &
                    & SOLVER_FINITE_ELASTICITY,err,error,*999)
                SOLVER_EQUATIONS_FINITE_ELASTICITY=>SOLVER_FINITE_ELASTICITY%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS_FINITE_ELASTICITY)) THEN
                  SOLVER_MAPPING_FINITE_ELASTICITY=>SOLVER_EQUATIONS_FINITE_ELASTICITY%SOLVER_MAPPING
                  IF(ASSOCIATED(SOLVER_MAPPING_FINITE_ELASTICITY)) THEN
                    EQUATIONS_SET_FINITE_ELASTICITY=>SOLVER_MAPPING_FINITE_ELASTICITY%EQUATIONS_SETS(1)%ptr
                    IF(ASSOCIATED(EQUATIONS_SET_FINITE_ELASTICITY)) THEN
                      DEPENDENT_FIELD_FINITE_ELASTICITY=>EQUATIONS_SET_FINITE_ELASTICITY%DEPENDENT%DEPENDENT_FIELD
                      IF(ASSOCIATED(DEPENDENT_FIELD_FINITE_ELASTICITY)) THEN
                          !No longer needed, since no more 'Field_ParametersToFieldParametersCopy'
!                         CALL FIELD_NUMBER_OF_COMPONENTS_GET(DEPENDENT_FIELD_FINITE_ELASTICITY, &
!                           & FIELD_U_VARIABLE_TYPE,NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_FINITE_ELASTICITY,err,error,*999)
                      ELSE
                        CALL FlagError("DEPENDENT_FIELD_FINITE_ELASTICITY is not associated.",err,error,*999)
                      END IF
                    ELSE
                      CALL FlagError("Finite elasticity equations set is not associated.",err,error,*999)
                    END IF
                  ELSE
                    CALL FlagError("Finite elasticity solver mapping is not associated.",err,error,*999)
                  END IF
                ELSE
                  CALL FlagError("Finite elasticity solver equations are not associated.",err,error,*999)
                END IF

                !--- Get the geometric field for the ALE Darcy equations
                CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,SOLVER_NUMBER_DARCY,SOLVER_DARCY,err,error,*999)
                SOLVER_EQUATIONS_DARCY=>SOLVER_DARCY%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS_DARCY)) THEN
                  SOLVER_MAPPING_DARCY=>SOLVER_EQUATIONS_DARCY%SOLVER_MAPPING
                  IF(ASSOCIATED(SOLVER_MAPPING_DARCY)) THEN
                    EQUATIONS_SET_DARCY=>SOLVER_MAPPING_DARCY%EQUATIONS_SETS(1)%ptr
                    IF(ASSOCIATED(EQUATIONS_SET_DARCY)) THEN
                      GEOMETRIC_FIELD_DARCY=>EQUATIONS_SET_DARCY%GEOMETRY%GEOMETRIC_FIELD
                      IF(ASSOCIATED(GEOMETRIC_FIELD_DARCY)) THEN
                          !No longer needed, since no more 'Field_ParametersToFieldParametersCopy'
!                         CALL FIELD_NUMBER_OF_COMPONENTS_GET(GEOMETRIC_FIELD_DARCY, &
!                           & FIELD_U_VARIABLE_TYPE,NUMBER_OF_COMPONENTS_GEOMETRIC_FIELD_DARCY,err,error,*999)
                      ELSE
                        CALL FlagError("GEOMETRIC_FIELD_DARCY is not associated.",err,error,*999)
                      END IF
                    ELSE
                      CALL FlagError("Darcy equations set is not associated.",err,error,*999)
                    END IF
                  ELSE
                    CALL FlagError("Darcy solver mapping is not associated.",err,error,*999)
                  END IF
                ELSE
                  CALL FlagError("Darcy solver equations are not associated.",err,error,*999)
                END IF

                !--- Copy the result from Finite-elasticity's dependent field to ALE Darcy's geometric field
                !--- First: FIELD_MESH_DISPLACEMENT_SET_TYPE = - FIELD_PREVIOUS_VALUES_SET_TYPE
                ALPHA=-1.0_DP
                CALL FIELD_PARAMETER_SETS_COPY(GEOMETRIC_FIELD_DARCY,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_PREVIOUS_VALUES_SET_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE,ALPHA,err,error,*999)

                ! Write 'FIELD_PREVIOUS_VALUES_SET_TYPE'
                IF(DIAGNOSTICS3) THEN
                  NULLIFY( DUMMY_VALUES2 )
                  CALL FIELD_PARAMETER_SET_DATA_GET(GEOMETRIC_FIELD_DARCY,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_PREVIOUS_VALUES_SET_TYPE,DUMMY_VALUES2,err,error,*999)
                  NDOFS_TO_PRINT = SIZE(DUMMY_VALUES2,1)
                  CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NDOFS_TO_PRINT,NDOFS_TO_PRINT,NDOFS_TO_PRINT,DUMMY_VALUES2, &
                    & '(" GEOMETRIC_FIELD_DARCY,FIELD_U_VARIABLE_TYPE,FIELD_PREVIOUS_VALUES_SET_TYPE = ",4(X,E13.6))',&
                    & '4(4(X,E13.6))',err,error,*999)
                  CALL FIELD_PARAMETER_SET_DATA_RESTORE(GEOMETRIC_FIELD_DARCY,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_PREVIOUS_VALUES_SET_TYPE,DUMMY_VALUES2,err,error,*999)
                ENDIF

                !--- Second: Get a pointer to the solution values of the solid
                !    (deformed absolute positions in x, y, z; possibly solid pressure)
                CALL FIELD_PARAMETER_SET_DATA_GET(DEPENDENT_FIELD_FINITE_ELASTICITY,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,SOLUTION_VALUES_SOLID,err,error,*999)
!                 CALL FIELD_PARAMETER_SET_DATA_RESTORE(DEPENDENT_FIELD_FINITE_ELASTICITY,FIELD_U_VARIABLE_TYPE, &
!                   & FIELD_VALUES_SET_TYPE,SOLUTION_VALUES_SOLID,err,error,*999) ! necessary ???

                ! Write 'DEPENDENT_FIELD_FINITE_ELASTICITY'
                IF(DIAGNOSTICS3) THEN
                  NULLIFY( DUMMY_VALUES2 )
                  CALL FIELD_PARAMETER_SET_DATA_GET(DEPENDENT_FIELD_FINITE_ELASTICITY,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,DUMMY_VALUES2,err,error,*999)
                  NDOFS_TO_PRINT = SIZE(DUMMY_VALUES2,1)
                  CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NDOFS_TO_PRINT,NDOFS_TO_PRINT,NDOFS_TO_PRINT,DUMMY_VALUES2, &
                    & '(" DEPENDENT_FIELD_FINITE_ELASTICITY,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE = ",4(X,E13.6))',&
                    & '4(4(X,E13.6))',err,error,*999)
                  CALL FIELD_PARAMETER_SET_DATA_RESTORE(DEPENDENT_FIELD_FINITE_ELASTICITY,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,DUMMY_VALUES2,err,error,*999)
                ENDIF

                !--- Third: FIELD_MESH_DISPLACEMENT_SET_TYPE += Deformed absolute position of solid
                NUMBER_OF_DOFS = GEOMETRIC_FIELD_DARCY%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%ptr%NUMBER_OF_DOFS
                DO dof_number=1,NUMBER_OF_DOFS
                  ! assumes fluid-geometry and solid-dependent mesh are identical \todo: introduce check
                  CALL FIELD_PARAMETER_SET_ADD_LOCAL_DOF(GEOMETRIC_FIELD_DARCY, &
                    & FIELD_U_VARIABLE_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE,dof_number, &
                    & SOLUTION_VALUES_SOLID(dof_number), &
                    & err,error,*999)

!---              !!! Why not directly do the mesh update here ??? !!!
                  CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(GEOMETRIC_FIELD_DARCY, &
                    & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,dof_number, &
                    & SOLUTION_VALUES_SOLID(dof_number), &
                    & err,error,*999)
!---

                END DO
                CALL FIELD_PARAMETER_SET_UPDATE_START(GEOMETRIC_FIELD_DARCY, &
                  & FIELD_U_VARIABLE_TYPE, FIELD_MESH_DISPLACEMENT_SET_TYPE,err,error,*999)
                CALL FIELD_PARAMETER_SET_UPDATE_FINISH(GEOMETRIC_FIELD_DARCY, &
                  & FIELD_U_VARIABLE_TYPE, FIELD_MESH_DISPLACEMENT_SET_TYPE,err,error,*999)
                !
                CALL FIELD_PARAMETER_SET_UPDATE_START(GEOMETRIC_FIELD_DARCY, &
                  & FIELD_U_VARIABLE_TYPE, FIELD_VALUES_SET_TYPE,err,error,*999)
                CALL FIELD_PARAMETER_SET_UPDATE_FINISH(GEOMETRIC_FIELD_DARCY, &
                  & FIELD_U_VARIABLE_TYPE, FIELD_VALUES_SET_TYPE,err,error,*999)

                ! Write 'FIELD_MESH_DISPLACEMENT_SET_TYPE'
                IF(DIAGNOSTICS3) THEN
                  NULLIFY( DUMMY_VALUES2 )
                  CALL FIELD_PARAMETER_SET_DATA_GET(GEOMETRIC_FIELD_DARCY,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_MESH_DISPLACEMENT_SET_TYPE,DUMMY_VALUES2,err,error,*999)
                  NDOFS_TO_PRINT = SIZE(DUMMY_VALUES2,1)
                  CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NDOFS_TO_PRINT,NDOFS_TO_PRINT,NDOFS_TO_PRINT,DUMMY_VALUES2, &
                    & '(" GEOMETRIC_FIELD_DARCY,FIELD_U_VARIABLE_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE = ",4(X,E13.6))',&
                    & '4(4(X,E13.6))',err,error,*999)
                  CALL FIELD_PARAMETER_SET_DATA_RESTORE(GEOMETRIC_FIELD_DARCY,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_MESH_DISPLACEMENT_SET_TYPE,DUMMY_VALUES2,err,error,*999)
                ENDIF
              ELSE
                ! do nothing
              END IF
            CASE DEFAULT
              localError="Problem subtype "//TRIM(NumberToVString(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",err,error))// &
                & " is not valid for a Darcy equation fluid type of a fluid mechanics problem class."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ELSE
          CALL FlagError("Problem is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Solver is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",err,error,*999)
    ENDIF

    EXITS("Darcy_PreSolveGetSolidDisplacement")
    RETURN
999 ERRORSEXITS("Darcy_PreSolveGetSolidDisplacement",err,error)
    RETURN 1
    
  END SUBROUTINE Darcy_PreSolveGetSolidDisplacement

  !
  !================================================================================================================================

  !>Store solution of previous subiteration iterate
  SUBROUTINE Darcy_PreSolveStorePreviousIterate(CONTROL_LOOP,SOLVER,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solvers
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS  !<A pointer to the solver equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equations_SET !<A pointer to the equations set
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: localError

    REAL(DP) :: ALPHA
    INTEGER(INTG) :: FIELD_VAR_TYPE,equations_set_idx


    ENTERS("Darcy_PreSolveStorePreviousIterate",err,error,*999)

    NULLIFY(DEPENDENT_FIELD)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(SOLVER_MAPPING)
    NULLIFY(EQUATIONS_SET)
    NULLIFY(vectorMapping)
    NULLIFY(FIELD_VARIABLE)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(SOLVER%GLOBAL_NUMBER==SOLVER_NUMBER_DARCY) THEN
          IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
            IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
              CALL FlagError("Problem specification is not allocated.",err,error,*999)
            ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
              CALL FlagError("Problem specification must have three entries for a Darcy equation problem.",err,error,*999)
            END IF
            SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
              CASE(PROBLEM_STANDARD_DARCY_SUBTYPE)
                ! do nothing
              CASE(PROBLEM_QUASISTATIC_DARCY_SUBTYPE)
                ! do nothing
              CASE(PROBLEM_ALE_DARCY_SUBTYPE,PROBLEM_PGM_DARCY_SUBTYPE,PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE, &
                & PROBLEM_PGM_ELASTICITY_DARCY_SUBTYPE,PROBLEM_PGM_TRANSIENT_DARCY_SUBTYPE, &
                & PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE,PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)
                SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                  SOLVER_MAPPING=>SOLVER_equations%SOLVER_MAPPING
                  IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                    !loop over the equations sets
                   DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                    EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%ptr
                    IF(ASSOCIATED(EQUATIONS_SET)) THEN
                      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
                        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
                      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
                        CALL FlagError("Equations set specification must have three entries for a Darcy type equations set.", &
                          & err,error,*999)
                      END IF
                      SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
                        CASE(EQUATIONS_SET_STANDARD_DARCY_SUBTYPE)
                          ! do nothing
                        CASE(EQUATIONS_SET_QUASISTATIC_DARCY_SUBTYPE)
                          ! do nothing
                        CASE(EQUATIONS_SET_ALE_DARCY_SUBTYPE,EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE, &
                          & EQUATIONS_SET_TRANSIENT_ALE_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
                          & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE, &
                          & EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
                          DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD

                          IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                            write(*,*)'-------------------------------------------------------'
                            write(*,*)'+++     Storing previous subiteration iterate       +++'
                            write(*,*)'-------------------------------------------------------'
                            !--- Store the DEPENDENT field values of the previous subiteration iterate
                            vectorMapping=>EQUATIONS_SET%equations%vectorEquations%vectorMapping
                            IF(ASSOCIATED(vectorMapping)) THEN
                              SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
                              CASE(EQUATIONS_SET_ALE_DARCY_SUBTYPE,EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE)
                              FIELD_VARIABLE=>vectorMapping%linearMapping%equationsMatrixToVarMaps(1)%VARIABLE
                              ! '1' associated with linear matrix
                              CASE(EQUATIONS_SET_TRANSIENT_ALE_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
                                  & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE, &
                                  & EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
                                FIELD_VARIABLE=>vectorMapping%dynamicMapping%equationsMatrixToVarMaps(1)%VARIABLE
                              END SELECT
                              IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                                FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
                                ALPHA = 1.0_DP
                                CALL FIELD_PARAMETER_SETS_COPY(DEPENDENT_FIELD,FIELD_VAR_TYPE, &
                                  & FIELD_VALUES_SET_TYPE,FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,ALPHA,err,error,*999)
                              ELSE
                                CALL FlagError("FIELD_VAR_TYPE is not associated.",err,error,*999)
                              ENDIF
                            ELSE
                              CALL FlagError("vectorMapping is not associated.",err,error,*999)
                            ENDIF
                          ELSE
                            CALL FlagError("Dependent field is not associated.",err,error,*999)
                          ENDIF
                        CASE DEFAULT
                          localError="Equations set subtype " &
                            & //TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
                            & " is not valid for a Darcy equation fluid type of a fluid mechanics problem class."
                          CALL FlagError(localError,err,error,*999)
                      END SELECT
                    ELSE
                      CALL FlagError("Equations set is not associated.",err,error,*999)
                    ENDIF
                   ENDDO
                  ELSE
                    CALL FlagError("Solver mapping is not associated.",err,error,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Solver equations is not associated.",err,error,*999)
                ENDIF
              CASE DEFAULT
                localError="Problem subtype "//TRIM(NumberToVString(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",err,error))// &
                  & " is not valid for a Darcy equation fluid type of a fluid mechanics problem class."
                CALL FlagError(localError,err,error,*999)
            END SELECT
          ELSE
            CALL FlagError("Problem is not associated.",err,error,*999)
          ENDIF
        ELSE
          ! do nothing
        ENDIF
      ELSE
        CALL FlagError("Solver is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",err,error,*999)
    ENDIF

    EXITS("Darcy_PreSolveStorePreviousIterate")
    RETURN
999 ERRORSEXITS("Darcy_PreSolveStorePreviousIterate",err,error)
    RETURN 1
    
  END SUBROUTINE Darcy_PreSolveStorePreviousIterate

  !
  !================================================================================================================================
  !
  !updates the boundary conditions etc to the required analytic values
  !for the case EquationsSetIncompElastDarcyAnalyticDarcy the pressure field obtained from the finite elasticity solve is overwritten
  !by the appropriate mass increase for that time step
  SUBROUTINE Darcy_PreSolveUpdateAnalyticValues(CONTROL_LOOP,SOLVER,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,GEOMETRIC_FIELD
!    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE,GEOMETRIC_VARIABLE
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS  !<A pointer to the solver equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equations_SET !<A pointer to the equations set
    TYPE(EquationsType), POINTER :: equations
!    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: DOMAIN_TOPOLOGY
    TYPE(VARYING_STRING) :: localError
!    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE
!    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
!    REAL(DP), POINTER :: BOUNDARY_VALUES(:)
    REAL(DP), POINTER :: GEOMETRIC_PARAMETERS(:)
    INTEGER(INTG) :: NUMBER_OF_DIMENSIONS
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_TIME_LOOP

    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT
!     REAL(DP) :: k_xx, k_yy, k_zz
    INTEGER(INTG) :: eqnset_idx,loop_idx
    INTEGER(INTG) :: VARIABLE_TYPE !<The field variable type to add \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    REAL(DP) :: A1,D1
!    INTEGER(INTG) :: FIELD_SET_TYPE !<The field parameter set identifier \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
!    INTEGER(INTG) :: DERIVATIVE_NUMBER !<The node derivative number
!    INTEGER(INTG) :: COMPONENT_NUMBER !<The field variable component number
!    INTEGER(INTG) :: TOTAL_NUMBER_OF_NODES !<The total number of (geometry) nodes
!    INTEGER(INTG) :: LOCAL_NODE_NUMBER
!    INTEGER(INTG) :: equations_SET_IDX
!    INTEGER(INTG) :: equations_row_number

    ENTERS("Darcy_PreSolveUpdateAnalyticValues",err,error,*999)


    A1 = 0.4_DP
    D1=1.0_DP

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      CONTROL_TIME_LOOP=>CONTROL_LOOP
      DO loop_idx=1,CONTROL_LOOP%CONTROL_LOOP_LEVEL
        IF(CONTROL_TIME_LOOP%LOOP_TYPE==PROBLEM_CONTROL_TIME_LOOP_TYPE) THEN
          CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_TIME_LOOP,CURRENT_TIME,TIME_INCREMENT,err,error,*999)
          EXIT
        ENDIF
        IF (ASSOCIATED(CONTROL_LOOP%PARENT_LOOP)) THEN
          CONTROL_TIME_LOOP=>CONTROL_TIME_LOOP%PARENT_LOOP
        ELSE
          CALL FlagError("Could not find a time control loop.",err,error,*999)
        ENDIF
      ENDDO

!     IF(ASSOCIATED(CONTROL_LOOP)) THEN
!       CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_LOOP,CURRENT_TIME,TIME_INCREMENT,err,error,*999)
       !write(*,*)'CURRENT_TIME = ',CURRENT_TIME
       !write(*,*)'TIME_INCREMENT = ',TIME_INCREMENT
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a Darcy equation problem.",err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE,PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)
                SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                !loop over all the equation sets and set the appropriate field variable type BCs and
                !the source field associated with each equation set
                DO eqnset_idx=1,SOLVER_equations%SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                  SOLVER_MAPPING=>SOLVER_equations%SOLVER_MAPPING
                  EQUATIONS=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(eqnset_idx)%EQUATIONS
                  IF(ASSOCIATED(EQUATIONS)) THEN
                    EQUATIONS_SET=>equations%equationsSet
                    IF(ASSOCIATED(EQUATIONS_SET)) THEN
                     IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
                      IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_INCOMP_ELAST_DARCY_ANALYTIC_DARCY)THEN
                      !for this analytic case we copy the mass variable to the pressure variable
                        DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                        IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                          GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                          IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN
                            CALL FIELD_NUMBER_OF_COMPONENTS_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,&
                              & NUMBER_OF_DIMENSIONS,err,error,*999)
                            NULLIFY(GEOMETRIC_VARIABLE)
                            NULLIFY(GEOMETRIC_PARAMETERS)
                            CALL Field_VariableGet(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,GEOMETRIC_VARIABLE,err,error,*999)
                            CALL FIELD_PARAMETER_SET_DATA_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,&
                              & GEOMETRIC_PARAMETERS,err,error,*999)
                             EQUATIONS_SET%ANALYTIC%ANALYTIC_USER_PARAMS(1)=CURRENT_TIME
!                              DO variable_idx=1,DEPENDENT_FIELD%NUMBER_OF_VARIABLES
                              !variable_type=DEPENDENT_FIELD%VARIABLES(2*eqnset_idx-1)%VARIABLE_TYPE
                              variable_type=FIELD_V_VARIABLE_TYPE
                              FIELD_VARIABLE=>DEPENDENT_FIELD%VARIABLE_TYPE_MAP(variable_type)%ptr
                              IF(ASSOCIATED(FIELD_VARIABLE)) THEN
!                                DO component_idx=4,FIELD_VARIABLE%NUMBER_OF_COMPONENTS


                               CALL Field_ParametersToFieldParametersCopy(DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                                 & FIELD_VALUES_SET_TYPE,4,DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                                 & FIELD_VALUES_SET_TYPE,4,err,error,*999)

!                                   IF(FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE== &
!                                     & FIELD_NODE_BASED_INTERPOLATION) THEN
!                                     DOMAIN=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN
!                                     IF(ASSOCIATED(DOMAIN)) THEN
!                                       IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
!                                         DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
!                                         IF(ASSOCIATED(DOMAIN_NODES)) THEN
!                                           !Loop over the local nodes excluding the ghosts.
!                                           DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
!                                           CALL FIELD_PARAMETER_SET_GET_NODE(DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
!                                              & FIELD_VALUES_SET_TYPE,1,node_idx,4,MASS_INCREASE,err,error,*999)
!                                           CALL FIELD_PARAMETER_SET_UPDATE_NODE(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
!                                              & FIELD_VALUES_SET_TYPE,1,node_idx,4,0.1*MASS_INCREASE,err,error,*999)
!                                           write(*,*) MASS_INCREASE
!
!                                             !!TODO \todo We should interpolate the geometric field here and the node position.
! !                                             DO dim_idx=1,NUMBER_OF_DIMENSIONS
! !                                               local_ny= &
! !                                           & GEOMETRIC_VARIABLE%COMPONENTS(dim_idx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP(1,node_idx)
! !                                               X(dim_idx)=GEOMETRIC_PARAMETERS(local_ny)
! !                                             ENDDO !dim_idx
! !                                             !Loop over the derivatives
! !                                             DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
! !                                               ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE
! !                                               GLOBAL_DERIV_INDEX=DOMAIN_NODES%NODES(node_idx)%GLOBAL_DERIVATIVE_INDEX(deriv_idx)
! ! !                                               CALL DIFFUSION_EQUATION_ANALYTIC_FUNCTIONS(VALUE,X, &
! ! !                                                 & CURRENT_TIME,variable_type,GLOBAL_DERIV_INDEX, &
! ! !                                                 & ANALYTIC_FUNCTION_TYPE,err,error,*999)
! !                                               local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
! !                                                 & NODE_PARAM2DOF_MAP(deriv_idx,node_idx)
! !                                               CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD,variable_type, &
! !                                                 & FIELD_ANALYTIC_VALUES_SET_TYPE,local_ny,VALUE,err,error,*999)
! ! !                                               BOUNDARY_CONDITION_CHECK_VARIABLE=SOLVER_equations%BOUNDARY_CONDITIONS% &
! ! !                                                 & BOUNDARY_CONDITIONS_VARIABLE_TYPE_MAP(variable_type)%ptr% &
! ! !                                                 & CONDITION_TYPES(local_ny)
! ! !                                               IF(BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_FIXED) THEN
! ! !                                                CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD, &
! ! !                                                  & variable_type,FIELD_VALUES_SET_TYPE,local_ny, &
! ! !                                                  & VALUE,err,error,*999)
! ! !                                               ENDIF
! !
! ! !                                              IF(variable_type==FIELD_U_VARIABLE_TYPE) THEN
! ! !                                                IF(DOMAIN_NODES%NODES(node_idx)%BOUNDARY_NODE) THEN
! !                                                   !If we are a boundary node then set the analytic value on the boundary
! ! !                                                  CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(BOUNDARY_CONDITIONS,variable_type,local_ny, &
! ! !                                                    & BOUNDARY_CONDITION_FIXED,VALUE,err,error,*999)
! ! !                                                ENDIF
! ! !                                              ENDIF
! !                                             ENDDO !deriv_idx
!                                           ENDDO !node_idx
!                                         ELSE
!                                           CALL FlagError("Domain topology nodes is not associated.",err,error,*999)
!                                         ENDIF
!                                       ELSE
!                                         CALL FlagError("Domain topology is not associated.",err,error,*999)
!                                       ENDIF
!                                     ELSE
!                                       CALL FlagError("Domain is not associated.",err,error,*999)
!                                     ENDIF
!                                   ELSE
!                                     CALL FlagError("Only node based interpolation is implemented.",err,error,*999)
!                                   ENDIF
!                                 ENDDO !component_idx
                                CALL FIELD_PARAMETER_SET_UPDATE_START(DEPENDENT_FIELD,variable_type, &
                                 & FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
                                CALL FIELD_PARAMETER_SET_UPDATE_FINISH(DEPENDENT_FIELD,variable_type, &
                                 & FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
                                CALL FIELD_PARAMETER_SET_UPDATE_START(DEPENDENT_FIELD,variable_type, &
                                 & FIELD_VALUES_SET_TYPE,err,error,*999)
                                CALL FIELD_PARAMETER_SET_UPDATE_FINISH(DEPENDENT_FIELD,variable_type, &
                                 & FIELD_VALUES_SET_TYPE,err,error,*999)
                              ELSE
                                CALL FlagError("Field variable is not associated.",err,error,*999)
                              ENDIF

!                              ENDDO !variable_idx
                             CALL FIELD_PARAMETER_SET_DATA_RESTORE(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,&
                              & FIELD_VALUES_SET_TYPE,GEOMETRIC_PARAMETERS,err,error,*999)
                          ELSE
                            CALL FlagError("Equations set geometric field is not associated.",err,error,*999)
                          ENDIF
                        ELSE
                          CALL FlagError("Equations set dependent field is not associated.",err,error,*999)
                        ENDIF
                       ENDIF
                      ELSE
                        !CALL FlagError("Equations set analytic is not associated.",err,error,*999)
                      ENDIF
                    ELSE
                      CALL FlagError("Equations set is not associated.",err,error,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Equations are not associated.",err,error,*999)
                  END IF
!                 ELSE
!                   CALL FlagError("Solver equations are not associated.",err,error,*999)
!                 END IF
                CALL FIELD_PARAMETER_SET_UPDATE_START(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,err,error,*999)
                CALL FIELD_PARAMETER_SET_UPDATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,err,error,*999)
!             IF(CONTROL_LOOP%PROBLEM%SPECIFICATION(3)==PROBLEM_LINEAR_SOURCE_DIFFUSION_SUBTYPE)THEN
!             !>Set the source field to a specified analytical function
!MAY NEED TO USE THIS ULTIMATELY - BUT WILL REQUIRE IMPLEMENTING SOURCE FIELD & VECTOR FUNCTIONALITY FOR DARCY EQUATION
!             IF(ASSOCIATED(EQUATIONS_SET)) THEN
!               IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
!                 sourceField=>EQUATIONS_SET%SOURCE%SOURCE_FIELD
!                 IF(ASSOCIATED(sourceField)) THEN
!                   GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
!                   IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN
!                     CALL FIELD_NUMBER_OF_COMPONENTS_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS,err,error,*999)
!                     NULLIFY(GEOMETRIC_VARIABLE)
!                     CALL Field_VariableGet(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,GEOMETRIC_VARIABLE,err,error,*999)
!                     CALL FIELD_PARAMETER_SET_DATA_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
!                       & GEOMETRIC_PARAMETERS,err,error,*999)
!                       variable_type=FIELD_U_VARIABLE_TYPE
!                       FIELD_VARIABLE=>sourceField%VARIABLE_TYPE_MAP(variable_type)%ptr
!                       IF(ASSOCIATED(FIELD_VARIABLE)) THEN
!                         DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
!                           IF(FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION) THEN
!                             DOMAIN=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN
!                             IF(ASSOCIATED(DOMAIN)) THEN
!                               IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
!                                 DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
!                                 IF(ASSOCIATED(DOMAIN_NODES)) THEN
!                                   !Loop over the local nodes excluding the ghosts.
!                                   DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
!                                     !!TODO \todo We should interpolate the geometric field here and the node position.
!                                     DO dim_idx=1,NUMBER_OF_DIMENSIONS
!                                       local_ny=GEOMETRIC_VARIABLE%COMPONENTS(dim_idx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP(1,&
!                                       & node_idx)
!                                       X(dim_idx)=GEOMETRIC_PARAMETERS(local_ny)
!                                     ENDDO !dim_idx
!                                     !Loop over the derivatives
!                                     DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
!                                       SELECT CASE(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE)
!                                       CASE(EQUATIONS_SET_DIFFUSION_EQUATION_THREE_DIM_1)
!                                           VALUE_SOURCE=-1*A1*EXP(-1*CURRENT_TIME)*(X(1)*X(1)+X(2)*X(2)+X(3)*X(3)+6)
!                                       CASE DEFAULT
!                                         localError="The analytic function type of "// &
!                                           & TRIM(NumberToVString(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE,"*",err,error))//&
!                                           & " is invalid."
!                                         CALL FlagError(localError,err,error,*999)
!                                       END SELECT
!                                       local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
!                                         & NODE_PARAM2DOF_MAP(deriv_idx,node_idx)
!                                       CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(sourceField,FIELD_U_VARIABLE_TYPE, &
!                                         & FIELD_VALUES_SET_TYPE,local_ny,VALUE_SOURCE,err,error,*999)
!                                     ENDDO !deriv_idx
!                                   ENDDO !node_idx
!                                 ELSE
!                                   CALL FlagError("Domain topology nodes is not associated.",err,error,*999)
!                                 ENDIF
!                               ELSE
!                                 CALL FlagError("Domain topology is not associated.",err,error,*999)
!                               ENDIF
!                             ELSE
!                               CALL FlagError("Domain is not associated.",err,error,*999)
!                             ENDIF
!                           ELSE
!                             CALL FlagError("Only node based interpolation is implemented.",err,error,*999)
!                           ENDIF
!                         ENDDO !component_idx
!                         CALL FIELD_PARAMETER_SET_UPDATE_START(sourceField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
!                           & err,error,*999)
!                         CALL FIELD_PARAMETER_SET_UPDATE_FINISH(sourceField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
!                           & err,error,*999)
!                       ELSE
!                         CALL FlagError("Field variable is not associated.",err,error,*999)
!                       ENDIF
!                     CALL FIELD_PARAMETER_SET_DATA_RESTORE(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
!                       & GEOMETRIC_PARAMETERS,err,error,*999)
!                   ELSE
!                     CALL FlagError("Equations set geometric field is not associated.",err,error,*999)
!                   ENDIF
!                 ELSE
!                   CALL FlagError("Equations set source field is not associated.",err,error,*999)
!                 ENDIF
!               ELSE
!                 CALL FlagError("Equations set analytic is not associated.",err,error,*999)
!               ENDIF
!             ELSE
!               CALL FlagError("Equations set is not associated.",err,error,*999)
!             ENDIF
!             ENDIF
            ENDDO !eqnset_idx
                ELSE
                  CALL FlagError("Solver equations are not associated.",err,error,*999)
                END IF
          CASE DEFAULT
              localError="Problem subtype "//TRIM(NumberToVString(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",err,error))// &
                & " is not valid for a Darcy equation type of a fluid mechanics problem class."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ELSE
          CALL FlagError("Problem is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Solver is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",err,error,*999)
    ENDIF
    
    EXITS("Darcy_PreSolveUpdateAnalyticValues")
    RETURN
999 ERRORSEXITS("Darcy_PreSolveUpdateAnalyticValues",err,error)
    RETURN 1
    
  END SUBROUTINE Darcy_PreSolveUpdateAnalyticValues

  !
  !================================================================================================================================
  !
  !> Monitor convergence of the Darcy solution
  SUBROUTINE DARCY_EQUATION_MONITOR_CONVERGENCE(CONTROL_LOOP,SOLVER,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS  !<A pointer to the solver equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equations_SET !<A pointer to the equations set
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: localError
    CHARACTER(25) :: FILENAME
    TYPE(VARYING_STRING) :: FILEPATH

    REAL(DP), POINTER :: ITERATION_VALUES_N(:),ITERATION_VALUES_N1(:)
    REAL(DP) :: RESIDUAL_NORM

    REAL(DP), PARAMETER :: RESIDUAL_TOLERANCE_RELATIVE=1.0E-05_DP
    REAL(DP), PARAMETER :: RESIDUAL_TOLERANCE_ABSOLUTE=1.0E-10_DP

    INTEGER(INTG) :: FIELD_VAR_TYPE
    INTEGER(INTG) :: dof_number,NUMBER_OF_DOFS,equations_set_idx
    INTEGER(INTG) :: COMPUTATIONAL_NODE_NUMBER
    INTEGER(INTG) :: FILEUNIT_N, FILEUNIT_N1

    ENTERS("DARCY_EQUATION_MONITOR_CONVERGENCE",err,error,*999)

    NULLIFY(DEPENDENT_FIELD)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(SOLVER_MAPPING)
    NULLIFY(EQUATIONS_SET)
    NULLIFY(EQUATIONS)
    NULLIFY(vectorMapping)
    NULLIFY(FIELD_VARIABLE)

    COMPUTATIONAL_NODE_NUMBER=ComputationalEnvironment_NodeNumberGet(err,error)
    WRITE(FILENAME,'("Darcy_",I3.3,".conv")') COMPUTATIONAL_NODE_NUMBER
    FILEPATH = "./output/"//FILENAME
    OPEN(UNIT=23, FILE=CHAR(FILEPATH),STATUS='unknown',ACCESS='append')

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(SOLVER%GLOBAL_NUMBER==SOLVER_NUMBER_DARCY) THEN
          IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
            IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
              CALL FlagError("Problem specification is not allocated.",err,error,*999)
            ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
              CALL FlagError("Problem specification must have three entries for a Darcy equation problem.",err,error,*999)
            END IF
            SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
              CASE(PROBLEM_STANDARD_DARCY_SUBTYPE,PROBLEM_QUASISTATIC_DARCY_SUBTYPE,PROBLEM_TRANSIENT_DARCY_SUBTYPE, &
                & PROBLEM_ALE_DARCY_SUBTYPE,PROBLEM_PGM_DARCY_SUBTYPE,PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE, &
                & PROBLEM_PGM_ELASTICITY_DARCY_SUBTYPE,PROBLEM_PGM_TRANSIENT_DARCY_SUBTYPE)
                ! do nothing
              CASE(PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE,PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)
                SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                  SOLVER_MAPPING=>SOLVER_equations%SOLVER_MAPPING
                  IF(ASSOCIATED(SOLVER_MAPPING)) THEN
!                     EQUATIONS=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(1)%EQUATIONS
!                     IF(ASSOCIATED(EQUATIONS)) THEN
!                       EQUATIONS_SET=>equations%equationsSet
                   DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                    EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%ptr
                      IF(ASSOCIATED(EQUATIONS_SET)) THEN
                        IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
                          CALL FlagError("Equations set specification is not allocated.",err,error,*999)
                        ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
                          CALL FlagError("Equations set specification must have three entries for a Darcy type equations set.", &
                            & err,error,*999)
                        END IF
                        SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
                          CASE(EQUATIONS_SET_STANDARD_DARCY_SUBTYPE,EQUATIONS_SET_QUASISTATIC_DARCY_SUBTYPE, &
                            & EQUATIONS_SET_ALE_DARCY_SUBTYPE,EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE, &
                            & EQUATIONS_SET_TRANSIENT_ALE_DARCY_SUBTYPE)
                            ! do nothing
                          CASE(EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
                              & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE, &
                              & EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
                            CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Darcy monitor convergence ... ",err,error,*999)
                            DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                            IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                              vectorMapping=>EQUATIONS_SET%equations%vectorEquations%vectorMapping
                              IF(ASSOCIATED(vectorMapping)) THEN
                                SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
                                CASE(EQUATIONS_SET_ALE_DARCY_SUBTYPE,EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE)
                                  FIELD_VARIABLE=>vectorMapping%linearMapping%equationsMatrixToVarMaps(1)%VARIABLE
                                  ! '1' associated with linear matrix
                                CASE(EQUATIONS_SET_TRANSIENT_ALE_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
                                    & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE, &
                                    & EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
                                  FIELD_VARIABLE=>vectorMapping%dynamicMapping%equationsMatrixToVarMaps(1)%VARIABLE
                                END SELECT
                                IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                                  FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE

                                  !iter 1
                                  NULLIFY(ITERATION_VALUES_N)
                                  CALL FIELD_PARAMETER_SET_DATA_GET(DEPENDENT_FIELD,FIELD_VAR_TYPE, &
                                    & FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,ITERATION_VALUES_N,err,error,*999)

                                  !iter 2
                                  NULLIFY(ITERATION_VALUES_N1)
                                  CALL FIELD_PARAMETER_SET_DATA_GET(DEPENDENT_FIELD,FIELD_VAR_TYPE, &
                                    & FIELD_VALUES_SET_TYPE,ITERATION_VALUES_N1,err,error,*999)

                                  RESIDUAL_NORM = 0.0_DP
                                  NUMBER_OF_DOFS = DEPENDENT_FIELD%VARIABLE_TYPE_MAP(FIELD_VAR_TYPE)%ptr%NUMBER_OF_DOFS
                                  DO dof_number=1,NUMBER_OF_DOFS
                                    RESIDUAL_NORM = RESIDUAL_NORM + &
                                      & ( ITERATION_VALUES_N1(dof_number) - ITERATION_VALUES_N(dof_number) )**2.0_DP
                                  END DO
                                  RESIDUAL_NORM = SQRT(RESIDUAL_NORM / NUMBER_OF_DOFS)

                                  IF(CONTROL_LOOP%LOOP_TYPE==PROBLEM_CONTROL_WHILE_LOOP_TYPE) THEN
                                    IF(CONTROL_LOOP%WHILE_LOOP%ITERATION_NUMBER>=2) THEN !Omit initialised solution
                                      IF(CONTROL_LOOP%WHILE_LOOP%ITERATION_NUMBER==2) THEN
                                      RESIDUAL_NORM_0 = RESIDUAL_NORM
                                      WRITE(23,*) 'RESIDUAL_NORM_0 = ',RESIDUAL_NORM_0
                                      WRITE(23,*) 'R / R0 :'
                                      ENDIF
                                      write(*,*)'-------------------------------------------------------'
                                      write(*,*)'+++     RESIDUAL_NORM   =        +++',RESIDUAL_NORM
                                      write(*,*)'+++     RESIDUAL_NORM_0 =        +++',RESIDUAL_NORM_0
                                      write(*,*)'+++     R / R_0         =        +++',RESIDUAL_NORM / RESIDUAL_NORM_0
                                      write(*,*)'-------------------------------------------------------'
                                      WRITE(23,*) RESIDUAL_NORM / RESIDUAL_NORM_0

                                      !End subiteration loop if residual is small relative to residual in first step
                                      IF((RESIDUAL_NORM/RESIDUAL_NORM_0)<=RESIDUAL_TOLERANCE_RELATIVE .OR. &
                                        & RESIDUAL_NORM<=RESIDUAL_TOLERANCE_ABSOLUTE ) THEN
                                        write(*,*)'++++++++++++++++++++++++++++++++++++'
                                        write(*,*)'+++    SUBITERATION CONVERGED    +++'
                                        write(*,*)'++++++++++++++++++++++++++++++++++++'
                                        CONTROL_LOOP%WHILE_LOOP%CONTINUE_LOOP=.FALSE.
                                      ELSE IF(CONTROL_LOOP%WHILE_LOOP%ITERATION_NUMBER== &
                                          & CONTROL_LOOP%WHILE_LOOP%MAXIMUM_NUMBER_OF_ITERATIONS) THEN
                                        CALL FLAG_WARNING("Subiterations between solid and fluid "// &
                                            & "equations did not converge.",err,error,*999)
                                      ENDIF
                                    ENDIF
                                  ELSE
                                    CALL FlagError("DARCY_EQUATION_MONITOR_CONVERGENCE must be called "// &
                                        & "with a while control loop",err,error,*999)
                                  ENDIF


!                                   SUBITERATION_NUMBER = CONTROL_LOOP%WHILE_LOOP%ITERATION_NUMBER
!
!                                   WRITE(FILENAME,'("Darcy_DOFs_N_",I2.2,".dat")') SUBITERATION_NUMBER
!                                   FILEPATH = "./output/"//FILENAME
!                                   FILEUNIT_N = 7777 + 2*SUBITERATION_NUMBER
!                                   OPEN(UNIT=FILEUNIT_N,FILE=CHAR(FILEPATH),STATUS='unknown',ACCESS='append')
!                                   DO dof_number=1,NUMBER_OF_DOFS
!                                     WRITE(FILEUNIT_N,*) ITERATION_VALUES_N(dof_number)
!                                   END DO
!
!
!                                   WRITE(FILENAME,'("Darcy_DOFs_N1_",I2.2,".dat")') SUBITERATION_NUMBER
!                                   FILEPATH = "./output/"//FILENAME
!                                   FILEUNIT_N1 = 7777 + 2*SUBITERATION_NUMBER+1
!                                   OPEN(UNIT=FILEUNIT_N1,FILE=CHAR(FILEPATH),STATUS='unknown',ACCESS='append')
!                                   DO dof_number=1,NUMBER_OF_DOFS
!                                     WRITE(FILEUNIT_N1,*) ITERATION_VALUES_N1(dof_number)
!                                   END DO


                                  CALL FIELD_PARAMETER_SET_DATA_RESTORE(DEPENDENT_FIELD,FIELD_VAR_TYPE, &
                                    & FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,ITERATION_VALUES_N,err,error,*999)
                                  CALL FIELD_PARAMETER_SET_DATA_RESTORE(DEPENDENT_FIELD,FIELD_VAR_TYPE, &
                                    & FIELD_VALUES_SET_TYPE,ITERATION_VALUES_N1,err,error,*999)

                                ELSE
                                  CALL FlagError("FIELD_VAR_TYPE is not associated.",err,error,*999)
                                ENDIF
                              ELSE
                                CALL FlagError("vectorMapping is not associated.",err,error,*999)
                              ENDIF
                            ELSE
                              CALL FlagError("Dependent field is not associated.",err,error,*999)
                            END IF
                          CASE DEFAULT
                            localError="Equations set subtype " &
                              & //TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
                              & " is not valid for a Darcy equation fluid type of a fluid mechanics problem class."
                            CALL FlagError(localError,err,error,*999)
                        END SELECT
                      ELSE
                        CALL FlagError("Equations set is not associated.",err,error,*999)
                      END IF
                     ENDDO
!                     ELSE
!                       CALL FlagError("Equations are not associated.",err,error,*999)
!                     END IF
                  ELSE
                    CALL FlagError("Solver mapping is not associated.",err,error,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Solver equations are not associated.",err,error,*999)
                END IF
              CASE DEFAULT
                localError="Problem subtype "//TRIM(NumberToVString(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",err,error))// &
                  & " is not valid for a Darcy equation fluid type of a fluid mechanics problem class."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          ELSE
            CALL FlagError("Problem is not associated.",err,error,*999)
          ENDIF
        ELSE
          ! do nothing
        ENDIF
      ELSE
        CALL FlagError("Solver is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",err,error,*999)
    ENDIF

    CLOSE(23)
    CLOSE(FILEUNIT_N)
    CLOSE(FILEUNIT_N1)

    EXITS("DARCY_EQUATION_MONITOR_CONVERGENCE")
    RETURN
999 ERRORSEXITS("DARCY_EQUATION_MONITOR_CONVERGENCE",err,error)
    RETURN 1
  END SUBROUTINE DARCY_EQUATION_MONITOR_CONVERGENCE

  !
  !================================================================================================================================
  !

  !> Accelerate convergence of the Darcy solution
  SUBROUTINE DARCY_EQUATION_ACCELERATE_CONVERGENCE(CONTROL_LOOP,SOLVER,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS  !<A pointer to the solver equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equations_SET !<A pointer to the equations set
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: localError

    REAL(DP), POINTER :: ITERATION_VALUES_N(:),ITERATION_VALUES_N1(:)
    REAL(DP) :: RELAXATION_PARAM,ACCELERATED_VALUE

    INTEGER(INTG) :: FIELD_VAR_TYPE
    INTEGER(INTG) :: dof_number,NUMBER_OF_DOFS,equations_set_idx


    ENTERS("DARCY_EQUATION_ACCELERATE_CONVERGENCE",err,error,*999)

    NULLIFY(DEPENDENT_FIELD)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(SOLVER_MAPPING)
    NULLIFY(EQUATIONS_SET)
    NULLIFY(EQUATIONS)
    NULLIFY(vectorMapping)
    NULLIFY(FIELD_VARIABLE)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(SOLVER%GLOBAL_NUMBER==SOLVER_NUMBER_DARCY) THEN
          IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
            IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
              CALL FlagError("Problem specification is not allocated.",err,error,*999)
            ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
              CALL FlagError("Problem specification must have three entries for a Darcy equation problem.",err,error,*999)
            END IF
            SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
              CASE(PROBLEM_STANDARD_DARCY_SUBTYPE,PROBLEM_QUASISTATIC_DARCY_SUBTYPE,PROBLEM_TRANSIENT_DARCY_SUBTYPE, &
                & PROBLEM_ALE_DARCY_SUBTYPE,PROBLEM_PGM_DARCY_SUBTYPE,PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE, &
                & PROBLEM_PGM_ELASTICITY_DARCY_SUBTYPE,PROBLEM_PGM_TRANSIENT_DARCY_SUBTYPE)
                ! do nothing
              CASE(PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE,PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)
                SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                  SOLVER_MAPPING=>SOLVER_equations%SOLVER_MAPPING
                  IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                   DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                    EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%ptr
!                     EQUATIONS=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(1)%EQUATIONS
!                     IF(ASSOCIATED(EQUATIONS)) THEN
!                       EQUATIONS_SET=>equations%equationsSet
                      IF(ASSOCIATED(EQUATIONS_SET)) THEN
                        IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
                          CALL FlagError("Equations set specification is not allocated.",err,error,*999)
                        ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
                          CALL FlagError("Equations set specification must have three entries for a Darcy type equations set.", &
                            & err,error,*999)
                        END IF
                        SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
                          CASE(EQUATIONS_SET_STANDARD_DARCY_SUBTYPE,EQUATIONS_SET_QUASISTATIC_DARCY_SUBTYPE, &
                            & EQUATIONS_SET_ALE_DARCY_SUBTYPE,EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE, &
                            & EQUATIONS_SET_TRANSIENT_ALE_DARCY_SUBTYPE)
                            ! do nothing
                          CASE(EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
                              & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE, &
                              & EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
!                             CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Darcy accelerate convergence ... ",err,error,*999)
                            DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                            IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                              vectorMapping=>EQUATIONS_SET%equations%vectorEquations%vectorMapping
                              IF(ASSOCIATED(vectorMapping)) THEN
                                SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
                                CASE(EQUATIONS_SET_ALE_DARCY_SUBTYPE,EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE)
                                  FIELD_VARIABLE=>vectorMapping%linearMapping%equationsMatrixToVarMaps(1)%VARIABLE
                                  ! '1' associated with linear matrix
                                CASE(EQUATIONS_SET_TRANSIENT_ALE_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
                                    & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE, &
                                    & EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
                                  FIELD_VARIABLE=>vectorMapping%dynamicMapping%equationsMatrixToVarMaps(1)%VARIABLE
                                END SELECT
                                IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                                  FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE

                                  !iter 1
                                  NULLIFY(ITERATION_VALUES_N)
                                  CALL FIELD_PARAMETER_SET_DATA_GET(DEPENDENT_FIELD,FIELD_VAR_TYPE, &
                                    & FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,ITERATION_VALUES_N,err,error,*999)

                                  !iter 2
                                  NULLIFY(ITERATION_VALUES_N1)
                                  CALL FIELD_PARAMETER_SET_DATA_GET(DEPENDENT_FIELD,FIELD_VAR_TYPE, &
                                    & FIELD_VALUES_SET_TYPE,ITERATION_VALUES_N1,err,error,*999)

!                                   RESIDUAL_NORM = 0.0_DP
                                  NUMBER_OF_DOFS = DEPENDENT_FIELD%VARIABLE_TYPE_MAP(FIELD_VAR_TYPE)%ptr%NUMBER_OF_DOFS

!                                   DO dof_number=1,NUMBER_OF_DOFS
!                                     RESIDUAL_NORM = RESIDUAL_NORM + &
!                                       & ( ITERATION_VALUES_N1(dof_number) - ITERATION_VALUES_N(dof_number) )**2.0_DP
!                                   END DO
!                                   RESIDUAL_NORM = SQRT(RESIDUAL_NORM / NUMBER_OF_DOFS)

                                  RELAXATION_PARAM = 2.0_DP  !\ToDo Devise better way of determining optimal Aitken parameter

                                  IF( CONTROL_LOOP%WHILE_LOOP%ITERATION_NUMBER>2 )THEN
                                    CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Darcy accelerate convergence ... ",err,error,*999)
                                    DO dof_number=1,NUMBER_OF_DOFS
                                      ACCELERATED_VALUE = ITERATION_VALUES_N(dof_number) &
                                        & + RELAXATION_PARAM * ( ITERATION_VALUES_N1(dof_number) - ITERATION_VALUES_N(dof_number) )
                                      CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD, &
                                        & FIELD_VAR_TYPE,FIELD_VALUES_SET_TYPE,dof_number, &
                                        & ACCELERATED_VALUE,err,error,*999)
                                    END DO
                                    CALL FIELD_PARAMETER_SET_UPDATE_START(DEPENDENT_FIELD, &
                                      & FIELD_VAR_TYPE, FIELD_VALUES_SET_TYPE,err,error,*999)
                                    CALL FIELD_PARAMETER_SET_UPDATE_FINISH(DEPENDENT_FIELD, &
                                      & FIELD_VAR_TYPE, FIELD_VALUES_SET_TYPE,err,error,*999)
                                  END IF
                                  CALL FIELD_PARAMETER_SET_DATA_RESTORE(DEPENDENT_FIELD,FIELD_VAR_TYPE, &
                                    & FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,ITERATION_VALUES_N,err,error,*999)
                                  CALL FIELD_PARAMETER_SET_DATA_RESTORE(DEPENDENT_FIELD,FIELD_VAR_TYPE, &
                                    & FIELD_VALUES_SET_TYPE,ITERATION_VALUES_N1,err,error,*999)

                                ELSE
                                  CALL FlagError("FIELD_VAR_TYPE is not associated.",err,error,*999)
                                ENDIF
                              ELSE
                                CALL FlagError("vectorMapping is not associated.",err,error,*999)
                              ENDIF
                            ELSE
                              CALL FlagError("Dependent field is not associated.",err,error,*999)
                            END IF
                          CASE DEFAULT
                            localError="Equations set subtype " &
                              & //TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
                              & " is not valid for a Darcy equation fluid type of a fluid mechanics problem class."
                            CALL FlagError(localError,err,error,*999)
                        END SELECT
                      ELSE
                        CALL FlagError("Equations set is not associated.",err,error,*999)
                      END IF
                    ENDDO
!                     ELSE
!                       CALL FlagError("Equations are not associated.",err,error,*999)
!                     END IF
                  ELSE
                    CALL FlagError("Solver mapping is not associated.",err,error,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Solver equations are not associated.",err,error,*999)
                END IF
              CASE DEFAULT
                localError="Problem subtype "//TRIM(NumberToVString(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",err,error))// &
                  & " is not valid for a Darcy equation fluid type of a fluid mechanics problem class."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          ELSE
            CALL FlagError("Problem is not associated.",err,error,*999)
          ENDIF
        ELSE
          ! do nothing
        ENDIF
      ELSE
        CALL FlagError("Solver is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",err,error,*999)
    ENDIF

    EXITS("DARCY_EQUATION_ACCELERATE_CONVERGENCE")
    RETURN
999 ERRORSEXITS("DARCY_EQUATION_ACCELERATE_CONVERGENCE",err,error)
    RETURN 1
  END SUBROUTINE DARCY_EQUATION_ACCELERATE_CONVERGENCE

  !
  !================================================================================================================================
  !



  !================================================================================================================================
  !

  !> Allows to set an explicit Darcy mass increase to test finite elasticity
  !> (and only then this function is called, but not for the coupled problem)
  SUBROUTINE DARCY_EQUATION_POST_SOLVE_SET_MASS_INCREASE(CONTROL_LOOP,SOLVER,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solvers
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER_FINITE_ELASTICITY, SOLVER_DARCY  !<A pointer to the solvers
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD_DARCY
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS_DARCY  !<A pointer to the solver equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING_DARCY !<A pointer to the solver mapping
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equations_SET_DARCY !<A pointer to the equations set
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_TIME_LOOP !<A pointer to the control time loop
    TYPE(CONTROL_LOOP_TYPE), POINTER :: ROOT_CONTROL_LOOP, CONTROL_LOOP_SOLID
    TYPE(VARYING_STRING) :: localError

    REAL(DP), POINTER :: MESH_DISPLACEMENT_VALUES(:),SOLUTION_VALUES_SOLID(:)
    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT,ALPHA

!     INTEGER(INTG) :: NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_FINITE_ELASTICITY,NUMBER_OF_COMPONENTS_GEOMETRIC_FIELD_DARCY
    INTEGER(INTG) :: dof_number,loop_idx,NUMBER_OF_DOFS


    ENTERS("DARCY_EQUATION_POST_SOLVE_SET_MASS_INCREASE",err,error,*999)

    NULLIFY(SOLVER_FINITE_ELASTICITY)
    NULLIFY(SOLVER_DARCY)
    NULLIFY(MESH_DISPLACEMENT_VALUES)
    NULLIFY(SOLUTION_VALUES_SOLID)
    NULLIFY(ROOT_CONTROL_LOOP)
    NULLIFY(CONTROL_LOOP_SOLID)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      CONTROL_TIME_LOOP=>CONTROL_LOOP
      DO loop_idx=1,CONTROL_LOOP%CONTROL_LOOP_LEVEL
        IF(CONTROL_TIME_LOOP%LOOP_TYPE==PROBLEM_CONTROL_TIME_LOOP_TYPE) THEN
          CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_TIME_LOOP,CURRENT_TIME,TIME_INCREMENT,err,error,*999)
          EXIT
        ENDIF
        IF (ASSOCIATED(CONTROL_LOOP%PARENT_LOOP)) THEN
          CONTROL_TIME_LOOP=>CONTROL_TIME_LOOP%PARENT_LOOP
        ELSE
          CALL FlagError("Could not find a time control loop.",err,error,*999)
        ENDIF
      ENDDO

      IF(DIAGNOSTICS1) THEN
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE, &
          & "*******************************************************************************************************", &
          & err,error,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"CURRENT_TIME   = ",CURRENT_TIME,err,error,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"TIME_INCREMENT = ",TIME_INCREMENT,err,error,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE, &
          & "*******************************************************************************************************", &
          & err,error,*999)
      ENDIF

      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          ROOT_CONTROL_LOOP=>CONTROL_LOOP%PROBLEM%CONTROL_LOOP
          IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a Darcy equation problem.",err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_STANDARD_DARCY_SUBTYPE)
              ! do nothing
            CASE(PROBLEM_QUASISTATIC_DARCY_SUBTYPE)
              ! do nothing
            CASE(PROBLEM_TRANSIENT_DARCY_SUBTYPE)
              ! do nothing
            CASE(PROBLEM_ALE_DARCY_SUBTYPE)
              ! do nothing
            CASE(PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE,PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE, &
              & PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)
              !--- Mass increase specified
              IF(SOLVER%GLOBAL_NUMBER==SOLVER_NUMBER_DARCY) THEN  !It is called with 'SOLVER%GLOBAL_NUMBER=SOLVER_NUMBER_DARCY', otherwise it doesn't work
                !--- Get the dependent field of the Darcy equations
                CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,SOLVER_NUMBER_DARCY,SOLVER_DARCY,err,error,*999)
                SOLVER_EQUATIONS_DARCY=>SOLVER_DARCY%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS_DARCY)) THEN
                  SOLVER_MAPPING_DARCY=>SOLVER_EQUATIONS_DARCY%SOLVER_MAPPING
                  IF(ASSOCIATED(SOLVER_MAPPING_DARCY)) THEN
                    EQUATIONS_SET_DARCY=>SOLVER_MAPPING_DARCY%EQUATIONS_SETS(1)%ptr
                    IF(ASSOCIATED(EQUATIONS_SET_DARCY)) THEN
                      DEPENDENT_FIELD_DARCY=>EQUATIONS_SET_DARCY%DEPENDENT%DEPENDENT_FIELD
                      IF(ASSOCIATED(DEPENDENT_FIELD_DARCY)) THEN
                        ! do nothing
                      ELSE
                        CALL FlagError("GEOMETRIC_FIELD_DARCY is not associated.",err,error,*999)
                      END IF
                    ELSE
                      CALL FlagError("Darcy equations set is not associated.",err,error,*999)
                    END IF
                  ELSE
                    CALL FlagError("Darcy solver mapping is not associated.",err,error,*999)
                  END IF
                ELSE
                  CALL FlagError("Darcy solver equations are not associated.",err,error,*999)
                END IF

                ! Set the mass increase for Darcy dependent field (u, v, w; m)

!                 ALPHA = 2.0E-03_DP

!                 ALPHA = 5.0E-04_DP * CURRENT_TIME / TIME_INCREMENT

                ALPHA = 5.0E-04_DP * SIN(2.0_DP * PI * CURRENT_TIME / TIME_INCREMENT / 20.0_DP)

                write(*,*)'ALPHA = ',ALPHA

                NUMBER_OF_DOFS = DEPENDENT_FIELD_DARCY%VARIABLE_TYPE_MAP(FIELD_V_VARIABLE_TYPE)%ptr%NUMBER_OF_DOFS

                DO dof_number = 3/4*NUMBER_OF_DOFS + 1, NUMBER_OF_DOFS
                  !'3/4' only works for equal order interpolation in (u,v,w) and p
                  CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD_DARCY, &
                    & FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,dof_number, &
                    & ALPHA,err,error,*999)
                END DO
                CALL FIELD_PARAMETER_SET_UPDATE_START(DEPENDENT_FIELD_DARCY, &
                  & FIELD_U_VARIABLE_TYPE, FIELD_VALUES_SET_TYPE,err,error,*999)
                CALL FIELD_PARAMETER_SET_UPDATE_FINISH(DEPENDENT_FIELD_DARCY, &
                  & FIELD_U_VARIABLE_TYPE, FIELD_VALUES_SET_TYPE,err,error,*999)

              ELSE
                ! do nothing
              END IF
            CASE DEFAULT
              localError="Problem subtype "//TRIM(NumberToVString(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",err,error))// &
                & " is not valid for a Darcy equation fluid type of a fluid mechanics problem class."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ELSE
          CALL FlagError("Problem is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Solver is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",err,error,*999)
    ENDIF

    EXITS("DARCY_EQUATION_POST_SOLVE_SET_MASS_INCREASE")
    RETURN
999 ERRORSEXITS("DARCY_EQUATION_POST_SOLVE_SET_MASS_INCREASE",err,error)
    RETURN 1
  END SUBROUTINE DARCY_EQUATION_POST_SOLVE_SET_MASS_INCREASE

  !
  !================================================================================================================================
  !

  !\ToDo: enable this penalty formulation also for (quasi-)static; as made available in solver_routines

  !Adds a penalty term to the equilibrium equations to enforce impermeability at certain boundaries
  ! derived from: "FINITE_ELASTICITY_SURFACE_PRESSURE_RESIDUAL_EVALUATE"; same restrictions apply
  SUBROUTINE DARCY_EQUATION_IMPERMEABLE_BC_VIA_PENALTY(EQUATIONS_SET,ELEMENT_NUMBER,err,error,*)
    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equations_SET !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FIELD_TYPE), POINTER :: dependentField, independentField
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DECOMPOSITION_ELEMENT_TYPE), POINTER :: DECOMP_ELEMENT
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatrixType), POINTER :: stiffnessMatrix
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: GEOMETRIC_INTERPOLATION_PARAMETERS
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: GEOMETRIC_INTERPOLATED_POINT
    TYPE(DECOMPOSITION_FACE_TYPE), POINTER :: DECOMP_FACE
    TYPE(DOMAIN_FACE_TYPE), POINTER :: DOMAIN_FACE
    TYPE(BASIS_TYPE), POINTER :: FACE_BASIS,DEPENDENT_BASIS
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: FACE_QUADRATURE_SCHEME
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: FACE_VELOCITY_INTERPOLATION_PARAMETERS
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: FACE_INTERPOLATED_POINT
    INTEGER(INTG) :: MESH_COMPONENT_NUMBER
    INTEGER(INTG) :: element_face_idx,face_number,normal_component_idx,gauss_idx
    INTEGER(INTG) :: FACE_NUMBER_OF_GAUSS_POINTS
    INTEGER(INTG) :: component_idx_1,element_base_dof_idx_1,face_node_idx_1
    INTEGER(INTG) :: element_node_derivative_idx_1,element_dof_idx_1,element_node_idx_1,parameter_idx_1
    INTEGER(INTG) :: face_parameter_idx_1,face_node_derivative_idx_1
    INTEGER(INTG) :: component_idx_2,element_base_dof_idx_2,face_node_idx_2
    INTEGER(INTG) :: element_node_derivative_idx_2,element_dof_idx_2,element_node_idx_2,parameter_idx_2
    INTEGER(INTG) :: face_parameter_idx_2,face_node_derivative_idx_2

    REAL(DP) :: GAUSS_WEIGHT,NORMAL_PROJECTION_1,NORMAL_PROJECTION_2, PENALTY_PARAM
    REAL(DP) :: DZDXI(3,3),DZDXIT(3,3),GIJL(3,3),GIJU(3,3),G,SQRT_G, PGM, PGN, SUM
    LOGICAL :: IMPERMEABLE_BC

    ENTERS("DARCY_EQUATION_IMPERMEABLE_BC_VIA_PENALTY",err,error,*999)

    !Make this routine conditional on (stiffnessMatrix%updateMatrix)

    NULLIFY(equations,dependentField,independentField)
    NULLIFY(dynamicMatrices,stiffnessMatrix,DECOMPOSITION)
    NULLIFY(DECOMP_ELEMENT)
    NULLIFY(GEOMETRIC_INTERPOLATION_PARAMETERS,GEOMETRIC_INTERPOLATED_POINT)
    NULLIFY(DECOMP_FACE,DOMAIN_FACE)
    NULLIFY(FACE_BASIS,DEPENDENT_BASIS,FACE_QUADRATURE_SCHEME,FACE_QUADRATURE_SCHEME)
    NULLIFY(FACE_VELOCITY_INTERPOLATION_PARAMETERS,FACE_INTERPOLATED_POINT)

    PENALTY_PARAM = 1.0e04_DP

    !Grab pointers of interest
    EQUATIONS=>EQUATIONS_SET%EQUATIONS
    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    dynamicMatrices=>vectorEquations%vectorMatrices%dynamicMatrices
    stiffnessMatrix=>dynamicMatrices%matrices(1)%ptr
    dependentField=>equations%interpolation%dependentField
    DECOMPOSITION  =>dependentField%DECOMPOSITION
    DECOMP_ELEMENT=>DECOMPOSITION%TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)

    independentField=>equations%interpolation%independentField

!     MESH_COMPONENT_NUMBER=DECOMPOSITION%MESH_COMPONENT_NUMBER
    MESH_COMPONENT_NUMBER = vectorEquations%vectorMapping%dynamicMapping%equationsMatrixToVarMaps(1)% &
      & VARIABLE%COMPONENTS(1)%MESH_COMPONENT_NUMBER

    DEPENDENT_BASIS=>DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%ptr%TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS

!     write(*,*)'ELEMENT_NUMBER = ',ELEMENT_NUMBER

    !Calculate penalty term to render surfaces impermeable: Loop over all faces
    DO element_face_idx=1,DEPENDENT_BASIS%NUMBER_OF_LOCAL_FACES
      face_number=DECOMP_ELEMENT%ELEMENT_FACES(element_face_idx)
      DECOMP_FACE=>DECOMPOSITION%TOPOLOGY%FACES%FACES(face_number)

      !Check if it's a boundary face
      IF(DECOMP_FACE%BOUNDARY_FACE) THEN !!temporary until MESH_FACE (or equivalent) is available (decomp face includes ghost faces?)

        !Grab normal xi direction of the face and the other two xi directions
        normal_component_idx=ABS(DECOMP_FACE%XI_DIRECTION)  ! if xi=0, this can be a negative number
!         FACE_COMPONENTS=OTHER_XI_DIRECTIONS3(normal_component_idx,2:3,1)  !Two xi directions for the current face
        !\todo: will FACE_COMPONENTS be a problem with sector elements? Check this.

        ! To find out which faces are set impermeable:
        FACE_VELOCITY_INTERPOLATION_PARAMETERS=>equations%interpolation%independentInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr
        CALL FIELD_INTERPOLATION_PARAMETERS_FACE_GET(FIELD_VALUES_SET_TYPE,face_number, &
          & FACE_VELOCITY_INTERPOLATION_PARAMETERS,err,error,*999)
        FACE_INTERPOLATED_POINT=>equations%interpolation%independentInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr


        !Check if impermeable boundary condition is defined on the face
        IMPERMEABLE_BC=.FALSE.
        IF(ANY(ABS(FACE_VELOCITY_INTERPOLATION_PARAMETERS%PARAMETERS(:,normal_component_idx))>ZERO_TOLERANCE)) THEN
          IMPERMEABLE_BC=.TRUE.
        ENDIF

        IF(IMPERMEABLE_BC) THEN

!           write(*,*)'element_face_idx = ',element_face_idx
!           write(*,*)'DECOMP_FACE%XI_DIRECTION = ',DECOMP_FACE%XI_DIRECTION

          !Grab some other pointers
          DOMAIN_FACE=>DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%ptr%TOPOLOGY%FACES%FACES(face_number)
          FACE_BASIS=>DOMAIN_FACE%BASIS
          FACE_QUADRATURE_SCHEME=>FACE_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
          FACE_NUMBER_OF_GAUSS_POINTS=FACE_QUADRATURE_SCHEME%NUMBER_OF_GAUSS

          !A single FACE_BASIS and DEPENDENT_BASIS should suffice, since we only deal with terms
          !  deriving from velocity test AND trial functions, and moreover use Galerkin,
          !  i.e. same basis functions for test and trial functions

          !Start integrating
!\todo: hopefully all quadrature stuff will always match up between face basis and local face stuff.
! Annoying issue here that p(appl) is interpolated using the face_basis, while dZdXI has to be evaluated
! using the 3D face interpolation... many variables are shared, probably supposed to be the same but I
! can't guarantee it and checking every single thing will be a fair bit of overhead
          DO gauss_idx=1,FACE_NUMBER_OF_GAUSS_POINTS
            GAUSS_WEIGHT=FACE_QUADRATURE_SCHEME%GAUSS_WEIGHTS(gauss_idx)
            !What happens with surface Jacobian ? SQRT_G ? - Apparently contained in normal calculation

            !Use (deformed) Geometric field to obtain delx_j/delxi_M = dZdxi at the face gauss point
            GEOMETRIC_INTERPOLATION_PARAMETERS=>equations%interpolation%geometricInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER, &
              & GEOMETRIC_INTERPOLATION_PARAMETERS,err,error,*999)
            GEOMETRIC_INTERPOLATED_POINT=>equations%interpolation%geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr
            CALL FIELD_INTERPOLATE_LOCAL_FACE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,element_face_idx,gauss_idx, &
              & GEOMETRIC_INTERPOLATED_POINT,err,error,*999)

            DZDXI=GEOMETRIC_INTERPOLATED_POINT%VALUES(1:3,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(1:3)) !(component,derivative)

!             write(*,*)'gauss_idx = ',gauss_idx
!             write(*,*)'GAUSS_COORDS = ',GEOMETRIC_INTERPOLATED_POINT%VALUES(1:3,NO_PART_DERIV) !(component,derivative)

            !Calculate covariant metric tensor
            CALL MatrixTranspose(DZDXI,DZDXIT,err,error,*999)
            CALL MatrixProduct(DZDXIT,DZDXI,GIJL,err,error,*999) !g_ij = dZdXI' * dZdXI
            CALL Invert(GIJL,GIJU,G,err,error,*999) !g^ij = inv(g_ij), G=DET(GIJL)
            SQRT_G=SQRT(G)

            !--- L o o p   1 : over element rows (3 velocity components) -----------------------------------
            DO component_idx_1=1,3
              !Calculate g^{normal_component_idx}M*dZ_j/dxi_M; this apparently includes the face Jacobian
              CALL DotProduct(GIJU(normal_component_idx,:),DZDXI(component_idx_1,:),NORMAL_PROJECTION_1,err,error,*999)

              IF(DECOMP_FACE%XI_DIRECTION<0) NORMAL_PROJECTION_1=-NORMAL_PROJECTION_1  !always outward normal

              IF(ABS(NORMAL_PROJECTION_1)<ZERO_TOLERANCE) CYCLE !Makes it a bit quicker

              element_base_dof_idx_1 = (component_idx_1-1) * DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS

              DO face_node_idx_1=1,FACE_BASIS%NUMBER_OF_NODES !nnf
                element_node_idx_1=DEPENDENT_BASIS%NODE_NUMBERS_IN_LOCAL_FACE(face_node_idx_1,element_face_idx) !nn

                DO face_node_derivative_idx_1=1,FACE_BASIS%NUMBER_OF_DERIVATIVES(face_node_idx_1) !nkf

                  element_node_derivative_idx_1=DEPENDENT_BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_FACE(face_node_derivative_idx_1, &
                    & face_node_idx_1,element_face_idx)

                  parameter_idx_1=DEPENDENT_BASIS%ELEMENT_PARAMETER_INDEX(element_node_derivative_idx_1,element_node_idx_1)

                  face_parameter_idx_1=FACE_BASIS%ELEMENT_PARAMETER_INDEX(face_node_derivative_idx_1,face_node_idx_1)

                  element_dof_idx_1=element_base_dof_idx_1+parameter_idx_1

                  !--- L o o p   2 : over element columns (3 velocity components) -----------------------------------
                  DO component_idx_2=1,3
                    !Calculate g^3M*dZ_j/dxi_M
                    CALL DotProduct(GIJU(normal_component_idx,:),DZDXI(component_idx_2,:),NORMAL_PROJECTION_2,err,error,*999)

                    IF(DECOMP_FACE%XI_DIRECTION<0) NORMAL_PROJECTION_2=-NORMAL_PROJECTION_2  !always outward normal

                    IF(ABS(NORMAL_PROJECTION_2)<ZERO_TOLERANCE) CYCLE !Makes it a bit quicker

                    element_base_dof_idx_2 = (component_idx_2-1) * DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS

                    DO face_node_idx_2=1,FACE_BASIS%NUMBER_OF_NODES !nnf
                      element_node_idx_2=DEPENDENT_BASIS%NODE_NUMBERS_IN_LOCAL_FACE(face_node_idx_2,element_face_idx) !nn

                      DO face_node_derivative_idx_2=1,FACE_BASIS%NUMBER_OF_DERIVATIVES(face_node_idx_2) !nkf

                        element_node_derivative_idx_2=DEPENDENT_BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_FACE(face_node_derivative_idx_2, &
                          & face_node_idx_2, element_face_idx)

                        parameter_idx_2=DEPENDENT_BASIS%ELEMENT_PARAMETER_INDEX(element_node_derivative_idx_2,element_node_idx_2)

                        face_parameter_idx_2=FACE_BASIS%ELEMENT_PARAMETER_INDEX(face_node_derivative_idx_2,face_node_idx_2)

                        element_dof_idx_2=element_base_dof_idx_2+parameter_idx_2

                        SUM = 0.0_DP

                        PGM=FACE_QUADRATURE_SCHEME%GAUSS_BASIS_FNS(face_parameter_idx_1,NO_PART_DERIV,gauss_idx)
                        PGN=FACE_QUADRATURE_SCHEME%GAUSS_BASIS_FNS(face_parameter_idx_2,NO_PART_DERIV,gauss_idx)

                        SUM = SUM + PENALTY_PARAM * PGM * NORMAL_PROJECTION_1 * SQRT_G * &
                                                  & PGN * NORMAL_PROJECTION_2 * SQRT_G

                        stiffnessMatrix%elementMatrix%matrix(element_dof_idx_1,element_dof_idx_2) = &
                          & stiffnessMatrix%elementMatrix%matrix(element_dof_idx_1,element_dof_idx_2) + &
                          & SUM * GAUSS_WEIGHT

                      ENDDO !element_node_derivative_idx_2
                    ENDDO !face_node_idx_2
                  ENDDO !component_idx_2

                ENDDO !element_node_derivative_idx_1
              ENDDO !face_node_idx_1

!               write(*,*)'component_idx_1 = ',component_idx_1
!               write(*,*)'NORMAL_PROJECTION_1 = ',NORMAL_PROJECTION_1
!               write(*,*)' '

            ENDDO !component_idx_1

          ENDDO !gauss_idx
        ENDIF !IMPERMEABLE_BC
      ENDIF !boundary face check
    ENDDO !element_face_idx

    EXITS("DARCY_EQUATION_IMPERMEABLE_BC_VIA_PENALTY")
    RETURN

999 ERRORSEXITS("DARCY_EQUATION_IMPERMEABLE_BC_VIA_PENALTY",err,error)
    RETURN 1
  END SUBROUTINE DARCY_EQUATION_IMPERMEABLE_BC_VIA_PENALTY

  !
  !================================================================================================================================
  !


END MODULE DARCY_EQUATIONS_ROUTINES

