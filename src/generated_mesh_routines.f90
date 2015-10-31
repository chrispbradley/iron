!> \file
!> \author Chris Bradley
!> \brief This module handles all generated mesh routines.
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
!> Contributor(s): Chris Bradley, Jack Lee
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

!> This module handles all generated mesh routines.
MODULE GeneratedMeshRoutines

  USE BASE_ROUTINES
  USE BASIS_ROUTINES
  USE COMP_ENVIRONMENT
  USE Constants
  USE COORDINATE_ROUTINES
  USE FIELD_ROUTINES
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE Kinds
  USE Maths
  USE MESH_ROUTINES
  USE NODE_ROUTINES
  USE Strings
  USE Types

#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !> \addtogroup GeneratedMeshRoutines_GeneratedMeshTypes GeneratedMeshRoutines::GeneratedMeshTypes
  !> \brief Generated mesh types.
  !> \see GeneratedMeshRoutines,OPENCMISS_GeneratedMeshTypes
  !>@{
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_REGULAR_MESH_TYPE=1 !<A regular generated mesh. \see GeneratedMeshRoutines_GeneratedMeshTypes,GeneratedMeshRoutines
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_POLAR_MESH_TYPE=2 !<A polar generated mesh. \see GeneratedMeshRoutines_GeneratedMeshTypes,GeneratedMeshRoutines
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_CYLINDER_MESH_TYPE=3 !<A cylinder generated mesh. \see GeneratedMeshRoutines_GeneratedMeshTypes,GeneratedMeshRoutines
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_ELLIPSOID_MESH_TYPE=4 !<An ellipsoid generated mesh. \see GeneratedMeshRoutines_GeneratedMeshTypes,GeneratedMeshRoutines
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_FRACTAL_TREE_MESH_TYPE=5 !<A fractal tree generated mesh. \see GeneratedMeshRoutines_GeneratedMeshTypes,GeneratedMeshRoutines
  !>@}

  !> \addtogroup GeneratedMeshRoutines_GeneratedMeshCylinderSurfaces GeneratedMeshRoutines::GeneratedMeshCylinderSurfaces
  !> \brief Generated mesh cylinder type surface types.
  !>@{
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_CYLINDER_INNER_SURFACE=1  !<Inner surface of the cylinder. \see GeneratedMeshRoutines_GeneratedMeshCylinderSurfaces,GeneratedMeshRoutines
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_CYLINDER_OUTER_SURFACE=2  !<Outer surface of the cylinder. \see GeneratedMeshRoutines_GeneratedMeshCylinderSurfaces,GeneratedMeshRoutines
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_CYLINDER_TOP_SURFACE=3    !<Top surface of the cylinder. \see GeneratedMeshRoutines_GeneratedMeshCylinderSurfaces,GeneratedMeshRoutines
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_CYLINDER_BOTTOM_SURFACE=4 !<Bottom surface of the cylinder. \see GeneratedMeshRoutines_GeneratedMeshCylinderSurfaces,GeneratedMeshRoutines
  !>@}

  !> \addtogroup GeneratedMeshRoutines_GeneratedMeshEllipsoidSurfaces GeneratedMeshRoutines::GeneratedMeshEllipsoidSurfaces
  !> \brief Generated mesh ellipsoid type surface types.
  !>@{
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_ELLIPSOID_INNER_SURFACE=5  !<Inner surface of the ellipsoid. \see GeneratedMeshRoutines_GeneratedMeshEllipsoidSurfaces,GeneratedMeshRoutines
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_ELLIPSOID_OUTER_SURFACE=6  !<Outer surface of the ellipsoid. \see GeneratedMeshRoutines_GeneratedMeshEllipsoidSurfaces,GeneratedMeshRoutines
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_ELLIPSOID_TOP_SURFACE=7    !<Top surface of the ellipsoid. \see GeneratedMeshRoutines_GeneratedMeshEllipsoidSurfaces,GeneratedMeshRoutines
  !>@}

  !> \addtogroup GeneratedMeshRoutines_GeneratedMeshRegularSurfaces GeneratedMeshRoutines::GeneratedMeshRegularSurfaces
  !> \brief Generated mesh regular type surface types.
  !>@{
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_REGULAR_LEFT_SURFACE=8    !<Left surface of the regular mesh. \see GeneratedMeshRoutines_GeneratedMeshRegularSurfaces,GeneratedMeshRoutines
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_REGULAR_RIGHT_SURFACE=9   !<Right surface of the regular mesh. \see GeneratedMeshRoutines_GeneratedMeshRegularSurfaces,GeneratedMeshRoutines
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_REGULAR_TOP_SURFACE=10    !<Top surface of the regular mesh. \see GeneratedMeshRoutines_GeneratedMeshRegularSurfaces,GeneratedMeshRoutines
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_REGULAR_BOTTOM_SURFACE=11 !<Bottom surface of the regular mesh. \see GeneratedMeshRoutines_GeneratedMeshRegularSurfaces,GeneratedMeshRoutines
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_REGULAR_FRONT_SURFACE=12  !<Front surface of the regular mesh. \see GeneratedMeshRoutines_GeneratedMeshRegularSurfaces,GeneratedMeshRoutines
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_REGULAR_BACK_SURFACE=13   !<Back surface of the regular mesh. \see GeneratedMeshRoutines_GeneratedMeshRegularSurfaces,GeneratedMeshRoutines
  !>@}

  !Module types

  !Interfaces

  !>Starts the process of creating a generated mesh
  INTERFACE GeneratedMesh_CreateStart
    MODULE PROCEDURE GeneratedMesh_CreateStartInterface
    MODULE PROCEDURE GeneratedMesh_CreateStartRegion
  END INTERFACE GeneratedMesh_CreateStart

  !>Initialises the generated meshes for a region or interface.
  INTERFACE GeneratedMeshes_Initialise
    MODULE PROCEDURE GeneratedMeshes_InitialiseInterface
    MODULE PROCEDURE GeneratedMeshes_InitialiseRegion 
  END INTERFACE GeneratedMeshesInitialise

  !>Finds a generated mesh in a list of generated meshes in a region or interface.
  INTERFACE GeneratedMesh_UserNumberFind
    MODULE PROCEDURE GeneratedMesh_UserNumberFindInterface
    MODULE PROCEDURE GeneratedMesh_UserNumberFindRegion 
  END INTERFACE GeneratedMesh_UserNumberFind

  PUBLIC GENERATED_MESH_REGULAR_MESH_TYPE,GENERATED_MESH_POLAR_MESH_TYPE,GENERATED_MESH_CYLINDER_MESH_TYPE, &
    & GENERATED_MESH_ELLIPSOID_MESH_TYPEGENERATED_MESH_FRACTAL_TREE_MESH_TYPE
  
  PUBLIC GENERATED_MESH_REGULAR_LEFT_SURFACE,GENERATED_MESH_REGULAR_RIGHT_SURFACE
  PUBLIC GENERATED_MESH_REGULAR_TOP_SURFACE,GENERATED_MESH_REGULAR_BOTTOM_SURFACE
  PUBLIC GENERATED_MESH_REGULAR_FRONT_SURFACE,GENERATED_MESH_REGULAR_BACK_SURFACE

  PUBLIC GENERATED_MESH_CYLINDER_INNER_SURFACE,GENERATED_MESH_CYLINDER_OUTER_SURFACE
  PUBLIC GENERATED_MESH_CYLINDER_TOP_SURFACE,GENERATED_MESH_CYLINDER_BOTTOM_SURFACE

  PUBLIC GENERATED_MESH_ELLIPSOID_INNER_SURFACE,GENERATED_MESH_ELLIPSOID_OUTER_SURFACE
  PUBLIC GENERATED_MESH_ELLIPSOID_TOP_SURFACE
  
  PUBLIC GeneratedMeshes_Initialise,GeneratedMeshes_Finalise

  PUBLIC GeneratedMesh_BaseVectorsSet

  PUBLIC GeneratedMesh_CoordinateSystemGet

  PUBLIC GeneratedMesh_CreateStart,GeneratedMesh_CreateFinish

  PUBLIC GeneratedMesh_Destroy

  PUBLIC GeneratedMesh_BasisSet,GeneratedMesh_BasisGet

  PUBLIC GeneratedMesh_NumberOfBasisGet

  PUBLIC GeneratedMesh_ExtentSet,GeneratedMesh_ExtentGet

  PUBLIC GeneratedMesh_NumberOfElementsSet,GeneratedMesh_NumberOfElementsGet

  PUBLIC GeneratedMesh_OriginSet,GeneratedMesh_OriginGet

  PUBLIC GeneratedMesh_TypeSet,GeneratedMesh_TypeGet

  PUBLIC GeneratedMesh_GeometricParametersCalculate

  PUBLIC GeneratedMesh_RegionGet

  PUBLIC GeneratedMesh_UserNumberFind
  
  PUBLIC GeneratedMesh_SurfaceGet

CONTAINS

  !
  !================================================================================================================================
  !

  !>Gets the basis of a generated mesh. \see OPENCMISS::CMISSGeneratedMeshBasisGet
  SUBROUTINE GeneratedMesh_BasisGet(generatedMesh,bases,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh to get the bases of
    TYPE(BASIS_PTR_TYPE) :: bases(:) !<bases(basisIdx). On return, the bases of mesh to generate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
<<<<<<< HEAD
    INTEGER(INTG) :: basisIdx
    TYPE(VARYING_STRING) :: localError

    ENTERS"GeneratedMesh_BasisGet",err,error,*999)

    IF(ASSOCIATED(generatedMesh)) THEN
      IF(generatedMesh%generatedMeshFinished) THEN
        IF(ALLOCATED(generatedMesh%bases)) THEN
          IF(SIZE(bases,1)>=SIZE(generatedMesh%bases,1)) THEN
            DO basisIdx=1,SIZE(generatedMesh%bases,1)
              IF(ASSOCIATED(bases(basisIdx)%ptr)) THEN
                localError="The pointer at location "//TRIM(NumberToVString(basisIdx,"*",err,error))// &
                  & " in the specified bases array is already associated."
                CALL FlagError(localError,err,error,*999)
              ELSE                    
                bases(basisIdx)%ptr=>generatedMesh%regularMesh%bases(basisIdx)%ptr
              ENDIF
            ENDDO !basisIdx
          ELSE
            localError="The size of the specified bases array, "// &
              & TRIM(NumberToVString(SIZE(bases,1),"*",err,error))// &
              & ", is too small to hold the number of bases in the generated mesh of "// &
              & TRIM(NumberToVString(SIZE(generatedMesh%bases,1),"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Generated mesh basis has not been set.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Generated mesh has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Generated mesh is not associated.",err,error,*999)
    ENDIF

    CALL Exits("GeneratedMeshBasisGet")
    RETURN
999 CALL Errors("GeneratedMeshBasisGet",err,error)
    CALL Exits("GeneratedMeshBasisGet")
=======
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    INTEGER(INTG) :: basis_idx,NUM_BASES

    ENTERS("GENERATED_MESH_BASIS_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      IF(GENERATED_MESH%GENERATED_MESH_FINISHED) THEN
        SELECT CASE(GENERATED_MESH%GENERATED_TYPE)
        CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
          IF(ASSOCIATED(GENERATED_MESH%REGULAR_MESH)) THEN
            IF(ALLOCATED(GENERATED_MESH%REGULAR_MESH%BASES)) THEN
              NUM_BASES=SIZE(GENERATED_MESH%REGULAR_MESH%BASES)
              ALLOCATE(BASES(NUM_BASES),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate bases.",ERR,ERROR,*999)
              DO basis_idx=1,NUM_BASES
                BASES(basis_idx)%PTR=>GENERATED_MESH%REGULAR_MESH%BASES(basis_idx)%PTR
              ENDDO
            ELSE
              CALL FlagError("Generated mesh bases are not allocated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Generated mesh regular mesh is not associated.",ERR,ERROR,*999)
          ENDIF
        CASE(GENERATED_MESH_POLAR_MESH_TYPE)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
          IF(ASSOCIATED(GENERATED_MESH%CYLINDER_MESH)) THEN
            IF(ALLOCATED(GENERATED_MESH%REGULAR_MESH%BASES)) THEN
              NUM_BASES=SIZE(GENERATED_MESH%CYLINDER_MESH%BASES)
              ALLOCATE(BASES(NUM_BASES),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate bases.",ERR,ERROR,*999)
              DO basis_idx=1,NUM_BASES
                BASES(basis_idx)%PTR=>GENERATED_MESH%CYLINDER_MESH%BASES(basis_idx)%PTR
              ENDDO
            ELSE
              CALL FlagError("Generated mesh bases are not allocated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Generated mesh cylinder mesh is not associated.",ERR,ERROR,*999)
          ENDIF
        CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
          IF(ASSOCIATED(GENERATED_MESH%ELLIPSOID_MESH)) THEN
            IF(ALLOCATED(GENERATED_MESH%REGULAR_MESH%BASES)) THEN
              NUM_BASES=SIZE(GENERATED_MESH%ELLIPSOID_MESH%BASES)
              ALLOCATE(BASES(NUM_BASES),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate bases.",ERR,ERROR,*999)
              DO basis_idx=1,NUM_BASES
                BASES(basis_idx)%PTR=>GENERATED_MESH%ELLIPSOID_MESH%BASES(basis_idx)%PTR
              ENDDO
            ELSE
              CALL FlagError("Generated mesh bases are not allocated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Generated mesh ellipsoid mesh is not associated.",ERR,ERROR,*999)
          ENDIF
        CASE DEFAULT
          LOCAL_ERROR="The generated mesh generated type of "// &
              & TRIM(NUMBER_TO_VSTRING(GENERATED_MESH%GENERATED_TYPE,"*",ERR,ERROR))//" is invalid."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FlagError("Generated mesh has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Generated mesh is already associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("GENERATED_MESH_BASIS_GET")
    RETURN
999 ERRORSEXITS("GENERATED_MESH_BASIS_GET",ERR,ERROR)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
    RETURN 1
    
  END SUBROUTINE GeneratedMeshBasisGet

  !
  !================================================================================================================================
  !

  !>Gets the number of bases in a generated mesh. \see OPENCMISS::CMISSGeneratedMeshNumberOfBasisGet
  SUBROUTINE GeneratedMeshNumberOfBasisGet(generatedMesh,numberOfBases,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh to get the bases of
    INTEGER(INTG), INTENT(OUT) :: numberOfBases !< On return, the number of bases in the generated mesh
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: basisIdx
    TYPE(VARYING_STRING) :: localError

<<<<<<< HEAD
    CALL Enters("GeneratedMeshNumberOfBasisGet",err,error,*999)

    IF(ASSOCIATED(generatedMesh)) THEN
      IF(generatedMesh%generatedMeshFinished) THEN
        IF(ALLOCATED(generatedMesh%bases)) THEN
          numberOfBases=SIZE(generatedMesh%bases,1)) THEN
        ELSE
          numberOfBases=0
        ENDIF
      ELSE
        CALL FlagError("Generated mesh has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Generated mesh is not associated.",err,error,*999)
    ENDIF

    CALL Exits("GeneratedMeshNumberOfBasisGet")
    RETURN
999 CALL Errors("GeneratedMeshNumberOfBasisGet",err,error)
    CALL Exits("GeneratedMeshNumberOfBasisGet")
=======
    ENTERS("GENERATED_MESH_BASIS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      IF(GENERATED_MESH%GENERATED_MESH_FINISHED) THEN
        CALL FlagError("Generated mesh has been finished.",ERR,ERROR,*999)
      ELSE
        NULLIFY(COORDINATE_SYSTEM)
        CALL GENERATED_MESH_COORDINATE_SYSTEM_GET(GENERATED_MESH,COORDINATE_SYSTEM,ERR,ERROR,*999)
        CALL COORDINATE_SYSTEM_DIMENSION_GET(COORDINATE_SYSTEM,COORDINATE_DIMENSION,ERR,ERROR,*999)
        NUM_BASES=SIZE(BASES)
        NUM_XI=BASES(1)%PTR%NUMBER_OF_XI
        BASIS_TYPE=BASES(1)%PTR%TYPE
        DO basis_idx=2,NUM_BASES
          IF(BASES(basis_idx)%PTR%NUMBER_OF_XI /= NUM_XI) THEN
            CALL FlagError("All bases must have the same number of xi.",ERR,ERROR,*999)
          ENDIF
          IF(BASES(basis_idx)%PTR%TYPE /= BASIS_TYPE) THEN
            CALL FlagError("Using different basis types is not supported for generated meshes.",ERR,ERROR,*999)
          ENDIF
        ENDDO
        SELECT CASE(GENERATED_MESH%GENERATED_TYPE)
        CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
          IF(ASSOCIATED(GENERATED_MESH%REGULAR_MESH)) THEN
            IF(ALLOCATED(GENERATED_MESH%REGULAR_MESH%BASE_VECTORS)) THEN
              CALL FlagError("Can not reset the basis if base vectors have been specified.",ERR,ERROR,*999)
            ELSE
              IF(ALLOCATED(GENERATED_MESH%REGULAR_MESH%BASES)) DEALLOCATE(GENERATED_MESH%REGULAR_MESH%BASES)
              ALLOCATE(GENERATED_MESH%REGULAR_MESH%BASES(NUM_BASES),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate bases.",ERR,ERROR,*999)
              DO basis_idx=1,NUM_BASES
                IF(ASSOCIATED(BASES(basis_idx)%PTR)) THEN
                  IF(BASES(basis_idx)%PTR%NUMBER_OF_XI<=COORDINATE_DIMENSION) THEN
                    GENERATED_MESH%REGULAR_MESH%BASES(basis_idx)%PTR=>BASES(basis_idx)%PTR
                  ELSE
                    LOCAL_ERROR="The basis number of xi dimensions of "// &
                      & TRIM(NUMBER_TO_VSTRING(BASES(basis_idx)%PTR%NUMBER_OF_XI,"*",ERR,ERROR))// &
                      & " is invalid. The number of xi dimensions must be <= the number of coordinate dimensions of "// &
                      & TRIM(NUMBER_TO_VSTRING(COORDINATE_DIMENSION,"*",ERR,ERROR))
                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The basis with index "//TRIM(NUMBER_TO_VSTRING(basis_idx,"*",ERR,ERROR))// &
                    & " is not associated."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ENDDO
            ENDIF
          ELSE
            CALL FlagError("Regular generated mesh is not associated.",ERR,ERROR,*999)
          ENDIF
        CASE(GENERATED_MESH_POLAR_MESH_TYPE)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
          IF(ASSOCIATED(GENERATED_MESH%CYLINDER_MESH)) THEN
            IF(ALLOCATED(GENERATED_MESH%CYLINDER_MESH%BASES)) DEALLOCATE(GENERATED_MESH%CYLINDER_MESH%BASES)
            ALLOCATE(GENERATED_MESH%CYLINDER_MESH%BASES(NUM_BASES),STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate bases.",ERR,ERROR,*999)
            DO basis_idx=1,NUM_BASES
              IF(ASSOCIATED(BASES(basis_idx)%PTR)) THEN
                GENERATED_MESH%CYLINDER_MESH%BASES(basis_idx)%PTR=>BASES(basis_idx)%PTR
              ELSE
                LOCAL_ERROR="The basis with index "//TRIM(NUMBER_TO_VSTRING(basis_idx,"*",ERR,ERROR))// &
                  & " is not associated."
                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ENDDO
          ELSE
            CALL FlagError("Cylinder generated mesh is not associated.",ERR,ERROR,*999)
          ENDIF
        CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
          IF(ASSOCIATED(GENERATED_MESH%ELLIPSOID_MESH)) THEN
            IF(ALLOCATED(GENERATED_MESH%ELLIPSOID_MESH%BASES)) DEALLOCATE(GENERATED_MESH%ELLIPSOID_MESH%BASES)
            ALLOCATE(GENERATED_MESH%ELLIPSOID_MESH%BASES(NUM_BASES),STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate bases.",ERR,ERROR,*999)
            DO basis_idx=1,NUM_BASES
              IF(ASSOCIATED(BASES(basis_idx)%PTR)) THEN
                GENERATED_MESH%ELLIPSOID_MESH%BASES(basis_idx)%PTR=>BASES(basis_idx)%PTR
              ELSE
                LOCAL_ERROR="The basis with index "//TRIM(NUMBER_TO_VSTRING(basis_idx,"*",ERR,ERROR))// &
                  & " is not associated."
                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ENDDO
          ELSE
            CALL FlagError("Ellpsoid generated mesh is not associated.",ERR,ERROR,*999)
          ENDIF
        CASE DEFAULT
          LOCAL_ERROR="The generated mesh type of "//TRIM(NUMBER_TO_VSTRING(GENERATED_MESH%GENERATED_TYPE,"*",ERR,ERROR))// &
            & " is invalid."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FlagError("Generated mesh is already associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("GENERATED_MESH_BASIS_SET")
    RETURN
999 ERRORSEXITS("GENERATED_MESH_BASIS_SET",ERR,ERROR)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
    RETURN 1
    
  END SUBROUTINE GeneratedMeshNumberOfBasisGet

  !
  !================================================================================================================================
  !

  !>Sets/changes the basis of a generated mesh. \see OPENCMISS::CMISSGeneratedMeshBasisSet
  SUBROUTINE GeneratedMeshBasisSet(generatedMesh,bases,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh to set the basis of
    TYPE(BASIS_PTR_TYPE) :: bases(:) !<bases(basisIdx). An array of pointers to the bases to generate the mesh with
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: basisIdx,basisType,coordinateDimension,coordinateIdx,count,numberOfBases,numberOfXi,xiIdx
    INTEGER(INTG), ALLOCATABLE :: newNumberOfElementsXi(:)
    REAL(DP), ALLOCATABLE :: newBaseVectors(:,:)
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: coordinateSystem
    TYPE(VARYING_STRING) :: localError

    NULLIFY(coordinateSystem)

<<<<<<< HEAD
    CALL Enters("GeneratedMeshBasisSet",err,error,*999)

    IF(ASSOCIATED(generatedMesh)) THEN
      IF(generatedMesh%generatedMeshFinished) THEN
        CALL FlagError("Generated mesh has been finished.",err,error,*999)
=======
    ENTERS("GENERATED_MESH_BASE_VECTORS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      IF(GENERATED_MESH%GENERATED_MESH_FINISHED) THEN
        CALL FlagError("Generated mesh has been finished.",ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
      ELSE
        numberOfBases=SIZE(bases,1)
        IF(numberOfBases >= 1) THEN
          !Check the supplied bases
          DO basisIdx=1,numberOfBases
            IF(ASSOCIATED(bases(basisIdx)%ptr)) THEN
              IF(bases(basisIdx)%ptr%NUMBER_OF_XI<=generatedMesh%coordinateDimension) THEN
                IF(basisIdx == 1) THEN
                  numberOfXi=bases(1)%ptr%NUMBER_OF_XI
                  basisType=bases(1)%ptr%type
                ELSE                  
                  IF(bases(basisIdx)%ptr%NUMBER_OF_XI /= numberOfXi) THEN
                    CALL FlagError("All bases must have the same number of xi directions.",err,error,*999)
                  ENDIF
                  IF(bases(basisIdx)%ptr%type /= basisType) THEN
                    CALL FlagError("Using different basis types is not supported for generated meshes.",err,error,*999)
                  ENDIF
                ENDIF
              ELSE
                localError="The basis number of xi dimensions of "// &
                  & TRIM(NumberToVString(bases(basisIdx)%ptr%NUMBER_OF_XI,"*",err,error))// &
                  & " is invalid. The number of xi dimensions must be <= the number of coordinate dimensions of "// &
                  & TRIM(NumberToVString(coordinateDimension,"*",err,error))
                CALL FlagError(localError,err,error,*999)
              ENDIF
            ELSE
              localError="The basis with index "//TRIM(NumberToVString(basisIdx,"*",err,error))// &
                & " is not associated."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ENDDO !basisIdx
          newMeshDimension=numberOfXi
          SELECT CASE(generatedMesh%generatedType)
          CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
<<<<<<< HEAD
            !Check the bases are not collapsed
            DO basisIdx=1,numberOfBases
              IF(.NOT.ALL(bases(basisIdx)%ptr%COLLAPSED_XI==BASIS_NOT_COLLAPSED)) THEN
                localError="Specified basis number "//TRIM(NumberToVString(basisIdx,"*",err,error))// &
                  & " is invalid. The bases must not be collapsed for a regular mesh."
                CALL FlagError(localError,err,error,*999)
              ENDIF
            ENDDO !basisIdx
            !All bases are OK
            IF(ASSOCIATED(generatedMesh%regularMesh)) THEN
              !Reset the number of elements in each xi direction
              ALLOCATE(newNumberOfElementsXi(newMeshDimension),STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate new number of elements xi.",err,error,*999)
              IF(generatedMesh%meshDimension==0) THEN
                !First time, default attributes
                newNumberOfElementsXi(1:newMeshDimension)=1
              ELSE IF(generatedMesh%meshDimension>newMeshDimension) THEN
                !New mesh dimension is less than the old mesh dimension
                newNumberOfElementsXi(1:numberOfXi)=generatedMesh%regularMesh%numberOfElementsXi(1:numberOfXi)
              ELSE
                !New mesh dimension is more than the old mesh dimension
                newNumberOfElementsXi(1:generatedMesh%meshDimension)= &
                  & generatedMesh%regularMesh%numberOfElementsXi(1:generatedMesh%meshDimension)
                newNumberOfElementsXi(generatedMesh%meshDimension+1:newMeshDimension)=1
              ENDIF
              CALL MOVE_ALLOC(newNumberOfElementsXi,generatedMesh%regularMesh%numberOfElementsXi)
            ELSE
              CALL FlagError("Regular generated mesh is not associated.",err,error,*999)
            ENDIF
          CASE(GENERATED_MESH_POLAR_MESH_TYPE)
            IF(ASSOCIATED(generatedMesh%polarMesh)) THEN
              !Check supplied bases
              IF(newMeshDimension==1) THEN
                !With one xi direction the mesh must be a circle.
                spherical=.FALSE.
              ELSE
                IF(newMeshDimension==2) THEN
                  !Need to determine if the mesh is circular of spherical. Check if there are any collapsed bases.
                  spherical=.FALSE.
                  DO basisIdx=1,numberOfBases
                    IF((bases(basisIdx)%ptr%COLLAPSED_XI(1)/=BASIS_NOT_COLLAPSED).OR.&
                      & (bases(basisIdx)%ptr%COLLAPSED_XI(2)/=BASIS_NOT_COLLAPSED)) THEN
                      spherical=.TRUE.
                      EXIT
                    ENDIF
                  ENDDO !basisIdx
                ELSE
                  spherical=.TRUE.                  
                ENDIF
              ENDIF
              IF(spherical) THEN
                !Mesh is spherical shells rather than circular
                IF(MOD(numberOfBases,3)/=0) THEN
                  localError="The number of supplied bases of "//TRIM(NumberToVString(numberOfBases,"*",err,error))// &
                    & " is invalid. The number of bases must be divisiable by 3 for wall, lower and upper apex bases."
                  CALL FlagError(localError,err,error,*999)
                ELSE
                  numberOfBasisSets=numberOfBases/3
                  !Check the collapse of the bases
                  DO basisSetIdx=1,numberOfBasisSets
                    !Wall basis
                    basisIdx=(basisSetIdx-1)*3+1
                    IF((bases(basisIdx)%ptr%COLLAPSED_XI(1)/=BASIS_NOT_COLLAPSED).OR.&
                      & (bases(basisIdx)%ptr%COLLAPSED_XI(2)/=BASIS_NOT_COLLAPSED)) THEN
                      localError=="Specified basis number "//TRIM(NumberToVString(basisIdx,"*",err,error))// &
                        & " is invalid. The wall bases must not be collapsed in a polar mesh."
                      CALL FlagError(localError,err,error,*999)
                    ENDIF
                    !Bottom apex basis
                    basisIdx=(basisSetIdx-1)*3+2
                    IF((bases(basisIdx)%ptr%COLLAPSED_XI(1)/=BASIS_XI_COLLAPSED).OR.&
                      & (bases(basisIdx)%ptr%COLLAPSED_XI(2)/=BASIS_COLLAPSED_AT_XI0)) THEN
                      localError="Specified basis number "//TRIM(NumberToVString(basisIdx,"*",err,error))// &
                        & " is invalid. The top apex bases must be collapsed correctly in a polar mesh."
                      CALL FlagError(localError,err,error,*999)
                    ENDIF
                    !Top apex basis
                    basisIdx=(basisSetIdx-1)*3+3
                    IF((bases(basisIdx)%ptr%COLLAPSED_XI(1)/=BASIS_XI_COLLAPSED).OR.&
                      & (bases(basisIdx)%ptr%COLLAPSED_XI(2)/=BASIS_COLLAPSED_AT_XI1)) THEN
                      localError="Specified basis number "//TRIM(NumberToVString(basisIdx,"*",err,error))// &
                        & " is invalid. The bottom apex bases must be collapsed correctly in a polar mesh."
                      CALL FlagError(localError,err,error,*999)
                    ENDIF
                  ENDDO !basisSetIdx
                  IF(newMeshDimension==3) THEN
                    DO basisIdx=1,numberOfBases
                      IF(bases(basisIdx)%ptr%COLLAPSED_XI(3)/=BASIS_NOT_COLLAPSED) THEN
                        localError="Specified basis number "//TRIM(NumberToVString(basisIdx,"*",err,error))// &
                          & " is invalid. The bases must not be collapsed in the radial direction in a polar mesh."
                        CALL FlagError(localError,err,error,*999)
                      ENDIF
                    ENDDO !basisIdx
                  ENDIF
                ENDIF
                numberOfMeshComponents=numberOfBasisSets
              ELSE
                numberOfMeshComponents=numberOfBases                  
              ENDIF
            ELSE
              CALL FlagError("Polar generated mesh is not associated.",err,error,*999)
            ENDIF
            !All bases are OK
            !Reset the number of elements in each xi direction
            ALLOCATE(newNumberOfElementsXi(newMeshDimension),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate new number of elements xi.",err,error,*999)
            IF(generatedMesh%meshDimension==0) THEN
              !First time, default attributes
              newNumberOfElementsXi(1:newMeshDimension)=1
            ELSE IF(generatedMesh%meshDimension>newMeshDimension) THEN
              !New mesh dimension is less than the old mesh dimension
              newNumberOfElementsXi(1:newMeshDimension)=generatedMesh%regularMesh%numberOfElementsXi(1:newMeshDimension)
            ELSE
              !New mesh dimension is more than the old mesh dimension
              newNumberOfElementsXi(1:generatedMesh%meshDimension)= &
                & generatedMesh%polarMesh%numberOfElementsXi(1:generatedMesh%meshDimension)
              newNumberOfElementsXi(generatedMesh%meshDimension+1:newMeshDimension)=1
            ENDIF
            CALL MOVE_ALLOC(newNumberOfElementsXi,generatedMesh%polarMesh%numberOfElementsXi)
            !Set the spherical nature of the mesh.
            generatedMesh%polarMesh%spherical=spherical
          CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
            !Check supplied bases
            IF(newMeshDimension<2) THEN
              CALL FlagError("The specified bases must have at least two xi directions for a cylindrical mesh.",err,error,*999)
            ELSE                
              IF(ASSOCIATED(generatedMesh%cylindricalMesh)) THEN
                !Need to determine if the mesh is open or closed. Check if there are any collapsed bases.
                closed=.FALSE.
                DO basisIdx=1,numberOfBases
                  IF((bases(basisIdx)%ptr%COLLAPSED_XI(1)/=BASIS_NOT_COLLAPSED).OR.&
                    & (bases(basisIdx)%ptr%COLLAPSED_XI(2)/=BASIS_NOT_COLLAPSED)) THEN
                    closed=.TRUE.
                    EXIT
                  ENDIF
                ENDDO !basisIdx
                IF(closed) THEN
                  !Mesh is a closed cylinder
                  IF(MOD(numberOfBases,3)/=0) THEN
                    localError="The number of supplied bases of "//TRIM(NumberToVString(numberOfBases,"*",err,error))// &
                      & " is invalid. The number of bases must be divisiable by 3 for wall, lower and upper apex bases."
                    CALL FlagError(localError,err,error,*999)
                  ELSE
                    numberOfBasisSets=numberOfBases/3
                    !Check the collapse of the bases
                    DO basisSetIdx=1,numberOfBasisSets
                      !Wall basis
                      basisIdx=(basisSetIdx-1)*3+1
                      IF((bases(basisIdx)%ptr%COLLAPSED_XI(1)/=BASIS_NOT_COLLAPSED).OR.&
                        & (bases(basisIdx)%ptr%COLLAPSED_XI(2)/=BASIS_NOT_COLLAPSED)) THEN
                        localError=="Specified basis number "//TRIM(NumberToVString(basisIdx,"*",err,error))// &
                          & " is invalid. The wall bases must not be collapsed in a cylindrical mesh."
                        CALL FlagError(localError,err,error,*999)
                      ENDIF
                      !Bottom apex basis
                      basisIdx=(basisSetIdx-1)*3+2
                      IF((bases(basisIdx)%ptr%COLLAPSED_XI(1)/=BASIS_XI_COLLAPSED).OR.&
                        & (bases(basisIdx)%ptr%COLLAPSED_XI(2)/=BASIS_COLLAPSED_AT_XI0)) THEN
                        localError="Specified basis number "//TRIM(NumberToVString(basisIdx,"*",err,error))// &
                          & " is invalid. The bottom apex bases must be collapsed correctly in a cylindrical mesh."
                        CALL FlagError(localError,err,error,*999)
                      ENDIF
                      !Top apex basis
                      basisIdx=(basisSetIdx-1)*3+3
                      IF((bases(basisIdx)%ptr%COLLAPSED_XI(1)/=BASIS_XI_COLLAPSED).OR.&
                        & (bases(basisIdx)%ptr%COLLAPSED_XI(2)/=BASIS_COLLAPSED_AT_XI1)) THEN
                        localError="Specified basis number "//TRIM(NumberToVString(basisIdx,"*",err,error))// &
                          & " is invalid. The top apex bases must be collapsed correctly in a cylindrical mesh."
                        CALL FlagError(localError,err,error,*999)
                      ENDIF
                    ENDDO !basisSetIdx
                    IF(newMeshDimension==3) THEN
                      DO basisIdx=1,numberOfBases
                        IF(bases(basisIdx)%ptr%COLLAPSED_XI(3)/=BASIS_NOT_COLLAPSED) THEN
                          localError="Specified basis number "//TRIM(NumberToVString(basisIdx,"*",err,error))// &
                            & " is invalid. The bases must not be collapsed in the radial direction in a cylindrical mesh."
                          CALL FlagError(localError,err,error,*999)
                        ENDIF
                      ENDDO !basisIdx
                    ENDIF
=======
            IF(ASSOCIATED(GENERATED_MESH%REGULAR_MESH)) THEN
              BASES=>GENERATED_MESH%REGULAR_MESH%BASES
              IF(ASSOCIATED(BASES)) THEN
                BASIS=>BASES(1)%PTR !Bases should all have same number of xi
                IF(ASSOCIATED(BASIS)) THEN
                  IF(SIZE(BASE_VECTORS,2)==BASIS%NUMBER_OF_XI) THEN
                    IF(ALLOCATED(GENERATED_MESH%REGULAR_MESH%BASE_VECTORS)) DEALLOCATE(GENERATED_MESH%REGULAR_MESH%BASE_VECTORS)
                    ALLOCATE(GENERATED_MESH%REGULAR_MESH%BASE_VECTORS(SIZE(BASE_VECTORS,1),SIZE(BASE_VECTORS,2)),STAT=ERR)
                    IF(ERR/=0) CALL FlagError("Could not allocate base vectors.",ERR,ERROR,*999)
                    GENERATED_MESH%REGULAR_MESH%BASE_VECTORS=BASE_VECTORS
                  ELSE
                    LOCAL_ERROR="The size of the second dimension of base vectors of "// &
                      & TRIM(NUMBER_TO_VSTRING(SIZE(BASE_VECTORS,2),"*",ERR,ERROR))// &
                      & " is invalid. The second dimension size must match the number of mesh dimensions of "// &
                      & TRIM(NUMBER_TO_VSTRING(BASIS%NUMBER_OF_XI,"*",ERR,ERROR))//"."
                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
                  ENDIF
                  numberOfMeshComponents=numberOfBasisSets
                ELSE
<<<<<<< HEAD
                  numberOfMeshComponents=numberOfBases                  
                ENDIF
              ELSE
                CALL FlagError("Cylindrical generated mesh is not associated.",err,error,*999)                
=======
                  CALL FlagError("Bases are not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("You must set the generated mesh basis before setting base vectors.",ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
              ENDIF
              !All bases are OK
              !Reset the number of elements in each xi direction
              ALLOCATE(newNumberOfElementsXi(newMeshDimension),STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate new number of elements xi.",err,error,*999)
              IF(generatedMesh%meshDimension==0) THEN
                !First time, default attributes
                newNumberOfElementsXi(1:newMeshDimension)=1
              ELSE IF(generatedMesh%meshDimension>newMeshDimension) THEN
                !New mesh dimension is less than the old mesh dimension
                newNumberOfElementsXi(1:newMeshDimension)=generatedMesh%cylinderMesh%numberOfElementsXi(1:newMeshDimension)
              ELSE
                !New mesh dimension is more than the old mesh dimension
                newNumberOfElementsXi(1:generatedMesh%meshDimension)= &
                  & generatedMesh%cylinderMesh%numberOfElementsXi(1:generatedMesh%meshDimension)
                newNumberOfElementsXi(generatedMesh%meshDimension+1:newMeshDimension)=1
              ENDIF
              CALL MOVE_ALLOC(newNumberOfElementsXi,generatedMesh%cylinderMesh%numberOfElementsXi)
              !Set the closed nature of the mesh.
              generatedMesh%cylinderMesh%closed=closed              
            CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
              !Check supplied bases
              IF(newMeshDimension<2) THEN
                CALL FlagError("The specified bases must have at least two xi directions for an ellipsoid mesh.",err,error,*999)
              ELSE                
                IF(ASSOCIATED(generatedMesh%ellipsoidMesh)) THEN
                  !Need to determine if the mesh is open or closed. Check if we have a top apex basis
                  closed=.FALSE.
                  numberOfBasesInASet=2
                  DO basisIdx=1,numberOfBases
                    IF((bases(basisIdx)%ptr%COLLAPSED_XI(1)/=BASIS_XI_COLLAPSED).OR.&
                      & (bases(basisIdx)%ptr%COLLAPSED_XI(2)/=BASIS_COLLAPSED_AT_XI1)) THEN
                      closed=.TRUE.
                      numberOfBasesInASet=3
                      EXIT
                    ENDIF
                  ENDDO !basisIdx
                  IF(closed) THEN
                    IF(MOD(numberOfBases,3)/=0) THEN
                      localError="The number of supplied bases of "//TRIM(NumberToVString(numberOfBases,"*",err,error))// &
                        & " is invalid. The number of bases must be divisiable by 3 for wall, lower and upper apex bases for "// &
                        & "a closed ellipsoid mesh."
                      CALL FlagError(localError,err,error,*999)
                    ENDIF
                  ELSE
                    IF(MOD(numberOfBases,2)/=0) THEN
                      localError="The number of supplied bases of "//TRIM(NumberToVString(numberOfBases,"*",err,error))// &
                        & " is invalid. The number of bases must be divisiable by 2 for wall and lower apex bases for "// &
                        & "an open ellipsoid mesh."
                      CALL FlagError(localError,err,error,*999)
                    ENDIF
                  ENDIF
                  numberOfBasisSets=numberOfBases/numberOfBasesInASet
                  !Check wall and bottom apex bases and if closed the top apex bases
                  DO basisSetIdx=1,numberOfBasisSets
                    !Wall basis
                    basisIdx=(basisSetIdx-1)*numberOfBasesInASet+1
                    IF((bases(basisIdx)%ptr%COLLAPSED_XI(1)/=BASIS_NOT_COLLAPSED).OR.&
                      & (bases(basisIdx)%ptr%COLLAPSED_XI(2)/=BASIS_NOT_COLLAPSED)) THEN
                      localError=="Specified basis number "//TRIM(NumberToVString(basisIdx,"*",err,error))// &
                        & " is invalid. The wall bases must not be collapsed in an ellipsoid mesh."
                      CALL FlagError(localError,err,error,*999)
                    ENDIF
                    !Bottom apex basis
                    basisIdx=(basisSetIdx-1)*numberOfBasesInASet+2
                    IF((bases(basisIdx)%ptr%COLLAPSED_XI(1)/=BASIS_XI_COLLAPSED).OR.&
                      & (bases(basisIdx)%ptr%COLLAPSED_XI(2)/=BASIS_COLLAPSED_AT_XI0)) THEN
                      localError="Specified basis number "//TRIM(NumberToVString(basisIdx,"*",err,error))// &
                        & " is invalid. The bottom apex bases must be collapsed correctly in an ellipsoid mesh."
                      CALL FlagError(localError,err,error,*999)
                    ENDIF
                    IF(closed) THEN
                      !Top apex basis
                      basisIdx=(basisSetIdx-1)*numberOfBasesInASet+3
                      IF((bases(basisIdx)%ptr%COLLAPSED_XI(1)/=BASIS_XI_COLLAPSED).OR.&
                        & (bases(basisIdx)%ptr%COLLAPSED_XI(2)/=BASIS_COLLAPSED_AT_XI1)) THEN
                        localError="Specified basis number "//TRIM(NumberToVString(basisIdx,"*",err,error))// &
                          & " is invalid. The top apex bases must be collapsed correctly in a closed ellipsoid mesh."
                        CALL FlagError(localError,err,error,*999)
                      ENDIF
                    ENDIF
                  ENDDO !basisSetIdx
                  IF(newMeshDimension==3) THEN
                    DO basisIdx=1,numberOfBases
                      IF(bases(basisIdx)%ptr%COLLAPSED_XI(3)/=BASIS_NOT_COLLAPSED) THEN
                        localError="Specified basis number "//TRIM(NumberToVString(basisIdx,"*",err,error))// &
                          & " is invalid. The bases must not be collapsed in the radial direction in an ellipsoid mesh."
                        CALL FlagError(localError,err,error,*999)
                      ENDIF
                    ENDDO !basisIdx
                  ENDIF
                ENDIF
                numberOfMeshComponents=numberOfBasisSets
              ELSE
                CALL FlagError("Ellipsoid generated mesh is not associated.",err,error,*999)                
              ENDIF
            ENDIF
            !All bases are OK
            !Reset the number of elements in each xi direction
            ALLOCATE(newNumberOfElementsXi(newMeshDimension),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate new number of elements xi.",err,error,*999)
            IF(generatedMesh%meshDimension==0) THEN
              !First time, default attributes
              newNumberOfElementsXi(1:newMeshDimension)=1
            ELSE IF(generatedMesh%meshDimension>newMeshDimension) THEN
              !New mesh dimension is less than the old mesh dimension
              newNumberOfElementsXi(1:newMeshDimension)=generatedMesh%ellipsoidMesh%numberOfElementsXi(1:newMeshDimension)
            ELSE
<<<<<<< HEAD
              !New mesh dimension is more than the old mesh dimension
              newNumberOfElementsXi(1:generatedMesh%meshDimension)= &
                & generatedMesh%ellipsoidMesh%numberOfElementsXi(1:generatedMesh%meshDimension)
              newNumberOfElementsXi(generatedMesh%meshDimension+1:newMeshDimension)=1
            ENDIF
            CALL MOVE_ALLOC(newNumberOfElementsXi,generatedMesh%ellipsoidMesh%numberOfElementsXi)
            !Set the closed nature of the mesh.
            generatedMesh%EllipsoidMesh%closed=closed              
          CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
            !Check supplied bases
            IF(newMeshDimension==1) THEN
              CALL FlagError("Not implemented.",err,error,*999)
            ELSE              
              CALL FlagError("The specified bases must have one xi direction for a fractal tree mesh.",err,error,*999)
            ENDIF
=======
              CALL FlagError("Regular generated mesh is not associated.",ERR,ERROR,*999)
            ENDIF
          CASE(GENERATED_MESH_POLAR_MESH_TYPE)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
          CASE DEFAULT
            localError="The generated mesh type of "//TRIM(NumberToVString(generatedMesh%generatedType,"*",err,error))// &
              & " is invalid."
<<<<<<< HEAD
            CALL FlagError(localError,err,error,*999)
=======
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
          END SELECT
          !Reset the base vectors given the new mesh dimension
          IF(newMeshDimension/=generatedMesh%meshDimension) THEN
            IF(ALLOCATED(generatedMesh%baseVectors)) THEN
              ALLOCATE(newBaseVectors(generatedMesh%coordinateDimension,newMeshDimension),STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate new base vectors.",err,error,*999)
              IF(generatedMesh%meshDimension>newMeshDimension) THEN
                !New number of xi is less than the old number of xi
                newBaseVectors(1:generatedMesh%coordinateDimension,1:newMeshDimension)= &
                  & generatedMesh%baseVectors(1:generatedMesh%coordinateDimension,1:newMeshDimension)
              ELSE
                !New number of xi is more than the old number of xi
                newBaseVectors(1:generatedMesh%coordinateDimension,1:generatedMesh%meshDimension)= &
                  & generatedMesh%baseVectors(1:generatedMesh%coordinateDimension,1:generatedMesh%meshDimension)
                newBaseVectors(1:generatedMesh%coordinateDimension,generatedMesh%meshDimension+1:newMeshDimension)=0.0_DP
              ENDIF
              CALL MOVE_ALLOC(newBaseVectors,generatedMesh%baseVectors)
            ENDIF
          ENDIF
          !Store the bases
          IF(ALLOCATED(generatedMesh%bases)) DEALLOCATE(generatedMesh%bases)
          ALLOCATE(generatedMesh%bases(numberOfBases),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate bases.",err,error,*999)
          DO basisIdx=1,numberOfBases
            generatedMesh%bases(basisIdx)%ptr=>bases(basisIdx)%ptr
          ENDDO !basisIdx
          !Set the number of mesh components
          generatedMesh%numberOfMeshComponents=numberOfMeshComponents
          !Reset the mesh dimension
          generatedMesh%meshDimension=newMeshDimension
        ELSE
<<<<<<< HEAD
          CALL FlagError("No bases where supplied.",err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Generated mesh is already associated.",err,error,*999)
    ENDIF

    CALL Exits("GeneratedMeshBasisSet")
    RETURN
999 IF(ALLOCATED(newNumberOfElementsXi)) DEALLOCATE(newNumberOfElementsXi)
    IF(ALLOCATED(newBaseVectors)) DEALLOCATE(newBaseVectors)
    CALL Errors("GeneratedMeshBasisSet",err,error)
    CALL Exits("GeneratedMeshBasisSet")
=======
          LOCAL_ERROR="The size of the first dimension of base vectors of "// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(BASE_VECTORS,1),"*",ERR,ERROR))// &
            & " is invalid. The first dimension size must match the coordinate system dimension of "// &
            & TRIM(NUMBER_TO_VSTRING(COORDINATE_DIMENSION,"*",ERR,ERROR))//"."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Generated mesh is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("GENERATED_MESH_BASE_VECTORS_SET")
    RETURN
999 ERRORSEXITS("GENERATED_MESH_BASE_VECTORS_SET",ERR,ERROR)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
    RETURN 1
    
  END SUBROUTINE GeneratedMeshBasisSet

  !
  !================================================================================================================================
  !

  !>Sets/changes the base vectors of a generated mesh. \see OPENCMISS::CMISSGeneratedMeshBaseVectorsSet
  SUBROUTINE GeneratedMeshBaseVectorsSet(generatedMesh,baseVectors,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh to set the base vectors fo
    REAL(DP), INTENT(IN) :: baseVectors(:,:) !<baseVectors(coordinateIdx,xiIdx). The base vectors for the generated mesh to set.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

<<<<<<< HEAD
    CALL Enters("GeneratedMeshBaseVectorsSet",err,error,*999)

    IF(ASSOCIATED(generatedMesh)) THEN
      IF(generatedMesh%generatedMeshFinished) THEN
        CALL FlagError("Generated mesh has been finished.",err,error,*999)
      ELSE
        IF(ALLOCATED(generatedMesh%bases)) THEN
          !Check the spatial dimension of the base vectors is the same as the coordinate space
          IF(SIZE(baseVectors,1)==generatedMesh%coordinateDimension) THEN
            !Check that the number of base vectors is equal to the mesh dimension
            IF(SIZE(baseVectors,2)==generatedMesh%meshDimension) THEN
              IF(ALLOCATED(generatedMesh%baseVectors)) DEALLOCATE(generatedMesh%baseVectors)
              ALLOCATE(generatedMesh%baseVectors(SIZE(baseVectors,1),SIZE(baseVectors,2)),STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate base vectors.",err,error,*999)
              generatedMesh%baseVectors(1:generatedMesh%coordinateDimension,1:generatedMesh%meshDimension)= &
                & baseVectors(1:generatedMesh%coordinateDimension,1:generatedMesh%meshDimension)
            ELSE
              localError="The size of the second dimension of base vectors of "// &
                & TRIM(NumberToVString(SIZE(baseVectors,2),"*",err,error))// &
                & " is invalid. The second dimension size must match the number of mesh dimensions of "// &
                & TRIM(NumberToVString(generatedMesh%meshDimension,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ELSE
            localError="The size of the first dimension of base vectors of "// &
              & TRIM(NumberToVString(SIZE(baseVectors,1),"*",err,error))// &
              & " is invalid. The first dimension size must match the coordinate system dimension of "// &
              & TRIM(NumberToVString(coordinateDimension,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("You must set the generated mesh basis before setting base vectors.",err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Generated mesh is not associated.",err,error,*999)
    ENDIF

    CALL Exits("GeneratedMeshBaseVectorsSet")
    RETURN
999 CALL Errors("GeneratedMeshBaseVectorsSet",err,error)
    CALL Exits("GeneratedMeshBaseVectorsSet")
    RETURN 1
    
  END SUBROUTINE GeneratedMeshBaseVectorsSet

  !
  !================================================================================================================================
  !

  !>Returns the coordinate system for a generated mesh accounting for regions and interfaces
  SUBROUTINE GeneratedMeshCoordinateSystemGet(generatedMesh,coordinateSystem,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh to get the coordinate system for
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: coordinateSystem !<On return, the generated meshes coordinate system. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(INTERFACE_TYPE), POINTER :: interface
    TYPE(REGION_TYPE), POINTER :: region
    TYPE(VARYING_STRING) :: localError

    NULLIFY(interface)
    NULLIFY(region)
    
    CALL Enters("GeneratedMeshCoordinateSystemGet",err,error,*999)

    IF(ASSOCIATED(generatedMesh)) THEN
      IF(ASSOCIATED(coordinateSystem)) THEN
        CALL FlagError("Coordinate system is already associated.",err,error,*999)
      ELSE
        NULLIFY(coordinateSystem)
        region=>generatedMesh%region
        IF(ASSOCIATED(region)) THEN
          coordinateSystem=>region%COORDINATE_SYSTEM
          IF(.NOT.ASSOCIATED(coordinateSystem)) THEN
            localError="The coordinate system is not associated for generated mesh number "// &
              & TRIM(NumberToVString(generatedMesh%userNumber,"*",err,error))//" of region number "// &
              & TRIM(NumberToVString(region%USER_NUMBER,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ELSE
          interface=>generatedMesh%interface
          IF(ASSOCIATED(interface)) THEN
            coordinateSystem=>interface%COORDINATE_SYSTEM
            IF(.NOT.ASSOCIATED(coordinateSystem)) THEN
              localError="The coordinate system is not associated for generated mesh number "// &
                & TRIM(NumberToVString(generatedMesh%userNumber,"*",err,error))//" of interface number "// &
                & TRIM(NumberToVString(interface%USER_NUMBER,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ELSE
            localError="The interface is not associated for generated mesh number "// &
              & TRIM(NumberToVString(generatedMesh%userNumber,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
=======
    ENTERS("GENERATED_MESH_COORDINATE_SYSTEM_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
        CALL FlagError("Coordinate system is already associated.",ERR,ERROR,*999)
      ELSE
        NULLIFY(COORDINATE_SYSTEM)
        REGION=>GENERATED_MESH%REGION
        IF(ASSOCIATED(REGION)) THEN
          COORDINATE_SYSTEM=>REGION%COORDINATE_SYSTEM
          IF(.NOT.ASSOCIATED(COORDINATE_SYSTEM)) THEN
            LOCAL_ERROR="The coordinate system is not associated for generated mesh number "// &
              & TRIM(NUMBER_TO_VSTRING(GENERATED_MESH%USER_NUMBER,"*",ERR,ERROR))//" of region number "// &
              & TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          NULLIFY(INTERFACE)
          INTERFACE=>GENERATED_MESH%INTERFACE
          IF(ASSOCIATED(INTERFACE)) THEN
            COORDINATE_SYSTEM=>INTERFACE%COORDINATE_SYSTEM
            IF(.NOT.ASSOCIATED(COORDINATE_SYSTEM)) THEN
              LOCAL_ERROR="The coordinate system is not associated for generated mesh number "// &
                & TRIM(NUMBER_TO_VSTRING(GENERATED_MESH%USER_NUMBER,"*",ERR,ERROR))//" of interface number "// &
                & TRIM(NUMBER_TO_VSTRING(INTERFACE%USER_NUMBER,"*",ERR,ERROR))//"."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The interface is not associated for generated mesh number "// &
              & TRIM(NUMBER_TO_VSTRING(GENERATED_MESH%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
          ENDIF
        ENDIF
      ENDIF
    ELSE
<<<<<<< HEAD
      CALL FlagError("Generated mesh is not associated.",err,error,*999)
    ENDIF

    CALL Exits("GeneratedMeshCoordinateSystemGet")
    RETURN
999 CALL Errors("GeneratedMeshCoordinateSystemGet",err,error)
    CALL Exits("GeneratedMeshCoordinateSystemGet")
=======
      CALL FlagError("Generated mesh is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("GENERATED_MESH_COORDINATE_SYSTEM_GET")
    RETURN
999 ERRORSEXITS("GENERATED_MESH_COORDINATE_SYSTEM_GET",ERR,ERROR)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
    RETURN 1
    
  END SUBROUTINE GeneratedMeshCoordinateSystemGet

  !
  !================================================================================================================================
  !

  !>Finishes the creation of a generated mesh. \see OPENCMISS::CMISSGeneratedMeshCreateFinish
  SUBROUTINE GeneratedMeshCreateFinish(generatedMesh,meshUserNumber,mesh,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh to finish the creation of
    INTEGER(INTG), INTENT(IN) :: meshUserNumber !<The mesh's user number
    TYPE(MESH_TYPE), POINTER :: mesh !<On exit, a pointer to the generated mesh. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

<<<<<<< HEAD
    CALL Enters("GeneratedMeshCreateFinish",err,error,*999)

    IF(ASSOCIATED(generatedMesh)) THEN
      IF(generatedMesh%generatedMeshFinished) THEN
        CALL FlagError("Generated mesh has already been finished.",err,error,*999)
      ELSE
        IF(ASSOCIATED(mesh)) THEN
          CALL FlagError("Mesh is already associated.",err,error,*999)
        ELSE
          IF(ALLOCATED(generatedMesh%bases)) THEN
            IF(.NOT.ALLOCATED(generatedMesh%baseVectors)) THEN
              !We don't have any base vectors defined so default the values to the extents
              ALLOCATE(newBaseVectors(generatedMesh%coordinateDimension,generatedMesh%meshDimension),STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate new base vectors.",err,error,*999)
              SELECT CASE(generatedMesh%generatedType)
              CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
                IF(ASSOCIATED(generatedMesh%regularMesh)) THEN
                  IF(generatedMesh%meshDimension==1) THEN
                    !The base vector is just the extent vector
                    newBaseVectors(1:generatedMesh%coordinateDimension,1)= &
                      & generatedMesh%regularMesh%maximumExtent(1:generatedMesh%coordinateDimension)
                  ELSE IF(generatedMesh%meshDimension<generatedMesh%coordinateDimension)             
                    !Find the first number of mesh dimensions for which the extent is non-zero.
                    count=0
                    coordinateIdx=1
                    DO xiIdx=1,generatedMesh%meshDimension
                      DO WHILE(ABS(generatedMesh%regularMesh%maximumExtent(coordinateIdx))<=ZERO_TOLERANCE)
                        coordinateIdx=coordinateIdx+1
                      ENDDO !While
                      newBaseVectors(coordinateIdx,xiIdx)=generatedMesh%regularMesh%maximumExtent(coordinateIdx)
                      coordinateIdx=coordinateIdx+1
                      count=count+1
                    ENDDO !xiIdx
                    IF(count/=generatedMesh%meshDimension)  &
                      & CALL FlagError("Invalid mesh extent. There number of non-zero components is < the mesh dimension.", &
                      & err,error,*999)
                  ELSE
                    !Number of xi is the same as the number of coordinates
                    !The default base vectors are aligned with the coordinate vectors
                    newBaseVectors=0.0_DP
                    DO coordinateIdx=1,coordinateDimension
                      newBaseVectors(coordinateIdx,coordinateIdx)=generatedMesh%regularMesh%maximumExtent(coordinateIdx)
                    ENDDO !coordinate_idx
                  ENDIF
                ELSE
                  CALL FlagError("Generated mesh regular mesh is not associated.",err,error,*999)
                ENDIF
              CASE(GENERATED_MESH_POLAR_MESH_TYPE,GENERATED_MESH_CYLINDER_MESH_TYPE,GENERATED_MESH_ELLIPSOID_MESH_TYPE, &
                & GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
                !The default base vectors are aligned with the coordinate vectors
                newBaseVectors=0.0_DP
                DO xiIdx=1,generatedMesh%meshDimension
                  newBaseVectors(xiIdx,xiIdx)=1.0_DP
                ENDDO !xiIdx
              CASE DEFAULT
                localError="The generated mesh type of "//TRIM(NumberToVString(generatedMesh%generatedType,"*",err,error))// &
                  & " is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ENDIF
            CALL MOVE_ALLOC(newBaseVectors,generatedMesh%baseVectors)
            !Finish the specific mesh type
            SELECT CASE(generatedMesh%generatedType)
            CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
              CALL GeneratedMeshRegularCreateFinish(generatedMesh%regularMesh,meshUserNumber,err,error,*999)
            CASE(GENERATED_MESH_POLAR_MESH_TYPE)
              CALL GeneratedMeshPolarCreateFinish(generatedMesh%polarMesh,meshUserNumber,err,error,*999)
            CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
              CALL GeneratedMeshCylinderCreateFinish(generatedMesh%cylinderMesh,meshUserNumber,err,error,*999)
            CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
              CALL GeneratedMeshEllipsoidCreateFinish(generatedMesh%ellipsoidMesh,meshUserNumber,err,error,*999)
            CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
              CALL GeneratedMeshFractalTreeCreateFinish(generatedMesh%fractalTreeMesh,meshUserNumber,err,error,*999)
            CASE DEFAULT
              localError="The generated mesh mesh type of "// &
                & TRIM(NumberToVString(generatedMesh%generatedType,"*",err,error))//" is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
            !Return the pointers
            mesh=>generatedMesh%mesh
            mesh%GENERATED_MESH=>generatedMesh
            generatedMesh%generatedMeshFinished=.TRUE.
          ELSE
            CALL FlagError("You must set the generated mesh basis before finishing the generated mesh.",err,error,*999)
          ENDIF
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Generated mesh is not associated.",err,error,*999)
    ENDIF

    CALL Exits("GeneratedMeshCreateFinish")
    RETURN
999 CALL Errors("GeneratedMeshCreateFinish",err,error)
    CALL Exits("GeneratedMeshCreateFinish")
=======
    ENTERS("GENERATED_MESH_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      IF(GENERATED_MESH%GENERATED_MESH_FINISHED) THEN
        CALL FlagError("Generated mesh has already been finished.",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(MESH)) THEN
          CALL FlagError("Mesh is already associated.",ERR,ERROR,*999)
        ELSE
          SELECT CASE(GENERATED_MESH%GENERATED_TYPE)
          CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
            CALL GENERATED_MESH_REGULAR_CREATE_FINISH(GENERATED_MESH,MESH_USER_NUMBER,ERR,ERROR,*999)
          CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
            CALL GENERATED_MESH_CYLINDER_CREATE_FINISH(GENERATED_MESH,MESH_USER_NUMBER,ERR,ERROR,*999)
          CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
            CALL GENERATED_MESH_ELLIPSOID_CREATE_FINISH(GENERATED_MESH,MESH_USER_NUMBER,ERR,ERROR,*999)
          CASE(GENERATED_MESH_POLAR_MESH_TYPE)
            CALL FlagError("Not implmented.",ERR,ERROR,*999)
          CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The generated mesh mesh type of "// &
              & TRIM(NUMBER_TO_VSTRING(GENERATED_MESH%GENERATED_TYPE,"*",ERR,ERROR))//" is invalid."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
          !Return the pointers
          MESH=>GENERATED_MESH%MESH
          MESH%GENERATED_MESH=>GENERATED_MESH
          GENERATED_MESH%GENERATED_MESH_FINISHED=.TRUE.
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Generated mesh is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("GENERATED_MESH_CREATE_FINISH")
    RETURN
999 ERRORSEXITS("GENERATED_MESH_CREATE_FINISH",ERR,ERROR)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
    RETURN 1
    
  END SUBROUTINE GeneratedMeshCreateFinish

  !
  !================================================================================================================================
  !

  !>Starts the creation of a generic generated mesh.
  SUBROUTINE GeneratedMeshCreateStartGeneric(generatedMeshes,userNumber,coordinateSystem,generatedMesh,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshesType), POINTER :: generatedMeshes !<A pointer to the generated meshes
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the generated mesh to create
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: coordinateSystem !<A pointer to the coordinate system to create the generated mesh in.
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<On exit, a pointer to the created generated mesh. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,generatedMeshIdx
    TYPE(GeneratedMeshType), POINTER :: newGeneratedMesh
    TYPE(GeneratedMeshPtrType), ALLOCATABLE :: newGeneratedMeshes(:)
    TYPE(VARYING_STRING) :: dummyError

    NULLIFY(newGeneratedMesh)

<<<<<<< HEAD
    CALL Enters("GeneratedMeshCreateStartGeneric",err,error,*997)

    IF(ASSOCIATED(generatedMeshes)) THEN
      IF(ASSOCIATED(generatedMesh)) THEN
        CALL FlagError("Generated mesh is already associated.",err,error,*997)
      ELSE
        IF(ASSOCIATED(coordianteSystem)) THEN
          !Initialise generated mesh
          CALL GeneratedMeshInitialise(coordinateSystem,newGeneratedMesh,err,error,*999)
          !Set default generated mesh values
          newGeneratedMesh%userNumber=userNumber
          newGeneratedMesh%globalNumber=generatedMeshes%numberOfGeneratedMeshes+1
          newGeneratedMesh%generatedMeshes=>generatedMeshes
          !Add new generated mesh into list of generated meshes
          ALLOCATE(newGeneratedMeshes(generatedMeshes%numberOfGeneratedMeshes+1),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate new generated meshes.",err,error,*999)
          DO generatedMeshIdx=1,generatedMeshes%numberOfGeneratedMeshes
            newGeneratedMeshes(generatedMeshIdx)%ptr=>generatedMeshes%generatedMeshes(generatedMeshIdx)%ptr
          ENDDO !generatedMeshIdx
          newGeneratedMeshes(generatedMeshes%numberOfGeneratedMeshes+1)%ptr=>newGeneratedMesh
          CALL MOVE_ALLOC(newGeneratedMeshes,generatedMeshes%generatedMeshes)
          generatedMeshes%numberOfGeneratedMeshes=generatedMeshes%numberOfGeneratedMeshes+1
          !Return the pointer
          generatedMesh=>newGeneratedMesh
        ELSE
          CALL FlagError("Coordinate system is not associated.",err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Generated meshes is not associated.",err,error,*997)
    ENDIF

    CALL Exits("GeneratedMeshCreateStartGeneric")
    RETURN
999 CALL GeneratedMeshFinalise(newGeneratedMesh,dummyErr,dummyError,*998)
998 IF(ALLOCATED(newGeneratedMeshes)) DEALLOCATE(newGeneratedMeshes)
    NULLIFY(generatedMesh)
997 CALL Errors("GeneratedMeshCreateStartGeneric",err,error)
    CALL Exits("GeneratedMeshCreateStartGeneric")
=======
    ENTERS("GENERATED_MESH_CREATE_START_GENERIC",ERR,ERROR,*997)

    IF(ASSOCIATED(GENERATED_MESHES)) THEN
      IF(ASSOCIATED(GENERATED_MESH)) THEN
        CALL FlagError("Generated mesh is already associated.",ERR,ERROR,*997)
      ELSE
        !Initialise generated mesh
        CALL GENERATED_MESH_INITIALISE(NEW_GENERATED_MESH,ERR,ERROR,*999)
        !Set default generated mesh values
        NEW_GENERATED_MESH%USER_NUMBER=USER_NUMBER
        NEW_GENERATED_MESH%GLOBAL_NUMBER=GENERATED_MESHES%NUMBER_OF_GENERATED_MESHES+1
        NEW_GENERATED_MESH%GENERATED_MESHES=>GENERATED_MESHES
        !Add new generated mesh into list of generated meshes
        ALLOCATE(NEW_GENERATED_MESHES(GENERATED_MESHES%NUMBER_OF_GENERATED_MESHES+1),STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate new generated meshes.",ERR,ERROR,*999)
        DO generated_mesh_idx=1,GENERATED_MESHES%NUMBER_OF_GENERATED_MESHES
          NEW_GENERATED_MESHES(generated_mesh_idx)%PTR=>GENERATED_MESHES%GENERATED_MESHES(generated_mesh_idx)%PTR
        ENDDO !generated_mesh_idx
        NEW_GENERATED_MESHES(GENERATED_MESHES%NUMBER_OF_GENERATED_MESHES+1)%PTR=>NEW_GENERATED_MESH
        IF(ASSOCIATED(GENERATED_MESHES%GENERATED_MESHES)) DEALLOCATE(GENERATED_MESHES%GENERATED_MESHES)
        GENERATED_MESHES%GENERATED_MESHES=>NEW_GENERATED_MESHES
        GENERATED_MESHES%NUMBER_OF_GENERATED_MESHES=GENERATED_MESHES%NUMBER_OF_GENERATED_MESHES+1
        !Return the pointer
        GENERATED_MESH=>NEW_GENERATED_MESH
      ENDIF
    ELSE
      CALL FlagError("Generated meshes is not associated.",ERR,ERROR,*997)
    ENDIF

    EXITS("GENERATED_MESH_CREATE_START_GENERIC")
    RETURN
999 CALL GENERATED_MESH_FINALISE(NEW_GENERATED_MESH,DUMMY_ERR,DUMMY_ERROR,*998)
998 IF(ASSOCIATED(NEW_GENERATED_MESHES)) DEALLOCATE(NEW_GENERATED_MESHES)
    NULLIFY(GENERATED_MESH)
997 ERRORSEXITS("GENERATED_MESH_CREATE_START_GENERIC",ERR,ERROR)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
    RETURN 1

  END SUBROUTINE GeneratedMeshCreateStartGeneric

  !
  !================================================================================================================================
  !

  !>Starts the creation of a generated mesh. \see OPENCMISS::CMISSGeneratedMeshCreateFinish
  SUBROUTINE GeneratedMeshCreateStartInterface(userNumber,interface,generatedMesh,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the generated mesh to create
    TYPE(INTERFACE_TYPE), POINTER :: interface !<A pointer to the interface to create the generated mesh on
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<On exit, a pointer to the created generated mesh. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(COORDINATE_SYSTEM), POINTER :: coordinateSystem
    TYPE(VARYING_STRING) :: localError

<<<<<<< HEAD
    NULLIFY(coordinateSystem)

    CALL Enters("GeneratedMeshCreateStartInterface",err,error,*999)

    IF(ASSOCIATED(interface)) THEN
      IF(ASSOCIATED(generatedMesh)) THEN
        CALL FlagError("Generated mesh is already associated.",err,error,*999)
=======
    ENTERS("GENERATED_MESH_CREATE_START_INTERFACE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE)) THEN
      IF(ASSOCIATED(GENERATED_MESH)) THEN
        CALL FlagError("Generated mesh is already associated.",ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
      ELSE
        NULLIFY(generatedMesh)
        CALL GeneratedMeshUserNumberFind(userNumber,interface,generatedMesh,err,error,*999)
        IF(ASSOCIATED(generatedMesh)) THEN
          localError="The specified user number of "//TRIM(NumberToVString(userNumber,"*",err,error))// &
            & " has already been used for a generated mesh."
<<<<<<< HEAD
          CALL FlagError(localError,err,error,*999)
=======
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
        ELSE
          IF(ASSOCIATED(INTERFACE%GENERATED_MESHES)) THEN
            CALL INTERFACE_COORDINATE_SYSTEM_GET(INTERFACE,coordinateSystem,err,error,*999)
            CALL GeneratedMeshCreateStartGeneric(INTERFACE%GENERATED_MESHES,userNumber,coordinateSystem,generatedMesh, &
              & err,error,*999)
            generatedMesh%interface=>interface
          ELSE
<<<<<<< HEAD
            CALL FlagError("Interface generated meshes is not associated.",err,error,*999)
=======
            CALL FlagError("Interface generated meshes is not associated.",ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
          ENDIF
        ENDIF
      ENDIF
    ELSE
<<<<<<< HEAD
      CALL FlagError("Interface is not associated.",err,error,*999)
    ENDIF

    CALL Exits("GeneratedMeshCreateStartInterface")
    RETURN
999 CALL Errors("GeneratedMeshCreateStartInterface",err,error)
    CALL Exits("GeneratedMeshCreateStartInterface")
=======
      CALL FlagError("Interface is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("GENERATED_MESH_CREATE_START_INTERFACE")
    RETURN
999 ERRORSEXITS("GENERATED_MESH_CREATE_START_INTERFACE",ERR,ERROR)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
    RETURN 1
    
  END SUBROUTINE GeneratedMeshCreateStartInterface

  !
  !================================================================================================================================
  !

  !>Starts the creation of a generated mesh. \see OPENCMISS::CMISSGeneratedMeshCreateFinish
  SUBROUTINE GeneratedMeshCreateStartRegion(userNumber,region,generatedMesh,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the generated mesh to create
    TYPE(REGION_TYPE), POINTER :: region !<A pointer to the region to create the generated mesh on
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<On exit, a pointer to the created generated mesh. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(COORDINATE_SYSTEM), POINTER :: coordinateSystem
    TYPE(VARYING_STRING) :: localError

<<<<<<< HEAD
    NULLIFY(coordinateSystem)
    
    CALL Enters("GeneratedMeshCreateStartRegion",err,error,*999)

    IF(ASSOCIATED(region)) THEN
      IF(ASSOCIATED(generatedMesh)) THEN
        CALL FlagError("Generated mesh is already associated.",err,error,*999)
=======
    ENTERS("GENERATED_MESH_CREATE_START_REGION",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(GENERATED_MESH)) THEN
        CALL FlagError("Generated mesh is already associated.",ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
      ELSE
        NULLIFY(generatedMesh)
        CALL GeneratedMeshUserNumberFind(userNumber,region,generatedMesh,err,error,*999)
        IF(ASSOCIATED(generatedMesh)) THEN
          localError="The specified user number of "//TRIM(NumberToVString(userNumber,"*",err,error))// &
            & " has already been used for a generated mesh."
<<<<<<< HEAD
          CALL FlagError(localError,err,error,*999)
=======
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
        ELSE
          IF(ASSOCIATED(region%GENERATED_MESHES)) THEN
            CALL REGION_COORDINATE_SYSTEM_GET(region,coordinateSystem,err,error,*999)
            CALL GeneratedMeshCreateStartGeneric(region%GENERATED_MESHES,userNumber,coordinateSystem,generatedMesh, &
              & err,error,*999)
            generatedMesh%region=>region
          ELSE
<<<<<<< HEAD
            CALL FlagError("Region generated meshes is not associated.",err,error,*999)
=======
            CALL FlagError("Region generated meshes is not associated.",ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
          ENDIF
        ENDIF
      ENDIF
    ELSE
<<<<<<< HEAD
      CALL FlagError("Region is not associated.",err,error,*999)
    ENDIF

    CALL Exits("GeneratedMeshCreateStartRegion")
    RETURN
999 CALL Errors("GeneratedMeshCreateStartRegion",err,error)
    CALL Exits("GeneratedMeshCreateStartRegion")
=======
      CALL FlagError("Region is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("GENERATED_MESH_CREATE_START_REGION")
    RETURN
999 ERRORSEXITS("GENERATED_MESH_CREATE_START_REGION",ERR,ERROR)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
    RETURN 1
    
  END SUBROUTINE GeneratedMeshCreateStartRegion

  !
  !================================================================================================================================
  !

  !>Destroys a generated mesh. \see OPENCMISS::CMISSGeneratedMeshDestroy
  SUBROUTINE GeneratedMeshDestroy(generatedMesh,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
<<<<<<< HEAD
    INTEGER(INTG) :: generatedMeshIdx,generatedMeshPosition
    TYPE(GeneratedMeshesType), POINTER :: generatedMeshes
    TYPE(GeneratedMeshPtrType), ALLOCATABLE :: newGeneratedMeshes(:)

    CALL Enters("GeneratedMeshDestroy",err,error,*998)

    IF(ASSOCIATED(generatedMesh)) THEN
      generatedMeshes=>generatedMesh%generatedMeshes
      IF(ASSOCIATED(generatedMeshes)) THEN
        IF(ALLOCATED(generatedMeshes%generatedMeshes)) THEN
          generatedMeshPosition=generatedMesh%globalNumber
          !Remove the generated mesh from the list of generated meshes
          IF(generatedMeshes%numberOfGeneratedMeshes>1) THEN
            ALLOCATE(newGeneratedMeshes(generatedMeshes%numberOfGeneratedMeshes-1),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate new generated meshes.",err,error,*999)
            DO generatedMeshIdx=1,generatedMeshes%numberOfGeneratedMeshes
              IF(generatedMeshIdx<generatedMeshPosition) THEN
                newGeneratedMeshes(generatedMeshIdx)%ptr=>generatedMeshes%generatedMeshes(generatedMeshIdx)%ptr
              ELSE IF(generatedMeshIdx>generatedMeshPosition) THEN
                generatedMeshes%generatedMeshes(generatedMeshIdx)%ptr%globalNumber=generatedMeshes% &
                  & generatedMeshes(generatedMeshIdx)%ptr%globalnumber-1
                newGeneratedMeshes(generatedMeshIdx-1)%ptr=>generatedMeshes%generatedMeshes(generatedMeshIdx)%ptr
=======
    INTEGER(INTG) :: generated_mesh_idx,generated_mesh_position
    TYPE(GENERATED_MESH_PTR_TYPE), POINTER :: NEW_GENERATED_MESHES(:)
    TYPE(GENERATED_MESHES_TYPE), POINTER :: GENERATED_MESHES

    ENTERS("GENERATED_MESH_DESTROY",ERR,ERROR,*998)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      GENERATED_MESHES=>GENERATED_MESH%GENERATED_MESHES
      IF(ASSOCIATED(GENERATED_MESHES)) THEN
        IF(ASSOCIATED(GENERATED_MESHES%GENERATED_MESHES)) THEN
          generated_mesh_position=GENERATED_MESH%GLOBAL_NUMBER
          CALL GENERATED_MESH_FINALISE(GENERATED_MESH,ERR,ERROR,*999)
          !Remove the generated mesh from the list of generated meshes
          IF(GENERATED_MESHES%NUMBER_OF_GENERATED_MESHES>1) THEN
            ALLOCATE(NEW_GENERATED_MESHES(GENERATED_MESHES%NUMBER_OF_GENERATED_MESHES-1),STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate new generated meshes.",ERR,ERROR,*999)
            DO generated_mesh_idx=1,GENERATED_MESHES%NUMBER_OF_GENERATED_MESHES
              IF(generated_mesh_idx<generated_mesh_position) THEN
                NEW_GENERATED_MESHES(generated_mesh_idx)%PTR=>GENERATED_MESHES%GENERATED_MESHES(generated_mesh_idx)%PTR
              ELSE IF(generated_mesh_idx>generated_mesh_position) THEN
                GENERATED_MESHES%GENERATED_MESHES(generated_mesh_idx)%PTR%GLOBAL_NUMBER=GENERATED_MESHES% &
                  & GENERATED_MESHES(generated_mesh_idx)%PTR%GLOBAL_NUMBER-1
                NEW_GENERATED_MESHES(generated_mesh_idx-1)%PTR=>GENERATED_MESHES%GENERATED_MESHES(generated_mesh_idx)%PTR
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
              ENDIF
            ENDDO !generatedMeshIdx
            CALL MOVE_ALLOC(newGeneratedMeshes,generatedMeshes%generatedMeshes)
            generatedMeshes%numberOfGeneratedMeshes=generatedMeshes%numberOfGeneratedMeshes-1
          ELSE
            DEALLOCATE(generatedMeshes%generatedMeshes)
            generatedMeshes%numberOfGeneratedMeshes=0
          ENDIF
          !Finalise the generated mesh
          CALL GeneratedMeshFinalise(generatedMesh,err,error,*999)
        ELSE
<<<<<<< HEAD
          CALL FlagError("Generated meshes have not been allocated.",err,error,*998)
        ENDIF
      ELSE
        CALL FlagError("Generated mesh generated meshes is not associated.",err,error,*998)
      ENDIF
    ELSE
      CALL FlagError("Generated mesh is not associated.",err,error,*998)
    END IF

    CALL Exits("GeneratedMeshDestroy")
    RETURN
999 IF(ALLOCATED(newGeneratedMeshes)) DEALLOCATE(newGeneratedMeshes)
998 CALL Errors("GeneratedMeshDestroy",err,error)
    CALL Exits("GeneratedMeshDestroy")
=======
          CALL FlagError("Generated meshes are not associated",ERR,ERROR,*998)
        ENDIF
      ELSE
        CALL FlagError("Generated mesh generated meshes is not associated.",ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FlagError("Generated mesh is not associated",ERR,ERROR,*998)
    END IF

    EXITS("GENERATED_MESH_DESTROY")
    RETURN
999 IF(ASSOCIATED(NEW_GENERATED_MESHES)) DEALLOCATE(NEW_GENERATED_MESHES)
998 ERRORSEXITS("GENERATED_MESH_DESTROY",ERR,ERROR)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
    RETURN 1
    
  END SUBROUTINE GeneratedMeshDestroy

  !
  !================================================================================================================================
  !

  
  !>Gets the extent of a generated mesh. \see OPENCMISS::CMISSGeneratedMeshExtentGet
  !>\todo Add in routine to return the number of extents
  SUBROUTINE GeneratedMeshExtentGet(generatedMesh,extent,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh to get the type of
    REAL(DP), INTENT(OUT) :: extent(:) !<On return, the mesh extent. For regular meshes this is the maximum extent per axis. For polar meshes it is the radius per shell. For cylindrical meshes this is the length and then radii of each cylinder. For ellipsoid meshes it is .....
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

<<<<<<< HEAD
    CALL Enters("GeneratedMeshExtentGet",err,error,*999)
=======
    ENTERS("GENERATED_MESH_EXTENT_GET",ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b

    IF(ASSOCIATED(generatedMesh)) THEN
      SELECT CASE(generatedMesh%generatedType)
      CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
        IF(ASSOCIATED(generatedMesh%regularMesh)) THEN
          numberOfExtents=SIZE(generatedMesh%polarMesh%polarExtent,1)
          IF(SIZE(extent,1)>=numberOfExtents) THEN
            extent(1:numberOfExtents)=generatedMesh%regularMesh%maximumExtent(1:numberOfExtents)
          ELSE
            localError="The size of extent is too small. The supplied size is "// &
              & TRIM(NumberToVString(SIZE(extent,1),"*",err,error))//" and it needs to be >= "// &
              & TRIM(NumberToVString(numberOfExtents,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ELSE
<<<<<<< HEAD
          CALL FlagError("Generated mesh regular mesh is not associated.",err,error,*999)
        ENDIF
      CASE(GENERATED_MESH_POLAR_MESH_TYPE)
        IF(ASSOCIATED(generatedMesh%polarMesh)) THEN
          numberOfExtents=SIZE(generatedMesh%polarMesh%polarExtent,1)
          IF(SIZE(extent,1)>=numberOfExtents) THEN
            extent(1:numberOfExtents)=generatedMesh%polarMesh%polarExtent(1:numberOfExtents)
          ELSE
            localError="The size of extent is too small. The supplied size is "// &
              & TRIM(NumberToVString(SIZE(extent,1),"*",err,error))//" and it needs to be >= "// &
              & TRIM(NumberToVString(numberOfExtents,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Generated mesh polar mesh is not associated.",err,error,*999)
        ENDIF
        CALL FlagError("Not implemented.",err,error,*999)
=======
          LOCAL_ERROR="The size of EXTENT is too small. The supplied size is "// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(EXTENT,1),"*",ERR,ERROR))//" and it needs to be >= "// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(GENERATED_MESH%REGULAR_MESH%MAXIMUM_EXTENT,1),"*",ERR,ERROR))//"."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      CASE(GENERATED_MESH_POLAR_MESH_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
      CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
        IF(ASSOCIATED(generatedMesh%cylinderMesh)) THEN
          numberOfExtents=SIZE(generatedMesh%cylinderMesh%cylinderExtent,1)
          IF(SIZE(extent,1)>=numberOfExtents) THEN
            extent(1:numberOfExtents)=generatedMesh%cylinderMesh%cylinderExtent(1:numberOfExtents)
          ELSE
            localError="The size of extent is too small. The supplied size is "// &
              & TRIM(NumberToVString(SIZE(extent,1),"*",err,error))//" and it needs to be >= "// &
              & TRIM(NumberToVString(numberOfExtents,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ELSE
<<<<<<< HEAD
          CALL FlagError("Generated mesh cylinder mesh is not associated.",err,error,*999)
=======
          LOCAL_ERROR="The size of EXTENT is too small. The supplied size is "// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(EXTENT,1),"*",ERR,ERROR))//" and it needs to be 3."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
        ENDIF
      CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
        IF(ASSOCIATED(generatedMesh%ellipsoidMesh)) THEN
          numberOfExtents=SIZE(generatedMesh%ellipsoidMesh%ellipsoidExtent,1)
          IF(SIZE(extent,1)>=numberOfExtents) THEN
            extent(1:numberOfExtents)=generatedMesh%ellipsoidMesh%ellipsoidExtent(1:numberOfExtents)
          ELSE
            localError="The size of extent is too small. The supplied size is "// &
              & TRIM(NumberToVString(SIZE(extent,1),"*",err,error))//" and it needs to be >= "// &
              & TRIM(NumberToVString(numberOfExtents,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Generated mesh cylinder mesh is not associated.",err,error,*999)
        ENDIF
      CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="The generated mesh mesh type of "//TRIM(NumberToVString(generatedMesh%generatedType,"*",err,error))// &
          & " is invalid."
<<<<<<< HEAD
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Generated mesh is not associated",err,error,*999)
    ENDIF

    CALL Exits("GeneratedMeshExtentGet")
    RETURN
999 CALL Errors("GeneratedMeshExtentGet",err,error)
    CALL Exits("GeneratedMeshExtentGet")
    RETURN 1
    
  END SUBROUTINE GeneratedMeshExtentGet
=======
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Generated mesh is not associated",ERR,ERROR,*999)
    ENDIF

    EXITS("GENERATED_MESH_EXTENT_GET")
    RETURN
999 ERRORSEXITS("GENERATED_MESH_EXTENT_GET",ERR,ERROR)
    RETURN
  END SUBROUTINE GENERATED_MESH_EXTENT_GET
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
  !
  !================================================================================================================================
  !

  !>Sets/changes the extent of a generated mesh. \see OPENCMISS::CMISSGeneratedMeshExtentSet
  SUBROUTINE GeneratedMeshExtentSet(generatedMesh,extent,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh to set the type of
    REAL(DP), INTENT(IN) :: extent(:) !<The extent of the generated mesh
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: extentIdx,minimumNumberOfExtents,shortAxisExtentStart
    REAL(DP) :: previousExtent
    REAL(DP), ALLOCATABLE :: newExtents(:)
    TYPE(VARYING_STRING) :: localError

<<<<<<< HEAD
    CALL Enters("GeneratedMeshExtentSet",err,error,*999)

    IF(ASSOCIATED(generatedMesh)) THEN
      IF(generatedMesh%generatedMeshFinished) THEN
        CALL FlagError("Generated mesh has been finished.",err,error,*999)
=======
    ENTERS("GENERATED_MESH_EXTENT_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      IF(GENERATED_MESH%GENERATED_MESH_FINISHED) THEN
        CALL FlagError("Generated mesh has been finished.",ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
      ELSE
        SELECT CASE(generatedMesh%generatedType)
        CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
<<<<<<< HEAD
          IF(ASSOCIATED(generatedMesh%regularMesh)) THEN
            IF(SIZE(extent,1)>=generatedMesh%coordinateDimension) THEN
              IF(L2Norm(extent(1:generatedMesh%coordinateDimension))>ZERO_TOLERANCE) THEN
                ALLOCATE(newExtents(1:generatedMesh%coordinateDimension),STAT=err)
                IF(err/=0) CALL FlagError("Could not allocate new extents.",err,error,*999)
                newExtents(1:generatedMesh%coordinateDimension)=extent(1:generatedMesh%coordinateDimension)
                CALL MOVE_ALLOC(newExtents,generatedMesh%regularMesh%maximumExtent)
              ELSE
                CALL FlagError("The norm of the mesh extent is zero.",err,error,*999)
              ENDIF
            ELSE
              localError="The specified number of extents of "//TRIM(NumberToVString(SIZE(extent,1),"*",err,error))// &
                & " is invalid. The number of extents must be >= the coordinate system dimension of "// &
                & TRIM(NumberToVString(generatedMesh%coordinateDimension,"*",err,error))// &
                & " for a regular generated mesh."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ELSE
            CALL FlagError("Regular generated mesh is not associated.",err,error,*999)
          ENDIF
        CASE(GENERATED_MESH_POLAR_MESH_TYPE)
          IF(ASSOCIATED(generatedMesh%polarMesh)) THEN
            IF(generatedMesh%polarMesh%sperhical) THEN
              IF(generatedMesh%meshDimension>2) THEN
                minimumNumberOfExtents=2
              ELSE
                minimumNumberOfExtents=1                
              ENDIF
            ELSE
              IF(generatedMesh%meshDimension>1) THEN
                minimumNumberOfExtents=2
              ELSE
                minimumNumberOfExtents=1                
              ENDIF              
            ENDIF
            numberOfExtents=SIZE(extent,1)
            IF(numberOfExtents>=minimumNumberOfExtents) THEN
              !Check the extents
              previousExtent=0.0_DP
              DO extentIdx=1,numberOfExtents
                IF(extent(extentIdx)>ZERO_TOLERANCE) THEN
                  IF(extent(extentIdx)>previousExtent) THEN
                    previousExtent=extent(extentIdx)+ZERO_TOLERANCE
                  ELSE
                    localError="Invalid extents. The specified extent index of "// &
                      & TRIM(NumberToVString(extentIdx,"*",err,error))// &
                      & " has an extent of "//TRIM(NumberToVString(extent(extentIdx),"*",err,error))// &
                      & ". The extent needs to be bigger than the previous extent of "// &
                      & TRIM(NumberToVString(previousExtent,"*",err,error))//"."
                    CALL FlagError(localError,err,error,*999)
                  ENDIF
                ELSE
                  localError="Invalid extents. The specified extent index of "// &
                    & TRIM(NumberToVString(extentIdx,"*",err,error))// &
                    & " has an extent of "//TRIM(NumberToVString(extent(extentIdx),"*",err,error))// &
                    & ". The extent needs to be > zero."
                  CALL FlagError(localError,err,error,*999)                  
                ENDIF
              ENDDO !extentIdx
              ALLOCATE(newExtents(1:numberOfExtents),STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate new extents.",err,error,*999)
              newExtents(1:numberOfExtents)=extent(1:numberOfExtents)
              CALL MOVE_ALLOC(newExtents,generatedMesh%polarMesh%polarExtent)
            ELSE
              localError="The specified number of extents of "//TRIM(NumberToVString(SIZE(extent,1),"*",err,error))// &
                & " is invalid. The number of extents must be >= "//TRIM(NumberToVString(minimumNumberOfExtents,"*",err,error))// &
                & " for a polar generated mesh."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ELSE
            CALL FlagError("Regular generated mesh is not associated.",err,error,*999)
          ENDIF
        CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
          IF(ASSOCIATED(generatedMesh%cylinderMesh)) THEN
            IF(generatedMesh%meshDimension==2) THEN
              minimumNumberOfExtents=2
            ELSE IF(generatedMesh%meshDimension==3) THEN
              minimumNumberOfExtents=3
            ELSE
              localError="The mesh dimension of "//TRIM(NumberToVString(generatedMesh%meshDimension,"*",err,error))// &
                & " is invalid for a cylindrical mesh. The mesh dimension must be either 2 or 3."
              CALL FlagError(localError,err,error,*999)             
            ENDIF
            numberOfExtents=SIZE(extent,1)
            IF(numberOfExtents>=minimumNumberOfExtents) THEN
              !Check the extents
              !Check the length of the cylinder
              IF(extent(1)<ZERO_TOLERANCE) THEN
                localError="The first specified extent of "//TRIM(NumberToVString(extent(1),"*",err,error))// &
                  & " is invalid as the length of a cylinder. The extent must be > 0."
                CALL FlagError(localError,err,error,*999)
              ENDIF
              !Check radii
              previousExtent=0.0_DP
              DO extentIdx=2,numberOfExtents
                IF(extent(extentIdx)>ZERO_TOLERANCE) THEN
                  IF(extent(extentIdx)>previousExtent) THEN
                    previousExtent=extent(extentIdx)+ZERO_TOLERANCE
                  ELSE
                    localError="The specified extent index of "// &
                      & TRIM(NumberToVString(extentIdx,"*",err,error))// &
                      & " which has an extent of "//TRIM(NumberToVString(extent(extentIdx),"*",err,error))// &
                      & " is invalid as a radius of a cylinder. The extent needs to be bigger than the previous extent of "// &
                      & TRIM(NumberToVString(previousExtent,"*",err,error))//"."
                    CALL FlagError(localError,err,error,*999)
                  ENDIF
                ELSE
                  localError="The specified extent index of "// &
                    & TRIM(NumberToVString(extentIdx,"*",err,error))// &
                    & " which has an extent of "//TRIM(NumberToVString(extent(extentIdx),"*",err,error))// &
                    & " is invalid as a radius of a cylinder. The extent needs to be > zero."
                  CALL FlagError(localError,err,error,*999)                  
                ENDIF
              ENDDO !extentIdx
              ALLOCATE(newExtents(1:numberOfExtents),STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate new extents.",err,error,*999)
              newExtents(1:numberOfExtents)=extent(1:numberOfExtents)
              CALL MOVE_ALLOC(newExtents,generatedMesh%cylinderMesh%polarExtent)
            ELSE
              localError="The extent size of "//TRIM(NumberToVString(SIZE(extent,1),"*",err,error))// &
                & " is invalid. The extent size must >= the mesh dimension of "// &
                & TRIM(NumberToVString(meshDimension,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ELSE
            CALL FlagError("Cylinder generated mesh is not associated.",err,error,*999)
          ENDIF
        CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
          IF(ASSOCIATED(generatedMesh%ellipsoidMesh)) THEN
            IF(generatedMesh%meshDimension==2) THEN
              IF(generatedMesh%ellipsoidMesh%closed) THEN
                minimumNumberOfExtents=2
                shortAxisExtentStart=2
              ELSE
                minimumNumberOfExtents=3
                shortAxisExtentStart=3
               ENDIF
            ELSE IF(generatedMesh%meshDimension==3) THEN
              IF(generatedMesh%ellipsoidMesh%closed) THEN
                minimumNumberOfExtents=3
                shortAxisExtentStart=2
              ELSE
                minimumNumberOfExtents=4
                shortAxisExtentStart=3
              ENDIF
            ELSE
              localError="The mesh dimension of "//TRIM(NumberToVString(generatedMesh%meshDimension,"*",err,error))// &
                & " is invalid for an ellipsoid mesh. The mesh dimension must be either 2 or 3."
              CALL FlagError(localError,err,error,*999)
              ENDIF
            ENDIF
            numberOfExtents=SIZE(extent,1)
            IF(numberOfExtents>=minimumNumberOfExtents) THEN
              !Check the extents
              !Check the size of the long axis of the ellipsoid
              IF(extent(1)<ZERO_TOLERANCE) THEN
                localError="The first specified extent of "//TRIM(NumberToVString(extent(1),"*",err,error))// &
                  & " is invalid as the size of the long axis of an ellipsoid. The extent must be > 0."
                CALL FlagError(localError,err,error,*999)
              ENDIF
              !If the ellipsoid is open check the opening angle
              IF(.NOT.generatedMesh%ellipsoidMesh%closed) THEN
                IF((extent(2)<ZERO_TOLERANCE).OR.(extent(2)>=PI)) THEN
                  localError="The second specified extent of "//TRIM(NumberToVString(extent(2),"*",err,error))// &
                    & " is invalid as the opening angle of an open ellipsoid. The extent must be > 0 and < PI."
                  CALL FlagError(localError,err,error,*999)
                ENDIF
              ENDIF
              !Check the ellipsoid short axis radii
              previousExtent=0.0_DP
              DO extentIdx=shortAxisExtentStart,numberOfExtents
                IF(extent(extentIdx)>ZERO_TOLERANCE) THEN
                  IF(extent(extentIdx)>previousExtent) THEN
                    previousExtent=extent(extentIdx)+ZERO_TOLERANCE
                  ELSE
                    localError="The specified extent index of "// &
                      & TRIM(NumberToVString(extentIdx,"*",err,error))// &
                      & " which has an extent of "//TRIM(NumberToVString(extent(extentIdx),"*",err,error))// &
                      & " is invalid as a short axis radius of an ellipsoid. The extent needs to be bigger "// &
                      & "than the previous extent of "//TRIM(NumberToVString(previousExtent,"*",err,error))//"."
                    CALL FlagError(localError,err,error,*999)
                  ENDIF
                ELSE
                  localError="The specified extent index of "// &
                    & TRIM(NumberToVString(extentIdx,"*",err,error))// &
                    & " which has an extent of "//TRIM(NumberToVString(extent(extentIdx),"*",err,error))// &
                    & " is invalid as a short axis radius of an ellipsoid. The extent needs to be > zero."
                  CALL FlagError(localError,err,error,*999)                  
                ENDIF
              ENDDO !extentIdx
              ALLOCATE(newExtents(1:numberOfExtents),STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate new extents.",err,error,*999)
              newExtents(1:numberOfExtents)=extent(1:numberOfExtents)
              CALL MOVE_ALLOC(newExtents,generatedMesh%ellipsoidMesh%polarExtent)
            ELSE
              localError="The specified number of extents of "//TRIM(NumberToVString(SIZE(extent,1),"*",err,error))// &
                & " is invalid. The number of extents must be >= "//TRIM(NumberToVString(minimumNumberOfExtents,"*",err,error))// &
                & " for an ellipsoid generated mesh."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ELSE
            CALL FlagError("Ellipsoid generated mesh is not associated.",err,error,*999)
=======
          IF(SIZE(EXTENT,1)==COORDINATE_DIMENSION) THEN
            IF(ASSOCIATED(GENERATED_MESH%REGULAR_MESH)) THEN
              IF(L2NORM(EXTENT)>ZERO_TOLERANCE) THEN
                IF(ALLOCATED(GENERATED_MESH%REGULAR_MESH%MAXIMUM_EXTENT)) &
                    & DEALLOCATE(GENERATED_MESH%REGULAR_MESH%MAXIMUM_EXTENT)
                ALLOCATE(GENERATED_MESH%REGULAR_MESH%MAXIMUM_EXTENT(COORDINATE_DIMENSION),STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate maximum extent.",ERR,ERROR,*999)
                GENERATED_MESH%REGULAR_MESH%MAXIMUM_EXTENT=EXTENT
              ELSE
                CALL FlagError("The norm of the mesh extent is zero.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Regular generated mesh is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The extent size of "//TRIM(NUMBER_TO_VSTRING(SIZE(EXTENT,1),"*",ERR,ERROR))// &
                & " is invalid. The extent size must match the coordinate system dimension of "// &
                & TRIM(NUMBER_TO_VSTRING(COORDINATE_DIMENSION,"*",ERR,ERROR))//"."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        CASE(GENERATED_MESH_POLAR_MESH_TYPE)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
          IF(SIZE(EXTENT,1)==COORDINATE_DIMENSION) THEN
            IF(ASSOCIATED(GENERATED_MESH%CYLINDER_MESH)) THEN
              ALLOCATE(GENERATED_MESH%CYLINDER_MESH%CYLINDER_EXTENT(SIZE(EXTENT)),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate maximum extent.",ERR,ERROR,*999)
              GENERATED_MESH%CYLINDER_MESH%CYLINDER_EXTENT=EXTENT
            ELSE
              CALL FlagError("Cylinder generated mesh is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The extent size of "//TRIM(NUMBER_TO_VSTRING(SIZE(EXTENT,1),"*",ERR,ERROR))// &
                & " is invalid. The extent size must match the coordinate system dimension of "// &
                & TRIM(NUMBER_TO_VSTRING(COORDINATE_DIMENSION,"*",ERR,ERROR))//"."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
          IF((SIZE(EXTENT,1)-1)==COORDINATE_DIMENSION) THEN
            IF(ASSOCIATED(GENERATED_MESH%ELLIPSOID_MESH)) THEN
              ALLOCATE(GENERATED_MESH%ELLIPSOID_MESH%ELLIPSOID_EXTENT(SIZE(EXTENT)),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate maximum extent.",ERR,ERROR,*999)
              GENERATED_MESH%ELLIPSOID_MESH%ELLIPSOID_EXTENT=EXTENT
            ELSE
              CALL FlagError("Ellipsoid generated mesh is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The extent size of "//TRIM(NUMBER_TO_VSTRING(SIZE(EXTENT,1),"*",ERR,ERROR))// &
                & " is invalid. The extent size must be equal one plus the coordinate system dimension of "// &
                & TRIM(NUMBER_TO_VSTRING(COORDINATE_DIMENSION,"*",ERR,ERROR))//"."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
          ENDIF
        CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE DEFAULT
<<<<<<< HEAD
          localError="The generated mesh mesh type of "//TRIM(NumberToVString(generatedMesh%generatedType,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FlagError("Generated mesh is not associated.",err,error,*999)
    ENDIF

    CALL Exits("GeneratedMeshExtentSet")
    RETURN
999 CALL Errors("GeneratedMeshExtentSet",err,error)
    CALL Exits("GeneratedMeshExtentSet")
=======
          LOCAL_ERROR="The generated mesh mesh type of "// &
              & TRIM(NUMBER_TO_VSTRING(GENERATED_MESH%GENERATED_TYPE,"*",ERR,ERROR))// &
              & " is invalid."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FlagError("Generated mesh is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("GENERATED_MESH_EXTENT_SET")
    RETURN
999 ERRORSEXITS("GENERATED_MESH_EXTENT_SET",ERR,ERROR)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
    RETURN 1
    
  END SUBROUTINE GeneratedMeshExtentSet

  !
  !================================================================================================================================
  !

  !>Get one of the surfaces of a generated mesh. \see OPENCMISS::CMISSGeneratedMeshSurfaceGet
  SUBROUTINE GeneratedMeshSurfaceGet(generatedMesh,meshComponent,surfaceType,surfaceNodes,normalXi,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh to get the type of
    INTEGER(INTG), INTENT(IN) :: meshComponent !<The mesh component to get the surface for.
    INTEGER(INTG), INTENT(IN) :: surfaceType !<The surface you are interested in
    INTEGER(INTG), ALLOCATABLE, INTENT(OUT) :: surfaceNodes (:) !<The nodes on the specified surface
    INTEGER(INTG), INTENT(OUT) :: normalXi !<The normal outward pointing xi direction
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

<<<<<<< HEAD
    CALL Enters("GeneratedMeshSurfaceGet",err,error,*999)
=======
    ENTERS("GENERATED_MESH_SURFACE_GET",ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b

    IF(ASSOCIATED(generatedMesh)) THEN
      SELECT CASE(generatedMesh%generatedType)
      CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
        CALL GeneratedMeshRegularSurfaceGet(generatedMesh%regularMesh,meshComponent,surfaceType,surfaceNodes,normalXi, &
          & err,error,*999)
      CASE(GENERATED_MESH_POLAR_MESH_TYPE)
<<<<<<< HEAD
        CALL GeneratedMeshPolarSurfaceGet(generatedMesh%polarMesh,meshComponent,surfaceType,surfaceNodes,normalXi, &
          & err,error,*999)
=======
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
      CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
        CALL GeneratedMeshCylinderSurfaceGet(generatedMesh%cylinderMesh,meshComponent,surfaceType,surfaceNodes,normalXi, &
          & err,error,*999)
      CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
        CALL GeneratedMeshEllipsoidSurfaceGet(generatedMesh%ellipsoidMesh,meshComponent,surfaceType,surfaceNodes,normalXi, &
          & err,error,*999)
      CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="The generated mesh mesh type of "//TRIM(NumberToVString(generatedMesh%generatedType,"*",err,error))// &
            & " is invalid."
<<<<<<< HEAD
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Generated mesh is not associated",err,error,*999)
    ENDIF

    CALL Exits("GeneratedMeshSurfaceGet")
    RETURN
999 CALL Errors("GeneratedMeshSurfaceGet",err,error)
    CALL Exits("GeneratedMeshSurfaceGet")
    RETURN 1
    
  END SUBROUTINE GeneratedMeshSurfaceGet
=======
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Generated mesh is not associated",ERR,ERROR,*999)
    ENDIF

    EXITS("GENERATED_MESH_SURFACE_GET")
    RETURN
999 ERRORSEXITS("GENERATED_MESH_SURFACE_GET",ERR,ERROR)
    RETURN
  END SUBROUTINE GENERATED_MESH_SURFACE_GET
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b

  !
  !================================================================================================================================
  !

  !>Finalises a generated mesh and dellocates all memory.
  SUBROUTINE GeneratedMeshFinalise(generatedMesh,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh to finalise.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

<<<<<<< HEAD
    CALL Enters("GeneratedMeshFinalise",err,error,*999)

    IF(ASSOCIATED(generatedMesh)) THEN
      CALL GeneratedMeshRegularFinalise(generatedMesh%regularMesh,err,error,*999)
      CALL GeneratedMeshPolarFinalise(generatedMesh%polarMesh,err,error,*999)
      CALL GeneratedMeshCylinderFinalise(generatedMesh%cylinderMesh,err,error,*999)
      CALL GeneratedMeshEllipsoidFinalise(generatedMesh%ellipsoidMesh,err,error,*999)
      CALL GeneratedMeshFractalTreeFinalise(generatedMesh%fractalTreeMesh,err,error,*999)
      IF(ALLOCATED(generatedMesh%bases)) DEALLOCATE(generatedMesh%bases)
      IF(ALLOCATED(generatedMesh%origin)) DEALLOCATE(generatedMesh%origin)
      IF(ALLOCATED(generatedMesh%baseVectors)) DEALLOCATE(generatedMesh%baseVectors)
      DEALLOCATE(generatedMesh)
    ENDIF

    CALL Exits("GeneratedMeshFinalise")
    RETURN
999 CALL Errors("GeneratedMeshFinalise",err,error)
    CALL Exits("GeneratedMeshFinalise")
=======
    ENTERS("GENERATED_MESH_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      CALL GENERATED_MESH_REGULAR_FINALISE(GENERATED_MESH%REGULAR_MESH,ERR,ERROR,*999)
      CALL GENERATED_MESH_CYLINDER_FINALISE(GENERATED_MESH%CYLINDER_MESH,ERR,ERROR,*999)
      CALL GENERATED_MESH_ELLIPSOID_FINALISE(GENERATED_MESH%ELLIPSOID_MESH,ERR,ERROR,*999)
      DEALLOCATE(GENERATED_MESH)
    ENDIF

    EXITS("GENERATED_MESH_FINALISE")
    RETURN
999 ERRORSEXITS("GENERATED_MESH_FINALISE",ERR,ERROR)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
    RETURN 1
    
  END SUBROUTINE GeneratedMeshFinalise

  !
  !================================================================================================================================
  !

  !>Initialises a generated mesh.
  SUBROUTINE GeneratedMeshInitialise(coordinateSystem,generatedMesh,err,error,*)

    !Argument variables
    TYPE(COORDINATE_SYSTEM), POINTER :: coordinateSystem !<A point to the coordinate system to initialise the generated mesh in.
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: coordinateDimension,coordinateIdx,dummyErr
    TYPE(VARYING_STRING) :: dummyError

<<<<<<< HEAD
    CALL Enters("GeneratedMeshInitialise",err,error,*998)

    IF(ASSOCIATED(coordinateSystem)) THEN
      IF(ASSOCIATED(generatedMesh)) THEN
        CALL FlagError("Generated mesh is already associated.",err,error,*998)
      ELSE
        NULLIFY(generatedMesh)
        CALL COORDINATE_SYSTEM_DIMENSION_GET(coordinateSystem,coordinateDimension,err,error,*999)
        ALLOCATE(generatedMesh,STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate generated mesh.",err,error,*999)
        generatedMesh%userNumber=0
        generatedMesh%globalNumber=0
        generatedMesh%generatedMeshFinished=.FALSE.
        NULLIFY(generatedMesh%region)
        NULLIFY(generatedMesh%INTERFACE)
        generatedMesh%generatedType=0
        NULLIFY(generatedMesh%regularMesh)
        NULLIFY(generatedMesh%polarMesh)
        NULLIFY(generatedMesh%cylinderMesh)
        NULLIFY(generatedMesh%ellipsoidMesh)
        NULLIFY(generatedMesh%fractalTreeMesh)
        ALLOCATE(generatedMesh%origin(coordinateDimension),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate origin.",err,error,*998)
        generatedMesh%origin=0.0_DP
        generatedMesh%coordinateDimension=coordinateDimension
        ALLOCATE(generatedMesh%baseVectors(coordinateDimension,coordinateDimension),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate base Vectors.",err,error,*998)
        generatedMesh%baseVectors=0.0_DP
        DO coordinateIdx=1,coordinateDimension
           generatedMesh%baseVectors(coordinateIdx,coordinateIdx)=1.0_DP
        ENDDO !coordinateidx
        generatedMesh%meshDimension=0
        generatedMesh%numberOfMeshComponents=0
        NULLIFY(generatedMesh%mesh)
        !Default to a regular mesh.
        CALL GeneratedMeshRegularInitialise(generatedMesh,err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Coordinate system is not associated.",err,error,*998)
    ENDIF
      
    CALL Exits("GeneratedMeshInitialise")
    RETURN
999 CALL GeneratedMeshFinalise(generatedMesh,dummyErr,dummyError,*998)
998 CALL Errors("GeneratedMeshInitialise",err,error)
    CALL Exits("GeneratedMeshInitialise")
=======
    ENTERS("GENERATED_MESH_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      CALL FlagError("Generated mesh is already associated.",ERR,ERROR,*999)
    ELSE
      ALLOCATE(GENERATED_MESH,STAT=ERR)
      IF(ERR/=0) CALL FlagError("Could not allocate generated mesh.",ERR,ERROR,*999)
      GENERATED_MESH%USER_NUMBER=0
      GENERATED_MESH%GLOBAL_NUMBER=0
      GENERATED_MESH%GENERATED_MESH_FINISHED=.FALSE.
      NULLIFY(GENERATED_MESH%REGION)
      NULLIFY(GENERATED_MESH%INTERFACE)
      GENERATED_MESH%GENERATED_TYPE=0
      NULLIFY(GENERATED_MESH%REGULAR_MESH)
      NULLIFY(GENERATED_MESH%CYLINDER_MESH)
      NULLIFY(GENERATED_MESH%ELLIPSOID_MESH)
      NULLIFY(GENERATED_MESH%MESH)
      !Default to a regular mesh.
      CALL GENERATED_MESH_REGULAR_INITIALISE(GENERATED_MESH,ERR,ERROR,*999)
    ENDIF

    EXITS("GENERATED_MESH_INITIALISE")
    RETURN
999 ERRORSEXITS("GENERATED_MESH_INITIALISE",ERR,ERROR)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
    RETURN 1
    
  END SUBROUTINE GeneratedMeshInitialise

  !
  !================================================================================================================================
  !

  !>Gets the number of elements in a generated mesh.  \see OPENCMISS::CMISSGeneratedMeshNumberOfElementsGet
  SUBROUTINE GeneratedMeshNumberOfElementsGet(generatedMesh,numberOfElements,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh to set the type of
    INTEGER(INTG), INTENT(OUT) :: numberOfElements(:) !<On return, number of elements per mesh dimension axis
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: meshDimension
    TYPE(VARYING_STRING) :: localError

<<<<<<< HEAD
    CALL Enters("GeneratedMeshNumberOfElementsGet",err,error,*999)
=======
    ENTERS("GENERATED_MESH_NUMBER_OF_ELEMENTS_GET",ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b

    IF(ASSOCIATED(generatedMesh)) THEN
      SELECT CASE(generatedMesh%generatedType)
      CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
        IF(ASSOCIATED(generatedMesh%regularMesh)) THEN
          IF(ALLOCATED(generatedMesh%regularMesh%numberOfElementsXi)) THEN
            IF(SIZE(numberOfElements,1)>=generatedMesh%meshDimension) THEN
              numberOfElements(1:generatedMesh%meshDimension)= &
                & generatedMesh%regularMesh%numberOfElementsXi(1:generatedMesh%meshDimension)
            ELSE
              localError="The size of number of elements is too small. The supplied size is "// &
                & TRIM(NumberToVString(SIZE(numberOfElements,1),"*",err,error))//" and it needs to be >= "// &
                & TRIM(NumberToVString(generatedMesh%meshDimension,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ELSE
            CALL FlagError("Number of elements has not yet been set.",err,error,*999)
          ENDIF
        ELSE
<<<<<<< HEAD
          CALL FlagError("Generated mesh regular mesh is not associated.",err,error,*999)
        ENDIF
      CASE(GENERATED_MESH_POLAR_MESH_TYPE)
        IF(ASSOCIATED(generatedMesh%polarMesh)) THEN
          IF(ALLOCATED(generatedMesh%polarMesh%numberOfElementsXi)) THEN
            IF(SIZE(numberOfElements,1)>=SIZE((generatedMesh%polarMesh%numberOfElementsXi,1)) THEN
              numberOfElements(1:SIZE((generatedMesh%polarMesh%numberOfElementsXi,1))= &
                & generatedMesh%polarMesh%numberOfElementsXi(1:SIZE((generatedMesh%polarMesh%numberOfElementsXi,1))
            ELSE
              localError="The size of number of elements is too small. The supplied size is "// &
                & TRIM(NumberToVString(SIZE(numberOfElements,1),"*",err,error))//" and it needs to be >= "// &
                & TRIM(NumberToVString(SIZE(generatedMesh%polarMesh%numberOfElementsXi,1),"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ELSE
            CALL FlagError("Number of elements has not yet been set.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Generated mesh polar mesh is not associated.",err,error,*999)
        ENDIF
=======
          LOCAL_ERROR="The size of NUMBER_OF_ELEMENTS is too small. The supplied size is "// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(NUMBER_OF_ELEMENTS,1),"*",ERR,ERROR))//" and it needs to be >= "// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(GENERATED_MESH%REGULAR_MESH%NUMBER_OF_ELEMENTS_XI,1),"*",ERR,ERROR))//"."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      CASE(GENERATED_MESH_POLAR_MESH_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
      CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
        IF(ASSOCIATED(generatedMesh%cylinderMesh)) THEN
          IF(ALLOCATED(generatedMesh%cylinderMesh%numberOfElementsXi)) THEN
            IF(SIZE(numberOfElements,1)>=SIZE((generatedMesh%cylinderMesh%numberOfElementsXi,1)) THEN
              numberOfElements(1:SIZE((generatedMesh%cylinderMesh%numberOfElementsXi,1))= &
                & generatedMesh%cylinderMesh%numberOfElementsXi(1:SIZE((generatedMesh%cylinderMesh%numberOfElementsXi,1))
            ELSE
              localError="The size of number of elements is too small. The supplied size is "// &
                & TRIM(NumberToVString(SIZE(numberOfElements,1),"*",err,error))//" and it needs to be >= "// &
                & TRIM(NumberToVString(SIZE(generatedMesh%cylinderMesh%numberOfElementsXi,1),"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ELSE
            CALL FlagError("Number of elements has not yet been set.",err,error,*999)
          ENDIF
        ELSE
<<<<<<< HEAD
          CALL FlagError("Generated mesh cylinder mesh is not associated.",err,error,*999)
=======
          LOCAL_ERROR="The size of NUMBER_OF_ELEMENTS is too small. The supplied size is "// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(NUMBER_OF_ELEMENTS,1),"*",ERR,ERROR))//" and it needs to be >= "// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(GENERATED_MESH%CYLINDER_MESH%NUMBER_OF_ELEMENTS_XI,1),"*",ERR,ERROR))//"."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
        ENDIF
      CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
        IF(ASSOCIATED(generatedMesh%ellipsoidMesh)) THEN
          IF(ALLOCATED(generatedMesh%ellipsoidMesh%numberOfElementsXi)) THEN
            IF(SIZE(numberOfElements,1)>=SIZE((generatedMesh%ellipsoidMesh%numberOfElementsXi,1)) THEN
              numberOfElements(1:SIZE((generatedMesh%ellipsoidMesh%numberOfElementsXi,1))= &
                & generatedMesh%ellipsoidMesh%numberOfElementsXi(1:SIZE((generatedMesh%ellipsoidMesh%numberOfElementsXi,1))
            ELSE
              localError="The size of number of elements is too small. The supplied size is "// &
                & TRIM(NumberToVString(SIZE(numberOfElements,1),"*",err,error))//" and it needs to be >= "// &
                & TRIM(NumberToVString(SIZE(generatedMesh%ellipsoidMesh%numberOfElementsXi,1),"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ELSE
            CALL FlagError("Number of elements has not yet been set.",err,error,*999)
          ENDIF
        ELSE
<<<<<<< HEAD
          CALL FlagError("Generated mesh ellipsoid mesh is not associated.",err,error,*999)
=======
          LOCAL_ERROR="The size of NUMBER_OF_ELEMENTS is too small. The supplied size is "// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(NUMBER_OF_ELEMENTS,1),"*",ERR,ERROR))//" and it needs to be >= "// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(GENERATED_MESH%ELLIPSOID_MESH%NUMBER_OF_ELEMENTS_XI,1),"*",ERR,ERROR))//"."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
        ENDIF
      CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="The generated mesh mesh type of "//TRIM(NumberToVString(generatedMesh%generatedType,"*",err,error))// &
          & " is invalid."
<<<<<<< HEAD
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Generated mesh is not associated.",err,error,*999)
    ENDIF

    CALL Exits("GeneratedMeshNumberOfElementsGet")
    RETURN
999 CALL Errors("GeneratedMeshNumberOfElementsGet",err,error)
    CALL Exits("GeneratedMeshNumberOfElementsGet")
=======
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Generated mesh is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("GENERATED_MESH_NUMBER_OF_ELEMENTS_GET")
    RETURN
999 ERRORSEXITS("GENERATED_MESH_NUMBER_OF_ELEMENTS_GET",ERR,ERROR)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
    RETURN 1
    
  END SUBROUTINE GeneratedMeshNumberOfElementsGet

  !
  !================================================================================================================================
  !

  !>Sets/changes the number of elements in a generated mesh. \see OPENCMISS::CMISSGeneratedMeshNumberOfElementsSet
  SUBROUTINE GeneratedMeshNumberOfElementsSet(generatedMesh,numberOfElements,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh to set the type of
    INTEGER(INTG), INTENT(IN) :: numberOfElements(:) !<numberOfElements(xiIdx). The number of elements in the xiIdx'th xi direction to set.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(GeneratedMeshCylinderType), POINTER :: cylinderMesh
    TYPE(GeneratedMeshEllipsoidType), POINTER :: ellipsoidMesh
    TYPE(GeneratedMeshPolarType), POINTER :: polarMesh
    TYPE(GeneratedMeshRegularType), POINTER :: regularMesh
    TYPE(VARYING_STRING) :: localError

<<<<<<< HEAD
    CALL Enters("GeneratedMeshNumberOfElementsSet",err,error,*999)

    IF(ASSOCIATED(generatedMesh)) THEN
      IF(generatedMesh%generatedMeshFinished) THEN
        CALL FlagError("Generated mesh has been finished.",err,error,*999)
      ELSE
        IF(ALLOCATED(generatedMesh%bases)) THEN
          SELECT CASE(generatedMesh%generatedType)
          CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
            regularMesh=>generatedMesh%regularMesh
            IF(ASSOCIATED(regularMesh)) THEN
              IF(SIZE(numberOfElements,1)>=generatedMesh%meshDimension) THEN
                IF(ALL(numberOfElements>0)) THEN
                  regularMesh%numberOfElementsXi(1:generatedMesh%meshDimension)=numberOfElements(1:generatedMesh%meshDimension)
                ELSE
                  CALL FlagError("Must have 1 or more elements in all directions.",err,error,*999)
                ENDIF
              ELSE
                localError="The number of elements xi size of "// &
                  & TRIM(NumberToVString(SIZE(numberOfElements,1),"*",err,error))// &
                  & " is invalid. The number of elements size must be >= the mesh dimension of "// &
                  & TRIM(NumberToVString(generatedMesh%meshDimension,"*",err,error))//"."
                CALL FlagError(localError,err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Regular generated mesh is not associated.",err,error,*999)
            ENDIF
          CASE(GENERATED_MESH_POLAR_MESH_TYPE)
            polarMesh=>generatedMesh%polarMesh
            IF(ASSOCIATED(polarMesh)) THEN
              IF(SIZE(numberOfElements,1)>=generatedMesh%meshDimension) THEN
                IF(ALL(numberOfElements>0)) THEN
                  polarMesh%numberOfElementsXi(1:generatedMesh%meshDimension)=numberOfElements(1:generatedMesh%meshDimension)
                ELSE
                  CALL FlagError("Must have 1 or more elements in all directions.",err,error,*999)
                ENDIF
              ELSE
                localError="The number of elements xi size of "// &
                  & TRIM(NumberToVString(SIZE(numberOfElements,1),"*",err,error))// &
                  & " is invalid. The number of elements size must be >= the mesh dimension of "// &
                  & TRIM(NumberToVString(generatedMesh%meshDimension,"*",err,error))//"."
                CALL FlagError(localError,err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Polar generated mesh is not associated.",err,error,*999)
            ENDIF
          CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
            cylinderMesh=>generatedMesh%cylinderMesh
            IF(ASSOCIATED(cylinderMesh)) THEN
              IF(SIZE(numberOfElements,1)>=generatedMesh%meshDimension) THEN
                IF(ALL(numberOfElements>0)) THEN
                  cylinderMesh%numberOfElementsXi(1:generatedMesh%meshDimension)=numberOfElements(1:generatedMesh%meshDimension)
                ELSE
                  CALL FlagError("Must have 1 or more elements in all directions.",err,error,*999)
                ENDIF
              ELSE
                localError="The number of elements xi size of "// &
                  & TRIM(NumberToVString(SIZE(numberOfElements,1),"*",err,error))// &
                  & " is invalid. The number of elements size must be >= the mesh dimension of "// &
                  & TRIM(NumberToVString(generatedMesh%meshDimension,"*",err,error))//"."
                CALL FlagError(localError,err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Cylinder generated mesh is not associated.",err,error,*999)
            ENDIF
          CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
            ellipsoidMesh=>generatedMesh%ellipsoidMesh
            IF(ASSOCIATED(ellipsoidMesh)) THEN
              IF(SIZE(numberOfElements,1)>=generatedMesh%meshDimension) THEN
                IF(ALL(numberOfElements>0)) THEN
                  ellipsoidMesh%numberOfElementsXi(1:generatedMesh%meshDimension)=numberOfElements(1:generatedMesh%meshDimension)
                ELSE
                  CALL FlagError("Must have 1 or more elements in all directions.",err,error,*999)
                ENDIF
              ELSE
                localError="The number of elements xi size of "// &
                  & TRIM(NumberToVString(SIZE(numberOfElements,1),"*",err,error))// &
                  & " is invalid. The number of elements size must be >= the mesh dimension of "// &
                  & TRIM(NumberToVString(generatedMesh%meshDimension,"*",err,error))//"."
                CALL FlagError(localError,err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Ellipsoid generated mesh is not associated.",err,error,*999)
            ENDIF
          CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The generated mesh mesh type of "//TRIM(NumberToVString(generatedMesh%generatedType,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ELSE
          CALL FlagError("Must set the generated mesh basis before setting the number of elements.",err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Generated mesh is not associated.",err,error,*999)
    ENDIF

    CALL Exits("GeneratedMeshNumberOfElementsSet")
    RETURN
999 CALL Errors("GeneratedMeshNumberOfElementsSet",err,error)
    CALL Exits("GeneratedMeshNumberOfElementsSet")
=======
    ENTERS("GENERATED_MESH_NUMBER_OF_ELEMENTS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      IF(GENERATED_MESH%GENERATED_MESH_FINISHED) THEN
        CALL FlagError("Generated mesh has been finished.",ERR,ERROR,*999)
      ELSE
        SELECT CASE(GENERATED_MESH%GENERATED_TYPE)
        CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
          REGULAR_MESH=>GENERATED_MESH%REGULAR_MESH
          IF(ASSOCIATED(REGULAR_MESH)) THEN
            IF(ALLOCATED(REGULAR_MESH%BASES)) THEN
              BASIS=>REGULAR_MESH%BASES(1)%PTR !Number of xi will be the same for all bases
              IF(ASSOCIATED(BASIS)) THEN
                IF(SIZE(NUMBER_OF_ELEMENTS_XI,1)==BASIS%NUMBER_OF_XI) THEN
                  IF(ALL(NUMBER_OF_ELEMENTS_XI>0)) THEN
                    IF(ALLOCATED(REGULAR_MESH%NUMBER_OF_ELEMENTS_XI)) DEALLOCATE(REGULAR_MESH%NUMBER_OF_ELEMENTS_XI)
                    ALLOCATE(REGULAR_MESH%NUMBER_OF_ELEMENTS_XI(BASIS%NUMBER_OF_XI),STAT=ERR)
                    IF(ERR/=0) CALL FlagError("Could not allocate number of elements xi.",ERR,ERROR,*999)
                    REGULAR_MESH%NUMBER_OF_ELEMENTS_XI(1:BASIS%NUMBER_OF_XI)=NUMBER_OF_ELEMENTS_XI(1:BASIS%NUMBER_OF_XI)
                  ELSE
                    CALL FlagError("Must have 1 or more elements in all directions.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The number of elements xi size of "// &
                    & TRIM(NUMBER_TO_VSTRING(SIZE(NUMBER_OF_ELEMENTS_XI,1),"*",ERR,ERROR))// &
                    & " is invalid. The number of elements xi size must match the basis number of xi dimensions of "// &
                    & TRIM(NUMBER_TO_VSTRING(BASIS%NUMBER_OF_XI,"*",ERR,ERROR))//"."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Must set the generated mesh basis before setting the number of elements.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Must set the generated mesh basis before setting the number of elements.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Regular generated mesh is not associated.",ERR,ERROR,*999)
          ENDIF
        CASE(GENERATED_MESH_POLAR_MESH_TYPE)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
          IF(ASSOCIATED(GENERATED_MESH%CYLINDER_MESH)) THEN
            ALLOCATE(GENERATED_MESH%CYLINDER_MESH%NUMBER_OF_ELEMENTS_XI(SIZE(NUMBER_OF_ELEMENTS_XI)),STAT=ERR)
            GENERATED_MESH%CYLINDER_MESH%NUMBER_OF_ELEMENTS_XI=NUMBER_OF_ELEMENTS_XI
          ELSE
            CALL FlagError("Cylinder generated mesh is not associated.",ERR,ERROR,*999)
          ENDIF
        CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
          IF(ASSOCIATED(GENERATED_MESH%ELLIPSOID_MESH)) THEN
            ALLOCATE(GENERATED_MESH%ELLIPSOID_MESH%NUMBER_OF_ELEMENTS_XI(SIZE(NUMBER_OF_ELEMENTS_XI)),STAT=ERR)
            GENERATED_MESH%ELLIPSOID_MESH%NUMBER_OF_ELEMENTS_XI=NUMBER_OF_ELEMENTS_XI
          ELSE
            CALL FlagError("Ellipsoid generated mesh is not associated.",ERR,ERROR,*999)
          ENDIF
        CASE DEFAULT
          LOCAL_ERROR="The generated mesh mesh type of "//TRIM(NUMBER_TO_VSTRING(GENERATED_MESH%GENERATED_TYPE,"*",ERR,ERROR))// &
            & " is invalid."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FlagError("Generated mesh is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("GENERATED_MESH_NUMBER_OF_ELEMENTS_SET")
    RETURN
999 ERRORSEXITS("GENERATED_MESH_NUMBER_OF_ELEMENTS_SET",ERR,ERROR)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
    RETURN 1
    
  END SUBROUTINE GeneratedMeshNumberOfElementsSet

  !
  !================================================================================================================================
  !

  !>Get the origin of a generated mesh. \see OPENCMISS::CMISSGeneratedMeshOriginGet
  SUBROUTINE GeneratedMeshOriginGet(generatedMesh,origin,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh to get the type of
    REAL(DP), INTENT(OUT) :: origin(:) !<On return, the origin coordinate for each axis
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

<<<<<<< HEAD
    CALL Enters("GeneratedMeshOriginGet",err,error,*999)

    IF(ASSOCIATED(generatedMesh)) THEN
      IF(SIZE(origin,1)>=generatedMesh%coordinateDimension) THEN
        origin(1:generatedMesh%coordinateDimension)=generatedMesh%origin(1:generatedMesh%coordinateDimension)
      ELSE
        localError="The size of origin is too small. The supplied size is "// &
          & TRIM(NumberToVString(SIZE(origin,1),"*",err,error))// &
          & " and it needs to be >= the coordinate dimension of "// &
          & TRIM(NumberToVString(generatedMesh%coordinateDimension,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Generated mesh is not associated",err,error,*999)
    ENDIF

    CALL Exits("GeneratedMeshOriginGet")
    RETURN
999 CALL Errors("GeneratedMeshOriginGet",err,error)
    CALL Exits("GeneratedMeshOriginGet")
=======
    ENTERS("GENERATED_MESH_ORIGIN_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      SELECT CASE(GENERATED_MESH%GENERATED_TYPE)
      CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
        IF(SIZE(ORIGIN,1)>=SIZE(GENERATED_MESH%REGULAR_MESH%ORIGIN,1)) THEN
          ORIGIN=GENERATED_MESH%REGULAR_MESH%ORIGIN
        ELSE
          LOCAL_ERROR="The size of ORIGIN is too small. The supplied size is "// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(ORIGIN,1),"*",ERR,ERROR))//" and it needs to be >= "// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(GENERATED_MESH%REGULAR_MESH%ORIGIN,1),"*",ERR,ERROR))//"."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      CASE(GENERATED_MESH_POLAR_MESH_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
        IF(SIZE(ORIGIN,1)>=SIZE(GENERATED_MESH%CYLINDER_MESH%ORIGIN,1)) THEN
          ORIGIN=GENERATED_MESH%CYLINDER_MESH%ORIGIN
        ELSE
          LOCAL_ERROR="The size of ORIGIN is too small. The supplied size is "// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(ORIGIN,1),"*",ERR,ERROR))//" and it needs to be 3."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
        IF(SIZE(ORIGIN,1)>=SIZE(GENERATED_MESH%ELLIPSOID_MESH%ORIGIN,1)) THEN
          ORIGIN=GENERATED_MESH%ELLIPSOID_MESH%ORIGIN
        ELSE
          LOCAL_ERROR="The size of ORIGIN is too small. The supplied size is "// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(ORIGIN,1),"*",ERR,ERROR))//" and it needs to be 3."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      CASE DEFAULT
        LOCAL_ERROR="The generated mesh mesh type of "//TRIM(NUMBER_TO_VSTRING(GENERATED_MESH%GENERATED_TYPE,"*",ERR,ERROR))// &
          & " is invalid."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Generated mesh is not associated",ERR,ERROR,*999)
    ENDIF

    EXITS("GENERATED_MESH_ORIGIN_GET")
    RETURN
999 ERRORSEXITS("GENERATED_MESH_ORIGIN_GET",ERR,ERROR)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
    RETURN 1
    
  END SUBROUTINE GeneratedMeshOriginGet

  !
  !================================================================================================================================
  !

  !>Sets/changes the origin of a generated mesh. \see OPENCMISS::CMISSGeneratedMeshOriginSet
  SUBROUTINE GeneratedMeshOriginSet(generatedMesh,origin,err,error,*)
    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh to set the type of
    REAL(DP), INTENT(IN) :: origin(:) !<origin(coordinateIdx). The coordinateIdx'th coordinate origin for the generated mesh
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

<<<<<<< HEAD
    CALL Enters("GeneratedMeshOriginSet",err,error,*999)

    IF(ASSOCIATED(generatedMesh)) THEN
      IF(generatedMesh%generatedMeshFinished) THEN
        CALL FlagError("Generated mesh has been finished.",err,error,*999)
      ELSE
        IF(SIZE(origin,1)>=generatedMesh%coordinateDimension) THEN
          generatedMesh%origin(1:generatedMesh%coordinateDimension)=origin(1:generatedMesh%coordinateDimension)
        ELSE
          localError="The origin size of "//TRIM(NumberToVString(SIZE(origin,1),"*",err,error))// &
            & " is invalid. The origin size must be >= the coordinate system dimension of "// &
            & TRIM(NumberToVString(generatedMesh%coordinateDimension,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Generated mesh is not associated.",err,error,*999)
    ENDIF

    CALL Exits("GeneratedMeshOriginSet")
    RETURN
999 CALL Errors("GeneratedMeshOriginSet",err,error)
    CALL Exits("GeneratedMeshOriginSet")
=======
    ENTERS("GENERATED_MESH_ORIGIN_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      IF(GENERATED_MESH%GENERATED_MESH_FINISHED) THEN
        CALL FlagError("Generated mesh has been finished.",ERR,ERROR,*999)
      ELSE
        NULLIFY(COORDINATE_SYSTEM)
        CALL GENERATED_MESH_COORDINATE_SYSTEM_GET(GENERATED_MESH,COORDINATE_SYSTEM,ERR,ERROR,*999)
        CALL COORDINATE_SYSTEM_DIMENSION_GET(COORDINATE_SYSTEM,COORDINATE_DIMENSION,ERR,ERROR,*999)
        IF(SIZE(ORIGIN,1)==COORDINATE_DIMENSION) THEN
          SELECT CASE(GENERATED_MESH%GENERATED_TYPE)
          CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
            IF(ASSOCIATED(GENERATED_MESH%REGULAR_MESH)) THEN
              IF(.NOT.ALLOCATED(GENERATED_MESH%REGULAR_MESH%ORIGIN)) THEN
                ALLOCATE(GENERATED_MESH%REGULAR_MESH%ORIGIN(SIZE(ORIGIN)),STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate origin.",ERR,ERROR,*999)
              ENDIF
              GENERATED_MESH%REGULAR_MESH%ORIGIN=ORIGIN
            ELSE
              CALL FlagError("Regular generated mesh is not associated.",ERR,ERROR,*999)
            ENDIF
          CASE(GENERATED_MESH_POLAR_MESH_TYPE)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
            IF(ASSOCIATED(GENERATED_MESH%CYLINDER_MESH)) THEN
              IF(SIZE(ORIGIN,1)==3) THEN
                IF(.NOT.ALLOCATED(GENERATED_MESH%CYLINDER_MESH%ORIGIN)) THEN
                  ALLOCATE(GENERATED_MESH%CYLINDER_MESH%ORIGIN(SIZE(ORIGIN)),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate origin.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Cylinder generated mesh is only supported for 3D.",ERR,ERROR,*999)
              ENDIF
              GENERATED_MESH%CYLINDER_MESH%ORIGIN=ORIGIN
            ELSE
              CALL FlagError("Cylinder generated mesh is not associated.",ERR,ERROR,*999)
            END IF
          CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
            IF(ASSOCIATED(GENERATED_MESH%ELLIPSOID_MESH)) THEN
              IF(SIZE(ORIGIN,1)==3) THEN
                IF(.NOT.ALLOCATED(GENERATED_MESH%ELLIPSOID_MESH%ORIGIN)) THEN
                  ALLOCATE(GENERATED_MESH%ELLIPSOID_MESH%ORIGIN(SIZE(ORIGIN)),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate origin.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Ellipsoid generated mesh is only supported for 3D.",ERR,ERROR,*999)
              ENDIF
              GENERATED_MESH%ELLIPSOID_MESH%ORIGIN=ORIGIN
            ELSE
              CALL FlagError("Ellipsoid generated mesh is not associated.",ERR,ERROR,*999)
            END IF
          CASE DEFAULT
            LOCAL_ERROR="The generated mesh mesh type of "//TRIM(NUMBER_TO_VSTRING(GENERATED_MESH%GENERATED_TYPE,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          LOCAL_ERROR="The origin size of "//TRIM(NUMBER_TO_VSTRING(SIZE(ORIGIN,1),"*",ERR,ERROR))// &
            & " is invalid. The extent size must match the coordinate system dimension of "// &
            & TRIM(NUMBER_TO_VSTRING(COORDINATE_DIMENSION,"*",ERR,ERROR))//"."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Generated mesh is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("GENERATED_MESH_ORIGIN_SET")
    RETURN
999 ERRORSEXITS("GENERATED_MESH_ORIGIN_SET",ERR,ERROR)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
    RETURN 1
    
  END SUBROUTINE GeneratedMeshOriginSet

  !
  !================================================================================================================================
  !

  !>Finishes the creation of a regular generated mesh.
  SUBROUTINE GeneratedMeshRegularCreateFinish(regularMesh,meshUserNumber,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshRegularType), POINTER :: regularMesh !<A pointer to the regular generated mesh to finish.
    INTEGER(INTG), INTENT(IN) :: meshUserNumber !<The user number for the mesh to generate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
<<<<<<< HEAD
    INTEGER(INTG) :: elementFactor,gridElementIdx,gridNumberOfElements,elementIdx,elementIdx1,elementIdx2, &
      & elementIdx3,nodeIdx,localNodeIdx,localNodeIdx1,localNodeIdx2,localNodeIdx3,coordinateType,coordinateDimension, &
      & numberOfElementsXi(3),totalNumberOfNodesXi(3),totalNumberOfNodes,numberOfCornerNodes, &
      & totalNumberOfElements,xiIdx,numberOfBases,basisIdx,basisNumberOfNodes
    INTEGER(INTG), ALLOCATABLE :: elementNodes(:),elementNodesUserNumbers(:)
    TYPE(BASIS_TYPE), POINTER :: basis
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: coordinateSystem
    TYPE(GeneratedMeshType), POINTER :: generatedMesh
    TYPE(INTERFACE_TYPE), POINTER :: interface
    TYPE(MeshComponentElementsType), POINTER :: meshElements
    TYPE(NODES_TYPE), POINTER :: nodes
    TYPE(REGION_TYPE), POINTER :: region
    TYPE(VARYING_STRING) :: localError

    NULLIFY(coordinateSystem)
    
    CALL Enters("GeneratedMeshRegularCreateFinish",err,error,*999)

    IF(ASSOCIATED(regularMesh)) THEN
      generatedMesh=>regularMesh%generatedMesh
      IF(ASSOCIATED(generatedMesh)) THEN
        CALL GeneratedMeshCoordinateSystemGet(generatedMesh,coordinateSystem,err,error,*999)
        CALL COORDINATE_SYSTEM_TYPE_GET(coordinateSystem,coordinateType,err,error,*999)
        region=>generatedMesh%region
        INTERFACE=>generatedMesh%INTERFACE
        SELECT CASE(coordinateType)
=======
    INTEGER(INTG) :: coordinate_idx,COUNT,ELEMENT_FACTOR,grid_ne,GRID_NUMBER_OF_ELEMENTS,ni,ne,ne1,ne2,ne3,nn,nn1,nn2,nn3,np, &
      & NUMBER_OF_ELEMENTS_XI(3),TOTAL_NUMBER_OF_NODES_XI(3),TOTAL_NUMBER_OF_NODES,NUMBER_OF_CORNER_NODES, &
      & TOTAL_NUMBER_OF_ELEMENTS,xi_idx,NUM_BASES,basis_idx,BASIS_NUMBER_OF_NODES(10)
    INTEGER(INTG), ALLOCATABLE :: ELEMENT_NODES(:),ELEMENT_NODES_USER_NUMBERS(:)
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
    TYPE(GENERATED_MESH_REGULAR_TYPE), POINTER :: REGULAR_MESH
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE
    TYPE(MeshElementsType), POINTER :: MESH_ELEMENTS
    TYPE(NODES_TYPE), POINTER :: NODES
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("GENERATED_MESH_REGULAR_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      NULLIFY(COORDINATE_SYSTEM)
      CALL GENERATED_MESH_COORDINATE_SYSTEM_GET(GENERATED_MESH,COORDINATE_SYSTEM,ERR,ERROR,*999)
      REGION=>GENERATED_MESH%REGION
      INTERFACE=>GENERATED_MESH%INTERFACE
      REGULAR_MESH=>GENERATED_MESH%REGULAR_MESH
      IF(ASSOCIATED(REGULAR_MESH)) THEN
        SELECT CASE(COORDINATE_SYSTEM%TYPE)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
        CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
          IF(ALLOCATED(generatedMesh%bases)) THEN
            !Use first basis to get the basis type
            basis=>regularMesh%bases(1)%ptr
            SELECT CASE(basis%TYPE)
            CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE,BASIS_SIMPLEX_TYPE)
<<<<<<< HEAD
=======
              IF(.NOT.ALL(BASIS%COLLAPSED_XI==BASIS_NOT_COLLAPSED)) &
                & CALL FlagError("Degenerate (collapsed) basis not implemented.",ERR,ERROR,*999)
              !Determine the coordinate system and create the regular mesh for that system
              REGULAR_MESH%COORDINATE_DIMENSION=COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS
              REGULAR_MESH%MESH_DIMENSION=BASIS%NUMBER_OF_XI
              IF(.NOT.ALLOCATED(REGULAR_MESH%ORIGIN)) THEN
                ALLOCATE(REGULAR_MESH%ORIGIN(REGULAR_MESH%COORDINATE_DIMENSION),STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate origin.",ERR,ERROR,*999)
                REGULAR_MESH%ORIGIN=0.0_DP
              ENDIF
              IF(.NOT.ALLOCATED(REGULAR_MESH%MAXIMUM_EXTENT)) THEN
                ALLOCATE(REGULAR_MESH%MAXIMUM_EXTENT(REGULAR_MESH%COORDINATE_DIMENSION),STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate maximum extent.",ERR,ERROR,*999)
                REGULAR_MESH%MAXIMUM_EXTENT=1.0_DP
              ENDIF
              IF(.NOT.ALLOCATED(REGULAR_MESH%NUMBER_OF_ELEMENTS_XI)) THEN
                ALLOCATE(REGULAR_MESH%NUMBER_OF_ELEMENTS_XI(REGULAR_MESH%MESH_DIMENSION),STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate number of elements xi.",ERR,ERROR,*999)
                REGULAR_MESH%NUMBER_OF_ELEMENTS_XI=1
              ENDIF
              IF(ALLOCATED(REGULAR_MESH%BASE_VECTORS)) THEN
                !!TODO: Check base vectors
              ELSE
                !Calculate base vectors
                ALLOCATE(REGULAR_MESH%BASE_VECTORS(REGULAR_MESH%COORDINATE_DIMENSION,REGULAR_MESH%MESH_DIMENSION),STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate number of elements xi.",ERR,ERROR,*999)
                REGULAR_MESH%BASE_VECTORS=0.0_DP
                IF(REGULAR_MESH%MESH_DIMENSION==1) THEN
                  !The base vector is just the extent vector
                  REGULAR_MESH%BASE_VECTORS(:,1)=REGULAR_MESH%MAXIMUM_EXTENT
                ELSE
                  IF(REGULAR_MESH%MESH_DIMENSION<REGULAR_MESH%COORDINATE_DIMENSION) THEN
                    !Find the first number of mesh dimensions for which the extent is non-zero.
                    COUNT=0
                    coordinate_idx=1
                    DO xi_idx=1,REGULAR_MESH%MESH_DIMENSION
                      DO WHILE(ABS(REGULAR_MESH%MAXIMUM_EXTENT(coordinate_idx))<=ZERO_TOLERANCE)
                        coordinate_idx=coordinate_idx+1
                      ENDDO
                      REGULAR_MESH%BASE_VECTORS(coordinate_idx,xi_idx)=REGULAR_MESH%MAXIMUM_EXTENT(coordinate_idx)
                      coordinate_idx=coordinate_idx+1
                      COUNT=COUNT+1
                    ENDDO !xi_idx
                    IF(COUNT/=REGULAR_MESH%MESH_DIMENSION)  &
                      & CALL FlagError("Invalid mesh extent. There number of non-zero components is < the mesh dimension.", &
                      & ERR,ERROR,*999)
                  ELSE IF(REGULAR_MESH%MESH_DIMENSION==REGULAR_MESH%COORDINATE_DIMENSION) THEN
                    !The default base vectors are aligned with the coordinate vectors
                    DO coordinate_idx=1,REGULAR_MESH%COORDINATE_DIMENSION
                      REGULAR_MESH%BASE_VECTORS(coordinate_idx,coordinate_idx)=REGULAR_MESH%MAXIMUM_EXTENT(coordinate_idx)
                    ENDDO !coordinate_idx
                  ELSE
                    CALL FlagError("The mesh dimension is greater than the coordinate dimension.",ERR,ERROR,*999)
                  ENDIF
                ENDIF
              ENDIF
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
              !Calculate the sizes of a regular grid of elements with the appropriate number of basis nodes in each dimension of
              !the grid element
              totalNumberOfNodes=1
              gridNumberOfElements=1
              numberOfElementsXi=1
              numberOfBases=SIZE(regularMesh%bases)
              DO xiIdx=1,generatedMesh%meshDimension
                !Set total number of nodes to corner nodes only
                totalNumberOfNodes=totalNumberOfNodes*(regularMesh%numberOfElementsXi(xiIdx)+1)
                numberOfElementsXi(xiIdx)=regularMesh%numberOfElementsXi(xiIdx)
                gridNumberOfElements=gridNumberOfElements*regularMesh%numberOfElementsXi(xiIdx)
              ENDDO !xiIdx
              numberOfCornerNodes=totalNumberOfNodes
              !Add extra nodes for each basis
              !Will end up with some duplicate nodes if bases have the same interpolation in one direction
<<<<<<< HEAD
              DO basisIdx=1,numberOfBases
                basis=>regularMesh%bases(basisIdx)%ptr
                basisNumberOfNodes=1
                DO xiIdx=1,generatedMesh%meshDimension
                  basisNumberofNodes=BasisNumberOfNodes*((basis%NUMBER_OF_NODES_XIC(xiIdx)-1)* &
                    & regularMesh%numberOfElementsXi(xiIdx)+1)
                ENDDO !xiIdx
                totalNumberOfNodes=totalNumberOfNodes+basisNumberOfNodes-numberOfCornerNodes
=======
              BASIS_NUMBER_OF_NODES=0
              DO basis_idx=1,NUM_BASES
                BASIS=>REGULAR_MESH%BASES(basis_idx)%PTR
                BASIS_NUMBER_OF_NODES(basis_idx)=1
                DO ni=1,REGULAR_MESH%MESH_DIMENSION
                  BASIS_NUMBER_OF_NODES(basis_idx)=BASIS_NUMBER_OF_NODES(basis_idx)*((BASIS%NUMBER_OF_NODES_XIC(ni)-1)* &
                      & REGULAR_MESH%NUMBER_OF_ELEMENTS_XI(ni)+1)
                ENDDO
                BASIS_NUMBER_OF_NODES(basis_idx)=TOTAL_NUMBER_OF_NODES+BASIS_NUMBER_OF_NODES(basis_idx)-NUMBER_OF_CORNER_NODES
                ! BASIS_NUMBER_OF_NODES=1
                ! DO ni=1,REGULAR_MESH%MESH_DIMENSION
                !   BASIS_NUMBER_OF_NODES=BASIS_NUMBER_OF_NODES*((BASIS%NUMBER_OF_NODES_XIC(ni)-1)* &
                !       & REGULAR_MESH%NUMBER_OF_ELEMENTS_XI(ni)+1)
                ! ENDDO
                ! TOTAL_NUMBER_OF_NODES=TOTAL_NUMBER_OF_NODES+BASIS_NUMBER_OF_NODES-NUMBER_OF_CORNER_NODES
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
              ENDDO
              TOTAL_NUMBER_OF_NODES=MAXVAL(BASIS_NUMBER_OF_NODES)
              !Compute the element factor i.e., the number of sub elements each grid element will be split into.
              IF(basis%TYPE==BASIS_LAGRANGE_HERMITE_TP_TYPE) THEN
                elementFactor=1
              ELSE
                SELECT CASE(generatedMesh%meshDimension)
                CASE(1)
                  elementFactor=1
                CASE(2)
                  elementFactor=2
                CASE(3)
                  elementFactor=6
                CASE DEFAULT
<<<<<<< HEAD
                  localError="The mesh dimension dimension of "// &
                    & TRIM(NumberToVString(generatedMesh%meshDimension,"*",err,error))//" is invalid."
                  CALL FlagError(localError,err,error,*999)
=======
                  LOCAL_ERROR="The mesh dimension dimension of "// &
                    & TRIM(NUMBER_TO_VSTRING(REGULAR_MESH%MESH_DIMENSION,"*",ERR,ERROR))//" is invalid."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
                END SELECT
              ENDIF
              totalNumberOfElements=elementFactor*gridNumberOfElements
              !Create the default node set
              NULLIFY(nodes)
              IF(ASSOCIATED(region)) THEN
                CALL NODES_CREATE_START(region,totalNumberOfNodes,nodes,err,error,*999)
              ELSE
                CALL NODES_CREATE_START(INTERFACE,totalNumberOfNodes,nodes,err,error,*999)
              ENDIF
              !Finish the nodes creation
              CALL NODES_CREATE_FINISH(nodes,err,error,*999)
              !Create the mesh
              IF(ASSOCIATED(region)) THEN
                CALL MESH_CREATE_START(meshUserNumber,region,generatedMesh%meshDimension,generatedMesh%mesh,err,error,*999)
              ELSE
                CALL MESH_CREATE_START(meshUserNumber,interface,generatedMesh%meshDimension,generatedMesh%mesh,err,error,*999)
              ENDIF
              !Set the number of mesh components
              CALL MESH_NUMBER_OF_COMPONENTS_SET(generatedMesh%mesh,generatedMesh%numberOfMeshComponents,err,error,*999)
              !Create the elements
              CALL MESH_NUMBER_OF_ELEMENTS_SET(generatedMesh%mesh,totalNumberOfElements,err,error,*999)
              DO meshComponentIdx=1,generatedMesh%numberOfMeshComponents
                basisIdx=meshComponentIdx
                basis=>regularMesh%bases(basisIdx)%ptr
                !Get number of nodes in each xi direction for this basis
                DO xiIdx=1,generatedMesh%meshDimension
                  totalNumberOfNodesXi(xiIdx)=(basis%NUMBER_OF_NODES_XIC(xiIdx)-1)*regularMesh%numberOfElementsXi(xiIdx)+1
                ENDDO
                NULLIFY(meshElements)
                CALL MESH_TOPOLOGY_ELEMENTS_CREATE_START(generatedMesh%mesh,meshComponentIdx,basis,meshElements,err,error,*999)
                !Set the elements for the regular mesh
<<<<<<< HEAD
                IF(ALLOCATED(elementNodes)) DEALLOCATE(elementNodes)
                ALLOCATE(elementNodes(basis%NUMBER_OF_NODES),STAT=err)
                IF (ALLOCATED(elementNodesUserNumbers)) DEALLOCATE(elementNodesUserNumbers)
                ALLOCATE(elementNodesUserNumbers(basis%NUMBER_OF_NODES),STAT=err)
                IF(err/=0) CALL FlagError("Could not allocate element nodes.",err,error,*999)
=======
                IF (ALLOCATED(ELEMENT_NODES)) DEALLOCATE(ELEMENT_NODES)
                ALLOCATE(ELEMENT_NODES(BASIS%NUMBER_OF_NODES),STAT=ERR)
                IF (ALLOCATED(ELEMENT_NODES_USER_NUMBERS)) DEALLOCATE(ELEMENT_NODES_USER_NUMBERS)
                ALLOCATE(ELEMENT_NODES_USER_NUMBERS(BASIS%NUMBER_OF_NODES),STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate element nodes.",ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
                !Step in the xi(3) direction
                DO elementIdx3=1,numberOfElementsXi(3)+1
                  DO elementIdx2=1,numberOfElementsXi(2)+1
                    DO elementIdx1=1,numberOfElementsXi(1)+1
                      IF(basis%NUMBER_OF_XI<3.OR.elementIdx3<=numberOfElementsXi(3)) THEN
                        IF(basis%NUMBER_OF_XI<2.OR.elementIdx2<=numberOfElementsXi(2)) THEN
                          IF(elementIdx1<=numberOfElementsXi(1)) THEN
                            gridElementIdx=elementIdx1
                            nodeIdx=1+(elementIdx1-1)*(basis%NUMBER_OF_NODES_XIC(1)-1)
                            IF(basis%NUMBER_OF_XI>1) THEN
                              gridElementIdx=gridElementIdx+(elementIdx2-1)*numberOfElementsXi(1)
                              nodeIdx=nodeIdx+(elementIdx2-1)*totalNumberOfNodesXi(1)*(basis%NUMBER_OF_NODES_XIC(2)-1)
                              IF(basis%NUMBER_OF_XI>2) THEN
                                gridElementIdx=gridElementIdx+(elementIdx3-1)*numberOfElementsXi(1)*numberOfElementsXi(2)
                                nodeIdx=nodeIdx+(elementIdx3-1)*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)* &
                                  & (basis%NUMBER_OF_NODES_XIC(3)-1)
                              ENDIF
                            ENDIF
                            IF(BASIS%TYPE==BASIS_LAGRANGE_HERMITE_TP_TYPE) THEN
                              !Lagrange Hermite TP elements
                              elementIdx=gridElementIdx
                              localNodeIdx=0
                              DO localNodeIdx1=1,basis%NUMBER_OF_NODES_XIC(1)
                                localNodeIdx=localNodeIdx+1
                                elementNodes(localNodeIdx)=nodeIdx+(localNodeIdx1-1)
                              ENDDO !localNodeIdx1
                              IF(basis%NUMBER_OF_XI>1) THEN
                                DO localNodeIdx2=2,basis%NUMBER_OF_NODES_XIC(2)
                                  DO localNodeIdx1=1,basis%NUMBER_OF_NODES_XIC(1)
                                    localNodeIdx=localNodeIdx+1
                                    elementNodes(localNodeIdx)=nodeIdx+(localNodeidx1-1)+(localNodeIdx2-1)*totalNumberOfNodesXi(1)
                                  ENDDO !localNodeIdx1
                                ENDDO !localNodeIdx2
                                IF(basis%NUMBER_OF_XI>2) THEN
                                  DO localNodeIdx3=2,basis%NUMBER_OF_NODES_XIC(3)
                                    DO localnodeIdx2=1,basis%NUMBER_OF_NODES_XIC(2)
                                      DO localNodeIdx1=1,basis%NUMBER_OF_NODES_XIC(1)
                                        localNodeIdx=localNodeidx+1
                                        elementNodes(localNodeIdx)=nodeIdx+(localNodeIdx1-1)+ &
                                          & (localnodeIdx2-1)* totalNumberOfNodesXi(1)+ &
                                          & (localNodeIdx3-1)*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                                      ENDDO !localNodeIdx1
                                    ENDDO !localnodeIdx2
                                  ENDDO !localNodeIdx3
                                ENDIF
                              ENDIF
<<<<<<< HEAD
                              CALL GeneratedMeshRegularComponentNodesToUserNumbers(regularMesh%generatedMesh, &
                                & basisIdx,elementNodes,elementNodesUserNumbers,err,error,*999)
                              CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(elementIdx,meshElements,elementNodesUserNumbers, &
                                & err,error,*999)
=======
                              CALL GeneratedMesh_RegularComponentNodesToUserNumbers(REGULAR_MESH%GENERATED_MESH, &
                                & basis_idx,ELEMENT_NODES,ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                              CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(ne,MESH_ELEMENTS, &
                                & ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
                            ELSE
                              !Simplex elements
                              SELECT CASE(basis%NUMBER_OF_XI)
                              CASE(1)
                                !Line element
                                elementIdx=gridElementIdx
                                localNodeIdx=0
                                DO localNodeIdx1=1,basis%NUMBER_OF_NODES_XIC(1)
                                  localNodeIdx=localNodeIdx+1
                                  elementNodes(localNodeIdx)=nodeIdx+(localNodeIdx1-1)
                                ENDDO !localNodeIdx1
                                CALL GeneratedMeshComponentNodesToUserNumbers(regularMesh%generatedMesh,basisIdx,elementNodes, &
                                  & elementNodesUserNumbers,err,error,*999)
                                CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(elementIdx,meshElements,elementNodesUserNumbers, &
                                  & err,error,*999)
                              CASE(2)
                                !Triangular element
                                !Break the grid square element into 2 triangles. The 2 triangles are
                                !Element 1: vertices {(0,0);(1,0);(1,1)}
                                !Element 2: vertices {(0,0);(1,1);(0,1)}
                                SELECT CASE(basis%INTERPOLATION_ORDER(1))
                                CASE(BASIS_LINEAR_INTERPOLATION_ORDER)
                                  !First sub-element
                                  elementIdx=(gridElementIdx-1)*elementFactor+1
                                  elementNodes(1)=nodeIdx
                                  elementNodes(2)=nodeIdx+1
                                  elementNodes(3)=nodeIdx+1+totalNumberOfNodesXi(1)
                                  CALL GeneratedMeshComponentNodesToUserNumbers(regularMesh%generatedMesh,basisIdx,elementNodes, &
                                    & elementNodesUserNumbers,err,error,*999)
                                  CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(elementIdx,meshElements,elementNodesUserNumbers, &
                                    & err,error,*999)
                                  !Second sub-element
                                  elementIdx=(gridElementIdx-1)*elementFactor+2
                                  elementNodes(1)=nodeIdx
                                  elementNodes(2)=nodeIdx+1+totalNumberOfNodesXi(1)
                                  elementNodes(3)=nodeIdx+totalNumberOfNodesXi(1)
                                  CALL GeneratedMeshComponentNodesToUserNumbers(regularMesh%generatedMesh,basisIdx,elementNodes, &
                                    & elementNodesUserNumbers,err,error,*999)
                                  CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(elementIdx,meshElements,elementNodesUserNumbers, &
                                    & err,error,*999)
                                CASE(BASIS_QUADRATIC_INTERPOLATION_ORDER)
                                  !First sub-element
                                  elementIdx=(gridElementIdx-1)*elementFactor+1
                                  elementNodes(1)=nodeIdx
                                  elementNodes(2)=nodeIdx+2
                                  elementNodes(3)=nodeIdx+2+2*totalNumberOfNodesXi(1)
                                  elementNodes(4)=nodeIdx+1
                                  elementNodes(5)=nodeIdx+2+totalNumberOfNodesXi(1)
                                  elementNodes(6)=nodeIdx+1+totalNumberOfNodesXi(1)
                                  CALL GeneratedMeshComponentNodesToUserNumbers(regularMesh%generatedMesh,basisIdx,elementNodes, &
                                    & elementNodesUserNumbers,err,error,*999)
                                  CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(elementIdx,meshElements,elementNodesUserNumbers, &
                                    & err,error,*999)
                                  !Second sub-element
                                  elementIdx=(gridElementIdx-1)*elementFactor+2
                                  elementNodes(1)=nodeIdx
                                  elementNodes(2)=nodeIdx+2+2*totalNumberOfNodesXi(1)
                                  elementNodes(3)=nodeIdx+2*totalNumberOfNodesXi(1)
                                  elementNodes(4)=nodeIdx+1+totalNumberOfNodesXi(1)
                                  elementNodes(5)=nodeIdx+1+2*totalNumberOfNodesXi(1)
                                  elementNodes(6)=nodeIdx+totalNumberOfNodesXi(1)
                                  CALL GeneratedMeshComponentNodesToUserNumbers(regularMesh%generatedMesh,basisIdx,elementNodes, &
                                    & elementNodesUserNumbers,err,error,*999)
                                  CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(elementIdx,meshElements,elementNodesUserNumbers, &
                                    & err,error,*999)
                                CASE(BASIS_CUBIC_INTERPOLATION_ORDER)
                                  !First sub-element
                                  elementIdx=(gridElementIdx-1)*elementFactor+1
                                  elementNodes(1)=nodeIdx
                                  elementNodes(2)=nodeIdx+3
                                  elementNodes(3)=nodeIdx+3+3*totalNumberOfNodesXi(1)
                                  elementNodes(4)=nodeIdx+1
                                  elementNodes(5)=nodeIdx+2
                                  elementNodes(6)=nodeIdx+3+totalNumberOfNodesXi(1)
                                  elementNodes(7)=nodeIdx+3+2*totalNumberOfNodesXi(1)
                                  elementNodes(8)=nodeIdx+2+2*totalNumberOfNodesXi(1)
                                  elementNodes(9)=nodeIdx+1+totalNumberOfNodesXi(1)
                                  elementNodes(10)=nodeIdx+2+totalNumberOfNodesXi(1)
                                  CALL GeneratedMeshComponentNodesToUserNumbers(regularMesh%generatedMesh,basisIdx,elementNodes, &
                                    & elementNodesUserNumbers,err,error,*999)
                                  CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(elementIdx,meshElements,elementNodesUserNumbers, &
                                    & err,error,*999)
                                  !Second sub-element
                                  elementIdx=(gridElementIdx-1)*elementFactor+2
                                  elementNodes(1)=nodeIdx
                                  elementNodes(2)=nodeIdx+3+3*totalNumberOfNodesXi(1)
                                  elementNodes(3)=nodeIdx+3*totalNumberOfNodesXi(1)
                                  elementNodes(4)=nodeIdx+2+2*totalNumberOfNodesXi(1)
                                  elementNodes(5)=nodeIdx+1+totalNumberOfNodesXi(1)
                                  elementNodes(6)=nodeIdx+2+3*totalNumberOfNodesXi(1)
                                  elementNodes(7)=nodeIdx+1+3*totalNumberOfNodesXi(1)
                                  elementNodes(8)=nodeIdx+totalNumberOfNodesXi(1)
                                  elementNodes(9)=nodeIdx+2*totalNumberOfNodesXi(1)
                                  elementNodes(10)=nodeIdx+1+2*totalNumberOfNodesXi(1)
                                  CALL GeneratedMeshComponentNodesToUserNumbers(regularMesh%generatedMesh,basisIdx,elementNodes, &
                                    & elementNodesUserNumbers,err,error,*999)
                                  CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(elementIdx,meshElements,elementNodesUserNumbers, &
                                    & err,error,*999)
                                CASE DEFAULT
                                  localError="The simplex basis interpolation order of "// &
                                    & TRIM(NumberToVString(basis%INTERPOLATION_ORDER(1),"*",err,error))// &
                                    & " is invalid."
<<<<<<< HEAD
                                  CALL FlagError(localError,err,error,*999)
=======
                                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
                                END SELECT
                              CASE(3)
                                !Tetrahedra element
                                !Break the grid cube element into 6 tetrahedra (so that we have a break down the main diagonal of
                                !the cube in order to allow for the middle node in quadratics to be included). The 6 tetrahedra are
                                !Element 1: vertices {(0,0,0);(1,0,0);(1,1,0);(1,1,1)}
                                !Element 2: vertices {(0,0,0);(1,1,0);(0,1,0);(1,1,1)}
                                !Element 3: vertices {(0,0,0);(1,0,1);(1,0,0);(1,1,1)}
                                !Element 4: vertices {(0,0,0);(0,0,1);(1,0,1);(1,1,1)}
                                !Element 5: vertices {(0,0,0);(0,1,0);(0,1,1);(1,1,1)}
                                !Element 6: vertices {(0,0,0);(0,1,1);(0,0,1);(1,1,1)}
                                SELECT CASE(basis%INTERPOLATION_ORDER(1))
                                CASE(BASIS_LINEAR_INTERPOLATION_ORDER)
                                  !First sub-element
                                  elementIdx=(gridElementIdx-1)*elementFactor+1
                                  elementNodes(1)=nodeIdx
                                  elementNodes(2)=nodeIdx+1
                                  elementNodes(3)=nodeIdx+1+totalNumberOfNodesXi(1)
                                  elementNodes(4)=nodeIdx+1+totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  CALL GeneratedMeshComponentNodesToUserNumbers(regularMesh%generatedMesh,basisIdx,elementNodes, &
                                    & elementNodesUserNumbers,err,error,*999)
                                  CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(elementIdx,meshElements,elementNodesUserNumbers, &
                                    & err,error,*999)
                                  !Second sub-element
                                  elementIdx=(gridElementIdx-1)*elementFactor+2
                                  elementNodes(1)=nodeIdx
                                  elementNodes(2)=nodeIdx+1+totalNumberOfNodesXi(1)
                                  elementNodes(3)=nodeIdx+totalNumberOfNodesXi(1)
                                  elementNodes(4)=nodeIdx+1+totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  CALL GeneratedMeshComponentNodesToUserNumbers(regularMesh%generatedMesh,basisIdx,elementNodes, &
                                    & elementNodesUserNumbers,err,error,*999)
                                  CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(elementIdx,meshElements,elementNodesUserNumbers, &
                                    & err,error,*999)
                                  !Third sub-element
                                  elementIdx=(gridElementIdx-1)*elementFactor+3
                                  elementNodes(1)=nodeIdx
                                  elementNodes(2)=nodeIdx+1+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                                  elementNodes(3)=nodeIdx+1
                                  elementNodes(4)=nodeIdx+1+totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  CALL GeneratedMeshComponentNodesToUserNumbers(regularMesh%generatedMesh,basisIdx,elementNodes, &
                                    & elementNodesUserNumbers,err,error,*999)
                                  CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(elementIdx,meshElements,elementNodesUserNumbers, &
                                    & err,error,*999)
                                  !Fourth sub-element
                                  elementIdx=(gridElementIdx-1)*elementFactor+4
                                  elementNodes(1)=nodeIdx
                                  elementNodes(2)=nodeIdx+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                                  elementNodes(3)=nodeIdx+1+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                                  elementNodes(4)=nodeIdx+1+totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  CALL GeneratedMeshComponentNodesToUserNumbers(regularMesh%generatedMesh,basisIdx,elementNodes, &
                                    & elementNodesUserNumbers,err,error,*999)
                                  CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(elementIdx,meshElements,elementNodesUserNumbers, &
                                    & err,error,*999)
                                  !Fifth sub-element
                                  elementIdx=(gridElementIdx-1)*elementFactor+5
                                  elementNodes(1)=nodeIdx
                                  elementNodes(2)=nodeIdx+totalNumberOfNodesXi(1)
                                  elementNodes(3)=nodeIdx+totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(4)=nodeIdx+1+totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  CALL GeneratedMeshComponentNodesToUserNumbers(regularMesh%generatedMesh,basisIdx,elementNodes, &
                                    & elementNodesUserNumbers,err,error,*999)
                                  CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(elementIdx,meshElements,elementNodesUserNumbers, &
                                    & err,error,*999)
                                  !Sixth sub-element
                                  elementIdx=(gridElementIdx-1)*elementFactor+6
                                  elementNodes(1)=nodeIdx
                                  elementNodes(2)=nodeIdx+totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(3)=nodeIdx+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                                  elementNodes(4)=nodeIdx+1+totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  CALL GeneratedMeshComponentNodesToUserNumbers(regularMesh%generatedMesh,basisIdx,elementNodes, &
                                    & elementNodesUserNumbers,err,error,*999)
                                  CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(elementIdx,meshElements,elementNodesUserNumbers, &
                                    & err,error,*999)
                                CASE(BASIS_QUADRATIC_INTERPOLATION_ORDER)
                                  !First sub-element
                                  elementIdx=(gridElementIdx-1)*elementFactor+1
                                  elementNodes(1)=nodeIdx
                                  elementNodes(2)=nodeIdx+2
                                  elementNodes(3)=nodeIdx+2+2*totalNumberOfNodesXi(1)
                                  elementNodes(4)=nodeIdx+2+2*totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(5)=nodeIdx+1
                                  elementNodes(6)=nodeIdx+1+totalNumberOfNodesXi(1)
                                  elementNodes(7)=nodeIdx+1+totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(8)=nodeIdx+2+totalNumberOfNodesXi(1)
                                  elementNodes(9)=nodeIdx+2+2*totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(10)=nodeIdx+2+totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  CALL GeneratedMeshComponentNodesToUserNumbers(regularMesh%generatedMesh,basisIdx,elementNodes, &
                                    & elementNodesUserNumbers,err,error,*999)
                                  CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(elementIdx,meshElements,elementNodesUserNumbers, &
                                    & err,error,*999)
                                  !Second sub-element
                                  elementIdx=(gridElementIdx-1)*elementFactor+2
                                  elementNodes(1)=nodeIdx
                                  elementNodes(2)=nodeIdx+2+2*totalNumberOfNodesXi(1)
                                  elementNodes(3)=nodeIdx+2*totalNumberOfNodesXi(1)
                                  elementNodes(4)=nodeIdx+2+2*totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(5)=nodeIdx+1+totalNumberOfNodesXi(1)
                                  elementNodes(6)=nodeIdx+totalNumberOfNodesXi(1)
                                  elementNodes(7)=nodeIdx+1+totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(8)=nodeIdx+1+2*totalNumberOfNodesXi(1)
                                  elementNodes(9)=nodeIdx+1+2*totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(10)=nodeIdx+2+2*totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  CALL GeneratedMeshComponentNodesToUserNumbers(regularMesh%generatedMesh,basisIdx,elementNodes, &
                                    & elementNodesUserNumbers,err,error,*999)
                                  CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(elementIdx,meshElements,elementNodesUserNumbers, &
                                    & err,error,*999)
                                  !Third sub-element
                                  elementIdx=(gridElementIdx-1)*elementFactor+3
                                  elementNodes(1)=nodeIdx
                                  elementNodes(2)=nodeIdx+2+2*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                                  elementNodes(3)=nodeIdx+2
                                  elementNodes(4)=nodeIdx+2+2*totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(5)=nodeIdx+1+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                                  elementNodes(6)=nodeIdx+1
                                  elementNodes(7)=nodeIdx+1+totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(8)=nodeIdx+2+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                                  elementNodes(9)=nodeIdx+2+totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(10)=nodeIdx+2+totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  CALL GeneratedMeshComponentNodesToUserNumbers(regularMesh%generatedMesh,basisIdx,elementNodes, &
                                    & elementNodesUserNumbers,err,error,*999)
                                  CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(elementIdx,meshElements,elementNodesUserNumbers, &
                                    & err,error,*999)
                                  !Fourth sub-element
                                  elementIdx=(gridElementIdx-1)*elementFactor+4
                                  elementNodes(1)=nodeIdx
                                  elementNodes(2)=nodeIdx+2*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                                  elementNodes(3)=nodeIdx+2+2*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                                  elementNodes(4)=nodeIdx+2+2*totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(5)=nodeIdx+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                                  elementNodes(6)=nodeIdx+1+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                                  elementNodes(7)=nodeIdx+1+totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(8)=nodeIdx+1+2*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                                  elementNodes(9)=nodeIdx+2+totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(10)=nodeIdx+1+totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  CALL GeneratedMeshComponentNodesToUserNumbers(regularMesh%generatedMesh,basisIdx,elementNodes, &
                                    & elementNodesUserNumbers,err,error,*999)
                                  CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(elementIdx,meshElements,elementNodesUserNumbers, &
                                    & err,error,*999)
                                  !Fifth sub-element
                                  elementIdx=(gridElementIdx-1)*elementFactor+5
                                  elementNodes(1)=nodeIdx
                                  elementNodes(2)=nodeIdx+2*totalNumberOfNodesXi(1)
                                  elementNodes(3)=nodeIdx+2*totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(4)=nodeIdx+2+2*totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(5)=nodeIdx+totalNumberOfNodesXi(1)
                                  elementNodes(6)=nodeIdx+totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(7)=nodeIdx+1+totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(8)=nodeIdx+2*totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(9)=nodeIdx+1+2*totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(10)=nodeIdx+1+2*totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  CALL GeneratedMeshComponentNodesToUserNumbers(regularMesh%generatedMesh,basisIdx,elementNodes, &
                                    & elementNodesUserNumbers,err,error,*999)
                                  CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(elementIdx,meshElements,elementNodesUserNumbers, &
                                    & err,error,*999)
                                  !Sixth sub-element
                                  elementIdx=(gridElementIdx-1)*elementFactor+6
                                  elementNodes(1)=nodeIdx
                                  elementNodes(2)=nodeIdx+2*totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(3)=nodeIdx+2*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                                  elementNodes(4)=nodeIdx+2+2*totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(5)=nodeIdx+totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(6)=nodeIdx+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                                  elementNodes(7)=nodeIdx+1+totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(8)=nodeIdx+totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(9)=nodeIdx+1+totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(10)=nodeIdx+1+2*totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  CALL GeneratedMeshComponentNodesToUserNumbers(regularMesh%generatedMesh,basisIdx,elementNodes, &
                                    & elementNodesUserNumbers,err,error,*999)
                                  CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(elementIdx,meshElements,elementNodesUserNumbers, &
                                    & err,error,*999)
                                CASE(BASIS_CUBIC_INTERPOLATION_ORDER)
                                  !First sub-element
                                  elementIdx=(gridElementIdx-1)*elementFactor+1
                                  elementNodes(1)=nodeIdx
                                  elementNodes(2)=nodeIdx+3
                                  elementNodes(3)=nodeIdx+3+3*totalNumberOfNodesXi(1)
                                  elementNodes(4)=nodeIdx+3+3*totalNumberOfNodesXi(1)+3*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(5)=nodeIdx+1
                                  elementNodes(6)=nodeIdx+2
                                  elementNodes(7)=nodeIdx+1+totalNumberOfNodesXi(1)
                                  elementNodes(8)=nodeIdx+2+2*totalNumberOfNodesXi(1)
                                  elementNodes(9)=nodeIdx+1+totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(10)=nodeIdx+2+2*totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(11)=nodeIdx+3+totalNumberOfNodesXi(1)
                                  elementNodes(12)=nodeIdx+3+2*totalNumberOfNodesXi(1)
                                  elementNodes(13)=nodeIdx+3+3*totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(14)=nodeIdx+3+3*totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(15)=nodeIdx+3+totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(16)=nodeIdx+3+2*totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(17)=nodeIdx+2+totalNumberOfNodesXi(1)
                                  elementNodes(18)=nodeIdx+2*totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(19)=nodeIdx+2+2*totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(20)=nodeIdx+3+2*totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  CALL GeneratedMeshComponentNodesToUserNumbers(regularMesh%generatedMesh,basisIdx,elementNodes, &
                                    & elementNodesUserNumbers,err,error,*999)
                                  CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(elementIdx,meshElements,elementNodesUserNumbers, &
                                    & err,error,*999)
                                  !Second sub-element
                                  elementIdx=(gridElementIdx-1)*elementFactor+2
                                  elementNodes(1)=nodeIdx
                                  elementNodes(2)=nodeIdx+3+3*totalNumberOfNodesXi(1)
                                  elementNodes(3)=nodeIdx+3*totalNumberOfNodesXi(1)
                                  elementNodes(4)=nodeIdx+3+3*totalNumberOfNodesXi(1)+3*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(5)=nodeIdx+1+totalNumberOfNodesXi(1)
                                  elementNodes(6)=nodeIdx+2+2*totalNumberOfNodesXi(1)
                                  elementNodes(7)=nodeIdx+totalNumberOfNodesXi(1)
                                  elementNodes(8)=nodeIdx+2*totalNumberOfNodesXi(1)
                                  elementNodes(9)=nodeIdx+1+totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(10)=nodeIdx+2+2*totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(11)=nodeIdx+2+3*totalNumberOfNodesXi(1)
                                  elementNodes(12)=nodeIdx+1+3*totalNumberOfNodesXi(1)
                                  elementNodes(13)=nodeIdx+1+3*totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(14)=nodeIdx+2+3*totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(15)=nodeIdx+3+3*totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(16)=nodeIdx+3+3*totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(17)=nodeIdx+1+2*totalNumberOfNodesXi(1)
                                  elementNodes(18)=nodeIdx+2+2*totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(19)=nodeIdx+1+2*totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(20)=nodeIdx+2+3*totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  CALL GeneratedMeshComponentNodesToUserNumbers(regularMesh%generatedMesh,basisIdx,elementNodes, &
                                    & elementNodesUserNumbers,err,error,*999)
                                  CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(elementIdx,meshElements,elementNodesUserNumbers, &
                                    & err,error,*999)
                                  !Third sub-element
                                  elementIdx=(gridElementIdx-1)*elementFactor+3
                                  elementNodes(1)=nodeIdx
                                  elementNodes(2)=nodeIdx+3+3*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                                  elementNodes(3)=nodeIdx+3
                                  elementNodes(4)=nodeIdx+3+3*totalNumberOfNodesXi(1)+3*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(5)=nodeIdx+1+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                                  elementNodes(6)=nodeIdx+2+2*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                                  elementNodes(7)=nodeIdx+1
                                  elementNodes(8)=nodeIdx+2
                                  elementNodes(9)=nodeIdx+1+totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(10)=nodeIdx+2+2*totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(11)=nodeIdx+3+2*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                                  elementNodes(12)=nodeIdx+3+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                                  elementNodes(13)=nodeIdx+3+totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(14)=nodeIdx+3+2*totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(15)=nodeIdx+3+totalNumberOfNodesXi(1)+3*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(16)=nodeIdx+3+2*totalNumberOfNodesXi(1)+3*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(17)=nodeIdx+2+2*totalNumberOfNodesXi(1)+3*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(18)=nodeIdx+2+totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(19)=nodeIdx+2+totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(20)=nodeIdx+3+totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  CALL GeneratedMeshComponentNodesToUserNumbers(regularMesh%generatedMesh,basisIdx,elementNodes, &
                                    & elementNodesUserNumbers,err,error,*999)
                                  CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(elementIdx,meshElements,elementNodesUserNumbers, &
                                    & err,error,*999)
                                  !Fourth sub-element
                                  elementIdx=(gridElementIdx-1)*elementFactor+4
                                  elementNodes(1)=nodeIdx
                                  elementNodes(2)=nodeIdx+3*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                                  elementNodes(3)=nodeIdx+3+3*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                                  elementNodes(4)=nodeIdx+3+3*totalNumberOfNodesXi(1)+3*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(5)=nodeIdx+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                                  elementNodes(6)=nodeIdx+2*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                                  elementNodes(7)=nodeIdx+1+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                                  elementNodes(8)=nodeIdx+2+2*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                                  elementNodes(9)=nodeIdx+1+totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(10)=nodeIdx+2+2*totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(11)=nodeIdx+1+3*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                                  elementNodes(12)=nodeIdx+2+3*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                                  elementNodes(13)=nodeIdx+3+totalNumberOfNodesXi(1)+3*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(14)=nodeIdx+3+2*totalNumberOfNodesXi(1)+3*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(15)=nodeIdx+1+totalNumberOfNodesXi(1)+3*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(16)=nodeIdx+2+2*totalNumberOfNodesXi(1)+3*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(17)=nodeIdx+1+2*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                                  elementNodes(18)=nodeIdx+1+totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(19)=nodeIdx+2+totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(20)=nodeIdx+2+totalNumberOfNodesXi(1)+3*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  CALL GeneratedMeshComponentNodesToUserNumbers(regularMesh%generatedMesh,basisIdx,elementNodes, &
                                    & elementNodesUserNumbers,err,error,*999)
                                  CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(elementIdx,meshElements,elementNodesUserNumbers, &
                                    & err,error,*999)
                                  !Fifth sub-element
                                  elementIdx=(gridElementIdx-1)*elementFactor+5
                                  elementNodes(1)=nodeIdx
                                  elementNodes(2)=nodeIdx+3*totalNumberOfNodesXi(1)
                                  elementNodes(3)=nodeIdx+3*totalNumberOfNodesXi(1)+3*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(4)=nodeIdx+3+3*totalNumberOfNodesXi(1)+3*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(5)=nodeIdx+totalNumberOfNodesXi(1)
                                  elementNodes(6)=nodeIdx+2*totalNumberOfNodesXi(1)
                                  elementNodes(7)=nodeIdx+totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(8)=nodeIdx+2*totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(9)=nodeIdx+1+totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(10)=nodeIdx+2+2*totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(11)=nodeIdx+3*totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(12)=nodeIdx+3*totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(13)=nodeIdx+1+3*totalNumberOfNodesXi(1)+3*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(14)=nodeIdx+2+3*totalNumberOfNodesXi(1)+3*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(15)=nodeIdx+1+3*totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(16)=nodeIdx+2+3*totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(17)=nodeIdx+2*totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(18)=nodeIdx+1+2*totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(19)=nodeIdx+1+2*totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(20)=nodeIdx+1+3*totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  CALL GeneratedMeshComponentNodesToUserNumbers(regularMesh%generatedMesh,basisIdx,elementNodes, &
                                    & elementNodesUserNumbers,err,error,*999)
                                  CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(elementIdx,meshElements,elementNodesUserNumbers, &
                                    & err,error,*999)
                                  !Sixth sub-element
                                  elementIdx=(gridElementIdx-1)*elementFactor+6
                                  elementNodes(1)=nodeIdx
                                  elementNodes(2)=nodeIdx+3*totalNumberOfNodesXi(1)+3*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(3)=nodeIdx+3*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                                  elementNodes(4)=nodeIdx+3+3*totalNumberOfNodesXi(1)+3*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(5)=nodeIdx+totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(6)=nodeIdx+2*totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(7)=nodeIdx+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                                  elementNodes(8)=nodeIdx+2*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                                  elementNodes(9)=nodeIdx+1+totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(10)=nodeIdx+2+2*totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(11)=nodeIdx+2*totalNumberOfNodesXi(1)+3*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(12)=nodeIdx+totalNumberOfNodesXi(1)+3*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(13)=nodeIdx+1+totalNumberOfNodesXi(1)+3*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(14)=nodeIdx+2+2*totalNumberOfNodesXi(1)+3*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(15)=nodeIdx+1+3*totalNumberOfNodesXi(1)+3*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(16)=nodeIdx+2+3*totalNumberOfNodesXi(1)+3*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(17)=nodeIdx+totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(18)=nodeIdx+1+2*totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(19)=nodeIdx+1+totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  elementNodes(20)=nodeIdx+1+2*totalNumberOfNodesXi(1)+3*totalNumberOfNodesXi(1)* &
                                    & totalNumberOfNodesXi(2)
                                  CALL GeneratedMeshComponentNodesToUserNumbers(regularMesh%generatedMesh,basisIdx,elementNodes, &
                                    & elementNodesUserNumbers,err,error,*999)
                                  CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(elementIdx,meshElements,elementNodesUserNumbers, &
                                    & err,error,*999)
                                CASE DEFAULT
                                  localError="The simplex basis interpolation order of "// &
                                    & TRIM(NumberToVString(basis%INTERPOLATION_ORDER(1),"*",err,error))// &
                                    & " is invalid."
<<<<<<< HEAD
                                  CALL FlagError(localError,err,error,*999)
=======
                                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
                                END SELECT
                              CASE DEFAULT
                                localError="The simplex number of xi directions of "// &
                                  & TRIM(NumberToVString(basis%NUMBER_OF_XI,"*",err,error))// &
                                  & " is invalid."
<<<<<<< HEAD
                                CALL FlagError(localError,err,error,*999)
=======
                                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
                              END SELECT
                            ENDIF
                          ENDIF
                        ENDIF
                      ENDIF
                    ENDDO !elementIdx1
                  ENDDO !elementIdx2
                ENDDO !elementIdx3
                CALL MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH(meshElements,err,error,*999)
              ENDDO !basisIdx
              !Finish the mesh
              CALL MESH_CREATE_FINISH(generatedMesh%mesh,err,error,*999)
            CASE DEFAULT
<<<<<<< HEAD
              CALL FlagError("Basis type is either invalid or not implemented.",err,error,*999)
            END SELECT
          ELSE
            CALL FlagError("Bases are not allocated.",err,error,*999)
          ENDIF
        CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
          CALL FlagError("Not implemented.",err,error,*999)
=======
              CALL FlagError("Basis type is either invalid or not implemented.",ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FlagError("Bases are not allocated.",ERR,ERROR,*999)
          ENDIF
        CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
        CASE DEFAULT
          localError="The coordinate system type of "//TRIM(NumberToVString(coordinateType,"*",err,error))// &
            & " is invalid."
<<<<<<< HEAD
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        CALL FlagError("Regular mesh generated Mesh is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Regular mesh is not associated.",err,error,*999)
    ENDIF
    
    IF(ALLOCATED(elementNodes)) DEALLOCATE(elementNodes)
    
    CALL Exits("GeneratedMeshRegularCreateFinish")
    RETURN
999 IF(ALLOCATED(elementNodes)) DEALLOCATE(elementNodes)
    CALL Errors("GeneratedMeshRegularCreateFinish",err,error)
    CALL Exits("GeneratedMeshRegularCreateFinish")
=======
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FlagError("Regular mesh is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Generated Mesh is not associated.",ERR,ERROR,*999)
    ENDIF

    IF(ALLOCATED(ELEMENT_NODES)) DEALLOCATE(ELEMENT_NODES)

    EXITS("GENERATED_MESH_REGULAR_CREATE_FINISH")
    RETURN
    ! TODO invalidate other associations
999 IF(ALLOCATED(ELEMENT_NODES)) DEALLOCATE(ELEMENT_NODES)
    ERRORSEXITS("GENERATED_MESH_REGULAR_CREATE_FINISH",ERR,ERROR)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
    RETURN 1
    
  END SUBROUTINE GeneratedMeshRegularCreateFinish

  !
  !================================================================================================================================
  !

  !>Finish the creation of a polar generated mesh.
  SUBROUTINE GeneratedMeshPolarCreateFinish(polarMesh,meshUserNumber,err,error,*)

    !Argument variables
    TYPE(GeneratedPolarMeshType), POINTER :: polarMesh !<A pointer to the polar generated mesh to finish.
    INTEGER(INTG), INTENT(IN) :: meshUserNumber !<The user number of the mesh to generate.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
<<<<<<< HEAD
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: coordinateSystem
    TYPE(GeneratedMeshType), POINTER :: generatedMesh
    TYPE(MeshComponentElementsType), POINTER :: meshElements
    TYPE(NODES_TYPE), POINTER :: nodes
    TYPE(REGION_TYPE), POINTER :: region

    NULLIFY(coordinateSystem)
    
    CALL Enters("GeneratedMeshPolarCreateFinish",err,error,*999)

    IF(ASSOCIATED(polarMesh)) THEN
      generatedMesh=>polarMesh%generatedMesh
      IF(ASSOCIATED(generatedMesh)) THEN
        IF(ALLOCATED(generatedMesh%bases)) THEN
          CALL GeneratedMeshCoordinateSystemGet(generatedMesh,coordinateSystem,err,error,*999)
          CALL COORDINATE_SYSTEM_TYPE_GET(coordinateSystem,coordinateType,err,error,*999)
          region=>generatedMesh%region
          INTERFACE=>generatedMesh%INTERFACE
          SELECT CASE(coordinateType)
          CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
            meshDimension=generatedMesh%meshDimension
            !Use first basis to get the basis type
            basis=>regularMesh%bases(1)%ptr
            SELECT CASE(basis%TYPE)
            CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE,BASIS_SIMPLEX_TYPE)
              !Calculate the number of elements
              gridNumberOfElements=1
              numberOfElementsXi
            CASE DEFAULT
              CALL FlagError("Basis type is either invalid or not implemented.",err,error,*999)
            END SELECT
          CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The coordinate system type of "//TRIM(NumberToVString(coordinateType,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ELSE
          CALL FlagError("Generated mesh bases are not allocated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Polar mesh generated Mesh is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Polar mesh is not associated.",err,error,*999)
    ENDIF

    CALL Exits("GeneratedMeshPolarCreateFinish")
    RETURN
999 CALL Errors("GeneratedMeshPolarCreateFinish",err,error)
    CALL Exits("GeneratedMeshPolarCreateFinish")
    RETURN 1
    
  END SUBROUTINE GeneratedMeshPolarCreateFinish

  !
  !================================================================================================================================
  !

  !>Start to create the ellipsoid generated mesh type
  SUBROUTINE GeneratedMeshEllipsoidCreateFinish(generatedMesh,meshUserNumber,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh
    INTEGER(INTG), INTENT(IN) :: meshUserNumber !<The user number for the mesh to generate.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: totalNumberOfNodes,totalNumberOfElements,meshDimension,coordinateDimension,coordinateType
    INTEGER(INTG) :: basisNumberOfNodes,cornerNumberOfNodes
    INTEGER(INTG) :: elementIdx1,elementIdx2,elementIdx3,localNodeIdx1,localNodeIdx2,localNodeIdx3,from1,from2,from3, &
      & localNodeIdx,elementIdx,mc
    INTEGER(INTG), ALLOCATABLE :: apexElementNodes(:), wallElementNodes(:)
    INTEGER(INTG), ALLOCATABLE :: apexElementNodesUserNumbers(:), wallElementNodesUserNumbers(:)
    INTEGER(INTG), ALLOCATABLE :: nodeIndices(:,:,:),cornerNodes(:,:,:),elementIndices(:,:,:)
    INTEGER(INTG), ALLOCATABLE :: numberOfElementsXi(:)
    REAL(DP) :: DELTA(3),DELTAi(3)
    TYPE(BASIS_TYPE), POINTER :: basis1,basis2
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: coordinateSystem
    TYPE(GeneratedMeshEllipsoidType), POINTER :: ellipsoidMesh
    TYPE(MeshComponentElementsType), POINTER :: meshElements
    TYPE(NODES_TYPE), POINTER :: nodes
    TYPE(REGION_TYPE), POINTER :: region

    NULLIFY(coordinateSystem)
    
    CALL Enters("GeneratedMeshEllipsoidCreateFinish",err,error,*999)

    IF(ASSOCIATED(generatedMesh)) THEN
      ellipsoidMesh=>generatedMesh%ellipsoidMesh
      IF(ASSOCIATED(ellipsoidMesh)) THEN
        region=>generatedMesh%region
        IF(ASSOCIATED(region)) THEN
          CALL GeneratedMeshCoordinateSystemGet(generatedMesh,coordinateSystem,err,error,*999)
          CALL COORDINATE_SYSTEM_DIMENSION_GET(coordinateSystem,coordinateDimension,err,error,*999)
          CALL COORDINATE_SYSTEM_TYPE_GET(coordinateSystem,coordinateType,err,error,*999)
          SELECT CASE(coordinateType)
          CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
            meshDimension=generatedMesh%meshDimension
            IF(meshDimension==3) THEN ! hard-coded for 3D only
              IF(SIZE(ellipsoidMesh%ellipsoidExtent)==4) THEN
                IF(ALLOCATED(ellipsoidMesh%bases)) THEN
                  IF(MOD(SIZE(ellipsoidMesh%bases),2)==0) THEN
                    ALLOCATE(numberOfElementsXi(SIZE(ellipsoidMesh%numberOfElementsXi)),STAT=err)
                    IF(err/=0) CALL FlagError("Could not allocate number of elements xi.",err,error,*999)
                    numberOfElementsXi(1:SIZE(ellipsoidMesh%numberOfElementsXi))= &
                      & ellipsoidMesh%numberOfElementsXi(1:SIZE(ellipsoidMesh%numberOfElementsXi))
                    !Calculate total number of nodes from all bases and start mesh
                    cornerNumberOfNodes=numberOfElementsXi(1)*(numberOfElementsXi(2)+1)*(numberOfElementsXi(3)+1)- &
                      & (numberOfElementsXi(1)-1)*(numberOfElementsXi(3)+1)
                    totalNumberOfNodes=cornerNumberOfNodes
                    DO mc=1,SIZE(ellipsoidMesh%bases),2
                      basis1=>ellipsoidMesh%bases(mc)%PTR
                      basisNumberOfNodes=numberOfElementsXi(1)*(basis1%NUMBER_OF_NODES_XIC(1)-1)* &
                        & (numberOfElementsXi(2)*(basis1%NUMBER_OF_NODES_XIC(2)-1)+1)* &
                        & (numberOfElementsXi(3)*(basis1%NUMBER_OF_NODES_XIC(3)-1)+1)- &
                        & (numberOfElementsXi(1)*(basis1%NUMBER_OF_NODES_XIC(1)-1)-1)* &
                        & (numberOfElementsXi(3)*(basis1%NUMBER_OF_NODES_XIC(3)-1)+1)
                      totalNumberOfNodes=totalNumberOfNodes+basisNumberOfNodes-cornerNumberOfNodes
                    ENDDO
                    NULLIFY(NODES)
                    CALL NODES_CREATE_START(region,totalNumberOfNodes,NODES,err,error,*999)
                    !Finish the nodes creation
                    CALL NODES_CREATE_FINISH(NODES,err,error,*999)
                    !Create the mesh
                    CALL MESH_CREATE_START(meshUserNumber,generatedMesh%region, &
                      & SIZE(numberOfElementsXi,1), generatedMesh%MESH,err,error,*999)
                    !Create the elements
                    CALL MESH_NUMBER_OF_COMPONENTS_SET(generatedMesh%MESH,SIZE(ellipsoidMesh%bases)/2,err,error,*999)
                    DO mc=1,SIZE(ellipsoidMesh%bases),2
                      IF((ellipsoidMesh%bases(mc)%PTR%NUMBER_OF_COLLAPSED_XI==0).AND. &
                        & (ellipsoidMesh%bases(mc+1)%PTR%NUMBER_OF_COLLAPSED_XI>0))THEN
                        !test for collapsed nodes and force non collapsed to wall elements and collapsed to apex elements
                        basis1=>ellipsoidMesh%bases(mc)%PTR
                        basis2=>ellipsoidMesh%bases(mc+1)%PTR
                      ELSE
                        CALL FlagError("For each basis, one non collapsed version (basis1) and one collapsed "// &
                          "version (basis2) is needed.",err,error,*999)
                      ENDIF
                      SELECT CASE(basis1%TYPE)
                        !should also test for basis2
                      CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
                        IF(basis1%NUMBER_OF_XI==SIZE(numberOfElementsXi,1).AND. &
                          & basis2%NUMBER_OF_XI==SIZE(numberOfElementsXi,1)) THEN
                          IF(.NOT.ALL(numberOfElementsXi>0)) &
                            & CALL FlagError("Must have 1 or more elements in all directions.",err,error,*999)
                          IF(numberOfElementsXi(1)<3) &
                            & CALL FlagError("Need >2 elements around the circumferential direction.", &
                            & err,error,*999)
                          !IF(.NOT.ALL(BASIS%COLLAPSED_XI==BASIS_NOT_COLLAPSED))  &
                          !    & CALL FlagError("Degenerate (collapsed) basis not implemented.",err,error,*999)
                          IF(ALLOCATED(nodeIndices)) DEALLOCATE(nodeIndices)
                          IF(ALLOCATED(elementIndices)) DEALLOCATE(elementIndices)
                          IF(ALLOCATED(cornerNodes)) DEALLOCATE(cornerNodes)
                          CALL GeneratedMeshEllipsoidBuildNodeIndices(numberOfElementsXi,basis1% &
                            NUMBER_OF_NODES_XIC, ellipsoidMesh%ellipsoidExtent, totalNumberOfNodes, &
                            totalNumberOfElements, nodeIndices,cornerNodes,elementIndices,DELTA,DELTAi,err,error,*999)
                          IF(mc==1) THEN
                            CALL MESH_NUMBER_OF_ELEMENTS_SET(generatedMesh%MESH,totalNumberOfElements, &
                              & err,error,*999)
                          ENDIF
                          
                          !Create the default node set
                          !TODO we finish create after the nodes are initialised?
                          
                          NULLIFY(meshElements)
                          CALL MESH_TOPOLOGY_ELEMENTS_CREATE_START(generatedMesh%mesh,mc/2+1,basis1,meshElements, &
                            err, error,*999)
                          !Set the elements for the ellipsoid mesh
                          IF(ALLOCATED(wallElementNodes)) DEALLOCATE(wallElementNodes)
                          IF(ALLOCATED(apexElementNodes)) DEALLOCATE(apexElementNodes)
                          IF(ALLOCATED(wallElementNodesUserNumbers)) DEALLOCATE(wallElementNodesUserNumbers)
                          IF(ALLOCATED(apexElementNodesUserNumbers)) DEALLOCATE(apexElementNodesUserNumbers)
                          ALLOCATE(wallElementNodes(basis1%NUMBER_OF_NODES),STAT=err)
                          IF(err/=0) CALL FlagError("Could not allocate wall element nodes.",err,error,*999)
                          ALLOCATE(apexElementNodes(basis2%NUMBER_OF_NODES),STAT=err)
                          IF(err/=0) CALL FlagError("Could not allocate apex element nodes.",err,error,*999)
                          ALLOCATE(wallElementNodesUserNumbers(basis1%NUMBER_OF_NODES),STAT=err)
                          IF(err/=0) CALL FlagError("Could not allocate wall element nodes.",err,error,*999)
                          ALLOCATE(apexElementNodesUserNumbers(basis2%NUMBER_OF_NODES),STAT=err)
                          IF(err/=0) CALL FlagError("Could not allocate apex element nodes.",err,error,*999)
                          ! calculate element topology (nodes per each element)
                          ! the idea is to translate given (r,theta,z) to nodeIndices equivalents, which include interior nodes
                          elementIdx=0
                          localNodeIdx=0
                          !fromJ=global J direction counting number of first node in element in J direction
                          DO elementIdx3=1,numberOfElementsXi(3)
                            from3=NINT(DELTA(3)*(elementIdx3-1)/DELTAi(3)+1)
                            elementIdx2=1
                            from2=NINT(DELTA(2)*(elementIdx2-1)/DELTAi(2)+1)
                            !apex elements
                            DO elementIdx1=1,numberOfElementsXi(1)
                              from1=NINT(DELTA(1)*(elementIdx1-1)/DELTAi(1)+1)
                              localNodeIdx=0
                              ! number of nodes in an element is dependent on basis used
                              DO localNodeIdx3=from3,from3+basis2%NUMBER_OF_NODES_XIC(3)-1
                                localNodeIdx2=1
                                localNodeIdx1=1
                                !central axis nodes
                                localNodeIdx=localNodeIdx+1
                                apexElementNodes(localNodeIdx)=nodeIndices(localNodeIdx1,localNodeIdx2,localNodeIdx3)
                                DO localNodeIdx2=from2+1,from2+basis2%NUMBER_OF_NODES_XIC(2)-1
                                  DO localNodeIdx1=from1,from1+basis2%NUMBER_OF_NODES_XIC(1)-1
                                    localNodeIdx=localNodeIdx+1
                                    ! circumferential loop-around
                                    IF(localNodeIdx1>SIZE(nodeIndices,1)) THEN
                                      apexElementNodes(localNodeIdx)=nodeIndices(1,localNodeIdx2,localNodeIdx3)
                                    ELSE
                                      apexElementNodes(localNodeIdx)=nodeIndices(localNodeIdx1,localNodeIdx2,localNodeIdx3)
                                    ENDIF
                                  ENDDO ! localNodeIdx1
                                ENDDO ! localNodeIdx2
                              ENDDO ! localNodeIdx3
                              elementIdx=elementIdx+1
                              CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_SET(elementIdx,meshElements,basis2,err,error,*999)
                              CALL GeneratedMeshComponentNodesToUserNumbers(ellipsoidMesh%generatedMesh,mc,apexElementNodes, &
                                & apexElementNodesUserNumbers,err,error,*999)
                              CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(elementIdx,meshElements,apexElementNodesUserNumbers, &
                                & err,error,*999)
                            ENDDO ! elementIdx1
                            !wall elements
                            DO elementIdx2=2,numberOfElementsXi(2)
                              from2=NINT(DELTA(2)*(elementIdx2-1)/DELTAi(2)+1)
                              DO elementIdx1=1,numberOfElementsXi(1)
                                from1=NINT(DELTA(1)*(elementIdx1-1)/DELTAi(1)+1)
                                localNodeIdx=0
                                ! number of nodes in an element is dependent on basis used
                                DO localNodeIdx3=from3,from3+basis1%NUMBER_OF_NODES_XIC(3)-1
                                  DO localNodeIdx2=from2,from2+basis1%NUMBER_OF_NODES_XIC(2)-1
                                    DO localNodeIdx1=from1,from1+basis1%NUMBER_OF_NODES_XIC(1)-1
                                      localNodeIdx=localNodeIdx+1
                                      ! circumferential loop-around
                                      IF(localNodeIdx1>SIZE(nodeIndices,1)) THEN
                                        wallElementNodes(localNodeIdx)=nodeIndices(1,localNodeIdx2,localNodeIdx3)
                                      ELSE
                                        wallElementNodes(localNodeIdx)=nodeIndices(localNodeIdx1,localNodeIdx2,localNodeIdx3)
                                      ENDIF
                                    ENDDO ! localNodeIdx1
                                  ENDDO ! localNodeIdx2
                                ENDDO ! localNodeIdx3
                                elementIdx=elementIdx+1
                                CALL GeneratedMeshComponentNodesToUserNumbers(ellipsoidMesh%generatedMesh,mc, &
                                  & wallElementNodes,wallElementNodesUserNumbers,err,error,*999)
                                CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(elementIdx,meshElements, &
                                  & wallElementNodesUserNumbers,err,error,*999)
                              ENDDO ! elementIdx1
                            ENDDO ! elementIdx2
                          ENDDO ! elementIdx3
                          CALL MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH(meshElements,err,error,*999)
                        ELSE
                          CALL FlagError("The number of xi directions of the given basis does not match the size of &
                            &the number of elements for the mesh.",err,error,*999)
                        ENDIF
                      CASE(BASIS_SIMPLEX_TYPE)
                        CALL FlagError("Ellipsoid meshes with simplex basis types is not implemented.",err,error,*999)
                      CASE DEFAULT
                        CALL FlagError("Basis type is either invalid or not implemented.",err,error,*999)
                      END SELECT
                    ENDDO
                    !Finish the mesh
                    CALL MESH_CREATE_FINISH(generatedMesh%MESH,err,error,*999)
                  ELSE
                    CALL FlagError("An ellipsoid mesh requires a collapsed basis for each basis,"// &
                      & " so there must be n*2 bases.",err,error,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Bases is not allocated.",err,error,*999)
                ENDIF
              ELSE
                CALL FlagError("For an ellipsoid mesh the following measures need to be given: &
                  & LONG_AXIS, SHORT_AXIS, WALL_THICKNESS and CUTOFF_ANGLE.",err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Ellipsoid mesh requires a 3 dimensional coordinate system.",err,error,*999)
            ENDIF
          CASE DEFAULT
            CALL FlagError("Coordinate type is either invalid or not implemented.",err,error,*999)
          END SELECT
        ELSE
          CALL FlagError("Region is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Ellipsoid mesh is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Generated Mesh is not associated.",err,error,*999)
    ENDIF

    CALL Exits("GeneratedMeshEllipsoidCreateFinish")
    RETURN
    ! TODO invalidate other associations
999 IF(ALLOCATED(nodeIndices)) DEALLOCATE(nodeIndices)
    IF(ALLOCATED(elementIndices)) DEALLOCATE(elementIndices)
    IF(ALLOCATED(cornerNodes)) DEALLOCATE(cornerNodes)
    IF(ALLOCATED(numberOfElementsXi)) DEALLOCATE(numberOfElementsXi)
    IF(ALLOCATED(wallElementNodes)) DEALLOCATE(wallElementNodes)
    IF(ALLOCATED(apexElementNodes)) DEALLOCATE(apexElementNodes)
    CALL Errors("GeneratedMeshEllipsoidCreateFinish",err,error)
    CALL Exits("GeneratedMeshEllipsoidCreateFinish")
=======
    TYPE(GENERATED_MESH_ELLIPSOID_TYPE), POINTER :: ELLIPSOID_MESH
    TYPE(BASIS_TYPE), POINTER :: BASIS1,BASIS2
    INTEGER(INTG), ALLOCATABLE :: NUMBER_ELEMENTS_XI(:)!,NUMBER_OF_NODES_XIC(:)
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(NODES_TYPE), POINTER :: NODES
    INTEGER(INTG) :: TOTAL_NUMBER_OF_NODES,TOTAL_NUMBER_OF_ELEMENTS,NUMBER_OF_DIMENSIONS
    INTEGER(INTG) :: BASIS_NUMBER_OF_NODES,CORNER_NUMBER_OF_NODES
    INTEGER(INTG) :: ne1,ne2,ne3,nn1,nn2,nn3,from1,from2,from3,nn,ne,mc
    INTEGER(INTG), ALLOCATABLE :: APEX_ELEMENT_NODES(:), WALL_ELEMENT_NODES(:)
    INTEGER(INTG), ALLOCATABLE :: APEX_ELEMENT_NODES_USER_NUMBERS(:), WALL_ELEMENT_NODES_USER_NUMBERS(:)
    INTEGER(INTG), ALLOCATABLE :: NIDX(:,:,:),CORNER_NODES(:,:,:),EIDX(:,:,:)
    REAL(DP) :: DELTA(3),DELTAi(3)
    TYPE(MeshElementsType), POINTER :: MESH_ELEMENTS

    ENTERS("GENERATED_MESH_ELLIPSOID_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      ELLIPSOID_MESH=>GENERATED_MESH%ELLIPSOID_MESH
      IF(ASSOCIATED(ELLIPSOID_MESH)) THEN
        REGION=>GENERATED_MESH%REGION
        IF(ASSOCIATED(REGION)) THEN
          IF(ASSOCIATED(REGION%COORDINATE_SYSTEM)) THEN
            SELECT CASE(REGION%COORDINATE_SYSTEM%TYPE)
            CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
              !Determine the coordinate system and create the regular mesh for that system
              ELLIPSOID_MESH%MESH_DIMENSION=REGION%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS
              NUMBER_OF_DIMENSIONS=ELLIPSOID_MESH%MESH_DIMENSION
              IF(NUMBER_OF_DIMENSIONS==3) THEN ! hard-coded for 3D only
                IF(.NOT.ALLOCATED(ELLIPSOID_MESH%ORIGIN)) THEN
                  ALLOCATE(ELLIPSOID_MESH%ORIGIN(NUMBER_OF_DIMENSIONS),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate origin.",ERR,ERROR,*999)
                  ELLIPSOID_MESH%ORIGIN=0.0_DP
                ENDIF
                IF(SIZE(ELLIPSOID_MESH%ORIGIN)==ELLIPSOID_MESH%MESH_DIMENSION) THEN
                  IF(SIZE(ELLIPSOID_MESH%ELLIPSOID_EXTENT)==4) THEN
                    IF(ALLOCATED(ELLIPSOID_MESH%BASES)) THEN
                      IF(MOD(SIZE(ELLIPSOID_MESH%BASES),2)==0) THEN
                        ALLOCATE(NUMBER_ELEMENTS_XI(SIZE(ELLIPSOID_MESH%NUMBER_OF_ELEMENTS_XI)),STAT=ERR)
                        IF(ERR/=0) CALL FlagError("Could not allocate number of elements xi.",ERR,ERROR,*999)
                        NUMBER_ELEMENTS_XI(1:SIZE(ELLIPSOID_MESH%NUMBER_OF_ELEMENTS_XI))= &
                          & ELLIPSOID_MESH%NUMBER_OF_ELEMENTS_XI(1:SIZE(ELLIPSOID_MESH%NUMBER_OF_ELEMENTS_XI))
                        !Calculate total number of nodes from all bases and start mesh
                        CORNER_NUMBER_OF_NODES=NUMBER_ELEMENTS_XI(1)*(NUMBER_ELEMENTS_XI(2)+1)*(NUMBER_ELEMENTS_XI(3)+1)- &
                          & (NUMBER_ELEMENTS_XI(1)-1)*(NUMBER_ELEMENTS_XI(3)+1)
                        TOTAL_NUMBER_OF_NODES=CORNER_NUMBER_OF_NODES
                        DO mc=1,SIZE(ELLIPSOID_MESH%BASES),2
                          BASIS1=>ELLIPSOID_MESH%BASES(mc)%PTR
                          BASIS_NUMBER_OF_NODES=NUMBER_ELEMENTS_XI(1)*(BASIS1%NUMBER_OF_NODES_XIC(1)-1)* &
                            & (NUMBER_ELEMENTS_XI(2)*(BASIS1%NUMBER_OF_NODES_XIC(2)-1)+1)* &
                            & (NUMBER_ELEMENTS_XI(3)*(BASIS1%NUMBER_OF_NODES_XIC(3)-1)+1)- &
                            & (NUMBER_ELEMENTS_XI(1)*(BASIS1%NUMBER_OF_NODES_XIC(1)-1)-1)* &
                            & (NUMBER_ELEMENTS_XI(3)*(BASIS1%NUMBER_OF_NODES_XIC(3)-1)+1)
                          TOTAL_NUMBER_OF_NODES=TOTAL_NUMBER_OF_NODES+BASIS_NUMBER_OF_NODES-CORNER_NUMBER_OF_NODES
                        ENDDO
                        NULLIFY(NODES)
                        CALL NODES_CREATE_START(REGION,TOTAL_NUMBER_OF_NODES,NODES,ERR,ERROR,*999)
                        !Finish the nodes creation
                        CALL NODES_CREATE_FINISH(NODES,ERR,ERROR,*999)
                        !Create the mesh
                        CALL MESH_CREATE_START(MESH_USER_NUMBER,GENERATED_MESH%REGION, &
                          & SIZE(NUMBER_ELEMENTS_XI,1), GENERATED_MESH%MESH,ERR,ERROR,*999)
                        !Create the elements
                        CALL MESH_NUMBER_OF_COMPONENTS_SET(GENERATED_MESH%MESH,SIZE(ELLIPSOID_MESH%BASES)/2,ERR,ERROR,*999)
                        DO mc=1,SIZE(ELLIPSOID_MESH%BASES),2
                          IF((ELLIPSOID_MESH%BASES(mc)%PTR%NUMBER_OF_COLLAPSED_XI==0).AND. &
                              & (ELLIPSOID_MESH%BASES(mc+1)%PTR%NUMBER_OF_COLLAPSED_XI>0))THEN
                            !test for collapsed nodes and force non collapsed to wall elements and collapsed to apex elements
                            BASIS1=>ELLIPSOID_MESH%BASES(mc)%PTR
                            BASIS2=>ELLIPSOID_MESH%BASES(mc+1)%PTR
                          ELSE
                            CALL FlagError("For each basis, one non collapsed version (basis1) and one collapsed "// &
                                "version (basis2) is needed.",ERR,ERROR,*999)
                          ENDIF
                          SELECT CASE(BASIS1%TYPE)
                          !should also test for basis2
                          CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
                            IF(BASIS1%NUMBER_OF_XI==SIZE(NUMBER_ELEMENTS_XI,1).AND. &
                                & BASIS2%NUMBER_OF_XI==SIZE(NUMBER_ELEMENTS_XI,1)) THEN
                              IF(.NOT.ALL(NUMBER_ELEMENTS_XI>0)) &
                                  & CALL FlagError("Must have 1 or more elements in all directions.",ERR,ERROR,*999)
                              IF(NUMBER_ELEMENTS_XI(1)<3) &
                                & CALL FlagError("Need >2 elements around the circumferential direction.", &
                                & ERR,ERROR,*999)
                              !IF(.NOT.ALL(BASIS%COLLAPSED_XI==BASIS_NOT_COLLAPSED))  &
                              !    & CALL FlagError("Degenerate (collapsed) basis not implemented.",ERR,ERROR,*999)
                              IF(ALLOCATED(NIDX)) DEALLOCATE(NIDX)
                              IF(ALLOCATED(EIDX)) DEALLOCATE(EIDX)
                              IF(ALLOCATED(CORNER_NODES)) DEALLOCATE(CORNER_NODES)
                              CALL GENERATED_MESH_ELLIPSOID_BUILD_NODE_INDICES(NUMBER_ELEMENTS_XI,BASIS1% &
                                  NUMBER_OF_NODES_XIC, ELLIPSOID_MESH%ELLIPSOID_EXTENT, TOTAL_NUMBER_OF_NODES, &
                                  TOTAL_NUMBER_OF_ELEMENTS, NIDX,CORNER_NODES,EIDX,DELTA,DELTAi,ERR,ERROR,*999)
                              IF(mc==1) THEN
                                CALL MESH_NUMBER_OF_ELEMENTS_SET(GENERATED_MESH%MESH,TOTAL_NUMBER_OF_ELEMENTS, &
                                  & ERR,ERROR,*999)
                              ENDIF

                              !Create the default node set
                              !TODO we finish create after the nodes are initialised?

                              NULLIFY(MESH_ELEMENTS)
                              CALL MESH_TOPOLOGY_ELEMENTS_CREATE_START(GENERATED_MESH%MESH,mc/2+1,BASIS1,MESH_ELEMENTS, &
                                  ERR, ERROR,*999)
                              !Set the elements for the ellipsoid mesh
                              IF(ALLOCATED(WALL_ELEMENT_NODES)) DEALLOCATE(WALL_ELEMENT_NODES)
                              IF(ALLOCATED(APEX_ELEMENT_NODES)) DEALLOCATE(APEX_ELEMENT_NODES)
                              IF(ALLOCATED(WALL_ELEMENT_NODES_USER_NUMBERS)) DEALLOCATE(WALL_ELEMENT_NODES_USER_NUMBERS)
                              IF(ALLOCATED(APEX_ELEMENT_NODES_USER_NUMBERS)) DEALLOCATE(APEX_ELEMENT_NODES_USER_NUMBERS)
                              ALLOCATE(WALL_ELEMENT_NODES(BASIS1%NUMBER_OF_NODES),STAT=ERR)
                              IF(ERR/=0) CALL FlagError("Could not allocate wall element nodes.",ERR,ERROR,*999)
                              ALLOCATE(APEX_ELEMENT_NODES(BASIS2%NUMBER_OF_NODES),STAT=ERR)
                              IF(ERR/=0) CALL FlagError("Could not allocate apex element nodes.",ERR,ERROR,*999)
                              ALLOCATE(WALL_ELEMENT_NODES_USER_NUMBERS(BASIS1%NUMBER_OF_NODES),STAT=ERR)
                              IF(ERR/=0) CALL FlagError("Could not allocate wall element nodes.",ERR,ERROR,*999)
                              ALLOCATE(APEX_ELEMENT_NODES_USER_NUMBERS(BASIS2%NUMBER_OF_NODES),STAT=ERR)
                              IF(ERR/=0) CALL FlagError("Could not allocate apex element nodes.",ERR,ERROR,*999)
                              ! calculate element topology (nodes per each element)
                              ! the idea is to translate given (r,theta,z) to NIDX equivalents, which include interior nodes
                              ne=0
                              nn=0
                              !fromJ=global J direction counting number of first node in element in J direction
                              DO ne3=1,NUMBER_ELEMENTS_XI(3)
                                from3=NINT(DELTA(3)*(ne3-1)/DELTAi(3)+1)
                                ne2=1
                                from2=NINT(DELTA(2)*(ne2-1)/DELTAi(2)+1)
                                !apex elements
                                DO ne1=1,NUMBER_ELEMENTS_XI(1)
                                  from1=NINT(DELTA(1)*(ne1-1)/DELTAi(1)+1)
                                  nn=0
                                  ! number of nodes in an element is dependent on basis used
                                  DO nn3=from3,from3+BASIS2%NUMBER_OF_NODES_XIC(3)-1
                                    nn2=1
                                    nn1=1
                                    !central axis nodes
                                    nn=nn+1
                                    APEX_ELEMENT_NODES(nn)=NIDX(nn1,nn2,nn3)
                                    DO nn2=from2+1,from2+BASIS2%NUMBER_OF_NODES_XIC(2)-1
                                      DO nn1=from1,from1+BASIS2%NUMBER_OF_NODES_XIC(1)-1
                                        nn=nn+1
                                        ! circumferential loop-around
                                        IF(nn1>SIZE(NIDX,1)) THEN
                                          APEX_ELEMENT_NODES(nn)=NIDX(1,nn2,nn3)
                                        ELSE
                                          APEX_ELEMENT_NODES(nn)=NIDX(nn1,nn2,nn3)
                                        ENDIF
                                      ENDDO ! nn1
                                    ENDDO ! nn2
                                  ENDDO ! nn3
                                  ne=ne+1
                                  CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_SET(ne,MESH_ELEMENTS,BASIS2,ERR,ERROR,*999)
                                  CALL COMPONENT_NODES_TO_USER_NUMBERS(ELLIPSOID_MESH%GENERATED_MESH,mc,APEX_ELEMENT_NODES, &
                                      & APEX_ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                                  CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(ne,MESH_ELEMENTS, &
                                    APEX_ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                                ENDDO ! ne1
                                !wall elements
                                DO ne2=2,NUMBER_ELEMENTS_XI(2)
                                  from2=NINT(DELTA(2)*(ne2-1)/DELTAi(2)+1)
                                  DO ne1=1,NUMBER_ELEMENTS_XI(1)
                                    from1=NINT(DELTA(1)*(ne1-1)/DELTAi(1)+1)
                                    nn=0
                                    ! number of nodes in an element is dependent on basis used
                                    DO nn3=from3,from3+BASIS1%NUMBER_OF_NODES_XIC(3)-1
                                      DO nn2=from2,from2+BASIS1%NUMBER_OF_NODES_XIC(2)-1
                                        DO nn1=from1,from1+BASIS1%NUMBER_OF_NODES_XIC(1)-1
                                          nn=nn+1
                                          ! circumferential loop-around
                                          IF(nn1>SIZE(NIDX,1)) THEN
                                            WALL_ELEMENT_NODES(nn)=NIDX(1,nn2,nn3)
                                          ELSE
                                            WALL_ELEMENT_NODES(nn)=NIDX(nn1,nn2,nn3)
                                          ENDIF
                                        ENDDO ! nn1
                                      ENDDO ! nn2
                                    ENDDO ! nn3
                                    ne=ne+1
                                    CALL COMPONENT_NODES_TO_USER_NUMBERS(ELLIPSOID_MESH%GENERATED_MESH,mc,WALL_ELEMENT_NODES, &
                                        & WALL_ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                                    CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(ne,MESH_ELEMENTS, &
                                        & WALL_ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                                  ENDDO ! ne1
                                ENDDO ! ne2
                              ENDDO ! ne3
                              CALL MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH(MESH_ELEMENTS,ERR,ERROR,*999)
                            ELSE
                              CALL FlagError("The number of xi directions of the given basis does not match the size of &
                                &the number of elements for the mesh.",ERR,ERROR,*999)
                            ENDIF
                          CASE(BASIS_SIMPLEX_TYPE)
                            CALL FlagError("Ellipsoid meshes with simplex basis types is not implemented.",ERR,ERROR,*999)
                          CASE DEFAULT
                            CALL FlagError("Basis type is either invalid or not implemented.",ERR,ERROR,*999)
                          END SELECT
                        ENDDO
                        !Finish the mesh
                        CALL MESH_CREATE_FINISH(GENERATED_MESH%MESH,ERR,ERROR,*999)
                      ELSE
                        CALL FlagError("An ellipsoid mesh requires a collapsed basis for each basis,"// &
                            & " so there must be n*2 bases.",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FlagError("Bases is not allocated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("For an ellipsoid mesh the following measures need to be given: &
                        & LONG_AXIS, SHORT_AXIS, WALL_THICKNESS and CUTOFF_ANGLE.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FlagError("The number of dimensions of the given regular mesh does not match the size of &
                      &the origin.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Ellipsoid mesh requires a 3 dimensional coordinate system.",ERR,ERROR,*999)
              ENDIF
            CASE DEFAULT
              CALL FlagError("Coordinate type is either invalid or not implemented.",ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FlagError("Coordiate System is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Region is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Ellipsoid mesh is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Generated Mesh is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("GENERATED_MESH_ELLIPSOID_CREATE_FINISH")
    RETURN
    ! TODO invalidate other associations
999 IF(ALLOCATED(NIDX)) DEALLOCATE(NIDX)
    IF(ALLOCATED(EIDX)) DEALLOCATE(EIDX)
    IF(ALLOCATED(CORNER_NODES)) DEALLOCATE(CORNER_NODES)
    IF(ALLOCATED(NUMBER_ELEMENTS_XI)) DEALLOCATE(NUMBER_ELEMENTS_XI)
    IF(ALLOCATED(WALL_ELEMENT_NODES)) DEALLOCATE(WALL_ELEMENT_NODES)
    IF(ALLOCATED(APEX_ELEMENT_NODES)) DEALLOCATE(APEX_ELEMENT_NODES)
    ERRORSEXITS("GENERATED_MESH_ELLIPSOID_CREATE_FINISH",ERR,ERROR)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
    RETURN 1
  END SUBROUTINE GeneratedMeshEllipsoidCreateFinish
  !
  !================================================================================================================================
  !

  !>Start to create the cylinder generated mesh type
  SUBROUTINE GeneratedMeshCylinderCreateFinish(generatedMesh,meshUserNumber,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh
    INTEGER(INTG), INTENT(IN) :: meshUserNumber !<The user number for the mesh to generate.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
<<<<<<< HEAD
    TYPE(GeneratedMeshCylinderType), POINTER :: cylinderMesh
    TYPE(BASIS_TYPE), POINTER :: basis
    INTEGER(INTG), ALLOCATABLE :: numberOfElementsXi(:)
    TYPE(REGION_TYPE), POINTER :: region
    TYPE(NODES_TYPE), POINTER :: nodes
    INTEGER(INTG) :: totalNumberOfNodes,totalNumberOfElements,numberOfDimensions
    INTEGER(INTG) :: cornerNumberOfNodes,basisNumberOfNodes
    INTEGER(INTG) :: elementIdx,elementIdx1,elementIdx2,elementIdx3,localNodeIdx1,localNodeIdx2,localNodeIdx3,from1,from2,from3, &
      & localNodeNumber,basisIdx
    INTEGER(INTG), ALLOCATABLE :: elementNodes(:),elementNodesUserNumbers(:)
    INTEGER(INTG), ALLOCATABLE :: nodeIndices(:,:,:),elementIndices(:,:,:)
    REAL(DP) :: deltaCoordinates(3),deltaCoordinatesXi(3)
    TYPE(MeshComponentElementsType), POINTER :: meshElements

    CALL Enters("GeneratedMeshCylinderCreateFinish",err,error,*999)

    IF(ASSOCIATED(generatedMesh)) THEN
      cylinderMesh=>generatedMesh%cylinderMesh
      IF(ASSOCIATED(cylinderMesh)) THEN
        region=>generatedMesh%region
        IF(ASSOCIATED(region)) THEN
          IF(ASSOCIATED(region%COORDINATE_SYSTEM)) THEN
=======
    TYPE(GENERATED_MESH_CYLINDER_TYPE), POINTER :: CYLINDER_MESH
    TYPE(BASIS_TYPE), POINTER :: BASIS
    INTEGER(INTG), ALLOCATABLE :: NUMBER_ELEMENTS_XI(:)
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(NODES_TYPE), POINTER :: NODES
    INTEGER(INTG) :: TOTAL_NUMBER_OF_NODES,TOTAL_NUMBER_OF_ELEMENTS,NUMBER_OF_DIMENSIONS
    INTEGER(INTG) :: CORNER_NUMBER_OF_NODES,BASIS_NUMBER_OF_NODES
    INTEGER(INTG) :: ne1,ne2,ne3,nn1,nn2,nn3,from1,from2,from3,nn,ne,basis_idx
    INTEGER(INTG), ALLOCATABLE :: ELEMENT_NODES(:),ELEMENT_NODES_USER_NUMBERS(:)
    INTEGER(INTG), ALLOCATABLE :: NIDX(:,:,:),EIDX(:,:,:)
    REAL(DP) :: DELTA(3),DELTAi(3)
    TYPE(MeshElementsType), POINTER :: MESH_ELEMENTS

    ENTERS("GENERATED_MESH_CYLINDER_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      CYLINDER_MESH=>GENERATED_MESH%CYLINDER_MESH
      IF(ASSOCIATED(CYLINDER_MESH)) THEN
        REGION=>GENERATED_MESH%REGION
        IF(ASSOCIATED(REGION)) THEN
          IF(ASSOCIATED(REGION%COORDINATE_SYSTEM)) THEN
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
            !TODO is regular type only for COORDINATE_RECTANGULAR_CARTESIAN_TYPE?
            !If that, should we use IF rather than select?
            SELECT CASE(region%COORDINATE_SYSTEM%TYPE)
            CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
<<<<<<< HEAD
              numberOfDimensions=generatedMesh%meshDimension
              IF(numberOfDimensions==3) THEN ! hard-coded for 3D only
                IF(.NOT.ALLOCATED(cylinderMesh%origin)) THEN
                  ALLOCATE(cylinderMesh%origin(numberOfDimensions),STAT=err)
                  IF(err/=0) CALL FlagError("Could not allocate origin.",err,error,*999)
                  cylinderMesh%origin=0.0_DP
                ENDIF
                IF(SIZE(generatedMesh%origin,1)==generatedMesh%coordinateDimension) THEN
                  IF(SIZE(cylinderMesh%cylinderExtent)==generatedMesh%meshDimension) THEN
                    IF(ALLOCATED(generatedMesh%bases)) THEN
                      ALLOCATE(numberOfElementsXi(SIZE(cylinderMesh%numberOfElementsXi)),STAT=err)
                      IF(err/=0) CALL FlagError("Could not allocate number of elements xi.",err,error,*999)
                      numberOfElementsXi(1:SIZE(cylinderMesh%numberOfElementsXi))= &
                        & cylinderMesh%numberOfElementsXi(1:SIZE(cylinderMesh%numberOfElementsXi))
                      CALL MESH_CREATE_START(meshUserNumber,generatedMesh%region,SIZE(numberOfElementsXi,1), &
                        & generatedMesh%mesh,err,error,*999)
                      CALL MESH_NUMBER_OF_COMPONENTS_SET(generatedMesh%mesh,SIZE(cylinderMesh%bases),err,error,*999)
=======
              !Determine the coordinate system and create the regular mesh for that system
              CYLINDER_MESH%MESH_DIMENSION=REGION%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS
              NUMBER_OF_DIMENSIONS=CYLINDER_MESH%MESH_DIMENSION
              IF(NUMBER_OF_DIMENSIONS==3) THEN ! hard-coded for 3D only
                IF(.NOT.ALLOCATED(CYLINDER_MESH%ORIGIN)) THEN
                  ALLOCATE(CYLINDER_MESH%ORIGIN(NUMBER_OF_DIMENSIONS),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate origin.",ERR,ERROR,*999)
                  CYLINDER_MESH%ORIGIN=0.0_DP
                ENDIF
                IF(SIZE(CYLINDER_MESH%ORIGIN)==CYLINDER_MESH%MESH_DIMENSION) THEN
                  IF(SIZE(CYLINDER_MESH%CYLINDER_EXTENT)==CYLINDER_MESH%MESH_DIMENSION) THEN
                    IF(ALLOCATED(CYLINDER_MESH%BASES)) THEN
                      ALLOCATE(NUMBER_ELEMENTS_XI(SIZE(CYLINDER_MESH%NUMBER_OF_ELEMENTS_XI)),STAT=ERR)
                      IF(ERR/=0) CALL FlagError("Could not allocate number of elements xi.",ERR,ERROR,*999)
                      NUMBER_ELEMENTS_XI(1:SIZE(CYLINDER_MESH%NUMBER_OF_ELEMENTS_XI))= &
                        & CYLINDER_MESH%NUMBER_OF_ELEMENTS_XI(1:SIZE(CYLINDER_MESH%NUMBER_OF_ELEMENTS_XI))
                      CALL MESH_CREATE_START(MESH_USER_NUMBER,GENERATED_MESH%REGION,SIZE(NUMBER_ELEMENTS_XI,1), &
                        & GENERATED_MESH%MESH,ERR,ERROR,*999)
                      CALL MESH_NUMBER_OF_COMPONENTS_SET(GENERATED_MESH%MESH,SIZE(CYLINDER_MESH%BASES),ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
                      !Calculate number of nodes
                      cornerNumberOfNodes=(numberOfElementsXi(3)+1)*numberOfElementsXi(2)*(numberOfElementsXi(1)+1)
                      totalNumberOfNodes=cornerNumberOfNodes
                      DO basisIdx=1,SIZE(generatedMesh%bases)
                        basis=>generatedMesh%bases(basisIdx)%ptr
                        IF(ASSOCIATED(basis)) THEN
                          basisNumberOfNodes=((basis%NUMBER_OF_NODES_XIC(3)-1)*numberOfElementsXi(3)+1)* &
                              & ((basis%NUMBER_OF_NODES_XIC(2)-1)*numberOfElementsXi(2))* &
                              & ((basis%NUMBER_OF_NODES_XIC(1)-1)*numberOfElementsXi(1)+1)
                          totalNumberOfNodes=totalNumberOfNodes+basisNumberOfNodes-cornerNumberOfNodes
                        ELSE
<<<<<<< HEAD
                          CALL FlagError("Basis is not associated.",err,error,*999)
=======
                          CALL FlagError("Basis is not associated.",ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
                        ENDIF
                      ENDDO
                      NULLIFY(nodes)
                      CALL NODES_CREATE_START(region,totalNumberOfNodes,nodes,err,error,*999)
                      !Finish the nodes creation
                      CALL NODES_CREATE_FINISH(nodes,err,error,*999)
                      !Set the total number of elements
                      totalNumberOfElements=numberOfElementsXi(1)*numberOfElementsXi(2)*numberOfElementsXi(3)
                      CALL MESH_NUMBER_OF_ELEMENTS_SET(generatedMesh%MESH,totalNumberOfElements,err,error,*999)
                      DO basisIdx=1,SIZE(cylinderMesh%bases)
                        basis=>cylinderMesh%bases(basisIdx)%ptr
                        IF(ASSOCIATED(basis)) THEN
                          SELECT CASE(basis%TYPE)
                          CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
<<<<<<< HEAD
                            IF(basis%NUMBER_OF_XI==SIZE(numberOfElementsXi,1)) THEN
                              IF(.NOT.ALL(numberOfElementsXi>0)) &
                                & CALL FlagError("Must have 1 or more elements in all directions.",err,error,*999)
                              IF(numberOfElementsXi(2)<3) &
                                CALL FlagError("Need >2 elements around the circumferential direction.",err,error,*999)
                              IF(.NOT.ALL(basis%COLLAPSED_XI==BASIS_NOT_COLLAPSED))  &
                                & CALL FlagError("Degenerate (collapsed) basis not implemented.",err,error,*999)
=======
                            IF(BASIS%NUMBER_OF_XI==SIZE(NUMBER_ELEMENTS_XI,1)) THEN
                              IF(.NOT.ALL(NUMBER_ELEMENTS_XI>0)) &
                                & CALL FlagError("Must have 1 or more elements in all directions.",ERR,ERROR,*999)
                              IF(NUMBER_ELEMENTS_XI(2)<3) &
                                CALL FlagError("Need >2 elements around the circumferential direction.",ERR,ERROR,*999)
                              IF(.NOT.ALL(BASIS%COLLAPSED_XI==BASIS_NOT_COLLAPSED))  &
                                & CALL FlagError("Degenerate (collapsed) basis not implemented.",ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
                              !Calculate nodes and element sizes
                              IF(ALLOCATED(nodeIndices)) DEALLOCATE(nodeIndices)
                              IF(ALLOCATED(elementIndices)) DEALLOCATE(elementIndices)
                              CALL GeneratedMeshCylinderBuildNodeIndices(numberOfElementsXi,basis%NUMBER_OF_NODES_XIC, &
                                & cylinderMesh%cylinderExtent, totalNumberOfNodes,totalNumberOfElements, &
                                & nodeIndices,elementIndices,deltaCoordinates,deltaCoordinatesXi,err,error,*999)
                              !Set the elements for the cylinder mesh
<<<<<<< HEAD
                              IF(ALLOCATED(elementNodes)) DEALLOCATE(elementNodes)
                              IF(ALLOCATED(elementNodesUserNumbers)) DEALLOCATE(elementNodesUserNumbers)
                              ALLOCATE(elementNodesUserNumbers(basis%NUMBER_OF_NODES),STAT=err)
                              ALLOCATE(elementNodes(basis%NUMBER_OF_NODES),STAT=err)
                              IF(err/=0) CALL FlagError("Could not allocate element nodes.",err,error,*999)
=======
                              IF(ALLOCATED(ELEMENT_NODES)) DEALLOCATE(ELEMENT_NODES)
                              IF(ALLOCATED(ELEMENT_NODES_USER_NUMBERS)) DEALLOCATE(ELEMENT_NODES_USER_NUMBERS)
                              ALLOCATE(ELEMENT_NODES_USER_NUMBERS(BASIS%NUMBER_OF_NODES),STAT=ERR)
                              ALLOCATE(ELEMENT_NODES(BASIS%NUMBER_OF_NODES),STAT=ERR)
                              IF(ERR/=0) CALL FlagError("Could not allocate element nodes.",ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
                              !Create the elements
                              NULLIFY(meshElements)
                              CALL MESH_TOPOLOGY_ELEMENTS_CREATE_START(generatedMesh%mesh,basisIdx,basis,meshElements, &
                                & err,error,*999)
                              ! calculate element topology (nodes per each element)
                              ! the idea is to translate given (r,theta,z) to nodeIndices equivalents, which include interior nodes
                              elementIdx=0
                              DO elementIdx3=1,numberOfElementsXi(3)
                                from3=NINT(deltaCoordinates(3)*(elementIdx3-1)/deltaCoordinatesXi(3)+1)
                                DO elementIdx2=1,numberOfElementsXi(2)
                                  from2=NINT(deltaCoordinates(2)*(elementIdx2-1)/deltaCoordinatesXi(2)+1)
                                  DO elementIdx1=1,numberOfElementsXi(1)
                                    from1=NINT(deltaCoordinates(1)*(elementIdx1-1)/deltaCoordinatesXi(1)+1)
                                    localNodeNumber=0
                                    ! number of nodes in an element is dependent on basis used
                                    DO localNodeIdx3=from3,from3+basis%NUMBER_OF_NODES_XIC(3)-1
                                      DO localNodeIdx2=from2,from2+basis%NUMBER_OF_NODES_XIC(2)-1
                                        DO localNodeIdx1=from1,from1+basis%NUMBER_OF_NODES_XIC(1)-1
                                          localNodeNumber=localNodeNumber+1
                                          ! compensate for circumferential loop-around
                                          IF(localNodeIdx2>SIZE(nodeIndices,2)) THEN
                                            ! DEBUG: little check here
<<<<<<< HEAD
                                            IF(localNodeIdx2>SIZE(nodeIndices,2)+1) &
                                              & CALL FlagError("nodeIndices needs debugging",err,error,*999)
                                            elementNodes(localNodeNumber)=nodeIndices(localNodeIdx1,1,localNodeIdx3)
=======
                                            IF(nn2>SIZE(NIDX,2)+1) CALL FlagError("NIDX needs debugging",ERR,ERROR,*999)
                                            ELEMENT_NODES(nn)=NIDX(nn1,1,nn3)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
                                          ELSE
                                            elementNodes(localNodeNumber)=nodeIndices(localNodeIdx1,localNodeIdx2,localNodeIdx3)
                                          ENDIF
                                        ENDDO ! localNodeIdx1
                                      ENDDO ! localNodeIdx2
                                    ENDDO ! localNodeIdx3
                                    elementIdx=elementIdx+1
                                    CALL GeneratedMeshComponentNodesToUserNumbers(cylinderMesh%generatedMesh,basisIdx, &
                                      & elementNodes,elementNodesUserNumbers,err,error,*999)
                                    CALL MESH_TOPOLOGY_ELEMENTS_elementNodes_SET(elementIdx,meshElements,elementNodesUserNumbers, &
                                      & err,error,*999)
                                  ENDDO ! elementIdx1
                                ENDDO ! elementIdx2
                              ENDDO ! elementIdx3
                              CALL MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH(meshElements,err,error,*999)
                            ELSE
                              CALL FlagError("The number of xi directions of the given basis does not match the size of &
<<<<<<< HEAD
                                &the number of elements for the mesh.",err,error,*999)
                            ENDIF
                          CASE(BASIS_SIMPLEX_TYPE)
                            CALL FlagError("Cylinder meshes with simplex basis types is not implemented.",err,error,*999)
                          CASE DEFAULT
                            CALL FlagError("Basis type is either invalid or not implemented.",err,error,*999)
                          END SELECT
                        ELSE
                          CALL FlagError("Basis is not associated.",err,error,*999)
=======
                                &the number of elements for the mesh.",ERR,ERROR,*999)
                            ENDIF
                          CASE(BASIS_SIMPLEX_TYPE)
                            CALL FlagError("Cylinder meshes with simplex basis types is not implemented.",ERR,ERROR,*999)
                          CASE DEFAULT
                            CALL FlagError("Basis type is either invalid or not implemented.",ERR,ERROR,*999)
                          END SELECT
                        ELSE
                          CALL FlagError("Basis is not associated.",ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
                        ENDIF
                      ENDDO
                      !Finish the mesh
                      CALL MESH_CREATE_FINISH(generatedMesh%MESH,err,error,*999)
                    ELSE
<<<<<<< HEAD
                      CALL FlagError("Bases are not allocated.",err,error,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("The number of dimensions of the given regular mesh does not match the size of &
                      &the maximum extent.",err,error,*999)
                  ENDIF
                ELSE
                  CALL FlagError("The number of dimensions of the given regular mesh does not match the size of &
                    &the origin.",err,error,*999)
                ENDIF
              ELSE
                CALL FlagError("Cylinder mesh requires a 3 dimensional coordinate system.",err,error,*999)
              ENDIF
            CASE DEFAULT
              CALL FlagError("Coordinate type is either invalid or not implemented.",err,error,*999)
            END SELECT
          ELSE
            CALL FlagError("Coordiate System is not associated.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Region is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Regular mesh is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Generated Mesh is not associated.",err,error,*999)
    ENDIF

    CALL Exits("GeneratedMeshCylinderCreateFinish")
    RETURN
    ! TODO invalidate other associations
999 IF(ALLOCATED(nodeIndices)) DEALLOCATE(nodeIndices)
    IF(ALLOCATED(elementIndices)) DEALLOCATE(elementIndices)
    IF(ALLOCATED(numberOfElementsXi)) DEALLOCATE(numberOfElementsXi)
    IF(ALLOCATED(elementNodes)) DEALLOCATE(elementNodes)
    CALL Errors("GeneratedMeshCylinderCreateFinish",err,error)
    CALL Exits("GeneratedMeshCylinderCreateFinish")
=======
                      CALL FlagError("Bases are not allocated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("The number of dimensions of the given regular mesh does not match the size of &
                      &the maximum extent.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FlagError("The number of dimensions of the given regular mesh does not match the size of &
                    &the origin.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Cylinder mesh requires a 3 dimensional coordinate system.",ERR,ERROR,*999)
              ENDIF
            CASE DEFAULT
              CALL FlagError("Coordinate type is either invalid or not implemented.",ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FlagError("Coordiate System is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Region is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Regular mesh is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Generated Mesh is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("GENERATED_MESH_CYLINDER_CREATE_FINISH")
    RETURN
    ! TODO invalidate other associations
999 IF(ALLOCATED(NIDX)) DEALLOCATE(NIDX)
    IF(ALLOCATED(EIDX)) DEALLOCATE(EIDX)
    IF(ALLOCATED(NUMBER_ELEMENTS_XI)) DEALLOCATE(NUMBER_ELEMENTS_XI)
    IF(ALLOCATED(ELEMENT_NODES)) DEALLOCATE(ELEMENT_NODES)
    ERRORSEXITS("GENERATED_MESH_CYLINDER_CREATE_FINISH",ERR,ERROR)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
    RETURN 1
  END SUBROUTINE GeneratedMeshCylinderCreateFinish

  !
  !================================================================================================================================
  !

  !>Finalise the cylinder mesh type
  SUBROUTINE GeneratedMeshCylinderFinalise(cylinderMesh,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshCylinderType), POINTER :: cylinderMesh !<A pointer to the cylinder generated mesh to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

<<<<<<< HEAD
    CALL Enters("GeneratedMeshCylinderFinalise",err,error,*999)
=======
    ENTERS("GENERATED_MESH_CYLINDER_FINALISE",ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b

    IF(ASSOCIATED(cylinderMesh)) THEN
      IF(ALLOCATED(cylinderMesh%cylinderExtent)) DEALLOCATE(cylinderMesh%cylinderExtent)
      IF(ALLOCATED(cylinderMesh%numberOfElementsXi)) DEALLOCATE(cylinderMesh%numberOfElementsXi)
      DEALLOCATE(cylinderMesh)
    ENDIF

<<<<<<< HEAD
    CALL Exits("GeneratedMeshCylinderFinalise")
    RETURN
    ! TODO invalidate other associations
999 CALL Errors("GeneratedMeshCylinderFinalise",err,error)
    CALL Exits("GeneratedMeshCylinderFinalise")
=======
    EXITS("GENERATED_MESH_CYLINDER_FINALISE")
    RETURN
    ! TODO invalidate other associations
999 ERRORSEXITS("GENERATED_MESH_CYLINDER_FINALISE",ERR,ERROR)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
    RETURN 1
    
  END SUBROUTINE GeneratedMeshCylinderFinalise

  !
  !================================================================================================================================
  !

  !>Initialise the cylinder generated mesh type
  SUBROUTINE GeneratedMeshCylinderInitialise(generatedMesh,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

<<<<<<< HEAD
    CALL Enters("GeneratedMeshCylinderInitialise",err,error,*998)

    IF(ASSOCIATED(generatedMesh)) THEN
      IF(ASSOCIATED(generatedMesh%cylinderMesh)) THEN
        CALL FlagError("Cylinder mesh is already associated for this generated mesh.",err,error,*998)
      ELSE
        IF(generatedMesh%coordinateDimension==3) THEN
          ALLOCATE(generatedMesh%cylinderMesh,STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate cylinder generated mesh.",err,error,*999)
          generatedMesh%cylinderMesh%generatedMesh=>generatedMesh
          generatedMesh%generatedType=GENERATED_MESH_CYLINDER_MESH_TYPE
          !Set defaults
          generatedMesh%cylinderMesh%closed=.TRUE.
          ALLOCATE(generatedMesh%cylinderMesh%cylinderExtent(1:generatedMesh%coordinateDimension),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate cylinder mesh cylinder extent.",err,error,*999)
          generatedMesh%cylinderMesh%cylinderExtent(1)=2.0_DP !Length
          generatedMesh%cylinderMesh%cylinderExtent(2)=1.0_DP !Inner radius
          generatedMesh%cylinderMesh%cylinderExtent(3)=1.2_DP !Outer radius          
        ELSE
          CALL FlagError("Invalid coordinate system for a cylinder generated mesh. Three coordinate dimensions "// &
            & "are required for a cylinder mesh.",err,error,*998)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Generated mesh is not associated.",err,error,*998)
    ENDIF

    CALL Exits("GeneratedMeshCylinderInitialise")
    RETURN
999 CALL GeneratedMeshCylinderFinalise(generatedMesh%cylinderMesh,dummyErr,dummyError,*998)
998 CALL Errors("GeneratedMeshCylinderInitialise",err,error)
    CALL Exits("GeneratedMeshCylinderInitialise")
=======
    ENTERS("GENERATED_MESH_CYLINDER_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      IF(ASSOCIATED(GENERATED_MESH%CYLINDER_MESH)) THEN
        CALL FlagError("Cylinder mesh is already associated for this generated mesh.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(GENERATED_MESH%CYLINDER_MESH,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate cylinder generated mesh.",ERR,ERROR,*999)
        GENERATED_MESH%CYLINDER_MESH%GENERATED_MESH=>GENERATED_MESH
        GENERATED_MESH%GENERATED_TYPE=GENERATED_MESH_CYLINDER_MESH_TYPE
      ENDIF
    ELSE
      CALL FlagError("Generated mesh is not associated.",ERR,ERROR,*998)
    ENDIF

    EXITS("GENERATED_MESH_CYLINDER_INITIALISE")
    RETURN
999 CALL GENERATED_MESH_CYLINDER_FINALISE(GENERATED_MESH%CYLINDER_MESH,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("GENERATED_MESH_CYLINDER_INITIALISE",ERR,ERROR)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
    RETURN 1
  END SUBROUTINE GeneratedMeshCylinderInitialise

  !
  !================================================================================================================================
  !

  !>Finalise the polar mesh type
  SUBROUTINE GeneratedMeshPolarFinalise(polarMesh,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshPolarType), POINTER :: polarMesh !<A pointer to the polar generated mesh to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

<<<<<<< HEAD
    CALL Enters("GeneratedMeshPolarFinalise",err,error,*999)
=======
    ENTERS("GENERATED_MESH_REGULAR_FINALISE",ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b

    IF(ASSOCIATED(polarMesh)) THEN
      IF(ALLOCATED(polarMesh%polarExtent)) DEALLOCATE(polarMesh%polarExtent)
      IF(ALLOCATED(polarMesh%numberOfElementsXi)) DEALLOCATE(polarMesh%numberOfElementsXi)
      DEALLOCATE(polarMesh)
    ENDIF

<<<<<<< HEAD
    CALL Exits("GeneratedMeshPolarFinalise")
    RETURN
    ! TODO invalidate other associations
999 CALL Errors("GeneratedMeshPolarFinalise",err,error)
    CALL Exits("GeneratedMeshPolarFinalise")
=======
    EXITS("GENERATED_MESH_REGULAR_FINALISE")
    RETURN
    ! TODO invalidate other associations
999 ERRORSEXITS("GENERATED_MESH_REGULAR_FINALISE",ERR,ERROR)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
    RETURN 1
    
  END SUBROUTINE GeneratedMeshPolarFinalise

  !
  !================================================================================================================================
  !

  !>Initialise the polar generated mesh type
  SUBROUTINE GeneratedMeshPolarInitialise(generatedMesh,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh to initialise the polar mesh for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    CALL Enters("GeneratedMeshPolarInitialise",err,error,*998)

    IF(ASSOCIATED(generatedMesh)) THEN
      IF(ASSOCIATED(generatedMesh%polarMesh)) THEN
        CALL FlagError("Polar mesh is already associated for this generated mesh.",err,error,*998)
      ELSE
        IF(generatedMesh%coordinateDimension>=2) THEN
          ALLOCATE(generatedMesh%polarMesh,STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate polar generated mesh.",err,error,*999)
          generatedMesh%polarMesh%generatedMesh=>generatedMesh
          generatedMesh%generatedType=GENERATED_MESH_POLAR_MESH_TYPE
          !Set defaults
          generatedMesh%polarMesh%spherical=(generatedMesh%coordinateDimension==3)
          ALLOCATE(generatedMesh%polarMesh%polarExtent(2),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate polar mesh polar extent.",err,error,*999)
          generatedMesh%polarMesh%polarExtent(1)=1.0_DP !Inner radius
          generatedMesh%polarMesh%polarExtent(2)=1.2_DP !Outer radius
        ELSE
          CALL FlagError("Invalid coordinate system for a polar generated mesh. Two or three coordinate dimensions "// &
            & "are required for a polar mesh.",err,error,*998)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Generated mesh is not associated.",err,error,*998)
    ENDIF

    CALL Exits("GeneratedMeshPolarInitialise")
    RETURN
999 CALL GeneratedMeshPolarFinalise(generatedMesh%polarMesh,dummyErr,dummyError,*998)
998 CALL Errors("GeneratedMeshPolarInitialise",err,error)
    CALL Exits("GeneratedMeshPolarInitialise")
    RETURN 1
  END SUBROUTINE GeneratedMeshPolarInitialise

  !
  !================================================================================================================================
  !

  !>Finalise the regular mesh type
  SUBROUTINE GeneratedMeshRegularFinalise(regularMesh,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshRegularType), POINTER :: regularMesh !<A pointer to the generated mesh
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

<<<<<<< HEAD
    CALL Enters("GeneratedMeshRegularFinalise",err,error,*999)

    IF(ASSOCIATED(regularMesh)) THEN
      IF(ALLOCATED(regularMesh%maximumExtent)) DEALLOCATE(regularMesh%maximumExtent)
      IF(ALLOCATED(regularMesh%numberOfElementsXi)) DEALLOCATE(regularMesh%numberOfElementsXi)
      DEALLOCATE(regularMesh)
    ENDIF

    CALL Exits("GeneratedMeshRegularFinalise")
    RETURN
    ! TODO invalidate other associations
999 CALL Errors("GeneratedMeshRegularFinalise",err,error)
    CALL Exits("GeneratedMeshRegularFinalise")
    RETURN 1
  END SUBROUTINE GeneratedMeshRegularFinalise

  !
  !================================================================================================================================
  !

  !>Initialise the regular generated mesh type
  SUBROUTINE GeneratedMeshRegularInitialise(generatedMesh,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
    
    CALL Enters("GeneratedMeshRegularInitialise",err,error,*998)

    IF(ASSOCIATED(generatedMesh)) THEN
      IF(ASSOCIATED(generatedMesh%regularMesh)) THEN
        CALL FlagError("Regular mesh is already associated for this generated mesh.",err,error,*998)
      ELSE
        ALLOCATE(generatedMesh%regularMesh,STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate regular generated mesh.",err,error,*999)
        generatedMesh%regularMesh%generatedMesh=>generatedMesh        
        generatedMesh%generatedType=GENERATED_MESH_REGULAR_MESH_TYPE
        !Set defaults
        ALLOCATE(generatedMesh%regularMesh%maximumExtent(gneratedMesh%coordinateDimension),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate maximum extent.",err,error,*999)
        generatedMesh%regularMesh%maximumExtent=1.0_DP
      ENDIF
    ELSE
      CALL FlagError("Generated mesh is not associated.",err,error,*998)
    ENDIF

    CALL Exits("GeneratedMeshRegularInitialise")
    RETURN
999 CALL GeneratedMeshRegularFinalise(generatedMesh%regularMesh,dummyErr,dummyError,*998)
998 CALL Errors("GeneratedMeshRegularInitialise",err,error)
    CALL Exits("GeneratedMeshRegularInitialise")
=======
    ENTERS("GENERATED_MESH_REGULAR_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      IF(ASSOCIATED(GENERATED_MESH%REGULAR_MESH)) THEN
        CALL FlagError("Regular mesh is already associated for this generated mesh.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(GENERATED_MESH%REGULAR_MESH,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate regular generated mesh.",ERR,ERROR,*999)
        GENERATED_MESH%REGULAR_MESH%GENERATED_MESH=>GENERATED_MESH
        GENERATED_MESH%GENERATED_TYPE=GENERATED_MESH_REGULAR_MESH_TYPE
      ENDIF
    ELSE
      CALL FlagError("Generated mesh is not associated.",ERR,ERROR,*998)
    ENDIF

    EXITS("GENERATED_MESH_REGULAR_INITIALISE")
    RETURN
999 CALL GENERATED_MESH_REGULAR_FINALISE(GENERATED_MESH%REGULAR_MESH,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("GENERATED_MESH_REGULAR_INITIALISE",ERR,ERROR)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
    RETURN 1
  END SUBROUTINE GeneratedMeshRegularInitialise

  !
  !================================================================================================================================
  !

  !>Finalise ellipsoid mesh type
  SUBROUTINE GeneratedMeshEllipsoidFinalise(ellipsoidMesh,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshEllipsoidType), POINTER :: ellipsoidMesh !<A pointer to the generated mesh
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

<<<<<<< HEAD
    CALL Enters("GeneratedMeshEllipsoidFinalise",err,error,*999)
=======
    ENTERS("GENERATED_MESH_ELLIPSOID_FINALISE",ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b

    IF(ASSOCIATED(ellipsoidMesh)) THEN
      IF(ALLOCATED(ellipsoidMesh%ellipsoidExtent)) DEALLOCATE(ellipsoidMesh%ellipsoidExtent)
      IF(ALLOCATED(ellipsoidMesh%numberOfElementsXi)) DEALLOCATE(ellipsoidMesh%numberOfElementsXi)
      DEALLOCATE(ellipsoidMesh)
    ENDIF

<<<<<<< HEAD
    CALL Exits("GeneratedMeshEllipsoidFinalise")
    RETURN
    ! TODO invalidate other associations
999 CALL Errors("GeneratedMeshEllipsoidFinalise",err,error)
    CALL Exits("GeneratedMeshEllipsoidFinalise")
=======
    EXITS("GENERATED_MESH_ELLIPSOID_FINALISE")
    RETURN
    ! TODO invalidate other associations
999 ERRORSEXITS("GENERATED_MESH_ELLIPSOID_FINALISE",ERR,ERROR)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
    RETURN 1
  END SUBROUTINE GeneratedMeshEllipsoidFinalise

  !
  !================================================================================================================================
  !

  !>Initialise the ellipsoid generated mesh type
  SUBROUTINE GeneratedMeshEllipsoidInitialise(generatedMesh,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

<<<<<<< HEAD
    CALL Enters("GeneratedMeshEllipsoidInitialise",err,error,*998)

    IF(ASSOCIATED(generatedMesh)) THEN
      IF(ASSOCIATED(generatedMesh%ellipsoidMesh)) THEN
        CALL FlagError("Ellipsoid mesh is already associated for this generated mesh.",err,error,*998)
      ELSE
        IF(generatedMesh%coordinateDimension==3) THEN
          ALLOCATE(generatedMesh%ellipsoidMesh,STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate ellipsoid generated mesh.",err,error,*999)
          generatedMesh%ellipsoidMesh%generatedMesh=>generatedMesh
          generatedMesh%generatedType=GENERATED_MESH_ELLIPSOID_MESH_TYPE
          !Set defaults
          generatedMesh%ellipsoidMesh%closed=.TRUE.
          ALLOCATE(generatedMesh%ellipsoidMesh%ellipsoidExtent(1:generatedMesh%coordinateDimension),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate cylinder mesh cylinder extent.",err,error,*999)
          generatedMesh%ellipsoidMesh%ellipsoidExtent(1)=2.0_DP !Long axis
          generatedMesh%ellipsoidMesh%ellipsoidExtent(1)=1.0_DP !Short axis inner radius
          generatedMesh%ellipsoidMesh%ellipsoidExtent(1)=1.2_DP !Short axis outer radius      
        ENDIF
      ELSE
         CALL FlagError("Invalid coordinate system for an ellipsoid generated mesh. Three coordinate dimensions "// &
           & "are required for an ellipsoid mesh.",err,error,*998)
       ENDIF
     ELSE
      CALL FlagError("Generated mesh is not associated.",err,error,*998)
    ENDIF

    CALL Exits("GeneratedMeshEllipsoidInitialise")
    RETURN
999 CALL GeneratedMeshEllipsoidFinalise(generatedMesh%ellipsoidMesh,dummyErr,dummyError,*998)
998 CALL Errors("GeneratedMeshEllipsoidInitialise",err,error)
    CALL Exits("GeneratedMeshEllipsoidInitialise")
    RETURN 1
  END SUBROUTINE GeneratedMeshEllipsoidInitialise

  !
  !================================================================================================================================
  !

  !>Finalise the fractal tree mesh type.
  SUBROUTINE GeneratedMeshFractalTreeFinalise(fractalTreeMesh,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshFractalTreeType), POINTER :: fractalTreeMesh !<A pointer to the fractal tree generated mesh to finalise.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    CALL Enters("GeneratedMeshFractalTreeFinalise",err,error,*999)

    IF(ASSOCIATED(fractalTreeMesh)) THEN
      IF(ALLOCATED(fractalTreeMesh%numberOfElementsGeneration)) DEALLOCATE(fractalTreeMesh%numberOfElementsGeneration)
      DEALLOCATE(factalTreeMesh)
    ENDIF

    CALL Exits("GeneratedMeshFractalTreeFinalise")
    RETURN
999 CALL Errors("GeneratedMeshFractalTreeFinalise",err,error)
    CALL Exits("GeneratedMeshFractalTreeFinalise")
    RETURN 1
  END SUBROUTINE GeneratedMeshFractalTreeFinalise

  !
  !================================================================================================================================
  !

  !>Initialise the fractal tree generated mesh type
  SUBROUTINE GeneratedMeshFractalTreeInitialise(generatedMesh,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh to initialise the fractal tree mesh for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
    
    CALL Enters("GeneratedMeshFractalTreeInitialise",err,error,*998)

    IF(ASSOCIATED(generatedMesh)) THEN
      IF(ASSOCIATED(generatedMesh%fractalTreeMesh)) THEN
        CALL FlagError("Fractal tree mesh is already associated for this generated mesh.",err,error,*998)
      ELSE
        ALLOCATE(generatedMesh%fractalTreeMesh,STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate fractal tree generated mesh.",err,error,*999)
        generatedMesh%fractalTreeMesh%generatedMesh=>generatedMesh        
        generatedMesh%generatedType=GENERATED_MESH_FRACTAL_TREE_MESH_TYPE
        !Set defaults
        generatedMesh%fractalTreeMesh%fractalTreeType=0
        generatedMesh%fractalTreeMesh%symmetric=.TRUE.
        generatedMesh%fractalTreeMesh%numberOfGenerations=0
      ENDIF
    ELSE
      CALL FlagError("Generated mesh is not associated.",err,error,*998)
    ENDIF

    CALL Exits("GeneratedMeshFractalTreeInitialise")
    RETURN
999 CALL GeneratedMeshFractalTreeFinalise(generatedMesh%fractalTreeMesh,dummyErr,dummyError,*998)
998 CALL Errors("GeneratedMeshFractalTreeInitialise",err,error)
    CALL Exits("GeneratedMeshFractalTreeInitialise")
=======
    ENTERS("GENERATED_MESH_ELLIPSOID_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      IF(ASSOCIATED(GENERATED_MESH%ELLIPSOID_MESH)) THEN
        CALL FlagError("Ellipsoid mesh is already associated for this generated mesh.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(GENERATED_MESH%ELLIPSOID_MESH,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate ellipsoid generated mesh.",ERR,ERROR,*999)
        GENERATED_MESH%ELLIPSOID_MESH%GENERATED_MESH=>GENERATED_MESH
        GENERATED_MESH%GENERATED_TYPE=GENERATED_MESH_ELLIPSOID_MESH_TYPE
      ENDIF
    ELSE
      CALL FlagError("Generated mesh is not associated.",ERR,ERROR,*998)
    ENDIF

    EXITS("GENERATED_MESH_ELLIPSOID_INITIALISE")
    RETURN
999 CALL GENERATED_MESH_ELLIPSOID_FINALISE(GENERATED_MESH%ELLIPSOID_MESH,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("GENERATED_MESH_ELLIPSOID_INITIALISE",ERR,ERROR)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
    RETURN 1
    
  END SUBROUTINE GeneratedMeshFractalTreeInitialise

  !
  !================================================================================================================================
  !

  !>Gets the type of a generated mesh. \see OPENCMISS::CMISSGeneratedMeshTypeGet
  SUBROUTINE GeneratedMeshTypeGet(generatedMesh,meshType,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh to set the type of
    INTEGER(INTG), INTENT(OUT) :: meshType !<On return, the type of mesh to generate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

<<<<<<< HEAD
    CALL Enters("GeneratedMeshTypeGet",err,error,*999)
=======
    ENTERS("GENERATED_MESH_TYPE_GET",ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b

    IF(ASSOCIATED(generatedMesh)) THEN
      meshType=generatedMesh%generatedType
    ELSE
<<<<<<< HEAD
      CALL FlagError("Generated mesh is not associated.",err,error,*999)
    ENDIF

    CALL Exits("GeneratedMeshTypeGet")
    RETURN
999 CALL Errors("GeneratedMeshTypeGet",err,error)
    CALL Exits("GeneratedMeshTypeGet")
=======
      CALL FlagError("Generated mesh is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("GENERATED_MESH_TYPE_GET")
    RETURN
999 ERRORSEXITS("GENERATED_MESH_TYPE_GET",ERR,ERROR)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
    RETURN 1
    
  END SUBROUTINE GeneratedMeshTypeGet

  !
  !================================================================================================================================
  !

  !>Sets/changes the type of a generated mesh. \see OPENCMISS::CMISSGeneratedMeshTypeSet
  SUBROUTINE GeneratedMeshTypeSet(generatedMesh,meshType,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh to set the type of
    INTEGER(INTG), INTENT(IN) :: meshType !<The type of mesh to generate \see GeneratedMeshRoutines_GeneratedMeshTypes,GeneratedMeshRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: oldMeshType
    TYPE(VARYING_STRING) :: localError

<<<<<<< HEAD
    CALL Enters("GeneratedMeshTypeSet",err,error,*999)

    IF(ASSOCIATED(generatedMesh)) THEN
      IF(generatedMesh%generatedMeshFinished) THEN
        CALL FlagError("Generated mesh has already been finished.",err,error,*999)
=======
    ENTERS("GENERATED_MESH_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      IF(GENERATED_MESH%GENERATED_MESH_FINISHED) THEN
        CALL FlagError("Generated mesh has already been finished.",ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
      ELSE
        oldMeshType=generatedMesh%generatedType
        IF(oldMeshType/=meshType) THEN
          !Initialise the new generated mesh type
          SELECT CASE(meshType)
          CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
            CALL GeneratedMeshRegularInitialise(generatedMesh,err,error,*999)
          CASE(GENERATED_MESH_POLAR_MESH_TYPE)
<<<<<<< HEAD
            CALL GeneratedMeshPolarInitialise(generatedMesh,err,error,*999)
=======
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
          CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
            CALL GeneratedMeshCylinderInitialise(generatedMesh,err,error,*999)
          CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
            CALL GeneratedMeshEllipsoidInitialise(generatedMesh,err,error,*999)
          CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
            CALL GeneratedMeshFractalTreeInitialise(generatedMesh,err,error,*999)
          CASE DEFAULT
            localError="The specified generated mesh mesh type of "//TRIM(NumberToVString(meshType,"*",err,error))// &
              & " is invalid."
<<<<<<< HEAD
            CALL FlagError(localError,err,error,*999)
=======
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
          END SELECT
          !Finalise the new generated mesh type
          SELECT CASE(oldMeshType)
          CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
            CALL GeneratedMeshRegularFinalise(generatedMesh%regularMesh,err,error,*999)
          CASE(GENERATED_MESH_POLAR_MESH_TYPE)
<<<<<<< HEAD
            CALL GeneratedMeshPolarFinalise(generatedMesh%polarMesh,err,error,*999)
=======
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
          CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
            CALL GeneratedMeshCylinderFinalise(generatedMesh%cylinderMesh,err,error,*999)
          CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
            CALL GeneratedMeshEllipsoidFinalise(generatedMesh%ellipsoidMesh,err,error,*999)
          CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
            CALL GeneratedMeshFractalTreeFinalise(generatedMesh%fractalTreeMesh,err,error,*999)
          CASE DEFAULT
            localError="The generated mesh mesh type of "//TRIM(NumberToVString(oldMeshType,"*",err,error))// &
              & " is invalid."
<<<<<<< HEAD
            CALL FlagError(localError,err,error,*999)
=======
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
          END SELECT
        ENDIF
      ENDIF
    ELSE
<<<<<<< HEAD
      CALL FlagError("Generated mesh is not associated.",err,error,*999)
    ENDIF

    CALL Exits("GeneratedMeshTypeSet")
    RETURN
999 CALL Errors("GeneratedMeshTypeSet",err,error)
    CALL Exits("GeneratedMeshTypeSet")
=======
      CALL FlagError("Generated mesh is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("GENERATED_MESH_TYPE_SET")
    RETURN
999 ERRORSEXITS("GENERATED_MESH_TYPE_SET",ERR,ERROR)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
    RETURN 1
  END SUBROUTINE GeneratedMeshTypeSet

  !
  !================================================================================================================================
  !

  !>Finds and returns in generated mesh a pointer to that identified by userNumber in the given list of generatedMeshes.
  !>If no generated mesh with that number exists generatedMesh is left nullified.
  SUBROUTINE GeneratedMeshUserNumberFindGeneric(userNumber,generatedMeshes,generatedMesh,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the mesh to find
    TYPE(GeneratedMeshesType), POINTER :: generatedMeshes !<A pointer to the generated meshes to find the user number in
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<On return, a pointer to the generated mesh of the specified user number. In no generated mesh with the specified user number exists the pointer is returned NULL.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: generatedMeshIdx

<<<<<<< HEAD
    CALL Enters("GeneratedMeshUserNumberFindGeneric",err,error,*999)

    IF(ASSOCIATED(generatedMeshes)) THEN
      IF(ASSOCIATED(generatedMesh)) THEN
        CALL FlagError("Generated mesh is already associated.",err,error,*999)
=======
    ENTERS("GENERATED_MESH_USER_NUMBER_FIND_GENERIC",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESHES)) THEN
      IF(ASSOCIATED(GENERATED_MESH)) THEN
        CALL FlagError("Generated mesh is already associated.",ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
      ELSE
        NULLIFY(generatedMesh)
        generatedMeshIdx=1
        DO WHILE(generatedMeshIdx<=generatedMeshes%numberOfGeneratedMeshes.AND..NOT.ASSOCIATED(generatedMesh))
          IF(generatedMeshes%generatedMeshes(generatedMeshIdx)%ptr%userNumber==userNumber) THEN
            generatedMesh=>generatedMeshes%generatedMeshes(generatedMeshIdx)%ptr
            EXIT
          ELSE
            generatedMeshIdx=generatedMeshIdx+1
          ENDIF
        ENDDO
      ENDIF
    ELSE
<<<<<<< HEAD
      CALL FlagError("Generated meshes is not associated.",err,error,*999)
    ENDIF

    CALL Exits("GeneratedMeshUserNumberFindGeneric")
    RETURN
999 CALL Errors("GeneratedMeshUserNumberFindGeneric",err,error)
    CALL Exits("GeneratedMeshUserNumberFindGeneric")
=======
      CALL FlagError("Generated meshes is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("GENERATED_MESH_USER_NUMBER_FIND_GENERIC")
    RETURN
999 ERRORSEXITS("GENERATED_MESH_USER_NUMBER_FIND_GENERIC",ERR,ERROR)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
    RETURN 1
  END SUBROUTINE GeneratedMeshUserNumberFindGeneric

  !
  !================================================================================================================================
  !

  !>Finds and returns in generated mesh a pointer to that identified by userNumber in the given interface.
  !>If no generated mesh with that number exists generatedMesh is left nullified.
  SUBROUTINE GeneratedMeshUserNumberFindInterface(userNumber,interface,generatedMesh,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the mesh to find
    TYPE(INTERFACE_TYPE), POINTER :: interface !<A pointer to the interface containing the generated mesh
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<On return, a pointer to the generated mesh of the specified user number. In no generated mesh with the specified user number exists the pointer is returned NULL.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

<<<<<<< HEAD
    CALL Enters("GeneratedMeshUserNumberFindInterface",err,error,*999)
=======
    ENTERS("GENERATED_MESH_USER_NUMBER_FIND_INTERFACE",ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b

    IF(ASSOCIATED(interface)) THEN
      CALL GeneratedMeshUserNumberFindGeneric(userNumber,interface%GENERATED_MESHES,generatedMesh,err,error,*999)
    ELSE
<<<<<<< HEAD
      CALL FlagError("Interface is not associated.",err,error,*999)
    ENDIF

    CALL Exits("GeneratedMeshUserNumberFindInterface")
    RETURN
999 CALL Errors("GeneratedMeshUserNumberFindInterface",err,error)
    CALL Exits("GeneratedMeshUserNumberFindInterface")
=======
      CALL FlagError("Interface is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("GENERATED_MESH_USER_NUMBER_FIND_INTERFACE")
    RETURN
999 ERRORSEXITS("GENERATED_MESH_USER_NUMBER_FIND_INTERFACE",ERR,ERROR)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
    RETURN 1
  END SUBROUTINE GeneratedMeshUserNumberFindInterface

  !
  !================================================================================================================================
  !

  !>Finds and returns in generated mesh a pointer to that identified by userNumber in the given region.
  !>If no generated mesh with that number exists generatedMesh is left nullified.
  SUBROUTINE GeneratedMeshUserNumberFindRegion(userNumber,region,generatedMesh,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the mesh to find
    TYPE(REGION_TYPE), POINTER :: region !<A pointer to the region containing the generated mesh
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<On return, a pointer to the generated mesh of the specified user number. In no generated mesh with the specified user number exists the pointer is returned NULL.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

<<<<<<< HEAD
    CALL Enters("GeneratedMeshUserNumberFindRegion",err,error,*999)
=======
    ENTERS("GENERATED_MESH_USER_NUMBER_FIND_REGION",ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b

    IF(ASSOCIATED(region)) THEN
      CALL GeneratedMeshUserNumberFindGeneric(userNumber,region%GENERATED_MESHES,generatedMesh,err,error,*999)
    ELSE
<<<<<<< HEAD
      CALL FlagError("Region is not associated.",err,error,*999)
    ENDIF

    CALL Exits("GeneratedMeshUserNumberFindRegion")
    RETURN
999 CALL Errors("GeneratedMeshUserNumberFindRegion",err,error)
    CALL Exits("GeneratedMeshUserNumberFindRegion")
=======
      CALL FlagError("Region is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("GENERATED_MESH_USER_NUMBER_FIND_REGION")
    RETURN
999 ERRORSEXITS("GENERATED_MESH_USER_NUMBER_FIND_REGION",ERR,ERROR)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
    RETURN 1

  END SUBROUTINE GeneratedMeshUserNumberFindRegion

  !
  !================================================================================================================================
  !

  !>Finalises all generated meshes and deallocates all memory.
  SUBROUTINE GeneratedMeshesFinalise(generatedMeshes,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshesType), POINTER :: generatedMeshes !<A pointer to the generated meshes to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh

<<<<<<< HEAD
    CALL Enters("GeneratedMeshesFinalise",err,error,*999)
=======
    ENTERS("GENERATED_MESHES_FINALISE",ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b

    IF(ASSOCIATED(generatedMeshes)) THEN
      DO WHILE(generatedMeshes%numberOfGeneratedMeshes>0)
        generatedMesh=>generatedMeshes%generatedMeshes(1)%ptr
        CALL GeneratedMeshDestroy(generatedMesh,err,error,*999)
      ENDDO !generatedMeshIdx
      DEALLOCATE(generatedMeshes)
    ENDIF

<<<<<<< HEAD
    CALL Exits("GeneratedMeshesFinalise")
    RETURN
999 CALL Errors("GeneratedMeshesFinalise",err,error)
    CALL Exits("GeneratedMeshesFinalise")
=======
    EXITS("GENERATED_MESHES_FINALISE")
    RETURN
999 ERRORSEXITS("GENERATED_MESHES_FINALISE",ERR,ERROR)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
    RETURN 1

  END SUBROUTINE GeneratedMeshesFinalise

  !
  !================================================================================================================================
  !

  !>Intialises all generated meshes.
  SUBROUTINE GeneratedMeshesInitialiseGeneric(generatedMeshes,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshesType), POINTER :: generatedMeshes !<A pointer to the generated meshes to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

<<<<<<< HEAD
    CALL Enters("GeneratedMeshesInitialiseGeneric",err,error,*998)

    IF(ASSOCIATED(generatedMeshes)) THEN
      CALL FlagError("Generated meshes is already associated.",err,error,*998)
    ELSE
      ALLOCATE(generatedMeshes,STAT=err)
      IF(err/=0) CALL FlagError("Generated meshes is not associated.",err,error,*999)
      generatedMeshes%numberOfGeneratedMeshes=0
      NULLIFY(generatedMeshes%region)
      NULLIFY(generatedMeshes%interface)
    ENDIF

    CALL Exits("GeneratedMeshesInitialiseGeneric")
    RETURN
999 CALL GeneratedMeshesFinalise(generatedMeshes,dummyErr,dummyError,*998)
998 CALL Errors("GeneratedMeshesInitialiseGeneric",err,error)
    CALL Exits("GeneratedMeshesInitialiseGeneric")
=======
    ENTERS("GENERATED_MESHES_INITIALISE_GENERIC",ERR,ERROR,*998)

    IF(ASSOCIATED(GENERATED_MESHES)) THEN
      CALL FlagError("Generated meshes is already associated.",ERR,ERROR,*998)
    ELSE
      ALLOCATE(GENERATED_MESHES,STAT=ERR)
      IF(ERR/=0) CALL FlagError("Generated meshes is not associated.",ERR,ERROR,*999)
      GENERATED_MESHES%NUMBER_OF_GENERATED_MESHES=0
      NULLIFY(GENERATED_MESHES%GENERATED_MESHES)
      NULLIFY(GENERATED_MESHES%REGION)
      NULLIFY(GENERATED_MESHES%INTERFACE)
    ENDIF

    EXITS("GENERATED_MESHES_INITIALISE_GENERIC")
    RETURN
999 CALL GENERATED_MESHES_FINALISE(GENERATED_MESHES,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("GENERATED_MESHES_INITIALISE_GENERIC",ERR,ERROR)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
    RETURN 1
  END SUBROUTINE GeneratedMeshesInitialiseGeneric

  !
  !================================================================================================================================
  !

  !>Intialises the generated meshes for an interface.
  SUBROUTINE GeneratedMeshesInitialiseInterface(interface,err,error,*)

    !Argument variables
    TYPE(INTERFACE_TYPE), POINTER :: interface !<A pointer to the interface to initialise the generated meshes for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

<<<<<<< HEAD
    CALL Enters("GeneratedMeshesInitialiseInterface",err,error,*999)

    IF(ASSOCIATED(interface)) THEN
      IF(ASSOCIATED(interface%GENERATED_MESHES)) THEN
        CALL FlagError("Interface generated meshes is already associated.",err,error,*999)
=======
    ENTERS("GENERATED_MESHES_INITIALISE_INTERFACE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE)) THEN
      IF(ASSOCIATED(INTERFACE%GENERATED_MESHES)) THEN
        CALL FlagError("Interface generated meshes is already associated.",ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
      ELSE
        CALL GeneratedMeshesInitialiseGeneric(interface%GENERATED_MESHES,err,error,*999)
        interface%GENERATED_MESHES%interface=>interface
      ENDIF
    ELSE
<<<<<<< HEAD
      CALL FlagError("Interface is not associated.",err,error,*999)
    ENDIF

    CALL Exits("GeneratedMeshesInitialiseInterface")
    RETURN
999 CALL Errors("GeneratedMeshesInitialiseInterface",err,error)
    CALL Exits("GeneratedMeshesInitialiseInterface")
=======
      CALL FlagError("Interface is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("GENERATED_MESHES_INITIALISE_INTERFACE")
    RETURN
999 ERRORSEXITS("GENERATED_MESHES_INITIALISE_INTERFACE",ERR,ERROR)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
    RETURN 1
  END SUBROUTINE GeneratedMeshesInitialiseInterface

  !
  !================================================================================================================================
  !

  !>Intialises the generated meshes for a region.
  SUBROUTINE GeneratedMeshesInitialiseRegion(region,err,error,*)

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: region !<A pointer to the region
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

<<<<<<< HEAD
    CALL Enters("GeneratedMeshesInitialiseRegion",err,error,*999)

    IF(ASSOCIATED(region)) THEN
      IF(ASSOCIATED(region%GENERATED_MESHES)) THEN
        CALL FlagError("Region generated meshes is already associated.",err,error,*999)
=======
    ENTERS("GENERATED_MESHES_INITIALISE_REGION",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(REGION%GENERATED_MESHES)) THEN
        CALL FlagError("Region generated meshes is already associated.",ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
      ELSE
        CALL GeneratedMeshesInitialiseGeneric(region%GENERATED_MESHES,err,error,*999)
        region%GENERATED_MESHES%region=>region
      ENDIF
    ELSE
<<<<<<< HEAD
      CALL FlagError("Region is not associated.",err,error,*999)
    ENDIF

    CALL Exits("GeneratedMeshesInitialiseRegion")
    RETURN
999 CALL Errors("GeneratedMeshesInitialiseRegion",err,error)
    CALL Exits("GeneratedMeshesInitialiseRegion")
=======
      CALL FlagError("Region is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("GENERATED_MESHES_INITIALISE_REGION")
    RETURN
999 ERRORSEXITS("GENERATED_MESHES_INITIALISE_REGION",ERR,ERROR)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
    RETURN 1
  END SUBROUTINE GeneratedMeshesInitialiseRegion

  !
  !================================================================================================================================
  !

  !>Updates the geometric field parameters from the initial nodal positions of the mesh. \see OPENCMISS::CMISSGeneratedMeshGeometricParametersCalculate
<<<<<<< HEAD
  SUBROUTINE GeneratedMeshGeometricParametersCalculate(generatedMesh,geometricField,err,error,*)
=======
  SUBROUTINE GeneratedMesh_GeometricParametersCalculate(FIELD,GENERATED_MESH,ERR,ERROR,*)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<The mesh which is generated by the generated mesh \todo is this necessary???
    TYPE(FIELD_TYPE), POINTER :: geometricField !<A pointer to the field to update the geometric parameters for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
<<<<<<< HEAD
    INTEGER(INTG) :: componentIdx
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable
    TYPE(REGION_TYPE), POINTER :: meshRegion,fieldRegion
    TYPE(VARYING_STRING) :: localError

    NULLIFY(meshRegion)
    NULLIFY(fieldRegion)
    NULLIFY(fieldVariable)

    CALL Enters("GeneratedMeshGeometricParametersCalculate",err,error,*999)

    IF(ASSOCIATED(generatedMesh)) THEN
      IF(generatedMesh%generatedMeshFinished) THEN
        IF(ASSOCIATED(geometricField)) THEN
          IF(geometricField%FIELD_FINISHED) THEN
            !Check that the field is a geometric field.
            IF(geometricField%TYPE==FIELD_GEOMETRIC_TYPE) THEN
              !Check the mesh and field have the same region
              CALL GeneratedMeshRegionGet(generatedMesh,meshRegion,err,error,*999)
              CALL FIELD_REGION_GET(geometricField,fieldRegion,err,error,*999)
              IF(ASSOCIATED(meshRegion,fieldRegion)) THEN
                !Check that the geometric interpolation types are nodally based.
                CALL FIELD_VARIABLE_GET(geometricField,FIELD_U_VARIABLE_TYPE,fieldVariable,err,error,*999)
                DO componentIdx=1,fieldVariable%NUMBER_OF_COMPONENTS
                  IF(fieldVariable%components(componentIdx)%INTERPOLATION_TYPE/=FIELD_NODE_BASED_INTERPOLATION) THEN
                    localError="Component number "//TRIM(NumberToVString(componentIdx,"*",err,error))// &
                      & " of geometric field number "//TRIM(NumberToVString(geometricField%USER_NUMBER,"*",err,error))// &
                      & " is not nodally interpolated."
                    CALL FlagError(localError,err,error,*999)
                  ENDIF
                ENDDO !componentIdx
                !End of checks, calculate the geometric parameters for each mesh type
                SELECT CASE(generatedMesh%generatedType)
                CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
                  CALL GeneratedMeshRegularGeometricParametersCalculate(generatedMesh%regularMesh,geometricField,err,error,*999)
                CASE(GENERATED_MESH_POLAR_MESH_TYPE)
                  CALL GeneratedMeshPolarGeometricParametersCalculate(generatedMesh%polarMesh,geometricField,err,error,*999)
                CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
                  CALL GeneratedMeshCylinderGeometricParametersCalculate(generatedMesh%cylinderMesh,geometricField,err,error,*999)
                CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
                  CALL GeneratedMeshEllipsoidGeometricParametersCalculate(generatedMesh%ellipsoidMesh,geometricField,err,error,*999)
                CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
                   CALL GeneratedMeshFractalTreeGeometricParametersCalculate(generatedMesh%fractalTreeMesh,geometricField, &
                     & err,error,*999)
                CASE DEFAULT
                  localError="The generated mesh mesh type of "// &
                    & TRIM(NumberToVString(generatedMesh%generatedType,"*",err,error))//" is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              ELSE
                localError="The generated mesh region user number of "// &
                  & TRIM(NumberToVString(meshRegion%USER_NUMBER,"*",err,error))// &
                  & " does not match the geometric field region user number of "//&
                  & TRIM(NumberToVString(fieldRegion%USER_NUMBER,"*",err,error))//"."
                CALL FlagError(localError,err,error,*999)
              ENDIF
            ELSE
              localError="Field number "//TRIM(NumberToVString(geometricField%USER_NUMBER,"*",err,error))// &
                & " is not a geometric field."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ELSE
            localError="Field number "//TRIM(NumberToVString(geometricField%USER_NUMBER,"*",err,error))//" has not been finished."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Field is not associated.",err,error,*999)
        ENDIF
      ELSE
          CALL FlagError("Generated mesh has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Generated mesh is not associated.",err,error,*999)
    ENDIF

    CALL Exits("GeneratedMeshGeometricParametersCalculate")
    RETURN
999 CALL Errors("GeneratedMeshGeometricParametersCalculate",err,error)
    CALL Exits("GeneratedMeshGeometricParametersCalculate")
    RETURN 1
    
  END SUBROUTINE GeneratedMeshGeometricParametersCalculate
=======
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("GeneratedMesh_GeometricParametersCalculate",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(ASSOCIATED(GENERATED_MESH)) THEN
          SELECT CASE(GENERATED_MESH%GENERATED_TYPE)
          CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
            CALL GeneratedMesh_RegularGeometricParametersCalculate(GENERATED_MESH%REGULAR_MESH, &
                & FIELD,ERR,ERROR,*999)
          CASE(GENERATED_MESH_POLAR_MESH_TYPE)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
            CALL GeneratedMesh_CylinderGeometricParametersCalculate(GENERATED_MESH%CYLINDER_MESH, &
                & FIELD,ERR,ERROR,*999)
          CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
            CALL GeneratedMesh_EllipsoidGeometricParametersCalculate(GENERATED_MESH%ELLIPSOID_MESH, &
                & FIELD,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The generated mesh mesh type of "// &
              & TRIM(NUMBER_TO_VSTRING(GENERATED_MESH%GENERATED_TYPE,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            CALL FlagError("Generated mesh type is either invalid or not implemented.",ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FlagError("Generated mesh is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" has not been finished."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("GeneratedMesh_GeometricParametersCalculate")
    RETURN
999 ERRORSEXITS("GeneratedMesh_GeometricParametersCalculate",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE GeneratedMesh_GeometricParametersCalculate
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b

  !
  !================================================================================================================================
  !

  !>Returns the region for a generated mesh accounting for regions and interfaces
  SUBROUTINE GeneratedMeshRegionGet(generatedMesh,region,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh to get the region for
    TYPE(REGION_TYPE), POINTER :: region !<On return, the generated meshes region. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(INTERFACE_TYPE), POINTER :: interface
    TYPE(REGION_TYPE), POINTER :: parentRegion
    TYPE(VARYING_STRING) :: localError

<<<<<<< HEAD
    CALL Enters("GeneratedMeshRegionGet",err,error,*999)

    IF(ASSOCIATED(generatedMesh)) THEN
      IF(ASSOCIATED(region)) THEN
        CALL FlagError("Region is already associated.",err,error,*999)
=======
    ENTERS("GENERATED_MESH_REGION_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      IF(ASSOCIATED(REGION)) THEN
        CALL FlagError("Region is already associated.",ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
      ELSE
        NULLIFY(region)
        NULLIFY(interface)
        region=>generatedMesh%region
        IF(.NOT.ASSOCIATED(region)) THEN
          interface=>generatedMesh%interface
          IF(ASSOCIATED(interface)) THEN
            parentRegion=>interface%PARENT_REGION
            IF(ASSOCIATED(parentRegion)) THEN
              region=>parentRegion
            ELSE
<<<<<<< HEAD
              localError="The parent region not associated for generated mesh number "// &
                & TRIM(NumberToVString(generatedMesh%userNumber,"*",err,error))//" of interface number "// &
                & TRIM(NumberToVString(interface%USER_NUMBER,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ELSE
            localError="The region or interface is not associated for generated mesh number "// &
              & TRIM(NumberToVString(generatedMesh%userNumber,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
=======
              LOCAL_ERROR="The parent region not associated for generated mesh number "// &
                & TRIM(NUMBER_TO_VSTRING(GENERATED_MESH%USER_NUMBER,"*",ERR,ERROR))//" of interface number "// &
                & TRIM(NUMBER_TO_VSTRING(INTERFACE%USER_NUMBER,"*",ERR,ERROR))//"."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The region or interface is not associated for generated mesh number "// &
              & TRIM(NUMBER_TO_VSTRING(GENERATED_MESH%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
          ENDIF
        ENDIF
      ENDIF
    ELSE
<<<<<<< HEAD
      CALL FlagError("Generated mesh is not associated.",err,error,*999)
    ENDIF

    CALL Exits("GeneratedMeshRegionGet")
    RETURN
999 CALL Errors("GeneratedMeshRegionGet",err,error)
    CALL Exits("GeneratedMeshRegionGet")
=======
      CALL FlagError("Generated mesh is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("GENERATED_MESH_REGION_GET")
    RETURN
999 ERRORSEXITS("GENERATED_MESH_REGION_GET",ERR,ERROR)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
    RETURN 1
    
  END SUBROUTINE GeneratedMeshRegionGet

  !
  !
  !================================================================================================================================
  !

  !>Updates the geometric field parameters from the initial nodal positions of the regular mesh.
<<<<<<< HEAD
  SUBROUTINE GeneratedMeshRegularGeometricParametersCalculate(regularMesh,geometricField,err,error,*)
=======
  SUBROUTINE GeneratedMesh_RegularGeometricParametersCalculate(REGULAR_MESH,FIELD,ERR,ERROR,*)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b

    !Argument variables
    TYPE(GeneratedMeshRegularType), POINTER :: regularMesh !<A pointer to the regular mesh object
    TYPE(FIELD_TYPE), POINTER :: geometricField !<A pointer to the field to update the geometric parameters for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    !Local variables
<<<<<<< HEAD
    INTEGER(INTG) :: componentIdx,derivativeIdx,componentNode,meshComponent,nodeIdx,nodePositionIdx(3), &
      & totalNumberOfNodesXi(3),xiIdx,nodeUserNumber,coordinateType,coordinateDimension,scalingType
    REAL(DP) :: deltaCoordinate(3,3),VALUE
    REAL(DP) :: derivativeValues(MAXIMUM_GLOBAL_DERIV_NUMBER)
    TYPE(BASIS_TYPE), POINTER :: basis
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: coordinateSystem
    TYPE(DOMAIN_TYPE), POINTER :: domain
    TYPE(DOMAIN_NODES_TYPE), POINTER :: domainNodes
    TYPE(FIELD_VARIABLE_COMPONENT_TYPE), POINTER :: fieldVariableComponent
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable
    TYPE(VARYING_STRING) :: localError
    LOGICAL :: nodeExists,ghostNode

    NULLIFY(coordinateSystem)
    NULLIFY(fieldVariable)

    CALL Enters("GeneratedMeshRegularGeometricParametersCalculate",err,error,*999)

    IF(ASSOCIATED(regularMesh)) THEN
      IF(ASSOCIATED(geometricField)) THEN
        generatedMesh=>regularMesh%generatedMesh
        IF(ASSOCIATED(generatedMesh)) THEN
          CALL GeneratedMeshCoordinateSystemGet(generatedMesh,coordinateSystem,err,error,*999)
          CALL COORDINATE_SYSTEM_TYPE_GET(coordinateSystem,coordinateType,err,error,*999)         
          IF(coordinateType==COORDINATE_RECTANGULAR_CARTESIAN_TYPE) THEN
            coordinateDimension=generatedMesh%coordinateDimension
            deltaCoordinate=0.0_DP
            totalNumberOfNodesXi=1
            CALL FIELD_VARAIBLE_GET(geometricField,FIELD_U_VARIABLE_TYPE,fieldVariable,err,error,*999)
            DO componentIdx=1,fieldVariable%NUMBER_OF_COMPONENTS
              fieldVariableComponent=>fieldVariable%components(componentIdx)
              meshComponent=fieldVariableComponent%MESH_COMPONENT_NUMBER
              basis=>regularMesh%bases(meshComponent)%ptr
              !Calculate the total number of nodes in each xi direction
              DO xiIdx=1,generatedMesh%meshDimension
                totalNumberOfNodesXi(xiIdx)=(basis%NUMBER_OF_NODES_XIC(xiIdx)-2)*regularMesh%numberOfElementsXi(xiIdx)+ &
                  & regularMesh%numberOfElementsXi(xiIdx)+1
              ENDDO !xiIdx
              !Calculate delta
              DO xiIdx=1,generatedMesh%meshDimension
                deltaCoordinate(1:coordinateDimension,xiIdx)=generatedMesh%baseVectors(1:coordinateDimension,xiIdx)/ &
                & REAL(totalNumberOfNodesXi(xiIdx)-1,DP)
            ENDDO !xiIdx
            derivativeValues=0.0_DP
            CALL FIELD_SCALING_TYPE_GET(geometricField,scalingType,err,error,*999)
            SELECT CASE(scalingType)
            CASE(FIELD_NO_SCALING,FIELD_UNIT_SCALING)
              IF(regularMesh%numberOfElementsXi(1)>0) THEN
                derivativeValues(GLOBAL_DERIV_S1)=generatedMesh%baseVectors(componentIdx,1)/regularMesh%numberOfElementsXi(1)
              END IF
              IF(generatedMesh%meshDimension>1) THEN
                IF(regularMesh%numberOfElementsXi(2)>0) THEN
                  derivativeValues(GLOBAL_DERIV_S2)=generatedMesh%baseVectors(componentIdx,2)/regularMesh%numberOfElementsXi(2)
                END IF
              ENDIF
              IF(generatedMesh%meshDimension>2) THEN
                IF(regularMesh%numberOfElementsXi(3)>0) THEN
                  derivativeValues(GLOBAL_DERIV_S3)=generatedMesh%baseVectors(componentIdx,3)/regularMesh%numberOfElementsXi(3)
                END IF
              ENDIF
            CASE(FIELD_ARC_LENGTH_SCALING,FIELD_ARITHMETIC_MEAN_SCALING, &
              & FIELD_GEOMETRIC_MEAN_SCALING,FIELD_HARMONIC_MEAN_SCALING)
              !Arc length scalings
              IF(regularMesh%numberOfElementsXi(1)>0) THEN
                derivativeValues(GLOBAL_DERIV_S1)=generatedMesh%baseVectors(componentIdx,1)/L2Norm(generatedMesh%baseVectors(:,1))
              END IF
              IF(generatedMesh%meshDimension>1) THEN
                IF(regularMesh%numberOfElementsXi(2)>0) THEN
                  derivativeValues(GLOBAL_DERIV_S2)=generatedMesh%baseVectors(componentIdx,2)/L2Norm(generatedMesh%baseVectors(:,2))
                END IF
              ENDIF
              IF(regularMesh%meshDimension>2) THEN
                IF(regularMesh%numberOfElementsXi(3)>0) THEN
                  derivativeValues(GLOBAL_DERIV_S3)=regularMesh%baseVectors(componentIdx,3)/L2Norm(regularMesh%baseVectors(:,3))
                END IF
              ENDIF
            CASE DEFAULT
              localError="The scaling type of "// &
                & TRIM(NumberToVString(geometricField%SCALINGS%SCALING_TYPE,"*",err,error))// &
                & " for geometric field number "//TRIM(NumberToVString(geometricField%USER_NUMBER,"*",err,error))// &
                & " is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
            !Update geometric parameters in this computational domain only
            domain=>fieldVariableComponent%domain
            domainNodes=>domain%TOPOLOGY%NODES
            DO componentNode=1,totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)*totalNumberOfNodesXi(3)
              !Regular meshes with Lagrange/Hermite elements use different node numberings to other mesh types
              IF(regularMesh%bases(meshComponent)%ptr%TYPE==BASIS_LAGRANGE_HERMITE_TP_TYPE) THEN
                CALL GeneratedMeshRegularComponentNodeToUserNumber(regularMesh%generatedMesh,meshComponent, &
                  & componentNode,nodeUserNumber,err,error,*999)
              ELSE
                nodeUserNumber=GeneratedMeshComponentNodeToUserNumber(regularMesh%generatedMesh,meshComponent,componentNode, &
                  & err,error)
              END IF
              CALL DOMAIN_TOPOLOGY_NODE_CHECK_EXISTS(fieldVariableComponent%domain%topology,nodeUserNumber, &
                nodeExists,nodeIdx,ghostNode,err,error,*999)
              IF(nodeExists.AND..NOT.ghostNode) THEN
                nodePositionIdx(3)=(componentNode-1)/(totalNumberOfNodesXi(2)*totalNumberOfNodesXi(1))+1
                nodePositionIdx(2)=MOD(componentNode-1,totalNumberOfNodesXi(2)*totalNumberOfNodesXi(1))/ &
                  & totalNumberOfNodesXi(1)+1
                nodePositionIdx(1)=MOD(MOD(componentNode-1,totalNumberOfNodesXi(2)*totalNumberOfNodesXi(1)), &
                  & totalNumberOfNodesXi(1))+1
                value=0.0_DP
                DO xiIdx=1,regularMesh%meshDimension
                  value=value+REAL(nodePositionIdx(xiIdx)-1,DP)*deltaCoordinate(componentIdx,xiIdx)
                ENDDO !xiIdx
                value=regularMesh%origin(componentIdx)+value
                CALL FIELD_PARAMETER_SET_UPDATE_NODE(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & 1,1,nodeUserNumber,componentIdx,VALUE,err,error,*999)
                !Set derivatives
                DO derivativeIdx=2,domainNodes%nodes(nodeIdx)%NUMBER_OF_DERIVATIVES
                  CALL FIELD_PARAMETER_SET_UPDATE_NODE(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                    & 1,derivativeIdx,nodeUserNumber,componentIdx,derivativeValues(derivativeIdx),err,error,*999)
                ENDDO !derivativeIdx
              ENDIF !node_exists
            ENDDO !nodeIdx
          ENDDO !componentIdx
!!TODO: do boundary nodes first then start the update to overlap computation and computation.
          CALL FIELD_PARAMETER_SET_UPDATE_START(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
          CALL FIELD_PARAMETER_SET_UPDATE_FINISH(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
        ELSE
          CALL FlagError("Non rectangular Cartesian coordinate systems are not implemented.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Field is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Regular mesh is not associated.",err,error,*999)
    ENDIF

    CALL Exits("GeneratedMeshRegularGeometricParametersCalculate")
    RETURN
999 CALL Errors("GeneratedMeshRegularGeometricParametersCalculate",err,error)
    CALL Exits("GeneratedMeshRegularGeometricParametersCalculate")
    RETURN 1
  END SUBROUTINE GeneratedMeshRegularGeometricParametersCalculate

  !
  !================================================================================================================================
  !

  !>Updates the geometric field parameters from the initial nodal positions of a polar generated mesh. 
  SUBROUTINE GeneratedMeshPolarGeometricParametersCalculate(polarMesh,geometricField,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshPolarType), POINTER :: polarMesh !<A pointer to the polar mesh object to calculate the geometric parameters for.
    TYPE(FIELD_TYPE), POINTER :: geometricField !<A pointer to the field to update the geometric parameters for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(BASIS_TYPE), POINTER :: basis
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: coordinateSystem
    TYPE(DOMAIN_TYPE),POINTER :: domain
    TYPE(DOMAIN_NODES_TYPE), POINTER :: domainNodes
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable
    TYPE(FIELD_VARIABLE_COMPONENT_TYPE), POINTER :: fieldVariableComponent
    INTEGER(INTG) :: totalNumberOfNodesXi(3)
    INTEGER(INTG) :: componentIdx,xiIdx
    INTEGER(INTG) :: nodeIdx,globalNodeIdx,componentNodeIdx,dofIdx,derivativeIdx
    INTEGER(INTG) :: numberOfPlanarNodes,scalingType,meshComponent,coordinateType,coordinateDimension
    INTEGER(INTG) :: nodePositionIdx(3) ! holds r,theta,z indices
    REAL(DP) :: deltaCoordinate(3),deltaCoordinateXi(3),polarCoordinates(3),rectangularCoordinates(3),deriv
    TYPE(VARYING_STRING) :: localError

    NULLIFY(coordinateSystem)
    NULLIFY(fieldVariable)

    CALL Enters("GeneratedMeshPolarGeometricParametersCalculate",err,error,*999)

    IF(ASSOCIATED(polarMesh)) THEN
      IF(ASSOCIATED(geometricField)) THEN
        generatedMesh=>polarMesh%generatedMesh
        IF(ASSOCIATED(generatedMesh)) THEN
          CALL GeneratedMeshCoordinateSystemGet(generatedMesh,coordinateSystem,err,error,*999)
          CALL COORDINATE_SYSTEM_TYPE_GET(coordinateSystem,coordinateType,err,error,*999)
          IF(coordinateType==COORDINATE_RECTANGULAR_CARTESIAN_TYPE) THEN
            CALL FIELD_SCALING_TYPE_GET(geometricField,scalingType,err,error,*999)
            CALL FIELD_VARIABLE_TYPE_GET(geometricField,FIELD_U_VARIABLE_TYPE,fieldVariable,err,error,*999)
            IF(fieldVariable%NUMBER_OF_COMPONENTS==generatedMesh%coordinateDimension) THEN
              !Get the first basis type
              basis=>generatedMesh%bases(1)%ptr
              basisType=basis%type
              DO componentIdx=1,generatedMesh%coordinateDimension
                fieldVariableComponent=>fieldVariable%components(componentIdx)
                meshComponent=fieldVariableComponent%MESH_COMPONENT_NUMBER
                IF(generatedMesh%coordinateDimension==2) THEN
                  basisIdx=meshComponent
                ELSE
                  
                ENDIF
                basis=>generatedMesh%bases(meshComponent)%ptr
                !Calculate the total number of nodes in each xi direction
                DO xiIdx=1,3
                  totalNumberOfNodesXi(xiIdx)=(basis%NUMBER_OF_NODES_XIC(xiIdx)-1)*cylinderMesh%numberOfElementsXi(xiIdx)+1
                ENDDO !xiIdx
                totalNumberOfNodesXi(2)=totalNumberOfNodesXi(2)-1 !Theta loops around so slightly different
                numberOfPlanarNodes=totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                !Calculate delta
                deltaCoordinate(1)=(cylinderMesh%cylinderExtent(2)-cylinderMesh%cylinderExtent(1))/ &
                  & cylinderMesh%numberOfElementsXi(1)
                deltaCoordinate(2)=TWOPI/cylinderMesh%numberOfElementsXi(2)
                deltaCoordinate(3)=cylinderMesh%cylinderExtent(3)/cylinderMesh%numberOfElementsXi(3)
                DO xiIdx=1,3
                  deltaCoordinateXi(xiIdx)=deltaCoordinate(xiIdx)/(basis%NUMBER_OF_NODES_XIC(xiIdx)-1)
                ENDDO !xiIdx
                !Update geometric parameters in this computational domain only
                domain=>fieldVariable%components(meshComponent)%domain
                domainNodes=>domain%topology%nodes
                DO nodeIdx=1,domainNodes%NUMBER_OF_NODES
                  globalNodeIdx=domainNodes%nodes(nodeIdx)%GLOBAL_NUMBER
                  componentNodeIdx=GeneratedMeshUserNumberToComponentNode(cylinderMesh%generatedMesh,meshComponent, &
                    & globalNodeIdx,err,error)
                  !Calculate nodePositionIdx which will be used to calculate (r,theta,z) then (x,y,z)
                  nodePositionIdx(3)=(componentNodeIdx-1)/numberOfPlanarNodes
                  nodePositionIdx(2)=(componentNodeIdx-1-(nodePositionIdx(3))*numberOfPlanarNodes)/totalNumberOfNodesXi(1)
                  nodePositionIdx(1)=MOD(componentNodeIdx-1-(nodePositionIdx(3))*numberOfPlanarNodes,totalNumberOfNodesXi(1))
                  DO xiIdx=1,3
                    polarCoordinates(xiIdx)=nodePositionIdx(xiIdx)*deltaCoordinateXi(xiIdx)
                  ENDDO !xiIdx
                  polarCoordinates(1)=nodePositionIdx(1)*deltaCoordinateXi(1)+cylinderMesh%cylinderExtent(1) !Add the inner radius
                  rectangularCoordinates(1)=polarCoordinates(1)*COS(polarCoordinates(2))
                  rectangularCoordinates(2)=polarCoordinates(1)*SIN(polarCoordinates(2))
                  rectangularCoordinates(3)=polarCoordinates(3)
                  rectangularCoordinates=rectangularCoordinates+cylinderMesh%origin
                  !Default to version 1 of each node derivative
                  dofIdx=fieldVariableComponent%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(nodeIdx)%DERIVATIVES(1)%VERSIONS(1)
                  CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,dofIdx, &
                    & rectangularCoordinates(componentIdx),err,error,*999)
                  IF(domainNodes%NODES(nodeIdx)%NUMBER_OF_DERIVATIVES>1) THEN
                    DO derivativeIdx=2,fieldVariableComponent%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(nodeIdx)% &
                      & NUMBER_OF_DERIVATIVES
                      SELECT CASE(scalingType)
                      CASE(FIELD_NO_SCALING,FIELD_UNIT_SCALING)
                        SELECT CASE(domainNodes%nodes(nodeIdx)%derivatives(derivativeIdx)%GLOBAL_DERIVATIVE_INDEX)
                        CASE(GLOBAL_DERIV_S1)
                          SELECT CASE(componentIdx)
                          CASE(1)
                            DERIV=COS(polarCoordinates(2))*deltaCoordinate(1)
                          CASE(2)
                            DERIV=SIN(polarCoordinates(2))*deltaCoordinate(1)
                          CASE DEFAULT
                            DERIV=0.0_DP
                          END SELECT
                        CASE(GLOBAL_DERIV_S2)
                          SELECT CASE(componentIdx)
                          CASE(1)
                            DERIV=-polarCoordinates(1)*SIN(polarCoordinates(2))*deltaCoordinate(2)
                          CASE(2)
                            DERIV=polarCoordinates(1)*COS(polarCoordinates(2))*deltaCoordinate(2)
                          CASE DEFAULT
                            DERIV=0.0_DP
                          END SELECT
                        CASE(GLOBAL_DERIV_S3)
                          IF(componentIdx==3) THEN
                            DERIV=deltaCoordinate(3)
                          ELSE
                            DERIV=0.0_DP
                          ENDIF
                        CASE(GLOBAL_DERIV_S1_S2)
                          SELECT CASE(componentIdx)
                          CASE(1)
                            DERIV=-SIN(polarCoordinates(2))*deltaCoordinate(1)*deltaCoordinate(2)
                          CASE(2)
                            DERIV=COS(polarCoordinates(2))*deltaCoordinate(1)*deltaCoordinate(2)
                          CASE DEFAULT
                            DERIV=0.0_DP
                          END SELECT
                        CASE DEFAULT  ! all other non-xy-planar cross derivatives
                          DERIV=0.0_DP
                        END SELECT
                      CASE(FIELD_ARC_LENGTH_SCALING,FIELD_ARITHMETIC_MEAN_SCALING, &
                        & FIELD_GEOMETRIC_MEAN_SCALING,FIELD_HARMONIC_MEAN_SCALING)
                        SELECT CASE(domainNodes%NODES(nodeIdx)%DERIVATIVES(derivativeIdx)%GLOBAL_DERIVATIVE_INDEX)
                        CASE(GLOBAL_DERIV_S1)
                          SELECT CASE(componentIdx)
                          CASE(1)
                            DERIV=COS(polarCoordinates(2))
                          CASE(2)
                            DERIV=SIN(polarCoordinates(2))
                          CASE DEFAULT
                            DERIV=0.0_DP
                          END SELECT
                        CASE(GLOBAL_DERIV_S2)
                          SELECT CASE(componentIdx)
                          CASE(1)
                            DERIV=-SIN(polarCoordinates(2))
                          CASE(2)
                            DERIV=COS(polarCoordinates(2))
                          CASE DEFAULT
                            DERIV=0.0_DP
                          END SELECT
                        CASE(GLOBAL_DERIV_S3)
                          IF(componentIdx==3) THEN
                            DERIV=1.0_DP
                          ELSE
                            DERIV=0.0_DP
                          ENDIF
                        CASE(GLOBAL_DERIV_S1_S2)
                          SELECT CASE(componentIdx)
                          CASE(1)
                            DERIV=-SIN(polarCoordinates(2))*deltaCoordinate(1)*deltaCoordinate(2)
                          CASE(2)
                            DERIV=COS(polarCoordinates(2))*deltaCoordinate(1)*deltaCoordinate(2)
                          CASE DEFAULT
                            DERIV=0.0_DP
                          END SELECT
                        CASE DEFAULT  ! all other non-xy-planar cross derivatives
                          DERIV=0.0_DP
                        END SELECT
                      CASE DEFAULT
                        localError="The field scaling type of "//TRIM(NumberToVString(scalingType,"*",err,error))// &
                          & " is invalid for field number "//TRIM(NumberToVString(geometricField%USER_NUMBER,"*",err,error))//"."
                        CALL FlagError(localError,err,error,*999)
                      END SELECT
                      !Assign derivative. Default to version 1 of each node derivative
                      dofIdx=fieldVariableComponent%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(nodeIdx)%DERIVATIVES(derivativeIdx)% &
                        & VERSIONS(1)
                      CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                        & dofIdx,DERIV,err,error,*999)
                    ENDDO !derivativeIdx
                  ENDIF !derivatives
                ENDDO !nodeIdx
              ELSE
                CALL FlagError("All field variable components must have node-based interpolation.",err,error,*999)
              ENDIF
            ENDDO !componentIdx
            CALL FIELD_PARAMETER_SET_UPDATE_START(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
            CALL FIELD_PARAMETER_SET_UPDATE_FINISH(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
          ELSE
            CALL FlagError("Geometric field must be three dimensional.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Only rectangular cartesian coordinates are implemented.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Geometric field is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Polar mesh is not associated.",err,error,*999)
    ENDIF
    
    CALL EXITS("GeneratedMeshPolarGeometricParametersCalculate")
    RETURN
999 CALL ERRORS("GeneratedMeshPolarGeometricParametersCalculate",err,error)
    CALL EXITS("GeneratedMeshPolarGeometricParametersCalculate")
    RETURN 1

  END SUBROUTINE GeneratedMeshPolarGeometricParametersCalculate
=======
    INTEGER(INTG) :: component_idx,derivative_idx, &
      & component_node,MESH_COMPONENT, &
      & node_idx,node_position_idx(3), &
      & TOTAL_NUMBER_OF_NODES_XI(3),xi_idx,NODE_USER_NUMBER
    REAL(DP) :: DELTA_COORD(3,3),MY_ORIGIN(3),VALUE
    REAL(DP) :: DERIVATIVE_VALUES(MAXIMUM_GLOBAL_DERIV_NUMBER)
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    TYPE(FIELD_VARIABLE_COMPONENT_TYPE), POINTER :: FIELD_VARIABLE_COMPONENT
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    LOGICAL :: NODE_EXISTS,GHOST_NODE

    ENTERS("GeneratedMesh_RegularGeometricParametersCalculate",ERR,ERROR,*999)

    IF(ASSOCIATED(REGULAR_MESH)) THEN
      IF(ASSOCIATED(FIELD)) THEN
        NULLIFY(COORDINATE_SYSTEM)
        CALL FIELD_COORDINATE_SYSTEM_GET(FIELD,COORDINATE_SYSTEM,ERR,ERROR,*999)
        IF(COORDINATE_SYSTEM%TYPE==COORDINATE_RECTANGULAR_CARTESIAN_TYPE) THEN

          MY_ORIGIN=0.0_DP
          MY_ORIGIN(1:REGULAR_MESH%COORDINATE_DIMENSION)=REGULAR_MESH%ORIGIN(1:REGULAR_MESH%COORDINATE_DIMENSION)
          DELTA_COORD=0.0_DP
          TOTAL_NUMBER_OF_NODES_XI=1
          IF(FIELD%TYPE==FIELD_GEOMETRIC_TYPE) THEN
            FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR
            IF(ASSOCIATED(FIELD_VARIABLE)) THEN
              DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                FIELD_VARIABLE_COMPONENT=>FIELD_VARIABLE%COMPONENTS(component_idx)
                MESH_COMPONENT=FIELD_VARIABLE_COMPONENT%MESH_COMPONENT_NUMBER
                IF(FIELD_VARIABLE_COMPONENT%INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION) THEN
                  DO xi_idx=1,REGULAR_MESH%MESH_DIMENSION
                    TOTAL_NUMBER_OF_NODES_XI(xi_idx)=(REGULAR_MESH%BASES(MESH_COMPONENT)% &
                        & PTR%NUMBER_OF_NODES_XIC(xi_idx)-2)*REGULAR_MESH%NUMBER_OF_ELEMENTS_XI(xi_idx)+ &
                        & REGULAR_MESH%NUMBER_OF_ELEMENTS_XI(xi_idx)+1
                  ENDDO !xi_idx
                  DO xi_idx=1,REGULAR_MESH%MESH_DIMENSION
                    DELTA_COORD(1:REGULAR_MESH%COORDINATE_DIMENSION,xi_idx)= &
                      & REGULAR_MESH%BASE_VECTORS(1:REGULAR_MESH%COORDINATE_DIMENSION,xi_idx)/ &
                      & REAL(TOTAL_NUMBER_OF_NODES_XI(xi_idx)-1,DP)
                  ENDDO !xi_idx
                  SELECT CASE(FIELD%SCALINGS%SCALING_TYPE)
                  CASE(FIELD_NO_SCALING,FIELD_UNIT_SCALING)
                    DERIVATIVE_VALUES=0.0_DP
                    IF(REGULAR_MESH%NUMBER_OF_ELEMENTS_XI(1)>0) THEN
                      DERIVATIVE_VALUES(GLOBAL_DERIV_S1)= &
                        & REGULAR_MESH%BASE_VECTORS(component_idx,1)/REGULAR_MESH%NUMBER_OF_ELEMENTS_XI(1)
                    END IF
                    IF(REGULAR_MESH%MESH_DIMENSION>1) THEN
                      IF(REGULAR_MESH%NUMBER_OF_ELEMENTS_XI(2)>0) THEN
                        DERIVATIVE_VALUES(GLOBAL_DERIV_S2)= &
                          & REGULAR_MESH%BASE_VECTORS(component_idx,2)/REGULAR_MESH%NUMBER_OF_ELEMENTS_XI(2)
                      END IF
                    ENDIF
                    IF(REGULAR_MESH%MESH_DIMENSION>2) THEN
                      IF(REGULAR_MESH%NUMBER_OF_ELEMENTS_XI(3)>0) THEN
                        DERIVATIVE_VALUES(GLOBAL_DERIV_S3)= &
                          & REGULAR_MESH%BASE_VECTORS(component_idx,3)/REGULAR_MESH%NUMBER_OF_ELEMENTS_XI(3)
                      END IF
                    ENDIF
                  CASE DEFAULT
                    !Arc length or arithmetic mean scaling
                    DERIVATIVE_VALUES=0.0_DP
                    IF(REGULAR_MESH%NUMBER_OF_ELEMENTS_XI(1)>0) THEN
                      DERIVATIVE_VALUES(GLOBAL_DERIV_S1)=REGULAR_MESH%BASE_VECTORS(component_idx,1)/ &
                        & L2NORM(REGULAR_MESH%BASE_VECTORS(:,1))
                    END IF
                    IF(REGULAR_MESH%MESH_DIMENSION>1) THEN
                      IF(REGULAR_MESH%NUMBER_OF_ELEMENTS_XI(2)>0) THEN
                        DERIVATIVE_VALUES(GLOBAL_DERIV_S2)=REGULAR_MESH%BASE_VECTORS(component_idx,2)/ &
                          & L2NORM(REGULAR_MESH%BASE_VECTORS(:,2))
                      END IF
                    ENDIF
                    IF(REGULAR_MESH%MESH_DIMENSION>2) THEN
                      IF(REGULAR_MESH%NUMBER_OF_ELEMENTS_XI(3)>0) THEN
                        DERIVATIVE_VALUES(GLOBAL_DERIV_S3)=REGULAR_MESH%BASE_VECTORS(component_idx,3)/ &
                          & L2NORM(REGULAR_MESH%BASE_VECTORS(:,3))
                      END IF
                    ENDIF
                  END SELECT
                  !Update geometric parameters in this computational domain only
                  DOMAIN=>FIELD_VARIABLE_COMPONENT%DOMAIN
                  DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
                  DO component_node=1,TOTAL_NUMBER_OF_NODES_XI(1)*TOTAL_NUMBER_OF_NODES_XI(2)*TOTAL_NUMBER_OF_NODES_XI(3)
                    !Regular meshes with Lagrange/Hermite elements use different node numberings to other mesh types
                    IF(REGULAR_MESH%BASES(MESH_COMPONENT)%PTR%TYPE==BASIS_LAGRANGE_HERMITE_TP_TYPE) THEN
                      CALL GeneratedMesh_RegularComponentNodeToUserNumber(REGULAR_MESH%GENERATED_MESH,MESH_COMPONENT, &
                        & component_node,NODE_USER_NUMBER,ERR,ERROR,*999)
                    ELSE
                      NODE_USER_NUMBER=COMPONENT_NODE_TO_USER_NUMBER(REGULAR_MESH%GENERATED_MESH,MESH_COMPONENT, &
                        & component_node,ERR,ERROR)
                    END IF
                    CALL DOMAIN_TOPOLOGY_NODE_CHECK_EXISTS(FIELD_VARIABLE_COMPONENT%DOMAIN%TOPOLOGY, &
                      & NODE_USER_NUMBER,NODE_EXISTS,node_idx,GHOST_NODE,ERR,ERROR,*999)
                    IF(NODE_EXISTS.AND..NOT.GHOST_NODE) THEN
                      node_position_idx(3)=(component_node-1)/(TOTAL_NUMBER_OF_NODES_XI(2)*TOTAL_NUMBER_OF_NODES_XI(1))+1
                      node_position_idx(2)=MOD(component_node-1,TOTAL_NUMBER_OF_NODES_XI(2)*TOTAL_NUMBER_OF_NODES_XI(1))/ &
                        & TOTAL_NUMBER_OF_NODES_XI(1)+1
                      node_position_idx(1)=MOD(MOD(component_node-1,TOTAL_NUMBER_OF_NODES_XI(2)*TOTAL_NUMBER_OF_NODES_XI(1)), &
                        & TOTAL_NUMBER_OF_NODES_XI(1))+1
                      VALUE=0.0_DP
                      DO xi_idx=1,REGULAR_MESH%MESH_DIMENSION
                        VALUE=VALUE+REAL(node_position_idx(xi_idx)-1,DP)*DELTA_COORD(component_idx,xi_idx)
                      ENDDO !xi_idx
                      VALUE=MY_ORIGIN(component_idx)+VALUE
                      CALL FIELD_PARAMETER_SET_UPDATE_NODE(FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                        & 1,1,NODE_USER_NUMBER,component_idx,VALUE,ERR,ERROR,*999)
                      !Set derivatives
                      DO derivative_idx=2,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
                        CALL FIELD_PARAMETER_SET_UPDATE_NODE(FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                          & 1,derivative_idx,NODE_USER_NUMBER,component_idx,DERIVATIVE_VALUES(derivative_idx),ERR,ERROR,*999)
                      END DO !derivative_idx
                    ENDIF !node_exists
                  ENDDO !node_idx
                ELSE
                  LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))// &
                    & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
                    & " does not have node based interpolation."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ENDDO !component_idx
!!TODO: do boundary nodes first then start the update to overlap computation and computation.
              CALL FIELD_PARAMETER_SET_UPDATE_START(FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
              CALL FIELD_PARAMETER_SET_UPDATE_FINISH(FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
            ELSE
              LOCAL_ERROR="The standard field variable is not associated for field number "// &
                & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" is not a geometric field."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Non rectangular Cartesian coordinate systems are not implemented.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Field is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Regular mesh is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("GeneratedMesh_RegularGeometricParametersCalculate")
    RETURN
999 ERRORS("GeneratedMesh_RegularGeometricParametersCalculate",ERR,ERROR)
    EXITS("GeneratedMesh_RegularGeometricParametersCalculate")
    RETURN 1
    
  END SUBROUTINE GeneratedMesh_RegularGeometricParametersCalculate
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b

  !
  !================================================================================================================================
  !

<<<<<<< HEAD
  !>Updates the geometric field parameters from the initial nodal positions of the mesh. Derivatives are averaged via straight
  !>line approximation, except for circumferential component
  SUBROUTINE GeneratedMeshCylinderGeometricParametersCalculate(cylinderMesh,geometricField,err,error,*)
=======
  !>Updates the geometric field parameters from the initial nodal positions of the mesh. Derivatives are averaged via straight line approximation, except for circumferential component
  SUBROUTINE GeneratedMesh_CylinderGeometricParametersCalculate(CYLINDER_MESH,FIELD,ERR,ERROR,*)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b

    ! Argument variables
    TYPE(GeneratedMeshCylinderType), POINTER :: cylinderMesh !<A pointer to the cylinder mesh object
    TYPE(FIELD_TYPE), POINTER :: geometricField !<A pointer to the field to update the geometric parameters for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    ! Local variables
<<<<<<< HEAD
    TYPE(BASIS_TYPE), POINTER :: basis
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: coordinateSystem
    TYPE(DOMAIN_TYPE),POINTER :: domain
    TYPE(DOMAIN_NODES_TYPE), POINTER :: domainNodes
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable
    TYPE(FIELD_VARIABLE_COMPONENT_TYPE), POINTER :: fieldVariableComponent
    INTEGER(INTG) :: totalNumberOfNodesXi(3)
    INTEGER(INTG) :: componentIdx,xiIdx
    INTEGER(INTG) :: nodeIdx,globalNodeIdx,componentNodeIdx,dofIdx,derivativeIdx
    INTEGER(INTG) :: numberOfPlanarNodes,scalingType,meshComponent,coordinateType,coordinateDimension
    INTEGER(INTG) :: nodePositionIdx(3) ! holds r,theta,z indices
    REAL(DP) :: deltaCoordinate(3),deltaCoordinateXi(3),polarCoordinates(3),rectangularCoordinates(3),deriv
    TYPE(VARYING_STRING) :: localError

    NULLIFY(coordinateSystem)
    NULLIFY(fieldVariable)

    CALL Enters("GeneratedMeshCylinderGeometricParametersCalculate",err,error,*999)

    IF(ASSOCIATED(cylinderMesh)) THEN
      IF(ASSOCIATED(geometricField)) THEN
        CALL GeneratedMeshCoordinateSystemGet(cylinderMesh%generatedMesh,coordinateSystem,err,error,*999)
        CALL COORDINATE_SYSTEM_DIMENSION_GET(coordinateSystem,coordinateDimension,err,error,*999)
        CALL COORDINATE_SYSTEM_TYPE_GET(coordinateSystem,coordinateType,err,error,*999)
        IF(coordinateType==COORDINATE_RECTANGULAR_CARTESIAN_TYPE) THEN
          CALL FIELD_SCALING_TYPE_GET(geometricField,scalingType,err,error,*999)
          CALL FIELD_VARIABLE_TYPE_GET(geometricField,FIELD_U_VARIABLE_TYPE,fieldVariable,err,error,*999)
          IF(fieldVariable%NUMBER_OF_COMPONENTS==3) THEN
            DO componentIdx=1,3
              fieldVariableComponent=>fieldVariable%components(componentIdx)
              IF(fieldVariableComponent%INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION) THEN
                meshComponent=fieldVariableComponent%MESH_COMPONENT_NUMBER
                basis=>cylinderMesh%bases(meshComponent)%ptr
                !Calculate the total number of nodes in each xi direction
                DO xiIdx=1,3
                  totalNumberOfNodesXi(xiIdx)=(basis%NUMBER_OF_NODES_XIC(xiIdx)-1)*cylinderMesh%numberOfElementsXi(xiIdx)+1
                ENDDO !xiIdx
                totalNumberOfNodesXi(2)=totalNumberOfNodesXi(2)-1 !Theta loops around so slightly different
                numberOfPlanarNodes=totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                !Calculate delta
                deltaCoordinate(1)=(cylinderMesh%cylinderExtent(2)-cylinderMesh%cylinderExtent(1))/ &
                  & cylinderMesh%numberOfElementsXi(1)
                deltaCoordinate(2)=TWOPI/cylinderMesh%numberOfElementsXi(2)
                deltaCoordinate(3)=cylinderMesh%cylinderExtent(3)/cylinderMesh%numberOfElementsXi(3)
                DO xiIdx=1,3
                  deltaCoordinateXi(xiIdx)=deltaCoordinate(xiIdx)/(basis%NUMBER_OF_NODES_XIC(xiIdx)-1)
                ENDDO !xiIdx
                !Update geometric parameters in this computational domain only
                domain=>fieldVariable%components(meshComponent)%domain
                domainNodes=>domain%topology%nodes
                DO nodeIdx=1,domainNodes%NUMBER_OF_NODES
                  globalNodeIdx=domainNodes%nodes(nodeIdx)%GLOBAL_NUMBER
                  componentNodeIdx=GeneratedMeshUserNumberToComponentNode(cylinderMesh%generatedMesh,meshComponent, &
                    & globalNodeIdx,err,error)
                  !Calculate nodePositionIdx which will be used to calculate (r,theta,z) then (x,y,z)
                  nodePositionIdx(3)=(componentNodeIdx-1)/numberOfPlanarNodes
                  nodePositionIdx(2)=(componentNodeIdx-1-(nodePositionIdx(3))*numberOfPlanarNodes)/totalNumberOfNodesXi(1)
                  nodePositionIdx(1)=MOD(componentNodeIdx-1-(nodePositionIdx(3))*numberOfPlanarNodes,totalNumberOfNodesXi(1))
                  DO xiIdx=1,3
                    polarCoordinates(xiIdx)=nodePositionIdx(xiIdx)*deltaCoordinateXi(xiIdx)
                  ENDDO !xiIdx
                  polarCoordinates(1)=nodePositionIdx(1)*deltaCoordinateXi(1)+cylinderMesh%cylinderExtent(1) !Add the inner radius
                  rectangularCoordinates(1)=polarCoordinates(1)*COS(polarCoordinates(2))
                  rectangularCoordinates(2)=polarCoordinates(1)*SIN(polarCoordinates(2))
                  rectangularCoordinates(3)=polarCoordinates(3)
                  rectangularCoordinates=rectangularCoordinates+cylinderMesh%origin
                  !Default to version 1 of each node derivative
                  dofIdx=fieldVariableComponent%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(nodeIdx)%DERIVATIVES(1)%VERSIONS(1)
                  CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,dofIdx, &
                    & rectangularCoordinates(componentIdx),err,error,*999)
                  IF(domainNodes%NODES(nodeIdx)%NUMBER_OF_DERIVATIVES>1) THEN
                    DO derivativeIdx=2,fieldVariableComponent%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(nodeIdx)% &
                      & NUMBER_OF_DERIVATIVES
                      SELECT CASE(scalingType)
                      CASE(FIELD_NO_SCALING,FIELD_UNIT_SCALING)
                        SELECT CASE(domainNodes%nodes(nodeIdx)%derivatives(derivativeIdx)%GLOBAL_DERIVATIVE_INDEX)
                        CASE(GLOBAL_DERIV_S1)
                          SELECT CASE(componentIdx)
                          CASE(1)
                            DERIV=COS(polarCoordinates(2))*deltaCoordinate(1)
                          CASE(2)
                            DERIV=SIN(polarCoordinates(2))*deltaCoordinate(1)
                          CASE DEFAULT
                            DERIV=0.0_DP
                          END SELECT
                        CASE(GLOBAL_DERIV_S2)
                          SELECT CASE(componentIdx)
                          CASE(1)
                            DERIV=-polarCoordinates(1)*SIN(polarCoordinates(2))*deltaCoordinate(2)
                          CASE(2)
                            DERIV=polarCoordinates(1)*COS(polarCoordinates(2))*deltaCoordinate(2)
                          CASE DEFAULT
                            DERIV=0.0_DP
                          END SELECT
                        CASE(GLOBAL_DERIV_S3)
                          IF(componentIdx==3) THEN
                            DERIV=deltaCoordinate(3)
                          ELSE
                            DERIV=0.0_DP
                          ENDIF
                        CASE(GLOBAL_DERIV_S1_S2)
                          SELECT CASE(componentIdx)
                          CASE(1)
                            DERIV=-SIN(polarCoordinates(2))*deltaCoordinate(1)*deltaCoordinate(2)
                          CASE(2)
                            DERIV=COS(polarCoordinates(2))*deltaCoordinate(1)*deltaCoordinate(2)
                          CASE DEFAULT
                            DERIV=0.0_DP
                          END SELECT
                        CASE DEFAULT  ! all other non-xy-planar cross derivatives
                          DERIV=0.0_DP
                        END SELECT
                      CASE(FIELD_ARC_LENGTH_SCALING,FIELD_ARITHMETIC_MEAN_SCALING, &
                        & FIELD_GEOMETRIC_MEAN_SCALING,FIELD_HARMONIC_MEAN_SCALING)
                        SELECT CASE(domainNodes%NODES(nodeIdx)%DERIVATIVES(derivativeIdx)%GLOBAL_DERIVATIVE_INDEX)
                        CASE(GLOBAL_DERIV_S1)
                          SELECT CASE(componentIdx)
                          CASE(1)
                            DERIV=COS(polarCoordinates(2))
                          CASE(2)
                            DERIV=SIN(polarCoordinates(2))
                          CASE DEFAULT
                            DERIV=0.0_DP
                          END SELECT
                        CASE(GLOBAL_DERIV_S2)
                          SELECT CASE(componentIdx)
                          CASE(1)
                            DERIV=-SIN(polarCoordinates(2))
                          CASE(2)
                            DERIV=COS(polarCoordinates(2))
                          CASE DEFAULT
                            DERIV=0.0_DP
                          END SELECT
                        CASE(GLOBAL_DERIV_S3)
                          IF(componentIdx==3) THEN
                            DERIV=1.0_DP
                          ELSE
                            DERIV=0.0_DP
                          ENDIF
                        CASE(GLOBAL_DERIV_S1_S2)
                          SELECT CASE(componentIdx)
                          CASE(1)
                            DERIV=-SIN(polarCoordinates(2))*deltaCoordinate(1)*deltaCoordinate(2)
                          CASE(2)
                            DERIV=COS(polarCoordinates(2))*deltaCoordinate(1)*deltaCoordinate(2)
                          CASE DEFAULT
                            DERIV=0.0_DP
                          END SELECT
                        CASE DEFAULT  ! all other non-xy-planar cross derivatives
                          DERIV=0.0_DP
                        END SELECT
=======
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(DOMAIN_TYPE),POINTER :: DOMAIN
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(FIELD_VARIABLE_COMPONENT_TYPE), POINTER :: FIELD_VARIABLE_COMPONENT
    INTEGER(INTG) :: NUMBER_ELEMENTS_XI(3),NUMBER_OF_NODES_XIC(3)
    INTEGER(INTG) :: TOTAL_NUMBER_NODES_XI(3),INTERPOLATION_TYPES(3)
    INTEGER(INTG) :: component_idx,xi_idx
    INTEGER(INTG) :: np,global_np,component_np,ny,nk
    INTEGER(INTG) :: NUMBER_OF_PLANAR_NODES,SCALING_TYPE,MESH_COMPONENT
    INTEGER(INTG), ALLOCATABLE :: NIDX(:,:,:),EIDX(:,:,:)
    INTEGER(INTG) :: node_idx(3) ! holds r,theta,z indices
    REAL(DP) :: DELTA(3),DELTAi(3),POLAR_COORDS(3),RECT_COORDS(3)
    REAL(DP) :: CYLINDER_EXTENT(3),DERIV
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("GeneratedMesh_CylinderGeometricParametersCalculate",ERR,ERROR,*999)

    ! assign to the field
    IF(FIELD%TYPE==FIELD_GEOMETRIC_TYPE) THEN
      FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR
      IF(ASSOCIATED(FIELD_VARIABLE)) THEN
        IF(FIELD_VARIABLE%NUMBER_OF_COMPONENTS==3) THEN
          CALL FIELD_SCALING_TYPE_GET(FIELD,SCALING_TYPE,ERR,ERROR,*999)
          IF(SCALING_TYPE/=FIELD_UNIT_SCALING) &
            & CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"  Note: If the cylinder looks wonky, set field scaling to&
            & unit scaling type.",ERR,ERROR,*999)
          DO component_idx=1,3
            INTERPOLATION_TYPES(component_idx)=FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE
          ENDDO
          IF(ALL(INTERPOLATION_TYPES==FIELD_NODE_BASED_INTERPOLATION)) THEN
            DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
              FIELD_VARIABLE_COMPONENT=>FIELD_VARIABLE%COMPONENTS(component_idx)
              MESH_COMPONENT=FIELD_VARIABLE_COMPONENT%MESH_COMPONENT_NUMBER
              ! calculate the total number of nodes in each xi direction
              BASIS=>CYLINDER_MESH%BASES(MESH_COMPONENT)%PTR
              NUMBER_ELEMENTS_XI=CYLINDER_MESH%NUMBER_OF_ELEMENTS_XI
              NUMBER_OF_NODES_XIC=BASIS%NUMBER_OF_NODES_XIC
              DO xi_idx=1,3
                TOTAL_NUMBER_NODES_XI(xi_idx)=(NUMBER_OF_NODES_XIC(xi_idx)-1)*NUMBER_ELEMENTS_XI(xi_idx)+1
              ENDDO
              TOTAL_NUMBER_NODES_XI(2)=TOTAL_NUMBER_NODES_XI(2)-1 ! theta loops around so slightly different
              NUMBER_OF_PLANAR_NODES=TOTAL_NUMBER_NODES_XI(1)*TOTAL_NUMBER_NODES_XI(2)
              DOMAIN=>FIELD_VARIABLE%COMPONENTS(MESH_COMPONENT)%DOMAIN
              DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
              ! calculate DELTAi now
              CYLINDER_EXTENT=CYLINDER_MESH%CYLINDER_EXTENT
              DELTA(1)=(CYLINDER_EXTENT(2)-CYLINDER_EXTENT(1))/NUMBER_ELEMENTS_XI(1)
              DELTA(2)=TWOPI/NUMBER_ELEMENTS_XI(2)
              DELTA(3)=CYLINDER_EXTENT(3)/NUMBER_ELEMENTS_XI(3)
              DO xi_idx=1,3
                DELTAi(xi_idx)=DELTA(xi_idx)/(NUMBER_OF_NODES_XIC(xi_idx)-1)
              ENDDO
              DO np=1,DOMAIN_NODES%NUMBER_OF_NODES
                global_np=DOMAIN_NODES%NODES(np)%GLOBAL_NUMBER
                component_np=USER_NUMBER_TO_COMPONENT_NODE(CYLINDER_MESH%GENERATED_MESH, &
                    & MESH_COMPONENT,global_np,ERR,ERROR)
                ! calculate node_idx which will be used to calculate (r,theta,z) then (x,y,z)
                component_np=component_np-1 ! let's go 0-based index for a bit
                node_idx(3)=component_np/NUMBER_OF_PLANAR_NODES
                node_idx(2)=(component_np-(node_idx(3))*NUMBER_OF_PLANAR_NODES)/TOTAL_NUMBER_NODES_XI(1)
                node_idx(1)=MOD(component_np-(node_idx(3))*NUMBER_OF_PLANAR_NODES,TOTAL_NUMBER_NODES_XI(1))
                DO xi_idx=1,3
                  POLAR_COORDS(xi_idx)=node_idx(xi_idx)*DELTAi(xi_idx)
                ENDDO
                POLAR_COORDS(1)=node_idx(1)*DELTAi(1)+CYLINDER_EXTENT(1) ! add the inner radius
                RECT_COORDS(1)=POLAR_COORDS(1)*COS(POLAR_COORDS(2))
                RECT_COORDS(2)=POLAR_COORDS(1)*SIN(POLAR_COORDS(2))
                RECT_COORDS(3)=POLAR_COORDS(3)
                RECT_COORDS=RECT_COORDS+CYLINDER_MESH%ORIGIN
                !Default to version 1 of each node derivative
                ny=FIELD_VARIABLE_COMPONENT%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(np)%DERIVATIVES(1)%VERSIONS(1)
                CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,ny, &
                  & RECT_COORDS(component_idx),ERR,ERROR,*999)
                ! Do derivatives: if there are derivatives, we can assume it's cubic hermite
                !   given that quadratic hermites are only used for collapsed hex elements,
                !   but NB mixed bases have to be handled (e.g. CH-CH-linear combinations)
                IF(DOMAIN_NODES%NODES(np)%NUMBER_OF_DERIVATIVES>1) THEN
                  ! Since I decided how xi 1,2,3 line up with the cylinder polar coordinates,
                  ! we know a priori that only some of the derivatives are nonzero (analytically).
                  ! NOTE: if hermite type used, should assign FIELD_UNIT_SCALING type for this to work
                  DO nk=2,FIELD_VARIABLE_COMPONENT%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(np)%NUMBER_OF_DERIVATIVES
                    SELECT CASE(DOMAIN_NODES%NODES(np)%DERIVATIVES(nk)%GLOBAL_DERIVATIVE_INDEX)
                    CASE(GLOBAL_DERIV_S1)
                      SELECT CASE(component_idx)
                      CASE(1)
                        DERIV=COS(POLAR_COORDS(2))*DELTA(1)
                      CASE(2)
                        DERIV=SIN(POLAR_COORDS(2))*DELTA(1)
                      CASE DEFAULT
                        DERIV=0.0_DP
                      END SELECT
                    CASE(GLOBAL_DERIV_S2)
                      SELECT CASE(component_idx)
                      CASE(1)
                        DERIV=-POLAR_COORDS(1)*SIN(POLAR_COORDS(2))*DELTA(2)
                      CASE(2)
                        DERIV=POLAR_COORDS(1)*COS(POLAR_COORDS(2))*DELTA(2)
                      CASE DEFAULT
                        DERIV=0.0_DP
                      END SELECT
                    CASE(GLOBAL_DERIV_S3)
                      IF(component_idx==3) THEN
                        DERIV=DELTA(3)
                      ELSE
                        DERIV=0.0_DP
                      ENDIF
                    CASE(GLOBAL_DERIV_S1_S2)
                      SELECT CASE(component_idx)
                      CASE(1)
                        DERIV=-SIN(POLAR_COORDS(2))*DELTA(1)*DELTA(2)
                      CASE(2)
                        DERIV=COS(POLAR_COORDS(2))*DELTA(1)*DELTA(2)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
                      CASE DEFAULT
                        localError="The field scaling type of "//TRIM(NumberToVString(scalingType,"*",err,error))// &
                          & " is invalid for field number "//TRIM(NumberToVString(geometricField%USER_NUMBER,"*",err,error))//"."
                        CALL FlagError(localError,err,error,*999)
                      END SELECT
                      !Assign derivative. Default to version 1 of each node derivative
                      dofIdx=fieldVariableComponent%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(nodeIdx)%DERIVATIVES(derivativeIdx)% &
                        & VERSIONS(1)
                      CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                        & dofIdx,DERIV,err,error,*999)
                    ENDDO !derivativeIdx
                  ENDIF !derivatives
                ENDDO !nodeIdx
              ELSE
                CALL FlagError("All field variable components must have node-based interpolation.",err,error,*999)
              ENDIF
            ENDDO !componentIdx
            CALL FIELD_PARAMETER_SET_UPDATE_START(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
            CALL FIELD_PARAMETER_SET_UPDATE_FINISH(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
          ELSE
<<<<<<< HEAD
            CALL FlagError("Geometric field must be three dimensional.",err,error,*999)
=======
            CALL FlagError("All field variable components must have node-based interpolation.",ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
          ENDIF
        ELSE
<<<<<<< HEAD
          CALL FlagError("Only rectangular cartesian coordinates are implemented.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Geometric field is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Cylinder mesh is not associated.",err,error,*999)
    ENDIF
    
    CALL EXITS("GeneratedMeshCylinderGeometricParametersCalculate")
    RETURN
999 CALL ERRORS("GeneratedMeshCylinderGeometricParametersCalculate",err,error)
    CALL EXITS("GeneratedMeshCylinderGeometricParametersCalculate")
    RETURN 1

  END SUBROUTINE GeneratedMeshCylinderGeometricParametersCalculate
=======
          CALL FlagError("Geometric field must be three dimensional.",ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="The standard field variable is not associated for field number "// &
          & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" is not a geometric field."
      CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
    ENDIF

    ! all done
    IF(ALLOCATED(NIDX)) DEALLOCATE(NIDX)
    IF(ALLOCATED(EIDX)) DEALLOCATE(EIDX)

    EXITS("GeneratedMesh_CylinderGeometricParametersCalculate")
    RETURN
999 IF(ALLOCATED(NIDX)) DEALLOCATE(NIDX)
    IF(ALLOCATED(EIDX)) DEALLOCATE(EIDX)
    ERRORS("GeneratedMesh_CylinderGeometricParametersCalculate",ERR,ERROR)
    EXITS("GeneratedMesh_CylinderGeometricParametersCalculate")
    RETURN 1

  END SUBROUTINE GeneratedMesh_CylinderGeometricParametersCalculate
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b

  !
  !===========================================================================================================================
  !

  !>Calculates the geometric field parameters from the initial nodal positions of the mesh.
  !>Derivatives are averaged via straight line approximation, except for circumferential component
<<<<<<< HEAD
  SUBROUTINE GeneratedMeshEllipsoidGeometricParametersCalculate(ellipsoidMesh,geometricField,err,error,*)
=======
  SUBROUTINE GeneratedMesh_EllipsoidGeometricParametersCalculate(ELLIPSOID_MESH,FIELD,ERR,ERROR,*)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b

    ! Argument variables
    TYPE(GeneratedMeshEllipsoidType), POINTER :: ellipsoidMesh !<A pointer to the ellipsoid mesh object
    TYPE(FIELD_TYPE), POINTER :: geometricField !<A pointer to the geometric field to calculate the geometric parameters for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    ! Local variables
    TYPE(BASIS_TYPE), POINTER :: basis
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: coordinateSystem
    TYPE(DOMAIN_TYPE),POINTER :: domain
    TYPE(DECOMPOSITION_TYPE), POINTER :: decomposition
    TYPE(DOMAIN_NODES_TYPE), POINTER :: domainNodes
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable
    TYPE(FIELD_VARIABLE_COMPONENT_TYPE), POINTER :: fieldVariableComponent
    INTEGER(INTG) :: myComputationalNode,domainNumber,meshComponent,basisIdx,coordinateDimension,coordinateType
    INTEGER(INTG) :: numberOfElementsXi(3),numberOfNodesXic(3)
    INTEGER(INTG) :: totalNumberOfNodesXi(3)
    INTEGER(INTG) :: componentIdx,xiIdx
    INTEGER(INTG) :: nodeIdx,globalNodeNumber,nodeIdx1,nodeIdx2,nodeIdx3,localNodeNumber
    INTEGER(INTG) :: scalingType!,NUMBER_OF_PLANAR_NODES
    INTEGER(INTG), ALLOCATABLE :: nodeIndices(:,:,:),elementIndices(:,:,:)
    !INTEGER(INTG) :: node_idx(3) ! holds r,theta,z indices
<<<<<<< HEAD
    REAL(DP) :: deltaCoordinate(3),DeltaCoordinateXi(3),rectangularCoordinates(3),t,phi,alpha,xi,nu,x,y,z
    TYPE(VARYING_STRING) :: localError

    NULLIFY(basis,domain,decomposition,domainNodes,fieldVariable,fieldVariableComponent)

    CALL Enters("GeneratedMeshEllipsoidGeometricParametersCalculate",err,error,*999)

    myComputationalNode=COMPUTATIONAL_NODE_NUMBER_GET(err,error)

    IF(ASSOCIATED(ellipsoidMesh)) THEN
      IF(ASSOCIATED(geometricField)) THEN
        IF(geometricField%TYPE==FIELD_GEOMETRIC_TYPE) THEN
          CALL GeneratedMeshCoordinateSystemGet(ellipsoidMesh%generatedMesh,coordinateSystem,err,error,*999)
          CALL COORDINATE_SYSTEM_DIMENSION_GET(coordinateSystem,coordinateDimension,err,error,*999)
          CALL COORDINATE_SYSTEM_TYPE_GET(coordinateSystem,coordinateType,err,error,*999)
          IF(coordinateType==COORDINATE_RECTANGULAR_CARTESIAN_TYPE) THEN
            CALL FIELD_SCALING_TYPE_GET(geometricField,scalingType,err,error,*999)
            CALL FIELD_VARIABLE_TYPE_GET(geometricField,FIELD_U_VARIABLE_TYPE,fieldVariable,err,error,*999)
            ! assign to the field
            nodeIdx=0
            
            IF(fieldVariable%NUMBER_OF_COMPONENTS==3) THEN
              DO componentIdx=1,3
                fieldVariableComponent=>fieldVariable%components(componentIdx)
                IF(fieldVariableComponent%INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION) THEN
                  meshComponent=fieldVariableComponent%MESH_COMPONENT_NUMBER
                  basisIdx=meshComponent*2-1
                  !< Ellipsoid_extent= inner long axis, inner short axis, wall thickness, top angle (from 0)
                  ! calculate the total number of nodes in each xi direction
                  !Check that the all geometric bases use the same mesh component
                  basis=>ellipsoidMesh%bases(basisIdx)%ptr
                  numberOfElementsXi=ellipsoidMesh%numberOfElementsXi
                  numberOfNodesXic=basis%NUMBER_OF_NODES_XIC
                  DO xiIdx=1,3
                    totalNumberOfNodesXi(xiIdx)=(numberOfNodesXic(xiIdx)-1)*numberOfElementsXi(xiIdx)+1
                  ENDDO !xiIdx
                  totalNumberOfNodesXi(1)=totalNumberOfNodesXi(1)-1 ! theta loops around so slightly different
                  ! calculate DeltaCoordinateXi now
                  deltaCoordinate(1)=TWOPI/numberOfElementsXi(1)
                  deltaCoordinate(2)=(PI-ellipsoidMesh%ellipsoidExtent(4))/numberOfElementsXi(2)
                  deltaCoordinate(3)=ellipsoidMesh%ellipsoidExtent(3)/numberOfElementsXi(3)
                  DO xiIdx=1,3
                    DeltaCoordinateXi(xiIdx)=deltaCoordinate(xiIdx)/(numberOfNodesXic(xiIdx)-1)
                  ENDDO !xiIdx
                  ! NUMBER_OF_PLANAR_NODES=totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                  domain=>fieldVariable%components(componentIdx)%domain
                  domainNodes=>domain%topology%nodes
                  decomposition=>geometricField%decomposition !\todo: test all these pointers
                  IF(ellipsoidMesh%ellipsoidExtent(1)>ellipsoidMesh%ellipsoidExtent(2)) THEN
                    !Prolate spheroid
                    nodeIdx3=1
                    !Inner surface
                    alpha=SQRT((ellipsoidMesh%ellipsoidExtent(1))**2-(ellipsoidMesh%ellipsoidExtent(2))**2)
                    !xi=log(ellipsoidMesh%ellipsoidExtent(1)/alpha+sqrt((ellipsoidMesh%ellipsoidExtent(1)/alpha)**2+1))
                    xi=acosh(ellipsoidMesh%ellipsoidExtent(1)/alpha)
                    
                    nodeIdx2=1
                    !Apex node
                    nodeIdx=1
                    globalNodeNumber=GeneratedMeshComponentNodeToUserNumber(ellipsoidMesh%generatedMesh,basisIdx,nodeIdx,err,error)
                    CALL DECOMPOSITION_NODE_DOMAIN_GET(decomposition,globalNodeNumber,meshComponent,domainNumber,err,error,*999)
                    IF(domainNumber==myComputationalNode) THEN
                      rectangularCoordinates(1)=0.0_DP
                      rectangularCoordinates(2)=0.0_DP
                      rectangularCoordinates(3)=-ellipsoidMesh%ellipsoidExtent(1)
                      !Default to version 1 of each node derivative
                      CALL FIELD_PARAMETER_SET_UPDATE_NODE(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1, &
                        & globalNodeNumber,componentIdx,rectangularCoordinates(componentIdx),err,error,*999)
                      localNodeNumber=domain%mappings%nodes%GLOBAL_TO_LOCAL_MAP(globalNodeNumber)%local_number(1)
                      IF(domainNodes%nodes(localNodeNumber)%NUMBER_OF_DERIVATIVES>1) THEN
                        CALL FlagError("Not generalized to hermittean elements.",err,error,*999)
                      ENDIF !derivatives
                    ENDIF
                    
                    DO nodeIdx2=2,totalNumberOfNodesXi(2)
                      !longitudinal loop
                      nu=PI-DeltaCoordinateXi(2)*(nodeIdx2-1)
                      DO nodeIdx1=1,totalNumberOfNodesXi(1)
                        !circumferential loop
                        phi=DeltaCoordinateXi(1)*(nodeIdx1-1)
                        rectangularCoordinates(1)=alpha*(SINH(xi)*SIN(nu)*COS(phi))
                        rectangularCoordinates(2)=alpha*(SINH(xi)*SIN(nu)*SIN(phi))
                        rectangularCoordinates(3)=alpha*(COSH(xi)*COS(nu))
                        nodeIdx=nodeIdx+1
                        globalNodeNumber=GeneratedMeshComponentNodeToUserNumber(ellipsoidMesh%generatedMesh,basisIdx,nodeIdx, &
                          & err,error)
                        CALL DECOMPOSITION_NODE_DOMAIN_GET(decomposition,globalNodeNumber,meshComponent,domainNumber,err,error,*999)
                        IF(domainNumber==myComputationalNode) THEN
                          CALL FIELD_PARAMETER_SET_UPDATE_NODE(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1, &
                            & globalNodeNumber,componentIdx,rectangularCoordinates(componentIdx),err,error,*999)
                          localNodeNumber=domain%mappings%nodes%GLOBAL_TO_LOCAL_MAP(globalNodeNumber)%local_number(1)
                          IF(domainNodes%nodes(localNodeNumber)%NUMBER_OF_DERIVATIVES>1) THEN
                            CALL FlagError("Not generalized to hermittean elements.",err,error,*999)
                          ENDIF !derivatives
                        ENDIF
                      ENDDO !nodeIdx1
                    ENDDO !nodeIdx2
                    
                    DO nodeIdx3=2,totalNumberOfNodesXi(3)
=======
    REAL(DP) :: DELTA(3),DELTAi(3),RECT_COORDS(3),t,phi,alpha,xi,nu,x,y,z
    REAL(DP) :: ELLIPSOID_EXTENT(4)
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    NULLIFY(BASIS,DOMAIN,DECOMPOSITION,DOMAIN_NODES,FIELD_VARIABLE,FIELD_VARIABLE_COMPONENT)

    ENTERS("GeneratedMesh_EllipsoidGeometricParametersCalculate",ERR,ERROR,*999)

    MY_COMPUTATIONAL_NODE=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)

    ! assign to the field
    np=0
    IF(FIELD%TYPE==FIELD_GEOMETRIC_TYPE) THEN
       FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR
       IF(ASSOCIATED(FIELD_VARIABLE)) THEN
          IF(FIELD_VARIABLE%NUMBER_OF_COMPONENTS==3) THEN
             MESH_COMPONENT=FIELD_VARIABLE%COMPONENTS(1)%MESH_COMPONENT_NUMBER
             DO component_idx=2,3
                IF(FIELD_VARIABLE%COMPONENTS(component_idx)%MESH_COMPONENT_NUMBER/=MESH_COMPONENT) THEN
                   CALL FlagError("Multiple mesh components for geometric components is not implemented.",ERR,ERROR,*999)
                ENDIF
             ENDDO
             basis_idx=MESH_COMPONENT*2-1
             !< Ellipsoid_extent= inner long axis, inner short axis, wall thickness, top angle (from 0)
             ! calculate the total number of nodes in each xi direction
             IF(ALLOCATED(ELLIPSOID_MESH%BASES)) THEN
                !Check that the all geometric bases use the same mesh component
                BASIS=>ELLIPSOID_MESH%BASES(basis_idx)%PTR
                NUMBER_ELEMENTS_XI=ELLIPSOID_MESH%NUMBER_OF_ELEMENTS_XI
                NUMBER_OF_NODES_XIC=BASIS%NUMBER_OF_NODES_XIC
                DO xi_idx=1,3
                   TOTAL_NUMBER_NODES_XI(xi_idx)=(NUMBER_OF_NODES_XIC(xi_idx)-1)*NUMBER_ELEMENTS_XI(xi_idx)+1
                ENDDO
                TOTAL_NUMBER_NODES_XI(1)=TOTAL_NUMBER_NODES_XI(1)-1 ! theta loops around so slightly different
                ! calculate DELTAi now
                ELLIPSOID_EXTENT=ELLIPSOID_MESH%ELLIPSOID_EXTENT
                DELTA(1)=TWOPI/NUMBER_ELEMENTS_XI(1)
                DELTA(2)=(PI-ELLIPSOID_EXTENT(4))/NUMBER_ELEMENTS_XI(2)
                DELTA(3)=ELLIPSOID_EXTENT(3)/NUMBER_ELEMENTS_XI(3)
                DO xi_idx=1,3
                   DELTAi(xi_idx)=DELTA(xi_idx)/(NUMBER_OF_NODES_XIC(xi_idx)-1)
                ENDDO
             ELSE
                CALL FlagError("Ellipsoid mesh does not have bases allocated.",ERR,ERROR,*999)
             ENDIF
             CALL FIELD_SCALING_TYPE_GET(FIELD,SCALING_TYPE,ERR,ERROR,*999)
             IF(SCALING_TYPE/=FIELD_UNIT_SCALING) &
                  & CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"  Note: If the ellipsoid looks wonky, set field scaling to &
                  & unit scaling type.",ERR,ERROR,*999)
             ! NUMBER_OF_PLANAR_NODES=TOTAL_NUMBER_NODES_XI(1)*TOTAL_NUMBER_NODES_XI(2)
             DO component_idx=1,3
                INTERPOLATION_TYPES(component_idx)=FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE
             ENDDO
             IF(ALL(INTERPOLATION_TYPES==FIELD_NODE_BASED_INTERPOLATION)) THEN
                DOMAIN=>FIELD_VARIABLE%COMPONENTS(1)%DOMAIN ! just grab the first one
                DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
                !DECOMPOSITION=>DOMAIN%DECOMPOSITION !\todo: test all these pointers
                DECOMPOSITION=>FIELD%DECOMPOSITION !\todo: test all these pointers
                IF (ELLIPSOID_EXTENT(1)>ELLIPSOID_EXTENT(2)) THEN
                   !Prolate spheroid
                   k=1
                   !inner surface
                   alpha=sqrt((ELLIPSOID_EXTENT(1))**2-(ELLIPSOID_EXTENT(2))**2)
                   !xi=log(ELLIPSOID_EXTENT(1)/alpha+sqrt((ELLIPSOID_EXTENT(1)/alpha)**2+1))
                   xi=acosh(ELLIPSOID_EXTENT(1)/alpha)

                   j=1
                   !apex node
                   np=1
                   npg=COMPONENT_NODE_TO_USER_NUMBER(ELLIPSOID_MESH%GENERATED_MESH,basis_idx,np,ERR,ERROR)
                   CALL DECOMPOSITION_NODE_DOMAIN_GET(DECOMPOSITION,npg,MESH_COMPONENT,DOMAIN_NUMBER,ERR,ERROR,*999)
                   IF(DOMAIN_NUMBER==MY_COMPUTATIONAL_NODE) THEN
                      RECT_COORDS(1)=0
                      RECT_COORDS(2)=0
                      RECT_COORDS(3)=-ELLIPSOID_EXTENT(1)
                      DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                         !Default to version 1 of each node derivative
                         CALL FIELD_PARAMETER_SET_UPDATE_NODE(FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,npg, &
                              & component_idx,RECT_COORDS(component_idx),ERR,ERROR,*999)
                         local_node=DOMAIN%MAPPINGS%NODES%GLOBAL_TO_LOCAL_MAP(npg)%local_number(1)
                         IF(DOMAIN_NODES%NODES(local_node)%NUMBER_OF_DERIVATIVES>1) THEN
                            CALL FlagError("Not generalized to hermittean elements.",ERR,ERROR,*999)
                         ENDIF !derivatives
                      ENDDO
                   ENDIF

                   DO j=2,TOTAL_NUMBER_NODES_XI(2)
                      !longitudinal loop
                      nu=PI-DELTAi(2)*(j-1)
                      DO i=1,TOTAL_NUMBER_NODES_XI(1)
                         !circumferential loop
                         phi=DELTAi(1)*(i-1)
                         RECT_COORDS(1)=alpha*(sinh(xi)*sin(nu)*cos(phi))
                         RECT_COORDS(2)=alpha*(sinh(xi)*sin(nu)*sin(phi))
                         RECT_COORDS(3)=alpha*(cosh(xi)*cos(nu))
                         np=np+1
                         npg=COMPONENT_NODE_TO_USER_NUMBER(ELLIPSOID_MESH%GENERATED_MESH,basis_idx,np,ERR,ERROR)
                         CALL DECOMPOSITION_NODE_DOMAIN_GET(DECOMPOSITION,npg,MESH_COMPONENT,DOMAIN_NUMBER,ERR,ERROR,*999)
                         IF(DOMAIN_NUMBER==MY_COMPUTATIONAL_NODE) THEN
                            DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                               CALL FIELD_PARAMETER_SET_UPDATE_NODE(FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,npg, &
                                    & component_idx,RECT_COORDS(component_idx),ERR,ERROR,*999)
                               local_node=DOMAIN%MAPPINGS%NODES%GLOBAL_TO_LOCAL_MAP(npg)%local_number(1)
                               IF(DOMAIN_NODES%NODES(local_node)%NUMBER_OF_DERIVATIVES>1) THEN
                                  CALL FlagError("Not generalized to hermittean elements.",ERR,ERROR,*999)
                               ENDIF !derivatives
                            ENDDO
                         ENDIF
                      ENDDO
                   ENDDO

                   DO k=2,TOTAL_NUMBER_NODES_XI(3)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
                      !transmural loop
                      nodeIdx2=1
                      !apex nodes
<<<<<<< HEAD
                      rectangularCoordinates(1)=0
                      rectangularCoordinates(2)=0
                      rectangularCoordinates(3)=-ellipsoidMesh%ellipsoidExtent(1)-(nodeIdx3-1)*(DeltaCoordinateXi(3))
                      nodeIdx=nodeIdx+1
                      globalNodeNumber=GeneratedMeshComponentNodeToUserNumber(ellipsoidMesh%generatedMesh,basisIdx,nodeIdx, &
                        & err,error)
                      CALL DECOMPOSITION_NODE_DOMAIN_GET(decomposition,globalNodeNumber,meshComponent,domainNumber,err,error,*999)
                      IF(domainNumber==myComputationalNode) THEN
                        CALL FIELD_PARAMETER_SET_UPDATE_NODE(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1, &
                          & globalNodeNumber,componentIdx,rectangularCoordinates(componentIdx),err,error,*999)
                        localNodeNumber=domain%mappings%nodes%GLOBAL_TO_LOCAL_MAP(globalNodeNumber)%local_number(1)
                        IF(domainNodes%nodes(localNodeNumber)%NUMBER_OF_DERIVATIVES>1) THEN
                          CALL FlagError("Not generalized to hermittean elements.",err,error,*999)
                        ENDIF !derivatives
                      ENDIF
                      
                      DO nodeIdx2=2,totalNumberOfNodesXi(2)
                        !longitudinal loop
                        nu=PI-DeltaCoordinateXi(2)*(nodeIdx2-1)
                        DO nodeIdx1=1,totalNumberOfNodesXi(1)
                          !circumferential loop
                          phi=DeltaCoordinateXi(1)*(nodeIdx1-1)
                          x=alpha*(SINH(xi)*SIN(nu)*COS(phi))
                          y=alpha*(SINH(xi)*SIN(nu)*SIN(phi))
                          z=alpha*(COSH(xi)*COS(nu))
                          !Normal vector from inner surface with length DeltaCoordinateXi(3)(nodeIdx3-1)
                          ! Finney&Thomas: Calculus, second edition, Addison-Wesley Publishing Company, 1994, page 847
                          !X=x(1+2t/a^2) Y=y(1+2t/a^2) Z=z(1+2t/c^2)
                          t=(DeltaCoordinateXi(3)*(nodeIdx3-1))/SQRT((4*x**2/(ellipsoidMesh%ellipsoidExtent(2))**4)+ &
                            & (4*y**2/(ellipsoidMesh%ellipsoidExtent(2))**4)+(4*z**2/(ellipsoidMesh%ellipsoidExtent(1))**4))
                          rectangularCoordinates(1)=x*(1+2*t/(ellipsoidMesh%ellipsoidExtent(2))**2)
                          rectangularCoordinates(2)=y*(1+2*t/(ellipsoidMesh%ellipsoidExtent(2))**2)
                          rectangularCoordinates(3)=z*(1+2*t/(ellipsoidMesh%ellipsoidExtent(1))**2)
                          nodeIdx=nodeIdx+1
                          globalNodeNumber=GeneratedMeshComponentNodeToUserNumber(ellipsoidMesh%generatedMesh,basisIdx,nodeIdx, &
                            & err,error)
                          CALL DECOMPOSITION_NODE_DOMAIN_GET(decomposition,globalNodeNumber,meshComponent,domainNumber, &
                            & err,error,*999)
                          IF(domainNumber==myComputationalNode) THEN
                            CALL FIELD_PARAMETER_SET_UPDATE_NODE(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1, &
                              & globalNodeNumber,componentIdx,rectangularCoordinates(componentIdx),err,error,*999)
                            localNodeNumber=domain%mappings%nodes%GLOBAL_TO_LOCAL_MAP(globalNodeNumber)%local_number(1)
                            IF(domainNodes%nodes(localNodeNumber)%NUMBER_OF_DERIVATIVES>1) THEN
                              CALL FlagError("Not generalized to hermittean elements.",err,error,*999)
                            ENDIF !derivatives
                          ENDIF
                        ENDDO !nodeIdx1
                      ENDDO !nodeIdx2
                    ENDDO !nodeIdx3
                  ELSE IF(ABS(ellipsoidMesh%ellipsoidExtent(1)-ellipsoidMesh%ellipsoidExtent(2))<ZERO_TOLERANCE) THEN
                    !Sphere
                    nodeIdx=0
                    DO nodeIdx3=1,totalNumberOfNodesXi(3)
=======
                      RECT_COORDS(1)=0
                      RECT_COORDS(2)=0
                      RECT_COORDS(3)=-ELLIPSOID_EXTENT(1)-(k-1)*(DELTAi(3))
                      np=np+1
                      npg=COMPONENT_NODE_TO_USER_NUMBER(ELLIPSOID_MESH%GENERATED_MESH,basis_idx,np,ERR,ERROR)
                      CALL DECOMPOSITION_NODE_DOMAIN_GET(DECOMPOSITION,npg,MESH_COMPONENT,DOMAIN_NUMBER,ERR,ERROR,*999)
                      IF(DOMAIN_NUMBER==MY_COMPUTATIONAL_NODE) THEN
                         DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                            CALL FIELD_PARAMETER_SET_UPDATE_NODE(FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,npg, &
                                 & component_idx,RECT_COORDS(component_idx),ERR,ERROR,*999)
                            local_node=DOMAIN%MAPPINGS%NODES%GLOBAL_TO_LOCAL_MAP(npg)%local_number(1)
                            IF(DOMAIN_NODES%NODES(local_node)%NUMBER_OF_DERIVATIVES>1) THEN
                               CALL FlagError("Not generalized to hermittean elements.",ERR,ERROR,*999)
                            ENDIF !derivatives
                         ENDDO
                      ENDIF

                      DO j=2,TOTAL_NUMBER_NODES_XI(2)
                         !longitudinal loop
                         nu=PI-DELTAi(2)*(j-1)
                         DO i=1,TOTAL_NUMBER_NODES_XI(1)
                            !circumferential loop
                            phi=DELTAi(1)*(i-1)
                            x=alpha*(sinh(xi)*sin(nu)*cos(phi))
                            y=alpha*(sinh(xi)*sin(nu)*sin(phi))
                            z=alpha*(cosh(xi)*cos(nu))
                            !Normal vector from inner surface with length DELTAi(3)(k-1)
                            ! Finney&Thomas: Calculus, second edition, Addison-Wesley Publishing Company, 1994, page 847
                            !X=x(1+2t/a^2) Y=y(1+2t/a^2) Z=z(1+2t/c^2)
                            t=(DELTAi(3)*(k-1))/sqrt((4*x**2/(ELLIPSOID_EXTENT(2))**4)+ &
                                 & (4*y**2/(ELLIPSOID_EXTENT(2))**4)+(4*z**2/(ELLIPSOID_EXTENT(1))**4))
                            RECT_COORDS(1)=x*(1+2*t/(ELLIPSOID_EXTENT(2))**2)
                            RECT_COORDS(2)=y*(1+2*t/(ELLIPSOID_EXTENT(2))**2)
                            RECT_COORDS(3)=z*(1+2*t/(ELLIPSOID_EXTENT(1))**2)
                            np=np+1
                            npg=COMPONENT_NODE_TO_USER_NUMBER(ELLIPSOID_MESH%GENERATED_MESH,basis_idx,np,ERR,ERROR)
                            CALL DECOMPOSITION_NODE_DOMAIN_GET(DECOMPOSITION,npg,MESH_COMPONENT,DOMAIN_NUMBER,ERR,ERROR,*999)
                            IF(DOMAIN_NUMBER==MY_COMPUTATIONAL_NODE) THEN
                               DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                                  CALL FIELD_PARAMETER_SET_UPDATE_NODE(FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,npg, &
                                       & component_idx,RECT_COORDS(component_idx),ERR,ERROR,*999)
                                  local_node=DOMAIN%MAPPINGS%NODES%GLOBAL_TO_LOCAL_MAP(npg)%local_number(1)
                                  IF(DOMAIN_NODES%NODES(local_node)%NUMBER_OF_DERIVATIVES>1) THEN
                                     CALL FlagError("Not generalized to hermittean elements.",ERR,ERROR,*999)
                                  ENDIF !derivatives
                               ENDDO
                            ENDIF
                         ENDDO
                      ENDDO
                   ENDDO
                ELSEIF (ABS(ELLIPSOID_EXTENT(1)-ELLIPSOID_EXTENT(2))<ZERO_TOLERANCE) THEN
                   !Sphere
                   np=0
                   DO k=1,TOTAL_NUMBER_NODES_XI(3)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
                      !transmural loop
                      alpha=ellipsoidMesh%ellipsoidExtent(1)+(nodeIdx3-1)*(DeltaCoordinateXi(3))
                      nodeIdx2=1
                      !apex nodes
<<<<<<< HEAD
                      rectangularCoordinates(1)=0
                      rectangularCoordinates(2)=0
                      rectangularCoordinates(3)=-alpha
                      nodeIdx=nodeIdx+1
                      globalNodeNumber=GeneratedMeshComponentNodeToUserNumber(ellipsoidMesh%generatedMesh,basisIdx,nodeIdx, &
                        & err,error)
                      CALL DECOMPOSITION_NODE_DOMAIN_GET(decomposition,globalNodeNumber,meshComponent,domainNumber,err,error,*999)
                      IF(domainNumber==myComputationalNode) THEN
                        CALL FIELD_PARAMETER_SET_UPDATE_NODE(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1, &
                          & globalNodeNumber,componentIdx,rectangularCoordinates(componentIdx),err,error,*999)
                        localNodeNumber=domain%mappings%nodes%GLOBAL_TO_LOCAL_MAP(globalNodeNumber)%local_number(1)
                        IF(domainNodes%nodes(localNodeNumber)%NUMBER_OF_DERIVATIVES>1) THEN
                          CALL FlagError("Not generalized to hermittean elements.",err,error,*999)
                        ENDIF !derivatives
                      ENDIF
                      
                      DO nodeIdx2=2,totalNumberOfNodesXi(2)
                        !longitudinal loop
                        nu=PI-DeltaCoordinateXi(2)*(nodeIdx2-1)
                        DO nodeIdx1=1,totalNumberOfNodesXi(1)
                          !circumferential loop
                          phi=DeltaCoordinateXi(1)*(nodeIdx1-1)
                          rectangularCoordinates(1)=alpha*SIN(nu)*COS(phi)
                          rectangularCoordinates(2)=alpha*SIN(nu)*SIN(phi)
                          rectangularCoordinates(3)=alpha*COS(nu)
                          nodeIdx=nodeIdx+1
                          globalNodeNumber=GeneratedMeshComponentNodeToUserNumber(ellipsoidMesh%generatedMesh,basisIdx,nodeIdx, &
                            & err,error)
                          CALL DECOMPOSITION_NODE_DOMAIN_GET(decomposition,globalNodeNumber,meshComponent,domainNumber, &
                            & err,error,*999)
                          IF(domainNumber==myComputationalNode) THEN
                            CALL FIELD_PARAMETER_SET_UPDATE_NODE(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1, &
                              & globalNodeNumber,componentIdx,rectangularCoordinates(componentIdx),err,error,*999)
                            localNodeNumber=domain%mappings%nodes%GLOBAL_TO_LOCAL_MAP(globalNodeNumber)%local_number(1)
                            IF(domainNodes%nodes(localNodeNumber)%NUMBER_OF_DERIVATIVES>1) THEN
                              CALL FlagError("Not generalized to hermittean elements.",err,error,*999)
                            ENDIF !derivatives
                          ENDIF
                        ENDDO !nodeIdx1
                      ENDDO !nodeIdx2
                    ENDDO !nodeIdx3
                    
                  ELSE IF(ellipsoidMesh%ellipsoidExtent(1)<ellipsoidMesh%ellipsoidExtent(2)) THEN
                    !Oblate spheroid
                    nodeIdx3=1
                    !inner surface
                    alpha=SQRT((ellipsoidMesh%ellipsoidExtent(2))**2-(ellipsoidMesh%ellipsoidExtent(1))**2)
                    !xi=log(ellipsoidMesh%ellipsoidExtent(1)/alpha+sqrt((ellipsoidMesh%ellipsoidExtent(1)/alpha)**2+1))
                    xi=acosh(ellipsoidMesh%ellipsoidExtent(2)/alpha)
                    
                    nodeIdx2=1
                    !apex node
                    nodeIdx=1
                    globalNodeNumber=GeneratedMeshComponentNodeToUserNumber(ellipsoidMesh%generatedMesh,basisIdx,nodeIdx,err,error)
                    CALL DECOMPOSITION_NODE_DOMAIN_GET(decomposition,globalNodeNumber,meshComponent,domainNumber,err,error,*999)
                    IF(domainNumber==myComputationalNode) THEN
                      rectangularCoordinates(1)=0
                      rectangularCoordinates(2)=0
                      rectangularCoordinates(3)=-ellipsoidMesh%ellipsoidExtent(1)
                      CALL FIELD_PARAMETER_SET_UPDATE_NODE(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1, &
                        & globalNodeNumber,componentIdx,rectangularCoordinates(componentIdx),err,error,*999)
                      localNodeNumber=domain%mappings%nodes%GLOBAL_TO_LOCAL_MAP(globalNodeNumber)%local_number(1)
                      IF(domainNodes%nodes(localNodeNumber)%NUMBER_OF_DERIVATIVES>1) THEN
                        CALL FlagError("Not generalized to hermittean elements.",err,error,*999)
                      ENDIF !derivatives
                    ENDIF
                    
                    DO nodeIdx2=2,totalNumberOfNodesXi(2)
                      !longitudinal loop
                      nu=-PI/2+DeltaCoordinateXi(2)*(nodeIdx2-1)
                      DO nodeIdx1=1,totalNumberOfNodesXi(1)
                        !circumferential loop
                        phi=DeltaCoordinateXi(1)*(nodeIdx1-1)
                        rectangularCoordinates(1)=alpha*(COSH(xi)*COS(nu)*COS(phi))
                        rectangularCoordinates(2)=alpha*(COSH(xi)*COS(nu)*SIN(phi))
                        rectangularCoordinates(3)=alpha*(SINH(xi)*SIN(nu))
                        nodeIdx=nodeIdx+1
                        globalNodeNumber=GeneratedMeshComponentNodeToUserNumber(ellipsoidMesh%generatedMesh,basisIdx,nodeIdx, &
                          & err,error)
                        CALL DECOMPOSITION_NODE_DOMAIN_GET(decomposition,globalNodeNumber,meshComponent,domainNumber,err,error,*999)
                        IF(domainNumber==myComputationalNode) THEN
                          CALL FIELD_PARAMETER_SET_UPDATE_NODE(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1, &
                            & globalNodeNumber,componentIdx,rectangularCoordinates(componentIdx),err,error,*999)
                          localNodeNumber=domain%mappings%nodes%GLOBAL_TO_LOCAL_MAP(globalNodeNumber)%local_number(1)
                          IF(domainNodes%nodes(localNodeNumber)%NUMBER_OF_DERIVATIVES>1) THEN
                            CALL FlagError("Not generalized to hermittean elements.",err,error,*999)
                          ENDIF !derivatives
                        ENDIF
                      ENDDO !nodeIdx1
                    ENDDO !nodeIdx2
                    
                    DO nodeIdx3=2,totalNumberOfNodesXi(3)
=======
                      RECT_COORDS(1)=0
                      RECT_COORDS(2)=0
                      RECT_COORDS(3)=-alpha
                      np=np+1
                      npg=COMPONENT_NODE_TO_USER_NUMBER(ELLIPSOID_MESH%GENERATED_MESH,basis_idx,np,ERR,ERROR)
                      CALL DECOMPOSITION_NODE_DOMAIN_GET(DECOMPOSITION,npg,MESH_COMPONENT,DOMAIN_NUMBER,ERR,ERROR,*999)
                      IF(DOMAIN_NUMBER==MY_COMPUTATIONAL_NODE) THEN
                         DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                            CALL FIELD_PARAMETER_SET_UPDATE_NODE(FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,npg, &
                                 & component_idx,RECT_COORDS(component_idx),ERR,ERROR,*999)
                            local_node=DOMAIN%MAPPINGS%NODES%GLOBAL_TO_LOCAL_MAP(npg)%local_number(1)
                            IF(DOMAIN_NODES%NODES(local_node)%NUMBER_OF_DERIVATIVES>1) THEN
                               CALL FlagError("Not generalized to hermittean elements.",ERR,ERROR,*999)
                            ENDIF !derivatives
                         ENDDO
                      ENDIF

                      DO j=2,TOTAL_NUMBER_NODES_XI(2)
                         !longitudinal loop
                         nu=PI-DELTAi(2)*(j-1)
                         DO i=1,TOTAL_NUMBER_NODES_XI(1)
                            !circumferential loop
                            phi=DELTAi(1)*(i-1)
                            RECT_COORDS(1)=alpha*sin(nu)*cos(phi)
                            RECT_COORDS(2)=alpha*sin(nu)*sin(phi)
                            RECT_COORDS(3)=alpha*cos(nu)
                            np=np+1
                            npg=COMPONENT_NODE_TO_USER_NUMBER(ELLIPSOID_MESH%GENERATED_MESH,basis_idx,np,ERR,ERROR)
                            CALL DECOMPOSITION_NODE_DOMAIN_GET(DECOMPOSITION,npg,MESH_COMPONENT,DOMAIN_NUMBER,ERR,ERROR,*999)
                            IF(DOMAIN_NUMBER==MY_COMPUTATIONAL_NODE) THEN
                               DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                                  CALL FIELD_PARAMETER_SET_UPDATE_NODE(FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,npg, &
                                       & component_idx,RECT_COORDS(component_idx),ERR,ERROR,*999)
                                  local_node=DOMAIN%MAPPINGS%NODES%GLOBAL_TO_LOCAL_MAP(npg)%local_number(1)
                                  IF(DOMAIN_NODES%NODES(local_node)%NUMBER_OF_DERIVATIVES>1) THEN
                                     CALL FlagError("Not generalized to hermittean elements.",ERR,ERROR,*999)
                                  ENDIF !derivatives
                               ENDDO
                            ENDIF
                         ENDDO
                      ENDDO
                   ENDDO

                ELSEIF (ELLIPSOID_EXTENT(1)<ELLIPSOID_EXTENT(2)) THEN
                   !Oblate spheroid
                   k=1
                   !inner surface
                   alpha=sqrt((ELLIPSOID_EXTENT(2))**2-(ELLIPSOID_EXTENT(1))**2)
                   !xi=log(ELLIPSOID_EXTENT(1)/alpha+sqrt((ELLIPSOID_EXTENT(1)/alpha)**2+1))
                   xi=acosh(ELLIPSOID_EXTENT(2)/alpha)

                   j=1
                   !apex node
                   np=1
                   npg=COMPONENT_NODE_TO_USER_NUMBER(ELLIPSOID_MESH%GENERATED_MESH,basis_idx,np,ERR,ERROR)
                   CALL DECOMPOSITION_NODE_DOMAIN_GET(DECOMPOSITION,npg,MESH_COMPONENT,DOMAIN_NUMBER,ERR,ERROR,*999)
                   IF(DOMAIN_NUMBER==MY_COMPUTATIONAL_NODE) THEN
                      RECT_COORDS(1)=0
                      RECT_COORDS(2)=0
                      RECT_COORDS(3)=-ELLIPSOID_EXTENT(1)
                      DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                         CALL FIELD_PARAMETER_SET_UPDATE_NODE(FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,npg, &
                              & component_idx,RECT_COORDS(component_idx),ERR,ERROR,*999)
                         local_node=DOMAIN%MAPPINGS%NODES%GLOBAL_TO_LOCAL_MAP(npg)%local_number(1)
                         IF(DOMAIN_NODES%NODES(local_node)%NUMBER_OF_DERIVATIVES>1) THEN
                            CALL FlagError("Not generalized to hermittean elements.",ERR,ERROR,*999)
                         ENDIF !derivatives
                      ENDDO
                   ENDIF

                   DO j=2,TOTAL_NUMBER_NODES_XI(2)
                      !longitudinal loop
                      nu=-PI/2+DELTAi(2)*(j-1)
                      DO i=1,TOTAL_NUMBER_NODES_XI(1)
                         !circumferential loop
                         phi=DELTAi(1)*(i-1)
                         RECT_COORDS(1)=alpha*(cosh(xi)*cos(nu)*cos(phi))
                         RECT_COORDS(2)=alpha*(cosh(xi)*cos(nu)*sin(phi))
                         RECT_COORDS(3)=alpha*(sinh(xi)*sin(nu))
                         np=np+1
                         npg=COMPONENT_NODE_TO_USER_NUMBER(ELLIPSOID_MESH%GENERATED_MESH,basis_idx,np,ERR,ERROR)
                         CALL DECOMPOSITION_NODE_DOMAIN_GET(DECOMPOSITION,npg,MESH_COMPONENT,DOMAIN_NUMBER,ERR,ERROR,*999)
                         IF(DOMAIN_NUMBER==MY_COMPUTATIONAL_NODE) THEN
                            DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                               CALL FIELD_PARAMETER_SET_UPDATE_NODE(FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,npg, &
                                    & component_idx,RECT_COORDS(component_idx),ERR,ERROR,*999)
                               local_node=DOMAIN%MAPPINGS%NODES%GLOBAL_TO_LOCAL_MAP(npg)%local_number(1)
                               IF(DOMAIN_NODES%NODES(local_node)%NUMBER_OF_DERIVATIVES>1) THEN
                                  CALL FlagError("Not generalized to hermittean elements.",ERR,ERROR,*999)
                               ENDIF !derivatives
                            ENDDO
                         ENDIF
                      ENDDO
                   ENDDO

                   DO k=2,TOTAL_NUMBER_NODES_XI(3)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
                      !transmural loop
                      nodeIdx2=1
                      !apex nodes
<<<<<<< HEAD
                      rectangularCoordinates(1)=0
                      rectangularCoordinates(2)=0
                      rectangularCoordinates(3)=-ellipsoidMesh%ellipsoidExtent(1)-(nodeIdx3-1)*(DeltaCoordinateXi(3))
                      nodeIdx=nodeIdx+1
                      globalNodeNumber=GeneratedMeshComponentNodeToUserNumber(ellipsoidMesh%generatedMesh,basisIdx,nodeIdx, &
                        & err,error)
                      CALL DECOMPOSITION_NODE_DOMAIN_GET(decomposition,globalNodeNumber,meshComponent,domainNumber,err,error,*999)
                      IF(domainNumber==myComputationalNode) THEN
                        CALL FIELD_PARAMETER_SET_UPDATE_NODE(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1, &
                          & globalNodeNumber,componentIdx,rectangularCoordinates(componentIdx),err,error,*999)
                        localNodeNumber=domain%mappings%nodes%GLOBAL_TO_LOCAL_MAP(globalNodeNumber)%local_number(1)
                        IF(domainNodes%nodes(localNodeNumber)%NUMBER_OF_DERIVATIVES>1) THEN
                          CALL FlagError("Not generalized to hermittean elements.",err,error,*999)
                        ENDIF !derivatives
                      ENDIF
                      
                      DO nodeIdx2=2,totalNumberOfNodesXi(2)
                        !longitudinal loop
                        nu=-PI/2+DeltaCoordinateXi(2)*(nodeIdx2-1)
                        DO nodeIdx1=1,totalNumberOfNodesXi(1)
                          !circumferential loop
                          phi=DeltaCoordinateXi(1)*(nodeIdx1-1)
                          x=alpha*(COSH(xi)*COS(nu)*COS(phi))
                          y=alpha*(COSH(xi)*COS(nu)*SIN(phi))
                          z=alpha*(SINH(xi)*SIN(nu))
                          !Normal vector from inner surface with length DeltaCoordinateXi(3)(nodeIdx3-1)
                          ! Finney&Thomas: Calculus, second edition, Addison-Wesley Publishing Company, 1994, page 847
                          !X=x(1+2t/a^2) Y=y(1+2t/a^2) Z=z(1+2t/c^2)
                          t=(DeltaCoordinateXi(3)*(nodeIdx3-1))/SQRT((4*x**2/(ellipsoidMesh%ellipsoidExtent(2))**4)+ &
                            & (4*y**2/(ellipsoidMesh%ellipsoidExtent(2))**4)+(4*z**2/(ellipsoidMesh%ellipsoidExtent(1))**4))
                          rectangularCoordinates(1)=x*(1+2*t/(ellipsoidMesh%ellipsoidExtent(2))**2)
                          rectangularCoordinates(2)=y*(1+2*t/(ellipsoidMesh%ellipsoidExtent(2))**2)
                          rectangularCoordinates(3)=z*(1+2*t/(ellipsoidMesh%ellipsoidExtent(1))**2)
                          nodeIdx=nodeIdx+1
                          globalNodeNumber=GeneratedMeshComponentNodeToUserNumber(ellipsoidMesh%generatedMesh,basisIdx,nodeIdx, &
                            & err,error)
                          CALL DECOMPOSITION_NODE_DOMAIN_GET(decomposition,globalNodeNumber,meshComponent,domainNumber, &
                            & err,error,*999)
                          IF(domainNumber==myComputationalNode) THEN
                            CALL FIELD_PARAMETER_SET_UPDATE_NODE(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,&
                              & globalNodeNumber,componentIdx,rectangularCoordinates(componentIdx),err,error,*999)
                            localNodeNumber=domain%mappings%nodes%GLOBAL_TO_LOCAL_MAP(globalNodeNumber)%local_number(1)
                            IF(domainNodes%nodes(localNodeNumber)%NUMBER_OF_DERIVATIVES>1) THEN
                              CALL FlagError("Not generalized to hermittean elements.",err,error,*999)
                            ENDIF !derivatives
                          ENDIF
                        ENDDO !nodeIdx1
                      ENDDO !nodeIdx2
                    ENDDO !nodeIdx3
                  ELSE
                    CALL FlagError("Not valid long axis - short axis relation",err,error,*999)
                  ENDIF
                ELSE
                  CALL FlagError("All field variable components must have node-based interpolation.",err,error,*999)
                ENDIF
              ENDDO !componentIdx
              CALL FIELD_PARAMETER_SET_UPDATE_START(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
              CALL FIELD_PARAMETER_SET_UPDATE_FINISH(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
            ELSE
              CALL FlagError("Geometric field must be three dimensional.",err,error,*999)
            ENDIF
          ELSE
            CALL FlagError("Non rectangular-cartesian coordinates are not implementd.",err,error,*999)
          ENDIF
        ELSE
          localError="Field number "//TRIM(NumberToVString(geometricField%USER_NUMBER,"*",err,error))// &
            & " is not a geometric field."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Geometric field is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Ellipsoid mesh is not associated.",err,error,*999)
=======
                      RECT_COORDS(1)=0
                      RECT_COORDS(2)=0
                      RECT_COORDS(3)=-ELLIPSOID_EXTENT(1)-(k-1)*(DELTAi(3))
                      np=np+1
                      npg=COMPONENT_NODE_TO_USER_NUMBER(ELLIPSOID_MESH%GENERATED_MESH,basis_idx,np,ERR,ERROR)
                      CALL DECOMPOSITION_NODE_DOMAIN_GET(DECOMPOSITION,npg,MESH_COMPONENT,DOMAIN_NUMBER,ERR,ERROR,*999)
                      IF(DOMAIN_NUMBER==MY_COMPUTATIONAL_NODE) THEN
                         DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                            CALL FIELD_PARAMETER_SET_UPDATE_NODE(FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,npg, &
                                 & component_idx,RECT_COORDS(component_idx),ERR,ERROR,*999)
                            local_node=DOMAIN%MAPPINGS%NODES%GLOBAL_TO_LOCAL_MAP(npg)%local_number(1)
                            IF(DOMAIN_NODES%NODES(local_node)%NUMBER_OF_DERIVATIVES>1) THEN
                               CALL FlagError("Not generalized to hermittean elements.",ERR,ERROR,*999)
                            ENDIF !derivatives
                         ENDDO
                      ENDIF

                      DO j=2,TOTAL_NUMBER_NODES_XI(2)
                         !longitudinal loop
                         nu=-PI/2+DELTAi(2)*(j-1)
                         DO i=1,TOTAL_NUMBER_NODES_XI(1)
                            !circumferential loop
                            phi=DELTAi(1)*(i-1)
                            x=alpha*(cosh(xi)*cos(nu)*cos(phi))
                            y=alpha*(cosh(xi)*cos(nu)*sin(phi))
                            z=alpha*(sinh(xi)*sin(nu))
                            !Normal vector from inner surface with length DELTAi(3)(k-1)
                            ! Finney&Thomas: Calculus, second edition, Addison-Wesley Publishing Company, 1994, page 847
                            !X=x(1+2t/a^2) Y=y(1+2t/a^2) Z=z(1+2t/c^2)
                            t=(DELTAi(3)*(k-1))/sqrt((4*x**2/(ELLIPSOID_EXTENT(2))**4)+ &
                                 & (4*y**2/(ELLIPSOID_EXTENT(2))**4)+(4*z**2/(ELLIPSOID_EXTENT(1))**4))
                            RECT_COORDS(1)=x*(1+2*t/(ELLIPSOID_EXTENT(2))**2)
                            RECT_COORDS(2)=y*(1+2*t/(ELLIPSOID_EXTENT(2))**2)
                            RECT_COORDS(3)=z*(1+2*t/(ELLIPSOID_EXTENT(1))**2)
                            np=np+1
                            npg=COMPONENT_NODE_TO_USER_NUMBER(ELLIPSOID_MESH%GENERATED_MESH,basis_idx,np,ERR,ERROR)
                            CALL DECOMPOSITION_NODE_DOMAIN_GET(DECOMPOSITION,npg,MESH_COMPONENT,DOMAIN_NUMBER,ERR,ERROR,*999)
                            IF(DOMAIN_NUMBER==MY_COMPUTATIONAL_NODE) THEN
                               DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                                  CALL FIELD_PARAMETER_SET_UPDATE_NODE(FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,npg, &
                                       & component_idx,RECT_COORDS(component_idx),ERR,ERROR,*999)
                                  local_node=DOMAIN%MAPPINGS%NODES%GLOBAL_TO_LOCAL_MAP(npg)%local_number(1)
                                  IF(DOMAIN_NODES%NODES(local_node)%NUMBER_OF_DERIVATIVES>1) THEN
                                     CALL FlagError("Not generalized to hermittean elements.",ERR,ERROR,*999)
                                  ENDIF !derivatives
                               ENDDO
                            ENDIF
                         ENDDO
                      ENDDO
                   ENDDO
                ELSE
                   CALL FlagError("Not valid long axis - short axis relation",ERR,ERROR,*999)
                ENDIF
             ELSE
                CALL FlagError("All field variable components must have node-based interpolation.",ERR,ERROR,*999)
             ENDIF
             CALL FIELD_PARAMETER_SET_UPDATE_START(FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
             CALL FIELD_PARAMETER_SET_UPDATE_FINISH(FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
          ELSE
             CALL FlagError("Geometric field must be three dimensional.",ERR,ERROR,*999)
          ENDIF
       ELSE
          LOCAL_ERROR="The standard field variable is not associated for field number "// &
               & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
       ENDIF
    ELSE
       LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" is not a geometric field."
       CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
    ENDIF
    
    ! all done
<<<<<<< HEAD
    IF(ALLOCATED(nodeIndices)) DEALLOCATE(nodeIndices)
    IF(ALLOCATED(elementIndices)) DEALLOCATE(elementIndices)
    
    CALL Exits("GeneratedMeshEllipsoidGeometricParametersCalculate")
    RETURN
999 IF(ALLOCATED(nodeIndices)) DEALLOCATE(nodeIndices)
    IF(ALLOCATED(elementIndices)) DEALLOCATE(elementIndices)
    CALL Errors("GeneratedMeshEllipsoidGeometricParametersCalculate",err,error)
    CALL Exits("GeneratedMeshEllipsoidGeometricParametersCalculate")
    RETURN 1

  END SUBROUTINE GeneratedMeshEllipsoidGeometricParametersCalculate
=======
    IF(ALLOCATED(NIDX)) DEALLOCATE(NIDX)
    IF(ALLOCATED(EIDX)) DEALLOCATE(EIDX)

    EXITS("GeneratedMesh_EllipsoidGeometricParametersCalculate")
    RETURN
999 IF(ALLOCATED(NIDX)) DEALLOCATE(NIDX)
    IF(ALLOCATED(EIDX)) DEALLOCATE(EIDX)
    ERRORS("GeneratedMesh_EllipsoidGeometricParametersCalculate",ERR,ERROR)
    EXITS("GeneratedMesh_EllipsoidGeometricParametersCalculate")
    RETURN 1

  END SUBROUTINE GeneratedMesh_EllipsoidGeometricParametersCalculate
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b

  !
  !================================================================================================================================
  !

  !>Provides an easy way to grab surfaces for boundary condition assignment
  SUBROUTINE GeneratedMeshRegularSurfaceGet(regularMesh,meshComponent,surfaceType,surfaceNodes,normalXi,err,error,*)

    ! Argument variables
    TYPE(GeneratedMeshRegularType), POINTER :: regularMesh !<A pointer to the regular mesh object
    INTEGER(INTG), INTENT(IN) :: meshComponent !<The mesh component to get the surface for.
    INTEGER(INTG), INTENT(IN) :: surfaceType !<A constant identifying the type of surface to get \see GeneratedMeshRoutines_GeneratedMeshRegularSurfaces,GeneratedMeshRoutines
    INTEGER(INTG), ALLOCATABLE, INTENT(OUT) :: surfaceNodes(:) !<On exit, contains the list of nodes belonging to the surface
    INTEGER(INTG), INTENT(OUT) :: normalXi !<On exit, contains the xi direction of the outward pointing normal of the surface
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    ! Local variables
<<<<<<< HEAD
    INTEGER(INTG) :: nodeCounter,nodeIdx1,nodeIdx2,nodeIdx3
    INTEGER(INTG),ALLOCATABLE :: nodeIndices(:,:,:),elementIndices(:,:,:)
    INTEGER(INTG) :: numberOfElementsXi(3) !Specified number of elements in each xi direction
    INTEGER(INTG) :: numberOfNodesXi(3) ! Number of nodes per element in each xi direction (basis property)
    INTEGER(INTG) :: meshDimension,totalNumberOfNodes,totalNumberOfElements,nodeUserNumber
    REAL(DP) :: deltaCoordinate(3),deltaCoordinateXi(3)
    TYPE(BASIS_TYPE), POINTER :: basis
    TYPE(VARYING_STRING) :: localError

    CALL Enters("GeneratedMeshRegularSurfaceGet",err,error,*999)

    IF(ALLOCATED(regularMesh%numberOfElementsXi)) THEN
      meshDimension=SIZE(regularMesh%numberOfElementsXi,1)
      numberOfElementsXi=1
      numberOfElementsXi(1:meshDimension)=regularMesh%numberOfElementsXi(1:meshDimension)
      basis=>regularMesh%bases(meshComponent)%ptr
      IF(ASSOCIATED(basis)) THEN
        IF(.NOT.ALLOCATED(surfaceNodes)) THEN
          !Node that simplex bases have an extra area coordinate so size of number_of_nodes_xic=meshDimension+1
          numberOfNodesXi(1:meshDimension)=basis%NUMBER_OF_NODES_XIC(1:meshDimension)
          numberOfNodesXi(meshDimension+1:3)=1
=======
    TYPE(BASIS_TYPE), POINTER :: BASIS
    INTEGER(INTG),ALLOCATABLE :: NIDX(:,:,:),EIDX(:,:,:)
    INTEGER(INTG) :: NUMBER_OF_ELEMENTS_XI(3) !Specified number of elements in each xi direction
    INTEGER(INTG) :: NUMBER_OF_NODES_XI(3) ! Number of nodes per element in each xi direction (basis property)
    INTEGER(INTG) :: num_dims,TOTAL_NUMBER_OF_NODES,TOTAL_NUMBER_OF_ELEMENTS,NODE_USER_NUMBER
    REAL(DP) :: DELTA(3),DELTAI(3)
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    INTEGER(INTG) :: node_counter,i,j,k

    ENTERS("GENERATED_MESH_REGULAR_SURFACE_GET",ERR,ERROR,*999)

    IF(ALLOCATED(REGULAR_MESH%NUMBER_OF_ELEMENTS_XI)) THEN
      num_dims=SIZE(REGULAR_MESH%NUMBER_OF_ELEMENTS_XI)
      IF(num_dims==2) THEN
        NUMBER_OF_ELEMENTS_XI(1:2)=REGULAR_MESH%NUMBER_OF_ELEMENTS_XI(1:2)
        NUMBER_OF_ELEMENTS_XI(3)=1
      ELSE IF (num_dims==1) THEN
        NUMBER_OF_ELEMENTS_XI(1)=REGULAR_MESH%NUMBER_OF_ELEMENTS_XI(1)
        NUMBER_OF_ELEMENTS_XI(2)=1
        NUMBER_OF_ELEMENTS_XI(3)=1
      ELSE
        NUMBER_OF_ELEMENTS_XI=REGULAR_MESH%NUMBER_OF_ELEMENTS_XI
      ENDIF
      IF(ASSOCIATED(REGULAR_MESH%BASES(MESH_COMPONENT)%PTR)) THEN
        BASIS=>REGULAR_MESH%BASES(MESH_COMPONENT)%PTR
        IF(.NOT.ALLOCATED(SURFACE_NODES)) THEN
          !Node that simplex bases have an extra area coordinate so size of number_of_nodes_xic=num_dims+1
          NUMBER_OF_NODES_XI(1:num_dims)=BASIS%NUMBER_OF_NODES_XIC(1:num_dims)
          NUMBER_OF_NODES_XI(num_dims+1:3)=1
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b

          ! build indices first (some of these are dummy arguments)
          CALL GeneratedMeshRegularBuildNodeIndices(numberOfElementsXi,numberOfNodesXi, &
              & regularMesh%maximumExtent,totalNumberOfNodes,totalNumberOfElements,nodeIndices,elementIndices,deltaCoordinate, &
              & deltaCoordinateXi,err,error,*999)
          nodeCounter=0
          SELECT CASE(surfaceType)
          CASE(GENERATED_MESH_REGULAR_LEFT_SURFACE)
<<<<<<< HEAD
            ALLOCATE(surfaceNodes((SIZE(nodeIndices,2))*(SIZE(nodeIndices,3))),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate NODES array.",err,error,*999)
            DO nodeIdx3=1,SIZE(nodeIndices,3)
              DO nodeIdx2=1,SIZE(nodeIndices,2)
                nodeCounter=nodeCounter+1
                surfaceNodes(nodeCounter)=nodeIndices(1,nodeIdx2,nodeIdx3)
              ENDDO !nodeIdx2
            ENDDO !nodeIdx3
            normalXi=-1
          CASE(GENERATED_MESH_REGULAR_RIGHT_SURFACE)
            ALLOCATE(surfaceNodes((SIZE(nodeIndices,2))*(SIZE(nodeIndices,3))),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate NODES array.",err,error,*999)
            DO nodeIdx3=1,SIZE(nodeIndices,3)
              DO nodeIdx2=1,SIZE(nodeIndices,2)
                nodeCounter=nodeCounter+1
                surfaceNodes(nodeCounter)=nodeIndices(SIZE(nodeIndices,1),nodeIdx2,nodeIdx3)
              ENDDO !nodeIdx3
            ENDDO !nodeIdx2
            normalXi=1
          CASE(GENERATED_MESH_REGULAR_TOP_SURFACE)
            ALLOCATE(surfaceNodes((SIZE(nodeIndices,1))*(SIZE(nodeIndices,2))),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate NODES array.",err,error,*999)
            DO nodeIdx2=1,SIZE(nodeIndices,2)
              DO nodeIdx1=1,SIZE(nodeIndices,1)
                nodeCounter=nodeCounter+1
                surfaceNodes(nodeCounter)=nodeIndices(nodeIdx1,nodeIdx2,SIZE(nodeIndices,3))
              ENDDO !nodeIdx1
            ENDDO !nodeIdx2
            normalXi=3
          CASE(GENERATED_MESH_REGULAR_BOTTOM_SURFACE)
            ALLOCATE(surfaceNodes((SIZE(nodeIndices,1))*(SIZE(nodeIndices,2))),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate NODES array.",err,error,*999)
            DO nodeIdx2=1,SIZE(nodeIndices,2)
              DO nodeIdx1=1,SIZE(nodeIndices,1)
                nodeCounter=nodeCounter+1
                surfaceNodes(nodeCounter)=nodeIndices(nodeIdx1,nodeIdx2,1)
              ENDDO !nodeIdx1
            ENDDO !nodeIdx2
            normalXi=-3
          CASE(GENERATED_MESH_REGULAR_FRONT_SURFACE)
            ALLOCATE(surfaceNodes((SIZE(nodeIndices,1))*(SIZE(nodeIndices,3))),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate NODES array.",err,error,*999)
            DO nodeIdx2=1,SIZE(nodeIndices,3)
              DO nodeIdx1=1,SIZE(nodeIndices,1)
                nodeCounter=nodeCounter+1
                surfaceNodes(nodeCounter)=nodeIndices(nodeIdx1,1,nodeIdx2)
              ENDDO !nodeIdx1
            ENDDO !nodeIdx2
            normalXi=-2
          CASE(GENERATED_MESH_REGULAR_BACK_SURFACE)
            ALLOCATE(surfaceNodes((SIZE(nodeIndices,1))*(SIZE(nodeIndices,3))),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate NODES array.",err,error,*999)
            DO nodeIdx2=1,SIZE(nodeIndices,3)
              DO nodeIdx1=1,SIZE(nodeIndices,1)
                nodeCounter=nodeCounter+1
                surfaceNodes(nodeCounter)=nodeIndices(nodeIdx1,SIZE(nodeIndices,2),nodeIdx2)
              ENDDO !nodeIdx1
            ENDDO !nodeIdx2
            normalXi=2
=======
            ALLOCATE(SURFACE_NODES((SIZE(NIDX,2))*(SIZE(NIDX,3))),STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate NODES array.",ERR,ERROR,*999)
            DO k=1,SIZE(NIDX,3)
              DO j=1,SIZE(NIDX,2)
                node_counter=node_counter+1
                SURFACE_NODES(node_counter)=NIDX(1,j,k)
              ENDDO
            ENDDO
            NORMAL_XI=-1
          CASE(GENERATED_MESH_REGULAR_RIGHT_SURFACE)
            ALLOCATE(SURFACE_NODES((SIZE(NIDX,2))*(SIZE(NIDX,3))),STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate NODES array.",ERR,ERROR,*999)
            DO k=1,SIZE(NIDX,3)
              DO j=1,SIZE(NIDX,2)
                node_counter=node_counter+1
                SURFACE_NODES(node_counter)=NIDX(SIZE(NIDX,1),j,k)
              ENDDO
            ENDDO
            NORMAL_XI=1
          CASE(GENERATED_MESH_REGULAR_TOP_SURFACE)
            ALLOCATE(SURFACE_NODES((SIZE(NIDX,1))*(SIZE(NIDX,2))),STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate NODES array.",ERR,ERROR,*999)
            DO j=1,SIZE(NIDX,2)
              DO i=1,SIZE(NIDX,1)
                node_counter=node_counter+1
                SURFACE_NODES(node_counter)=NIDX(i,j,SIZE(NIDX,3))
              ENDDO
            ENDDO
            NORMAL_XI=3
          CASE(GENERATED_MESH_REGULAR_BOTTOM_SURFACE)
            ALLOCATE(SURFACE_NODES((SIZE(NIDX,1))*(SIZE(NIDX,2))),STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate NODES array.",ERR,ERROR,*999)
            DO j=1,SIZE(NIDX,2)
              DO i=1,SIZE(NIDX,1)
                node_counter=node_counter+1
                SURFACE_NODES(node_counter)=NIDX(i,j,1)
              ENDDO
            ENDDO
            NORMAL_XI=-3
          CASE(GENERATED_MESH_REGULAR_FRONT_SURFACE)
            ALLOCATE(SURFACE_NODES((SIZE(NIDX,1))*(SIZE(NIDX,3))),STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate NODES array.",ERR,ERROR,*999)
            DO j=1,SIZE(NIDX,3)
              DO i=1,SIZE(NIDX,1)
                node_counter=node_counter+1
                SURFACE_NODES(node_counter)=NIDX(i,1,j)
              ENDDO
            ENDDO
            NORMAL_XI=-2
          CASE(GENERATED_MESH_REGULAR_BACK_SURFACE)
            ALLOCATE(SURFACE_NODES((SIZE(NIDX,1))*(SIZE(NIDX,3))),STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate NODES array.",ERR,ERROR,*999)
            DO j=1,SIZE(NIDX,3)
              DO i=1,SIZE(NIDX,1)
                node_counter=node_counter+1
                SURFACE_NODES(node_counter)=NIDX(i,SIZE(NIDX,2),j)
              ENDDO
            ENDDO
            NORMAL_XI=2
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
          CASE DEFAULT
            localError="The specified surface type of "//TRIM(NumberToVString(surfaceType,"*",err,error))// &
              & " is invalid for a regular mesh."
<<<<<<< HEAD
            CALL FlagError(localError,err,error,*999)
=======
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
          END SELECT
          !Now convert the component node numbering to user numbers if a mesh has multiple components
          DO nodeCounter=1,SIZE(surfaceNodes,1)
            SELECT CASE(regularMesh%bases(meshComponent)%ptr%TYPE)
            CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
<<<<<<< HEAD
              CALL GeneratedMeshRegularComponentNodeToUserNumber(regularMesh%generatedMesh,meshComponent, &
                & surfaceNodes(nodeCounter),nodeUserNumber,err,error,*999)
              surfaceNodes(nodeCounter)=nodeUserNumber
=======
              CALL GeneratedMesh_RegularComponentNodeToUserNumber(REGULAR_MESH%GENERATED_MESH,MESH_COMPONENT, &
                  & SURFACE_NODES(node_counter),NODE_USER_NUMBER,ERR,ERROR,*999)
              SURFACE_NODES(node_counter)=NODE_USER_NUMBER
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
            CASE(BASIS_SIMPLEX_TYPE)
              surfaceNodes(nodeCounter)=GeneratedMeshComponentNodeToUserNumber(regularMesh%generatedMesh,meshComponent, &
                & surfaceNodes(nodeCounter),err,error)
              IF(err/=0) GOTO 999
            CASE DEFAULT
<<<<<<< HEAD
              CALL FlagError("The basis type of "//TRIM(NumberToVString(regularMesh%bases(meshComponent)%ptr%type, &
                & "*",err,error))//" is not implemented when getting a regular mesh surface.",err,error,*999)
            END SELECT
          END DO
        ELSE
          CALL FlagError("Output surface nodes array is already allocated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Regular mesh object does not have a basis associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Regular mesh object does not have number of elements property specified.",err,error,*999)
    ENDIF

    CALL Exits("GeneratedMeshRegularSurfaceGet")
    RETURN
999 CALL Errors("GeneratedMeshRegularSurfaceGet",err,error)
    CALL Exits("GeneratedMeshRegularSurfaceGet")
=======
              CALL FlagError("The basis type of "//TRIM(NUMBER_TO_VSTRING(REGULAR_MESH%BASES(MESH_COMPONENT)%PTR%TYPE, &
                & "*",ERR,ERROR))//" is not implemented when getting a regular mesh surface.",ERR,ERROR,*999)
            END SELECT
          END DO
        ELSE
          CALL FlagError("Output SURFACE_NODES array is already allocated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Regular mesh object does not have a basis associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Regular mesh object does not have number of elements property specified.",ERR,ERROR,*999)
    ENDIF

    EXITS("GENERATED_MESH_REGULAR_SURFACE_GET")
    RETURN
999 ERRORSEXITS("GENERATED_MESH_REGULAR_SURFACE_GET",ERR,ERROR)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
    RETURN 1
  END SUBROUTINE GeneratedMeshRegularSurfaceGet

  !
  !================================================================================================================================
  !

  !>Provides an easy way to grab surfaces for boundary condition assignment
  SUBROUTINE GeneratedMeshCylinderSurfaceGet(cylinderMesh,meshComponent,surfaceType,surfaceNodes,normalXi,err,error,*)

    ! Argument variables
    TYPE(GeneratedMeshCylinderType), POINTER :: cylinderMesh !<A pointer to the cylinder mesh object
    INTEGER(INTG), INTENT(IN) :: meshComponent !<The mesh component to get the surface for.
    INTEGER(INTG), INTENT(IN) :: surfaceType !<A constant identifying the type of surface to get \see GeneratedMeshRoutines_GeneratedMeshCylinderSurfaces,GeneratedMeshRoutines
    INTEGER(INTG), ALLOCATABLE, INTENT(OUT) :: surfaceNodes(:) !<On exit, contains the list of nodes belonging to the surface
    INTEGER(INTG), INTENT(OUT) :: normalXi !<On exit, contains the xi direction of the outward pointing normal of the surface
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    ! Local variables
    TYPE(BASIS_TYPE), POINTER :: basis
    INTEGER(INTG),ALLOCATABLE :: nodeIndices(:,:,:),elementIndices(:,:,:)
    INTEGER(INTG) :: numberOfElementsXi(3) !Specified number of elements in each xi direction
    INTEGER(INTG) :: numberOfNodesXi(3) ! Number of nodes per element in each xi direction (basis property)
    INTEGER(INTG) :: total_number_of_nodes,total_number_of_elements
    REAL(DP) :: delta(3),deltaCoordinateXi(3)
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: nodeCounter,nodeIdx1,nodeIdx2,nodeIdx3

<<<<<<< HEAD
    CALL Enters("GeneratedMeshCylinderSurfaceGet",err,error,*999)
=======
    ENTERS("GENERATED_MESH_CYLINDER_SURFACE_GET",ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b

    ! let's go
    IF(ALLOCATED(cylinderMesh%numberOfElementsXi)) THEN
      numberOfElementsXi=cylinderMesh%numberOfElementsXi
      IF(ASSOCIATED(cylinderMesh%bases(meshComponent)%ptr)) THEN
        basis=>cylinderMesh%bases(meshComponent)%ptr
        IF(.NOT.ALLOCATED(surfaceNodes)) THEN
          numberOfNodesXi=basis%NUMBER_OF_NODES_XIC
          ! build indices first (some of these are dummy arguments)
          CALL GeneratedMeshCylinderBuildNodeIndices(numberOfElementsXi,numberOfNodesXi,cylinderMesh%cylinderExtent, &
            & total_number_of_nodes,total_number_of_elements,nodeIndices,elementIndices,delta,deltaCoordinateXi,err,error,*999)
          nodeCounter=0
          SELECT CASE(surfaceType)
          CASE(GENERATED_MESH_CYLINDER_INNER_SURFACE)
<<<<<<< HEAD
            ALLOCATE(surfaceNodes((SIZE(nodeIndices,2))*(SIZE(nodeIndices,3))),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate NODES array.",err,error,*999)
            DO nodeIdx3=1,SIZE(nodeIndices,3)
              DO nodeIdx2=1,SIZE(nodeIndices,2)
                nodeCounter=nodeCounter+1
                surfaceNodes(nodeCounter)=GeneratedMeshComponentNodeToUserNumber(cylinderMesh%generatedMesh,meshComponent, &
                  & nodeIndices(1,nodeIdx2,nodeIdx3),err,error)
              ENDDO !nodeIdx2
            ENDDO !nodeIdx3
            normalXi=-1
          CASE(GENERATED_MESH_CYLINDER_OUTER_SURFACE)
            ALLOCATE(surfaceNodes((SIZE(nodeIndices,2))*(SIZE(nodeIndices,3))),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate NODES array.",err,error,*999)
            DO nodeIdx3=1,SIZE(nodeIndices,3)
              DO nodeIdx2=1,SIZE(nodeIndices,2)
                nodeCounter=nodeCounter+1
                surfaceNodes(nodeCounter)=GeneratedMeshComponentNodeToUserNumber(cylinderMesh%generatedMesh,meshComponent, &
                  & nodeIndices(SIZE(nodeIndices,1),nodeIdx2,nodeIdx3),err,error)
              ENDDO !nodeIdx2
            ENDDO !nodeIdx3
            normalXi=1
          CASE(GENERATED_MESH_CYLINDER_TOP_SURFACE)
            ALLOCATE(surfaceNodes((SIZE(nodeIndices,1))*(SIZE(nodeIndices,2))),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate NODES array.",err,error,*999)
            DO nodeIdx2=1,SIZE(nodeIndices,2)
              DO nodeIdx1=1,SIZE(nodeIndices,1)
                nodeCounter=nodeCounter+1
                surfaceNodes(nodeCounter)=GeneratedMeshComponentNodeToUserNumber(cylinderMesh%generatedMesh,meshComponent, &
                  & nodeIndices(nodeIdx1,nodeIdx2,SIZE(nodeIndices,3)),err,error)
              ENDDO !nodeIdx1
            ENDDO !nodeIdx2
            normalXi=3
          CASE(GENERATED_MESH_CYLINDER_BOTTOM_SURFACE)
            ALLOCATE(surfaceNodes((SIZE(nodeIndices,1))*(SIZE(nodeIndices,2))),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate NODES array.",err,error,*999)
            DO nodeIdx2=1,SIZE(nodeIndices,2)
              DO nodeIdx1=1,SIZE(nodeIndices,1)
                nodeCounter=nodeCounter+1
                surfaceNodes(nodeCounter)=GeneratedMeshComponentNodeToUserNumber(cylinderMesh%generatedMesh,meshComponent, &
                  & nodeIndices(nodeIdx1,nodeIdx2,1),err,error)
              ENDDO !nodeIdx1
            ENDDO !nodeIdx2
            normalXi=-3
          CASE DEFAULT
            localError="The specified surface type of "//TRIM(NumberToVString(surfaceType,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ELSE
          CALL FlagError("Output surfaceNodes array is already allocated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Cylinder mesh object does not have a basis associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Cylinder mesh object does not have number of elements property specified.",err,error,*999)
    ENDIF

    CALL Exits("GeneratedMeshCylinderSurfaceGet")
    RETURN
999 CALL Errors("GeneratedMeshCylinderSurfaceGet",err,error)
    CALL Exits("GeneratedMeshCylinderSurfaceGet")
=======
            ALLOCATE(SURFACE_NODES((SIZE(NIDX,2))*(SIZE(NIDX,3))),STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate NODES array.",ERR,ERROR,*999)
            DO k=1,SIZE(NIDX,3)
              DO j=1,SIZE(NIDX,2)
                node_counter=node_counter+1
                SURFACE_NODES(node_counter)=COMPONENT_NODE_TO_USER_NUMBER(CYLINDER_MESH%GENERATED_MESH,MESH_COMPONENT, &
                    & NIDX(1,j,k),ERR,ERROR)
              ENDDO
            ENDDO
            NORMAL_XI=-1
          CASE(GENERATED_MESH_CYLINDER_OUTER_SURFACE)
            ALLOCATE(SURFACE_NODES((SIZE(NIDX,2))*(SIZE(NIDX,3))),STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate NODES array.",ERR,ERROR,*999)
            DO k=1,SIZE(NIDX,3)
              DO j=1,SIZE(NIDX,2)
                node_counter=node_counter+1
                SURFACE_NODES(node_counter)=COMPONENT_NODE_TO_USER_NUMBER(CYLINDER_MESH%GENERATED_MESH,MESH_COMPONENT, &
                    & NIDX(SIZE(NIDX,1),j,k),ERR,ERROR)
              ENDDO
            ENDDO
            NORMAL_XI=1
          CASE(GENERATED_MESH_CYLINDER_TOP_SURFACE)
            ALLOCATE(SURFACE_NODES((SIZE(NIDX,1))*(SIZE(NIDX,2))),STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate NODES array.",ERR,ERROR,*999)
            DO j=1,SIZE(NIDX,2)
              DO i=1,SIZE(NIDX,1)
                node_counter=node_counter+1
                SURFACE_NODES(node_counter)=COMPONENT_NODE_TO_USER_NUMBER(CYLINDER_MESH%GENERATED_MESH,MESH_COMPONENT, &
                    & NIDX(i,j,SIZE(NIDX,3)),ERR,ERROR)
              ENDDO
            ENDDO
            NORMAL_XI=3
          CASE(GENERATED_MESH_CYLINDER_BOTTOM_SURFACE)
            ALLOCATE(SURFACE_NODES((SIZE(NIDX,1))*(SIZE(NIDX,2))),STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate NODES array.",ERR,ERROR,*999)
            DO j=1,SIZE(NIDX,2)
              DO i=1,SIZE(NIDX,1)
                node_counter=node_counter+1
                SURFACE_NODES(node_counter)=COMPONENT_NODE_TO_USER_NUMBER(CYLINDER_MESH%GENERATED_MESH,MESH_COMPONENT, &
                    & NIDX(i,j,1),ERR,ERROR)
              ENDDO
            ENDDO
            NORMAL_XI=-3
          CASE DEFAULT
            LOCAL_ERROR="The specified surface type of "//TRIM(NUMBER_TO_VSTRING(SURFACE_TYPE,"*",ERR,ERROR))//" is invalid."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FlagError("Output SURFACE_NODES array is already allocated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Cylinder mesh object does not have a basis associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Cylinder mesh object does not have number of elements property specified.",ERR,ERROR,*999)
    ENDIF

    EXITS("GENERATED_MESH_CYLINDER_SURFACE_GET")
    RETURN
999 ERRORSEXITS("GENERATED_MESH_CYLINDER_SURFACE_GET",ERR,ERROR)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
    RETURN 1
  END SUBROUTINE GeneratedMeshCylinderSurfaceGet
  !
  !================================================================================================================================
  !

  !>Provides an easy way to grab surfaces for boundary condition assignment
  SUBROUTINE GeneratedMeshEllipsoidSurfaceGet(ellipsoidMesh,meshComponent,surfaceType,surfaceNodes,normalXi,err,error,*)

    ! Argument variables
    TYPE(GeneratedMeshEllipsoidType), POINTER :: ellipsoidMesh !<A pointer to the ellipsoid mesh object
    INTEGER(INTG), INTENT(IN) :: meshComponent !<The mesh component to get the surface for.
    INTEGER(INTG), INTENT(IN) :: surfaceType !<A constant identifying the type of surface to get \see GeneratedMeshRoutines_GeneratedMeshEllipsoidSurfaces,GeneratedMeshRoutines
    INTEGER(INTG), ALLOCATABLE, INTENT(OUT) :: surfaceNodes(:) !<On exit, contains the list of nodes belonging to the surface
    INTEGER(INTG), INTENT(OUT) :: normalXi !<On exit, contains the xi direction of the outward pointing normal of the surface
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    ! Local variables
    TYPE(BASIS_TYPE), POINTER :: basis
    INTEGER(INTG),ALLOCATABLE :: nodeIndices(:,:,:),elementIndices(:,:,:),cornerNodes(:,:,:)
    INTEGER(INTG) :: numberOfElementsXi(3) !Specified number of elements in each xi direction
    INTEGER(INTG) :: numberOfNodesXi(3) ! Number of nodes per element in each xi direction (basis property)
    INTEGER(INTG) :: totalNumberOfNodes,totalNumberOfElements
    REAL(DP) :: deltaCoordinate(3),deltaCoordinateXi(3)
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: nodeCounter,nodeIdx1,nodeIdx2,nodeIdx3

<<<<<<< HEAD
    CALL Enters("GeneratedMeshEllipsoidSurfaceGet",err,error,*999)
=======
    ENTERS("GENERATED_MESH_ELLIPSOID_SURFACE_GET",ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b

    ! let's go
    IF(ALLOCATED(ellipsoidMesh%numberOfElementsXi)) THEN
      numberOfElementsXi=ellipsoidMesh%numberOfElementsXi
      IF(ALLOCATED(ellipsoidMesh%bases)) THEN

!         !Below, there is an issue:
!         !  basis=>ellipsoidMesh%bases(meshComponent)%ptr does not account for the fact that:
!         !  in 'GeneratedMeshEllipsoidCreateFinish' the following is done:
!         !  CALL MESH_NUMBER_OF_COMPONENTS_SET(generatedMesh%MESH,SIZE(ellipsoidMesh%bases)/2,err,error,*999)
!         !A temporary work around is the following (although this bug may need to be fixed in several places):
!
!         IF(meshComponent==2) THEN
!           basis_COMPONENT = meshComponent + 1
!         ELSE
!           basis_COMPONENT = meshComponent
!         ENDIF
!
!         IF(ASSOCIATED(ellipsoidMesh%bases(BASIS_COMPONENT)%ptr)) THEN
!           basis=>ellipsoidMesh%bases(BASIS_COMPONENT)%ptr

        IF(ASSOCIATED(ellipsoidMesh%bases(meshComponent)%ptr)) THEN
          basis=>ellipsoidMesh%bases(meshComponent)%ptr
          IF(.NOT.ALLOCATED(surfaceNodes)) THEN
            numberOfNodesXi=basis%NUMBER_OF_NODES_XIC
            ! build indices first (some of these are dummy arguments)
            CALL GeneratedMeshEllipsoidBuildNodeIndices(numberOfElementsXi,numberOfNodesXi, &
              & ellipsoidMesh%ellipsoidExtent,totalNumberOfNodes,totalNumberOfElements,nodeIndices, &
              & cornerNodes,elementIndices,deltaCoordinate,deltaCoordinateXi,err,error,*999)
            nodeCounter=0

            SELECT CASE(surfaceType)

            CASE(GENERATED_MESH_ELLIPSOID_INNER_SURFACE)
<<<<<<< HEAD
              ALLOCATE(surfaceNodes((SIZE(nodeIndices,1))*(SIZE(nodeIndices,2)-1)+1),STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate NODES array.",err,error,*999)
              nodeIdx2=1
              nodeIdx1=1
              nodeCounter=nodeCounter+1
              surfaceNodes(nodeCounter)=GeneratedMeshComponentNodeToUserNumber(ellipsoidMesh%generatedMesh,meshComponent, &
                  nodeIndices(nodeIdx1,nodeIdx2,1),err,error)
              DO nodeIdx2=2,SIZE(nodeIndices,2)
                DO nodeIdx1=1, SIZE(nodeIndices,1)
                  nodeCounter=nodeCounter+1
                  IF (nodeIndices(nodeIdx1,nodeIdx2,1)/=0) THEN
                    surfaceNodes(nodeCounter)=GeneratedMeshComponentNodeToUserNumber(ellipsoidMesh%generatedMesh,meshComponent, &
                        & nodeIndices(nodeIdx1,nodeIdx2,1),err,error)
=======
              ALLOCATE(SURFACE_NODES((SIZE(NIDX,1))*(SIZE(NIDX,2)-1)+1),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate NODES array.",ERR,ERROR,*999)
              j=1
              i=1
              node_counter=node_counter+1
              SURFACE_NODES(node_counter)=COMPONENT_NODE_TO_USER_NUMBER(ELLIPSOID_MESH%GENERATED_MESH,MESH_COMPONENT, &
                  NIDX(i,j,1),ERR,ERROR)
              DO j=2,SIZE(NIDX,2)
                DO i=1, SIZE(NIDX,1)
                  node_counter=node_counter+1
                  IF (NIDX(i,j,1)/=0) THEN
                    SURFACE_NODES(node_counter)=COMPONENT_NODE_TO_USER_NUMBER(ELLIPSOID_MESH%GENERATED_MESH,MESH_COMPONENT, &
                        & NIDX(i,j,1),ERR,ERROR)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
                  ELSE
                    nodeCounter=nodeCounter-1
                  ENDIF
                ENDDO
              ENDDO
              normalXi=-3

            CASE(GENERATED_MESH_ELLIPSOID_OUTER_SURFACE)
<<<<<<< HEAD
              ALLOCATE(surfaceNodes((SIZE(nodeIndices,1))*(SIZE(nodeIndices,2)-1)+1),STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate NODES array.",err,error,*999)
              nodeIdx2=1
              nodeIdx1=1
              nodeCounter=nodeCounter+1
              surfaceNodes(nodeCounter)=GeneratedMeshComponentNodeToUserNumber(ellipsoidMesh%generatedMesh,meshComponent, &
                  & nodeIndices(nodeIdx1,nodeIdx2,SIZE(nodeIndices,3)),err,error)
              DO nodeIdx2=2,SIZE(nodeIndices,2)
                DO nodeIdx1=1, SIZE(nodeIndices,1)
                  nodeCounter=nodeCounter+1
                  IF (nodeIndices(nodeIdx1,nodeIdx2,SIZE(nodeIndices,3))/=0) THEN
                    surfaceNodes(nodeCounter)=GeneratedMeshComponentNodeToUserNumber(ellipsoidMesh%generatedMesh,meshComponent, &
                        & nodeIndices(nodeIdx1,nodeIdx2,SIZE(nodeIndices,3)),err,error)
=======
              ALLOCATE(SURFACE_NODES((SIZE(NIDX,1))*(SIZE(NIDX,2)-1)+1),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate NODES array.",ERR,ERROR,*999)
              j=1
              i=1
              node_counter=node_counter+1
              SURFACE_NODES(node_counter)=COMPONENT_NODE_TO_USER_NUMBER(ELLIPSOID_MESH%GENERATED_MESH,MESH_COMPONENT, &
                  & NIDX(i,j,SIZE(NIDX,3)),ERR,ERROR)
              DO j=2,SIZE(NIDX,2)
                DO i=1, SIZE(NIDX,1)
                  node_counter=node_counter+1
                  IF (NIDX(i,j,SIZE(NIDX,3))/=0) THEN
                    SURFACE_NODES(node_counter)=COMPONENT_NODE_TO_USER_NUMBER(ELLIPSOID_MESH%GENERATED_MESH,MESH_COMPONENT, &
                        & NIDX(i,j,SIZE(NIDX,3)),ERR,ERROR)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
                  ELSE
                    nodeCounter=nodeCounter-1
                  ENDIF
                ENDDO
              ENDDO
              normalXi=3

            CASE(GENERATED_MESH_ELLIPSOID_TOP_SURFACE)
<<<<<<< HEAD
              ALLOCATE(surfaceNodes((SIZE(nodeIndices,1))*(SIZE(nodeIndices,3))),STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate NODES array.",err,error,*999)
              DO nodeIdx3=1,SIZE(nodeIndices,3)
                DO nodeIdx1=1, SIZE(nodeIndices,1)
                  nodeCounter=nodeCounter+1
                  IF (nodeIndices(nodeIdx1,SIZE(nodeIndices,2),nodeIdx3)/=0) THEN
                    surfaceNodes(nodeCounter)=GeneratedMeshComponentNodeToUserNumber(ellipsoidMesh%generatedMesh,meshComponent, &
                        & nodeIndices(nodeIdx1,SIZE(nodeIndices,2),nodeIdx3),err,error)
=======
              ALLOCATE(SURFACE_NODES((SIZE(NIDX,1))*(SIZE(NIDX,3))),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate NODES array.",ERR,ERROR,*999)
              DO k=1,SIZE(NIDX,3)
                DO i=1, SIZE(NIDX,1)
                  node_counter=node_counter+1
                  IF (NIDX(i,SIZE(NIDX,2),k)/=0) THEN
                    SURFACE_NODES(node_counter)=COMPONENT_NODE_TO_USER_NUMBER(ELLIPSOID_MESH%GENERATED_MESH,MESH_COMPONENT, &
                        & NIDX(i,SIZE(NIDX,2),k),ERR,ERROR)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
                  ELSE
                    nodeCounter=nodeCounter-1
                  ENDIF
                ENDDO
              ENDDO
              normalXi=2
            CASE DEFAULT
<<<<<<< HEAD
              localError="The specified surface type of "//TRIM(NumberToVString(surfaceType,"*",err,error))//" is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          ELSE
            CALL FlagError("Output surface nodes array is already allocated.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Ellipsoid mesh object does not have the first basis associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Ellipsoid mesh object does not have bases allocated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Ellipsoid mesh object does not have number of elements property specified.",err,error,*999)
    ENDIF

    CALL Exits("GeneratedMeshEllipsoidSurfaceGet")
    RETURN
999 CALL Errors("GeneratedMeshEllipsoidSurfaceGet",err,error)
    CALL Exits("GeneratedMeshEllipsoidSurfaceGet")
=======
              LOCAL_ERROR="The specified surface type of "//TRIM(NUMBER_TO_VSTRING(SURFACE_TYPE,"*",ERR,ERROR))//" is invalid."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FlagError("Output SURFACE_NODES array is already allocated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Ellipsoid mesh object does not have the first basis associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Ellipsoid mesh object does not have bases allocated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Ellipsoid mesh object does not have number of elements property specified.",ERR,ERROR,*999)
    ENDIF

    EXITS("GENERATED_MESH_ELLIPSOID_SURFACE_GET")
    RETURN
999 ERRORSEXITS("GENERATED_MESH_ELLIPSOID_SURFACE_GET",ERR,ERROR)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
    RETURN 1
  END SUBROUTINE GeneratedMeshEllipsoidSurfaceGet

  !
  !================================================================================================================================
  !

  !>Calculates the mesh topology information for a given regular mesh (Not to be called by user)
  SUBROUTINE GeneratedMeshRegularBuildNodeIndices(numberOfElementsXi,numberOfNodesXic,maximumExtent,totalNumberOfNodes, &
    & totalNumberOfElements,nodeIndices,elementIndices,deltaCoordinate,DeltaCoordinateXi,err,error,*)

    ! Argument variables
    INTEGER(INTG),INTENT(IN) :: numberOfElementsXi(3) !<Specified number of elements in each xi direction
    INTEGER(INTG),INTENT(IN) :: numberOfNodesXic(3) !<Number of nodes per element in each xi direction for this basis
    REAL(DP),INTENT(IN) :: maximumExtent(3) !<width, length and height of regular mesh
    INTEGER(INTG),INTENT(OUT) :: totalNumberOfNodes    !<On exit, contains total number of nodes in regular mesh component
    INTEGER(INTG),INTENT(OUT) :: totalNumberOfElements !<On exit, contains total number of elements in regular mesh
    INTEGER(INTG),ALLOCATABLE,INTENT(OUT) :: nodeIndices(:,:,:)  !<Mapping array to find a node number for a given (x,y,z)
    INTEGER(INTG),ALLOCATABLE,INTENT(OUT) :: elementIndices(:,:,:)  !<Mapping array to find an element number for a given (x,y,z)
    REAL(DP),INTENT(OUT) :: deltaCoordinate(3)  !<Step sizes in each of (x,y,z) for elements
    REAL(DP),INTENT(OUT) :: DeltaCoordinateXi(3) !<Step sizes in each of (x,y,z) for node (identical to deltaCoordinate if 2 nodes per xi direction)
    INTEGER(INTG) :: ERR !<The error code
    TYPE(VARYING_STRING) :: ERROR !<The error string

    ! Local variables
<<<<<<< HEAD
    INTEGER(INTG) :: xiIdx,elementIdx1,elementIdx2,elementIdx3,localNodeIdx1,localNodeIdx2,localNodeIdx3,localNode,elementNumber
    INTEGER(INTG) :: totalNumberOfNodesXi(3) !<Total number of nodes per element in each xi direction for this basis

    CALL Enters("GeneratedMeshRegularBuildNodeIndices",err,error,*999)

    IF(.NOT.ALLOCATED(nodeIndices)) THEN
      IF(.NOT.ALLOCATED(elementIndices)) THEN
        !Calculate deltaCoordinate and DeltaCoordinateXi
        deltaCoordinate(1)=maximumExtent(1)/numberOfElementsXi(1)
        deltaCoordinate(2)=maximumExtent(2)/numberOfElementsXi(2)
        deltaCoordinate(3)=maximumExtent(3)/numberOfElementsXi(3)
        DO xiIdx=1,3
          IF(numberOfNodesXic(xiIdx)>1) DeltaCoordinateXi(xiIdx)=deltaCoordinate(xiIdx)/(numberOfNodesXic(xiIdx)-1)
        ENDDO !xiIdx

        !Calculate total elements and nodes
        DO xiIdx=1,3
          totalNumberOfNodesXi(xiIdx)=(numberOfNodesXic(xiIdx)-1)*numberOfElementsXi(xiIdx)+1
        ENDDO !xiIdx
        totalNumberOfElements=PRODUCT(numberOfElementsXi)

        ! calculate nodeIndices first
        ALLOCATE(nodeIndices(totalNumberOfNodesXi(1),totalNumberOfNodesXi(2),totalNumberOfNodesXi(3)),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate node indices array.",err,error,*999)
        localNode=0
        DO localNodeIdx3=1,totalNumberOfNodesXi(3)
          DO localNodeIdx2=1,totalNumberOfNodesXi(2)
            DO localNodeIdx1=1,totalNumberOfNodesXi(1)
              localNode=localNode+1
              nodeIndices(localNodeIdx1,localNodeIdx2,localNodeIdx3)=localNode
            ENDDO ! localNodeIdx1
          ENDDO ! localNodeIdx2
        ENDDO ! localNodeIdx3
        totalNumberOfNodes=localNode

        ! now do elementIndices
        ALLOCATE(elementIndices(numberOfElementsXi(1),numberOfElementsXi(2),numberOfElementsXi(3)),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate element indices array.",err,error,*999)
        elementNumber=0
        DO elementIdx3=1,numberOfElementsXi(3)
          DO elementIdx2=1,numberOfElementsXi(2)
            DO elementIdx1=1,numberOfElementsXi(1)
              elementNumber=elementNumber+1
              elementIndices(elementIdx1,elementIdx2,elementIdx3)=elementNumber
=======
    INTEGER(INTG) :: xi_idx,ne1,ne2,ne3,nn1,nn2,nn3,NN,NE
    INTEGER(INTG) :: TOTAL_NUMBER_NODES_XI(3) !<Total number of nodes per element in each xi direction for this basis

    ENTERS("GENERATED_MESH_REGULAR_BUILD_NODE_INDICES",ERR,ERROR,*999)

    IF(.NOT.ALLOCATED(NIDX)) THEN
      IF(.NOT.ALLOCATED(EIDX)) THEN
        ! calculate DELTA and DELTAi
        DELTA(1)=MAXIMUM_EXTENT(1)/NUMBER_ELEMENTS_XI(1)
        DELTA(2)=MAXIMUM_EXTENT(2)/NUMBER_ELEMENTS_XI(2)
        DELTA(3)=MAXIMUM_EXTENT(3)/NUMBER_ELEMENTS_XI(3)
        DO xi_idx=1,3
          IF(NUMBER_OF_NODES_XIC(xi_idx)>1) DELTAi(xi_idx)=DELTA(xi_idx)/(NUMBER_OF_NODES_XIC(xi_idx)-1)
        ENDDO

        ! calculate total elements and nodes
        DO xi_idx=1,3
          TOTAL_NUMBER_NODES_XI(xi_idx)=(NUMBER_OF_NODES_XIC(xi_idx)-1)*NUMBER_ELEMENTS_XI(xi_idx)+1
        ENDDO
        TOTAL_NUMBER_OF_ELEMENTS=PRODUCT(NUMBER_ELEMENTS_XI)

        ! calculate NIDX first
        ALLOCATE(NIDX(TOTAL_NUMBER_NODES_XI(1),TOTAL_NUMBER_NODES_XI(2),TOTAL_NUMBER_NODES_XI(3)),STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate NIDX array.",ERR,ERROR,*999)
        NN=0
        DO nn3=1,TOTAL_NUMBER_NODES_XI(3)
          DO nn2=1,TOTAL_NUMBER_NODES_XI(2)
            DO nn1=1,TOTAL_NUMBER_NODES_XI(1)
              NN=NN+1
              NIDX(nn1,nn2,nn3)=NN
            ENDDO ! nn1
          ENDDO ! nn2
        ENDDO ! nn3
        TOTAL_NUMBER_OF_NODES=NN

        ! now do EIDX
        ALLOCATE(EIDX(NUMBER_ELEMENTS_XI(1),NUMBER_ELEMENTS_XI(2),NUMBER_ELEMENTS_XI(3)),STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate EIDX array.",ERR,ERROR,*999)
        NE=0
        DO ne3=1,NUMBER_ELEMENTS_XI(3)
          DO ne2=1,NUMBER_ELEMENTS_XI(2)
            DO ne1=1,NUMBER_ELEMENTS_XI(1)
              NE=NE+1
              EIDX(ne1,ne2,ne3)=NE
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
            ENDDO
          ENDDO
        ENDDO
        totalNumberOfElements=elementNumber

      ELSE
<<<<<<< HEAD
        CALL FlagError("nodeIndices array is already allocated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("elementIndices array is already allocated.",ERR,error,*999)
    ENDIF

    CALL Exits("GeneratedMeshRegularBuildNodeIndices")
    RETURN
999 CALL Errors("GeneratedMeshRegularBuildNodeIndices",err,error)
    CALL Exits("GeneratedMeshRegularBuildNodeIndices")
=======
        CALL FlagError("NIDX array is already allocated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("EIDX array is already allocated.",ERR,error,*999)
    ENDIF

    EXITS("GENERATED_MESH_REGULAR_BUILD_NODE_INDICES")
    RETURN
999 ERRORSEXITS("GENERATED_MESH_REGULAR_BUILD_NODE_INDICES",ERR,ERROR)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
    RETURN 1
  END SUBROUTINE GeneratedMeshRegularBuildNodeIndices

  !
  !================================================================================================================================
  !

  !>Calculates the mesh topology information for a given cylinder (Not to be called by user)
  SUBROUTINE GeneratedMeshCylinderBuildNodeIndices(numberOfElementsXi,numberOfNodesXic,cylinderExtent,totalNumberOfNodes, &
    & totalNumberOfElements,nodeIndices,elementIndices,deltaCoordinate,deltaCoordinateXi,err,error,*)

    ! Argument variables
    INTEGER(INTG),INTENT(IN) :: numberOfElementsXi(3) !<Specified number of elements in each xi direction
    INTEGER(INTG),INTENT(IN) :: numberOfNodesXic(3) !<Number of nodes per element in each xi direction (basis property)
    REAL(DP),INTENT(IN) :: cylinderExtent(3) !<inner & outer radii and height of cylinder
    INTEGER(INTG),INTENT(OUT) :: totalNumberOfNodes !<On exit, contains total number of nodes in cylinder mesh
    INTEGER(INTG),INTENT(OUT) :: totalNumberOfElements !<On exit, contains total number of elements in cylinder mesh
    INTEGER(INTG),ALLOCATABLE,INTENT(OUT) :: nodeIndices(:,:,:)  !<Mapping array to find a node number for a given (r,theta,z)
    INTEGER(INTG),ALLOCATABLE,INTENT(OUT) :: elementIndices(:,:,:)  !<Mapping array to find an element number for a given (r,theta,z)
    REAL(DP),INTENT(OUT) :: deltaCoordinate(3)  !<Step sizes in each of (r,theta,z) for elements
    REAL(DP),INTENT(OUT) :: deltaCoordinateXi(3) !<Step sizes in each of (r,theta,z) for node (identical to deltaCoordinate if 2 nodes per xi direction)
    INTEGER(INTG) :: err !<The error code
    TYPE(VARYING_STRING) :: error !<The error string

    ! Local variables
    INTEGER(INTG) :: xiIdx,elementIdx1,elementIdx2,elementIdx3,localNodeIdx1,localNodeIdx2,localNodeIdx3,localNode,elementNumber
    INTEGER(INTG) :: totalNumberOfNodesXi(3) ! total number of nodes in each xi direction

<<<<<<< HEAD
    CALL Enters("GeneratedMeshCylinderBuildNodeIndices",err,error,*999)
=======
    ENTERS("GENERATED_MESH_CYLINDER_BUILD_NODE_INDICES",ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b

    ! Can skip most of the testing as this subroutine is only to be called by
    ! GeneratedMeshCylinderCreateFinish, which tests the input params.
    IF(.NOT.ALLOCATED(nodeIndices)) THEN
      IF(.NOT.ALLOCATED(elementIndices)) THEN
        ! calculate deltaCoordinate and deltaCoordinateXi
        deltaCoordinate(1)=(cylinderExtent(2)-cylinderExtent(1))/numberOfElementsXi(1)
        deltaCoordinate(2)=TWOPI/numberOfElementsXi(2)
        deltaCoordinate(3)=cylinderExtent(3)/numberOfElementsXi(3)
        DO xiIdx=1,3
          deltaCoordinateXi(xiIdx)=deltaCoordinate(xiIdx)/(numberOfNodesXic(xiIdx)-1)
        ENDDO !xiIdx

        ! calculate total elements and nodes
<<<<<<< HEAD
        DO xiIdx=1,3
          totalNumberOfNodesXi(xiIdx)=(numberOfNodesXic(xiIdx)-1)*numberOfElementsXi(xiIdx)+1
        ENDDO !xiIdx
        totalNumberOfNodesXi(2)=totalNumberOfNodesXi(2)-1 ! theta loops around so slightly different
        !totalNumberOfElements=PRODUCT(numberOfElementsXi)

        ! calculate nodeIndices first
        ALLOCATE(nodeIndices(totalNumberOfNodesXi(1),totalNumberOfNodesXi(2),totalNumberOfNodesXi(3)),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate nodeIndices array.",err,error,*999)
        localNode=0
        DO localNodeIdx3=1,totalNumberOfNodesXi(3)
          DO localNodeIdx2=1,totalNumberOfNodesXi(2)
            DO localNodeIdx1=1,totalNumberOfNodesXi(1)
              localNode=localNode+1
              nodeIndices(localNodeIdx1,localNodeIdx2,localNodeIdx3)=localNode
            ENDDO ! localNodeIdx1
          ENDDO ! localNodeIdx2
        ENDDO ! localNodeIdx3
        totalNumberOfNodes=localNode

        ! now do elementIndices
        ALLOCATE(elementIndices(numberOfElementsXi(1),numberOfElementsXi(2),numberOfElementsXi(3)),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate elementIndices array.",err,error,*999)
        elementNumber=0
        DO elementIdx3=1,numberOfElementsXi(3)
          DO elementIdx2=1,numberOfElementsXi(2)
            DO elementIdx1=1,numberOfElementsXi(1)
              elementNumber=elementNumber+1
              elementIndices(elementIdx1,elementIdx2,elementIdx3)=elementNumber
            ENDDO !elementIdx1
          ENDDO !elementIdx2
        ENDDO !elementIdx3
        totalNumberOfElements=elementNumber

      ELSE
        CALL FlagError("nodeIndices array is already allocated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("elementIndices array is already allocated.",ERR,error,*999)
    ENDIF

    CALL Exits("GeneratedMeshCylinderBuildNodeIndices")
    RETURN
999 CALL Errors("GeneratedMeshCylinderBuildNodeIndices",err,error)
    CALL Exits("GeneratedMeshCylinderBuildNodeIndices")
=======
        DO xi_idx=1,3
          TOTAL_NUMBER_NODES_XI(xi_idx)=(NUMBER_OF_NODES_XIC(xi_idx)-1)*NUMBER_ELEMENTS_XI(xi_idx)+1
        ENDDO
        TOTAL_NUMBER_NODES_XI(2)=TOTAL_NUMBER_NODES_XI(2)-1 ! theta loops around so slightly different
        !TOTAL_NUMBER_OF_ELEMENTS=PRODUCT(NUMBER_ELEMENTS_XI)

        ! calculate NIDX first
        ALLOCATE(NIDX(TOTAL_NUMBER_NODES_XI(1),TOTAL_NUMBER_NODES_XI(2),TOTAL_NUMBER_NODES_XI(3)),STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate NIDX array.",ERR,ERROR,*999)
        NN=0
        DO nn3=1,TOTAL_NUMBER_NODES_XI(3)
          DO nn2=1,TOTAL_NUMBER_NODES_XI(2)
            DO nn1=1,TOTAL_NUMBER_NODES_XI(1)
              NN=NN+1
              NIDX(nn1,nn2,nn3)=NN
            ENDDO ! nn1
          ENDDO ! nn2
        ENDDO ! nn3
        TOTAL_NUMBER_OF_NODES=NN

        ! now do EIDX
        ALLOCATE(EIDX(NUMBER_ELEMENTS_XI(1),NUMBER_ELEMENTS_XI(2),NUMBER_ELEMENTS_XI(3)),STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate EIDX array.",ERR,ERROR,*999)
        NE=0
        DO ne3=1,NUMBER_ELEMENTS_XI(3)
          DO ne2=1,NUMBER_ELEMENTS_XI(2)
            DO ne1=1,NUMBER_ELEMENTS_XI(1)
              NE=NE+1
              EIDX(ne1,ne2,ne3)=NE
            ENDDO
          ENDDO
        ENDDO
        TOTAL_NUMBER_OF_ELEMENTS=NE

      ELSE
        CALL FlagError("NIDX array is already allocated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("EIDX array is already allocated.",ERR,error,*999)
    ENDIF

    EXITS("GENERATED_MESH_CYLINDER_BUILD_NODE_INDICES")
    RETURN
999 ERRORSEXITS("GENERATED_MESH_CYLINDER_BUILD_NODE_INDICES",ERR,ERROR)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
    RETURN 1
  END SUBROUTINE GeneratedMeshCylinderBuildNodeIndices

  !
  !================================================================================================================================
  !

  !>Calculate the mesh topology information for a given ellipsoid (Not to be called by user)
  SUBROUTINE GeneratedMeshEllipsoidBuildNodeIndices(numberOfElementsXi,numberOfNodesXi,ellipsoidExtent,totalNumberOfNodes, &
    & totalNumberOfElements,nodeIndices,cornerNodes,elementIndices,deltaCoordinate,deltaCoordinateXi,err,error,*)

    ! Argument variables
    INTEGER(INTG),INTENT(IN) :: numberOfElementsXi(3) !<Specified number of elements in each xi direction
    INTEGER(INTG),INTENT(IN) :: numberOfNodesXi(3) !<Number of nodes per element in each xi direction (basis property)
    REAL(DP),INTENT(IN) :: ellipsoidExtent(4) !< long axis, short axis, wall thickness, top angle
    INTEGER(INTG),INTENT(OUT) :: totalNumberOfNodes !<On exit, contains total number of nodes in ellipsoid mesh
    INTEGER(INTG),INTENT(OUT) :: totalNumberOfElements !<On exit, contains total number of elements in ellipsoid mesh
    INTEGER(INTG),ALLOCATABLE,INTENT(OUT) :: nodeIndices(:,:,:) !<Mapping array to find a node number for a given (r,theta,z)
    INTEGER(INTG),ALLOCATABLE,INTENT(OUT) :: cornerNodes(:,:,:) ! Returns the array of corner nodes numbered
    INTEGER(INTG),ALLOCATABLE,INTENT(OUT) :: elementIndices(:,:,:)  !<Mapping array to find an element number for a given (r,theta,z)
    REAL(DP),INTENT(OUT) :: deltaCoordinate(3)  !<Step sizes in each of (r,theta,z) for elements
    REAL(DP),INTENT(OUT) :: deltaCoordinateXi(3) !<Step sizes in each of (r,theta,z) for node (identical to deltaCoordinate if 2 nodes per xi direction)
    INTEGER(INTG) :: err !<The error code
    TYPE(VARYING_STRING) :: error !<The error string

    ! Local variables
    INTEGER(INTG) :: xiIdx,elementIdx1,elementIdx2,elementIdx3,localNodeIdx1,localNodeIdx2,localNodeIdx3,nodeIdx1,nodeIdx2, &
      & nodeIdx3,localNode,elementNumber
    INTEGER(INTG) :: totalNumberOfNodesXi(3) ! total number of nodes in each xi direction

<<<<<<< HEAD
    CALL Enters("GeneratedMeshEllipsoidBuildNodeIndices",err,error,*999)
=======
    ENTERS("GENERATED_MESH_ELLIPSOID_BUILD_NODE_INDICES",ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b

    ! Can skip most of the testing as this subroutine is only to be called by
    ! GeneratedMeshEllipsoidCreateFinish, which tests the input params.
    IF(.NOT.ALLOCATED(nodeIndices)) THEN
      IF(.NOT.ALLOCATED(elementIndices)) THEN
        ! calculate deltaCoordinate and deltaCoordinateXi
        deltaCoordinate(1)=TWOPI/numberOfElementsXi(1)
        deltaCoordinate(2)=(PI-ellipsoidExtent(4))/numberOfElementsXi(2)
        deltaCoordinate(3)=ellipsoidExtent(3)/numberOfElementsXi(3)
        DO xiIdx=1,3
          deltaCoordinateXi(xiIdx)=deltaCoordinate(xiIdx)/(numberOfNodesXi(xiIdx)-1)
        ENDDO !xiIdx

        ! calculate total elements and nodes
<<<<<<< HEAD
        DO xiIdx=1,3
          totalNumberOfNodesXi(xiIdx)=(numberOfNodesXi(xiIdx)-1)*numberOfElementsXi(xiIdx)+1
        ENDDO !xiIdx
        totalNumberOfNodesXi(1)=totalNumberOfNodesXi(1)-1 ! circumferential loops around so slightly different
        totalNumberOfElements=PRODUCT(numberOfElementsXi)

        ! calculate nodeIndices first
        ALLOCATE(nodeIndices(totalNumberOfNodesXi(1),totalNumberOfNodesXi(2),totalNumberOfNodesXi(3)),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate nodeIndices array.",err,error,*999)
        ALLOCATE(cornerNodes(numberOfElementsXi(1),numberOfElementsXi(2)+1,numberOfElementsXi(3)+1),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate nodeIndices array.",err,error,*999)

        !localNode: node number inside element in certain direction
        !elementNumber: element number in certain direction
        !nodIdx: global node number in certain direction
        !localNode: Node counter
=======
        DO xi_idx=1,3
          TOTAL_NUMBER_NODES_XI(xi_idx)=(NUMBER_OF_NODES_XI(xi_idx)-1)*NUMBER_ELEMENTS_XI(xi_idx)+1
        ENDDO
        TOTAL_NUMBER_NODES_XI(1)=TOTAL_NUMBER_NODES_XI(1)-1 ! circumferential loops around so slightly different
        TOTAL_NUMBER_OF_ELEMENTS=PRODUCT(NUMBER_ELEMENTS_XI)

        ! calculate NIDX first
        ALLOCATE(NIDX(TOTAL_NUMBER_NODES_XI(1),TOTAL_NUMBER_NODES_XI(2),TOTAL_NUMBER_NODES_XI(3)),STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate NIDX array.",ERR,ERROR,*999)
        ALLOCATE(CORNER_NODES(NUMBER_ELEMENTS_XI(1),NUMBER_ELEMENTS_XI(2)+1,NUMBER_ELEMENTS_XI(3)+1),STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate NIDX array.",ERR,ERROR,*999)

        !nn: node number inside element in certain direction
        !ne: element number in certain direction
        !tn: global node number in certain direction
        !NN: Node counter
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
        !Due to one more corner node than elements in transmural direction, first shell is taken separatly
        localNode=0
        elementIdx3=1
        localNodeIdx3=1
        !Due to one more corner node than elements in longitudinal direction, apex elements are taken separatly
        elementIdx2=1
        localNodeIdx2=1
        elementIdx1=1
        localNodeIdx1=1
        !apex nodes
        localNode=localNode+1
        nodeIdx3=1
        nodeIdx2=1
        nodeIdx1=1
        nodeIndices(nodeIdx1,nodeIdx2,nodeIdx3)=localNode
        cornerNodes(elementIdx1,elementIdx2,elementIdx3)=localNode
        DO elementIdx2=1,numberOfElementsXi(2)
          DO localNodeIdx2=2,(numberOfNodesXi(2))
            nodeIdx2=nodeIdx2+1
            nodeIdx1=0
            DO elementIdx1=1,numberOfElementsXi(1)
              DO localNodeIdx1=1,(numberOfNodesXi(1)-1)
                nodeIdx1=nodeIdx1+1
                localNode=localNode+1
                nodeIndices(nodeIdx1,nodeIdx2,nodeIdx3)=localNode
                IF ((localNodeIdx1==1).AND.(localNodeIdx2==numberOfNodesXi(2))) THEN
                  cornerNodes(elementIdx1,elementIdx2+1,elementIdx3)=localNode
                ENDIF
              ENDDO !localNodeIdx1
            ENDDO !elementIdx1
          ENDDO !localNodeIdx2
        ENDDO !elementIdx2
        DO elementIdx3=1,numberOfElementsXi(3)
          DO localNodeIdx3=2,numberOfNodesXi(3)
            elementIdx2=1
            localNodeIdx2=1
            elementIdx1=1
            localNodeIdx1=1
            !apex nodes
            localNode=localNode+1
            nodeIdx3=nodeIdx3+1
            nodeIdx2=1
            nodeIdx1=1
            nodeIndices(nodeIdx1,nodeIdx2,nodeIdx3)=localNode
            IF (localNodeIdx3==numberOfNodesXi(3)) THEN
              cornerNodes(elementIdx1,elementIdx2,elementIdx3+1)=localNode
            ENDIF
            DO elementIdx2=1,numberOfElementsXi(2)
              DO localNodeIdx2=2,(numberOfNodesXi(2))
                nodeIdx2=nodeIdx2+1
                nodeIdx1=0
                DO elementIdx1=1,numberOfElementsXi(1)
                  DO localNodeIdx1=1,(numberOfNodesXi(1)-1)
                    nodeIdx1=nodeIdx1+1
                    localNode=localNode+1
                    nodeIndices(nodeIdx1,nodeIdx2,nodeIdx3)=localNode
                    IF ((localNodeIdx1==1).AND.(localNodeIdx3==numberOfNodesXi(3)).AND.(localNodeIdx2==numberOfNodesXi(2))) THEN
                      cornerNodes(elementIdx1,elementIdx2+1,elementIdx3+1)=localNode
                    ENDIF
<<<<<<< HEAD
                  ENDDO !localNodeIdx1
                ENDDO !elementIdx1
              ENDDO !localNodeIdx2
            ENDDO !elementIdx2
          ENDDO !localNodeIdx3
        ENDDO !elementIdx3
        totalNumberOfNodes=localNode

        ! now do elementIndices
        ALLOCATE(elementIndices(numberOfElementsXi(1),numberOfElementsXi(2),numberOfElementsXi(3)),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate elementIndices array.",err,error,*999)
        elementNumber=0
        DO elementIdx3=1,numberOfElementsXi(3)
          DO elementIdx2=1,numberOfElementsXi(2)
            DO elementIdx1=1,numberOfElementsXi(1)
              elementNumber=elementNumber+1
              elementIndices(elementIdx1,elementIdx2,elementIdx3)=elementNumber
            ENDDO !elementIdx1
          ENDDO !elementIdx2
        ENDDO !elementIdx3
        totalNumberOfElements=elementNumber

      ELSE
        CALL FlagError("nodeIndices array is already allocated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("elementIndices array is already allocated.",ERR,error,*999)
    ENDIF

    CALL Exits("GeneratedMeshEllipsoidBuildNodeIndices")
    RETURN
999 CALL Errors("GeneratedMeshEllipsoidBuildNodeIndices",err,error)
    CALL Exits("GeneratedMeshEllipsoidBuildNodeIndices")
=======
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        TOTAL_NUMBER_OF_NODES=NN

        ! now do EIDX
        ALLOCATE(EIDX(NUMBER_ELEMENTS_XI(1),NUMBER_ELEMENTS_XI(2),NUMBER_ELEMENTS_XI(3)),STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate EIDX array.",ERR,ERROR,*999)
        NE=0
        DO ne3=1,NUMBER_ELEMENTS_XI(3)
          DO ne2=1,NUMBER_ELEMENTS_XI(2)
            DO ne1=1,NUMBER_ELEMENTS_XI(1)
              NE=NE+1
              EIDX(ne1,ne2,ne3)=NE
            ENDDO
          ENDDO
        ENDDO
        TOTAL_NUMBER_OF_ELEMENTS=NE

      ELSE
        CALL FlagError("NIDX array is already allocated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("EIDX array is already allocated.",ERR,error,*999)
    ENDIF

    EXITS("GENERATED_MESH_ELLIPSOID_BUILD_NODE_INDICES")
    RETURN
999 ERRORSEXITS("GENERATED_MESH_ELLIPSOID_BUILD_NODE_INDICES",ERR,ERROR)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
    RETURN 1
  END SUBROUTINE GeneratedMeshEllipsoidBuildNodeIndices

  !
  !================================================================================================================================
  !

  !>Calculates the user node numbers for an array of nodes numbered using one basis
  SUBROUTINE GeneratedMeshComponentNodesToUserNumbers(generatedMesh,basisIndex,nodeComponentNumbers,nodeUserNumbers,err,error,*)
    ! Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh object
    INTEGER(INTG),INTENT(IN) :: basisIndex !<The number of the basis being used
    INTEGER(INTG),INTENT(IN) :: nodeComponentNumbers(:) !<The node numbers for this component basis
    INTEGER(INTG),INTENT(INOUT) :: nodeUserNumbers(:) !<On return, the corresponding user numbers
    INTEGER(INTG) :: err !<The error code
    TYPE(VARYING_STRING) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: nodeIdx

<<<<<<< HEAD
    CALL Enters("GeneratedMeshComponentNodesToUserNumbers",err,error,*999)
=======
    ENTERS("COMPONENT_NODES_TO_USER_NUMBERS",ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b

    IF(SIZE(nodeUserNumbers)==SIZE(nodeComponentNumbers)) THEN
      DO nodeIdx=1,SIZE(nodeComponentNumbers)
        nodeUserNumbers(nodeIdx)=GeneratedMeshComponentNodeToUserNumber(generatedMesh,basisIndex,nodeComponentNumbers(nodeIdx), &
          & err,error)
      ENDDO !nodeIdx
    ELSE
<<<<<<< HEAD
      CALL FlagError("nodeComponentNumbers and nodeUserNumbers arrays have different sizes.",err,error,*999)
    ENDIF

    CALL Exits("GeneratedMeshComponentNodesToUserNumbers")
    RETURN
999 CALL Errors("GeneratedMeshComponentNodesToUserNumbers",err,error)
    CALL Exits("GeneratedMeshComponentNodesToUserNumbers")
=======
      CALL FlagError("NODE_COMPONENT_NUMBERS and NODE_USER_NUMBERS arrays have different sizes.",ERR,ERROR,*999)
    ENDIF

    EXITS("COMPONENT_NODES_TO_USER_NUMBERS")
    RETURN
999 ERRORSEXITS("COMPONENT_NODES_TO_USER_NUMBERS",ERR,ERROR)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
    RETURN 1
    
  END SUBROUTINE GeneratedMeshComponentNodesToUserNumbers

  !
  !================================================================================================================================
  !

  !>Calculates the user node number for a node numbered using one basis
<<<<<<< HEAD
  FUNCTION GeneratedMeshComponentNodeToUserNumber(generatedMesh,basisIndex,nodeComponentNumber,err,error)
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh object
    INTEGER(INTG),INTENT(IN) :: basisIndex !<The number of the basis being used
    INTEGER(INTG),INTENT(IN) :: nodeComponentNumber !<The node number for this component basis
    INTEGER(INTG) :: err !<The error code
    TYPE(VARYING_STRING) :: error !<The error string
    !Function variable
    INTEGER(INTG) :: GeneratedMeshComponentNodeToUserNumber !<On return, the corresponding user node number

    !Local variables
    INTEGER(INTG) :: numberOfBases,meshDimension,basisIdx,xiIdx,remainder,remainder2,temporaryTerm,numberCornerNodes,nodeOffset, &
      & basisNumberOfNodes,position(3),position2(3),cornerNodeFactor(3),basisNodeFactor(3),basisElementFactor(3), &
      & numberOfPreviousCorners,step,numberOfElementsXi(3)
    TYPE(BASIS_TYPE), POINTER :: basis
    TYPE(BASIS_PTR_TYPE), POINTER :: bases(:)
    LOGICAL :: cornerNode,finishedCount
    TYPE(VARYING_STRING) :: localError

    NULLIFY(basis)
    NULLIFY(bases)
    
    CALL Enters("GeneratedMeshComponentNodeToUserNumber",err,error,*999)

    numberCornerNodes=1
    remainder=nodeComponentNumber-1 !use zero based numbering
    remainder2=nodeComponentNumber-1
    GeneratedMeshComponentNodeToUserNumber=0
    position=0
    position2=0

    IF(ASSOCIATED(generatedMesh)) THEN
      SELECT CASE(generatedMesh%generatedType)
=======
  FUNCTION COMPONENT_NODE_TO_USER_NUMBER(GENERATED_MESH,BASIS_INDEX,NODE_COMPONENT_NUMBER,ERR,ERROR)
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH  !<A pointer to the generated mesh object
    INTEGER(INTG),INTENT(IN) :: BASIS_INDEX                  !<The number of the basis being used
    INTEGER(INTG),INTENT(IN) :: NODE_COMPONENT_NUMBER           !<The node number for this component basis
    INTEGER(INTG) :: ERR          !<The error code
    TYPE(VARYING_STRING) :: ERROR !<The error string
    !function variable
    INTEGER(INTG) :: COMPONENT_NODE_TO_USER_NUMBER !<On return, the corresponding user node number

    !local variables
    INTEGER(INTG) :: NUM_BASES,NUM_DIMS,basis_idx,ni,REMAINDER,REMAINDER2,TEMP_TERM,NUM_CORNER_NODES,NODE_OFFSET,BASIS_NUM_NODES
    INTEGER(INTG) :: POS(3),POS2(3),CORNER_NODE_FACTOR(3),BASIS_NODE_FACTOR(3),BASIS_ELEMENT_FACTOR(3),NUM_PREVIOUS_CORNERS,STEP
    INTEGER(INTG), POINTER :: NUMBER_OF_ELEMENTS_XI(:)
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(BASIS_PTR_TYPE), POINTER :: BASES(:)
    LOGICAL :: CORNER_NODE,FINISHED_COUNT
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("COMPONENT_NODE_TO_USER_NUMBER",ERR,ERROR,*999)

    NULLIFY(BASIS)
    NULLIFY(BASES)
    NUM_CORNER_NODES=1
    REMAINDER=NODE_COMPONENT_NUMBER-1 !use zero based numbering
    REMAINDER2=NODE_COMPONENT_NUMBER-1
    COMPONENT_NODE_TO_USER_NUMBER=0
    POS=0
    POS2=0

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      SELECT CASE(GENERATED_MESH%GENERATED_TYPE)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
      CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
        IF(ASSOCIATED(generatedMesh%regularMesh)) THEN
          numberOfBases=SIZE(generatedMesh%regularMesh%bases)
          meshDimension=generatedMesh%regularMesh%meshDimension
          bases=>generatedMesh%regularMesh%bases
          numberOfElementsXi(1:meshDimension)=generatedMesh%regularMesh%numberOfElementsXi(1:meshDimension)
        ELSE
<<<<<<< HEAD
          CALL FlagError("The regular mesh for this generated mesh is not associated.",err,error,*999)
        ENDIF
      CASE(GENERATED_MESH_POLAR_MESH_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
=======
          CALL FlagError("The regular mesh for this generated mesh is not associated.",ERR,ERROR,*999)
        ENDIF
      CASE(GENERATED_MESH_POLAR_MESH_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
      CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
        IF(ASSOCIATED(generatedMesh%cylinderMesh)) THEN
          numberOfBases=SIZE(generatedMesh%cylinderMesh%bases)
          meshDimension=generatedMesh%cylinderMesh%meshDimension
          bases=>generatedMesh%cylinderMesh%bases
          numberOfElementsXi(1:meshDimension)=generatedMesh%cylinderMesh%numberOfElementsXi(1:meshDimension)
        ELSE
<<<<<<< HEAD
          CALL FlagError("The cylinder mesh for this generated mesh is not associated.",err,error,*999)
=======
          CALL FlagError("The cylinder mesh for this generated mesh is not associated.",ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
        ENDIF
      CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
        IF(ASSOCIATED(generatedMesh%ellipsoidMesh)) THEN
          numberOfBases=SIZE(generatedMesh%ellipsoidMesh%bases)
          meshDimension=generatedMesh%ellipsoidMesh%meshDimension
          bases=>generatedMesh%ellipsoidMesh%bases
          numberOfElementsXi(1:meshDimension)=generatedMesh%ellipsoidMesh%numberOfElementsXi(1:meshDimension)
        ELSE
<<<<<<< HEAD
          CALL FlagError("The ellipsoid mesh for this generated mesh is not associated.",err,error,*999)
        ENDIF
      CASE DEFAULT
        localError="The generated mesh generated type of "//TRIM(NumberToVString(generatedMesh%generatedType,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
=======
        CALL FlagError("The ellipsoid mesh for this generated mesh is not associated.",ERR,ERROR,*999)
        ENDIF
      CASE DEFAULT
        LOCAL_ERROR="The generated mesh generated type of "// &
            & TRIM(NUMBER_TO_VSTRING(GENERATED_MESH%GENERATED_TYPE,"*",ERR,ERROR))//" is invalid."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
      END SELECT
      IF(basisIndex<=numberOfBases) THEN
        IF(numberOfBases==1) THEN
          !If is the only basis, don't do anything
          GeneratedMeshComponentNodeToUserNumber=nodeComponentNumber
        ELSE
          temporaryTerm=1
          numberCornerNodes=1
          DO xiIdx=1,meshDimension
            numberCornerNodes=numberCornerNodes*(numberOfElementsXi(xiIdx)+1)
            cornerNodeFactor(xiIdx)=1
            IF(xiIdx>1) THEN
              temporaryTerm=temporaryTerm*(numberOfElementsXi(xiIdx-1)+1)
              cornerNodeFactor(xiIdx)=cornerNodeFactor(xiIdx)*temporaryTerm
            ENDIF
          ENDDO !xiIdx
          !Adjust for other mesh types
          IF(generatedMesh%generatedType==GENERATED_MESH_CYLINDER_MESH_TYPE) THEN
            cornerNodeFactor(3)=cornerNodeFactor(3)-numberOfElementsXi(1)-1
            numberCornerNodes=numberCornerNodes-(numberOfElementsXi(1)+1)*(numberOfElementsXi(3)+1)
          ELSE IF(generatedMesh%generatedType==GENERATED_MESH_ELLIPSOID_MESH_TYPE) THEN
            cornerNodeFactor(3)=cornerNodeFactor(3)-numberOfElementsXi(1)-numberOfElementsXi(2)
            cornerNodeFactor(2)=cornerNodeFactor(2)-1
            numberCornerNodes=numberCornerNodes-(numberOfElementsXi(2)+1)*(numberOfElementsXi(3)+1)- &
                & (numberOfElementsXi(1)-1)*(numberOfElementsXi(3)+1)
          ENDIF
          nodeOffset=numberCornerNodes
          IF(generatedMesh%generatedType==GENERATED_MESH_ELLIPSOID_MESH_TYPE) THEN
            !Every second mesh component is the collapsed node version
            step=2
          ELSE
            step=1
          ENDIF
          DO basisIdx=1,basisIndex-1,step
            basis=>bases(basisIdx)%ptr
            basisNumberOfNodes=1
            DO xiIdx=1,meshDimension
              basisNumberOfNodes=basisNumberOfNodes*(numberOfElementsXi(xiIdx)*(basis%NUMBER_OF_NODES_XIC(xiIdx)-1)+1)
            ENDDO !xiIdx
            !Adjust for other mesh types
            IF(generatedMesh%generatedType==GENERATED_MESH_CYLINDER_MESH_TYPE) THEN
              basisNumberOfNodes=basisNumberOfNodes-(numberOfElementsXi(1)+1)*(basis%NUMBER_OF_nodes_xic(1)-1)* &
                  & (numberOfElementsXi(3)+1)*(basis%NUMBER_OF_nodes_xic(3)-1)
            ELSE IF(generatedMesh%generatedType==GENERATED_MESH_ELLIPSOID_MESH_TYPE) THEN
              basisNumberOfNodes=basisNumberOfNodes-(numberOfElementsXi(2)*(basis%NUMBER_OF_NODES_XIC(2)-1)+1)* &
                & (numberOfElementsXi(3)*(basis%NUMBER_OF_NODES_XIC(3)-1)+1)- &
                & (numberOfElementsXi(1)*(basis%NUMBER_OF_NODES_XIC(1)-1)-1)* &
                & (numberOfElementsXi(3)*(basis%NUMBER_OF_NODES_XIC(3)-1)+1)
            ENDIF
            nodeOffset=nodeOffset+basisNumberOfNodes-numberCornerNodes
          ENDDO !basisIdx
          basis=>bases(basisIndex)%ptr
          temporaryTerm=1
          DO xiIdx=1,meshDimension
            basisNodeFactor(xiIdx)=1
            basisElementFactor(xiIdx)=basis%NUMBER_OF_NODES_XIC(xiIdx)-1
            IF(xiIdx>1) THEN
              temporaryTerm=temporaryTerm*((basis%NUMBER_OF_NODES_XIC(xiIdx-1)-1)*numberOfElementsXi(xiIdx-1)+1)
              basisNodeFactor(xiIdx)=basisNodeFactor(xiIdx)*temporaryTerm
              basisElementFactor(xiIdx)=basisElementFactor(xiIdx)*temporaryTerm
            ENDIF
          ENDDO !xiIdx
          !Adjust for other mesh types
          IF(generatedMesh%generatedType==GENERATED_MESH_CYLINDER_MESH_TYPE) THEN
            !subtract nodes along line where y wraps around
            basisNodeFactor(3)=basisNodeFactor(3)-numberOfElementsXi(1)*(basis%NUMBER_OF_NODES_XIC(1)-1)-1
            basisElementFactor(3)=basisElementFactor(3)-(numberOfElementsXi(1)* &
              & (basis%NUMBER_OF_NODES_XIC(1)-1)+1)*(basis%NUMBER_OF_NODES_XIC(3)-1)
          ELSE IF(generatedMesh%generatedType==GENERATED_MESH_ELLIPSOID_MESH_TYPE) THEN
            !subtract missing nodes at apex
            basisNodeFactor(3)=basisNodeFactor(3)-numberOfElementsXi(1)*(basis%NUMBER_OF_NODES_XIC(1)-1)+1
            basisElementFactor(3)=basisElementFactor(3)-(numberOfElementsXi(1)* &
              & (basis%NUMBER_OF_NODES_XIC(1)-1)+1)*(basis%NUMBER_OF_NODES_XIC(3)-1)
            !subtract nodes along line where x wraps around
            basisNodeFactor(3)=basisNodeFactor(3)-numberOfElementsXi(2)*(basis%NUMBER_OF_NODES_XIC(2)-1)-1
            basisElementFactor(3)=basisElementFactor(3)-(numberOfElementsXi(2)*(basis%NUMBER_OF_NODES_XIC(2)-1)-1)* &
              & (basis%NUMBER_OF_NODES_XIC(3)-1)
            basisNodeFactor(2)=basisNodeFactor(2)-1
            basisElementFactor(2)=basisElementFactor(2)-(basis%NUMBER_OF_NODES_XIC(2)-1)
          ENDIF
          !Work out if we have a corner node, otherwise add node numbers used by corners and
          !previous basis interpolations and subtract number of corner nodes used before the
          !given component node number to get the user number
          cornerNode=.TRUE.
          IF(meshDimension>2) THEN
            position(3)=remainder/basisNodeFactor(3)
            position2(3)=remainder2/basisElementFactor(3)
            remainder=MOD(remainder,basisNodeFactor(3))
            remainder2=MOD(remainder2,basisElementFactor(3))
            IF(MOD(position(3),basis%NUMBER_OF_NODES_XIC(3)-1)/=0) cornerNode=.FALSE.
          ENDIF
          IF(meshDimension>1) THEN
            IF(generatedMesh%generatedType==GENERATED_MESH_ELLIPSOID_MESH_TYPE) THEN
              !Need to account for missing nodes at apex
              IF(remainder>0) THEN
                remainder=remainder+numberOfElementsXi(1)*(basis%NUMBER_OF_NODES_XIC(1)-1)-1
                remainder2=remainder2+numberOfElementsXi(1)*(basis%NUMBER_OF_NODES_XIC(1)-1)-1
              ENDIF
            ENDIF
            position(2)=remainder/basisNodeFactor(2)
            position2(2)=remainder2/basisElementFactor(2)
            remainder=MOD(remainder,basisNodeFactor(2))
            remainder2=MOD(remainder2,basisElementFactor(2))
            IF(MOD(position(2),basis%NUMBER_OF_NODES_XIC(2)-1)/=0) cornerNode=.FALSE.
          ENDIF
          position(1)=remainder/basisNodeFactor(1)
          position2(1)=remainder2/basisElementFactor(1)
          IF(MOD(position(1),basis%NUMBER_OF_NODES_XIC(1)-1)/=0) cornerNode=.FALSE.
          IF(cornerNode) THEN
            GeneratedMeshComponentNodeToUserNumber=position2(1)*cornerNodeFactor(1)+position2(2)*cornerNodeFactor(2)+ &
                & position2(3)*cornerNodeFactor(3)
            IF(generatedMesh%generatedType==GENERATED_MESH_ELLIPSOID_MESH_TYPE.AND.position2(2)/=0) THEN
              !Subtract off non-existent nodes at apex
              GeneratedMeshComponentNodeToUserNumber=GeneratedMeshComponentNodeToUserNumber-(numberOfElementsXi(1)-1)
            ENDIF
            GeneratedMeshComponentNodeToUserNumber=GeneratedMeshComponentNodeToUserNumber+1
          ELSE
            !subtract previous corner nodes from node offset
            numberOfPreviousCorners=0
            finishedCount=.FALSE.
            IF(meshDimension>2) THEN
              IF(MOD(position(3),basis%NUMBER_OF_NODES_XIC(3)-1)/=0) THEN
                numberOfPreviousCorners=numberOfPreviousCorners+cornerNodeFactor(3)*(position2(3)+1)
                finishedCount=.TRUE.
              ELSE
                numberOfPreviousCorners=numberOfPreviousCorners+cornerNodeFactor(3)*position2(3)
              ENDIF
            ENDIF
            IF((meshDimension>1) .AND. (finishedCount.NEQV..TRUE.)) THEN
              IF(MOD(position(2),basis%NUMBER_OF_NODES_XIC(2)-1)/=0) THEN
                numberOfPreviousCorners=numberOfPreviousCorners+cornerNodeFactor(2)*(position2(2)+1)
                finishedCount=.TRUE.
              ELSE
                numberOfPreviousCorners=numberOfPreviousCorners+cornerNodeFactor(2)*position2(2)
              ENDIF
              IF(generatedMesh%generatedType==GENERATED_MESH_ELLIPSOID_MESH_TYPE) THEN
                numberOfPreviousCorners=numberOfPreviousCorners-(numberOfElementsXi(1)-1)
              ENDIF
            ENDIF
            IF(finishedCount.NEQV..TRUE.) THEN
              IF(MOD(position(1),basis%NUMBER_OF_NODES_XIC(1)-1)/=0) THEN
                numberOfPreviousCorners=numberOfPreviousCorners+cornerNodeFactor(1)*(position2(1)+1)
              ELSE
                numberOfPreviousCorners=numberOfPreviousCorners+cornerNodeFactor(1)*position2(1)
              ENDIF
            ENDIF
            nodeOffset=nodeOffset-numberOfPreviousCorners
            GeneratedMeshComponentNodeToUserNumber=nodeOffset+nodeComponentNumber
          ENDIF
        ENDIF
      ELSE
<<<<<<< HEAD
        localError="Mesh component must be less than or equal to "//(NumberToVString(numberOfBases,"*",err,error))// &
            & " but it is "//(NumberToVString(basisIndex,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Generated mesh is not associated",err,error,*999)
    ENDIF

    CALL Exits("GeneratedMeshComponentNodeToUserNumber")
    RETURN
999 CALL Errors("GeneratedMeshComponentNodeToUserNumber",err,error)
    CALL Exits("GeneratedMeshComponentNodeToUserNumber")
=======
        LOCAL_ERROR="Mesh component must be less than or equal to "//(NUMBER_TO_VSTRING(NUM_BASES,"*",ERR,ERROR))// &
            & " but it is "//(NUMBER_TO_VSTRING(BASIS_INDEX,"*",ERR,ERROR))//"."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Generated mesh is not associated",ERR,ERROR,*999)
    ENDIF

    EXITS("COMPONENT_NODE_TO_USER_NUMBER")
    RETURN
999 ERRORSEXITS("COMPONENT_NODE_TO_USER_NUMBER",ERR,ERROR)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
    RETURN
  END FUNCTION GeneratedMeshComponentNodeToUserNumber

  !
  !================================================================================================================================
  !
  !>Calculates the user node numbers for an array of nodes numbered using one basis for regular mesh type

  !1. For the current mesh component/basis, search previous basis to see if the
  !current basis has occurred.
  !2(1). If occurred, reuse user node number (i.e. same mesh topology)--> finish.
  !2(2). If not occurred (i.e. different mesh topology), reuse corner nodes
  !3. Search previous basis to see if current interpolation scheme in xi1/2/3
  !direction has occurred in the same xi direction if previous basis.
  !4(1). If occurred in xi1/2/3 direction, reuse node user numbers on
  !corresponding edges/faces. e.g. linear-quadratic scheme v.s. biquadratic
  !scheme, then node user numbers on edges alone xi2 direction can be reused.
  !4(2). If never occurred (i.e. completely different basis. e.g. biquadratic v.s.
  !bicubic), do nothing.
  !5. Search previous basis to find the largest node user number, any new node
  !user number will increment based on the current largest.
  !6. Give node user numbers to nodes that have never appeared in previous
  !basis.--> finish.

<<<<<<< HEAD
  SUBROUTINE GeneratedMeshRegularComponentNodesToUserNumbers(generatedMesh,basisIndex,nodeComponentNumbers,nodeUserNumbers, &
    & err,error,*)
=======
  SUBROUTINE GeneratedMesh_RegularComponentNodesToUserNumbers(GENERATED_MESH,BASIS_INDEX, &
      & NODE_COMPONENT_NUMBERS,NODE_USER_NUMBERS,ERR,ERROR,*)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh object
    INTEGER(INTG),INTENT(IN) :: basisIndex !<The number of the basis being used
    INTEGER(INTG),INTENT(IN) :: nodeComponentNumbers(:) !<The node numbers for this component basis
    INTEGER(INTG),INTENT(INOUT) :: nodeUserNumbers(:) !<On return, the corresponding user numbers
    INTEGER(INTG) :: err !<The error code
    TYPE(VARYING_STRING) :: error !<The error string
    !Local variables

<<<<<<< HEAD
    TYPE(BASIS_PTR_TYPE), POINTER :: bases(:)
    TYPE(BASIS_TYPE), POINTER :: basisFirstComponent,basisPrevious
    INTEGER(INTG) :: numberOfBases,meshDimension,nodeOffsetLastBasis,lastElementNumber,nodeOffsetElement,offsetUnit,elementNumber
    INTEGER(INTG) :: nodeOffsetXi2Accumulated,nodeOffsetXi2,nodeOffset,nodeOffsetXi3Accumulated
    INTEGER(INTG) :: nodeIdxCurrent,nodeIdxFirst,nodeIdxPrevious
    INTEGER(INTG) :: nodeIdx,localNodeIdx1,localNodeIdx2,localNodeIdx3,xiIdx,basisIdx
    INTEGER(INTG) :: elementIdx(3),sameBasis(3),numberOfNodesXic(3),numberOfElementsXi(3),remainderTemp
    INTEGER(INTG) :: nodeCount,indexCount,zeroCountXi1(16)
    INTEGER(INTG) :: zeroCountXi12(4),edgeNode(16),totalNumberZeroNodes,nodeOffsetElementXi12
    INTEGER(INTG) :: numberofNodesLayer,numberLocalNodes,localNodeIdx
    LOGICAL::basisAppeared
    TYPE(VARYING_STRING) :: localError

    CALL Enters("GeneratedMeshRegularComponentNodesToUserNumbers",err,error,*999)

    IF(SIZE(nodeUserNumbers)==SIZE(nodeComponentNumbers)) THEN
      nodeUserNumbers=0
      IF(ASSOCIATED(generatedMesh)) THEN
        IF(ASSOCIATED(generatedMesh%regularMesh)) THEN
          numberOfBases=SIZE(generatedMesh%regularMesh%bases)
          meshDimension=generatedMesh%regularMesh%meshDimension
          bases=>generatedMesh%regularMesh%bases
          numberOfElementsXi=1
          DO xiIdx=1,meshDimension
            numberOfElementsXi(xiIdx)=generatedMesh%regularMesh%numberOfElementsXi(xiIdx)
          ENDDO !xiIdx
        ELSE
          CALL FlagError("The regular mesh for this generated mesh is not associated.",err,error,*999)
=======
    TYPE(BASIS_PTR_TYPE), POINTER :: BASES(:)
    TYPE(BASIS_TYPE), POINTER :: BASIS_FIRST_COMP,BASIS_PRE
    INTEGER(INTG) :: NUM_BASES,NUM_DIMS,NODE_OFFSET_LAST_BASIS,LAST_ELEM_NO,NODE_OFFSET_ELEM,OFFSET_UNIT,ELEMENT_NO
    INTEGER(INTG) :: NODE_OFFSET_XI2_ACCUM,NODE_OFFSET_XI2,NODE_OFFSET,NODE_OFFSET_XI3_ACCUM
    INTEGER(INTG) :: NODE_IDX_CUR,NODE_IDX_FIRST,NODE_IDX_PRE
    INTEGER(INTG) :: node_idx,nn1,nn2,nn3,xi_idx,basis_idx
    INTEGER(INTG) :: ELEM_IDX(3),SAME_BASIS(3),NUMBER_OF_NODES_XIC(3),NUMBER_OF_ELEMENTS_XI(3),REMINDER_TEMP
    INTEGER(INTG) :: number_of_nodes_temp,node_index_temp,NODE_COUNT,INDEX_COUNT,ZERO_COUNT_XI1(16)
    INTEGER(INTG) :: ZERO_COUNT_XI12(4),EDGE_NODE(16),TOTAL_ZERO_NODE,NODE_OFFSET_ELEM_XI12
    INTEGER(INTG) :: NUMBER_OF_NODES_LAYER
    LOGICAL::BASIS_APPEARED

    ENTERS("GeneratedMesh_RegularComponentNodesToUserNumbers",ERR,ERROR,*999)

    IF(SIZE(NODE_USER_NUMBERS)==SIZE(NODE_COMPONENT_NUMBERS)) THEN
      NODE_USER_NUMBERS=0
      IF(ASSOCIATED(GENERATED_MESH)) THEN
        IF(ASSOCIATED(GENERATED_MESH%REGULAR_MESH)) THEN
          NUM_BASES=SIZE(GENERATED_MESH%REGULAR_MESH%BASES)
          NUM_DIMS=GENERATED_MESH%REGULAR_MESH%MESH_DIMENSION
          BASES=>GENERATED_MESH%REGULAR_MESH%BASES
          NUMBER_OF_ELEMENTS_XI=1
          DO xi_idx=1,NUM_DIMS
            NUMBER_OF_ELEMENTS_XI(xi_idx)=GENERATED_MESH%REGULAR_MESH%NUMBER_OF_ELEMENTS_XI(xi_idx)
          ENDDO
        ELSE
        CALL FlagError("The regular mesh for this generated mesh is not associated.",ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
        ENDIF

        !Number of nodes in each xi direction
        numberOfNodesXic=1
        DO xiIdx=1,meshDimension
          numberOfNodesXic(xiIdx)=bases(basisIndex)%ptr%NUMBER_OF_NODES_XIC(xiIdx)
        ENDDO !xiIdx

        !Calculate current element indices and number
        remainderTemp=0;
        elementIdx=1;
        SELECT CASE(meshDimension)
        CASE(1)
          !Calculate xi1 element index
          elementIdx(1)=(nodeComponentNumbers(1)-1)/(numberOfNodesXic(1)-1)+1
          !Calculate element number
          elementNumber=elementIdx(1)
        CASE(2)
          !Calculate xi2 element index
          numberofNodesLayer=((numberOfNodesXic(1)-1)*numberOfElementsXi(1)+1)*(numberOfNodesXic(2)-1)
          elementIdx(2)=nodeComponentNumbers(1)/numberofNodesLayer+1
          remainderTemp=MOD(nodeComponentNumbers(1),numberofNodesLayer)
          !Calculate xi1 element index
          elementIdx(1)=(remainderTemp-1)/(numberOfNodesXic(1)-1)+1
          !Calculate element number
          elementNumber=(elementIdx(2)-1)*numberOfElementsXi(1)+elementIdx(1)
        CASE(3)
          !Calculate xi3 element index
          numberofNodesLayer=((numberOfNodesXic(1)-1)*numberOfElementsXi(1)+1)*((numberOfNodesXic(2)-1)* &
            & numberOfElementsXi(2)+1)*(numberOfNodesXic(3)-1)
          elementIdx(3)=nodeComponentNumbers(1)/numberofNodesLayer+1
          remainderTemp=MOD(nodeComponentNumbers(1),numberofNodesLayer)
          !Calculate xi2 element index
          numberofNodesLayer=((numberOfNodesXic(1)-1)*numberOfElementsXi(1)+1)*(numberOfNodesXic(2)-1)
          elementIdx(2)=remainderTemp/numberofNodesLayer+1
          remainderTemp=MOD(remainderTemp,numberofNodesLayer)
          !Calculate xi1 element index
          elementIdx(1)=(remainderTemp-1)/(numberOfNodesXic(1)-1)+1
          !Calculate element number
          elementNumber=(elementIdx(3)-1)*numberOfElementsXi(1)*numberOfElementsXi(2)+ &
            & (elementIdx(2)-1)*numberOfElementsXi(1)+elementIdx(1)
        CASE DEFAULT
          localError="A mesh dimension of "//TRIM(NumberToVString(meshDimension,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        
        !If not the first basis, check if previous basis have same interpolation order in each xi direction
        !sameBasis(3) is initialised to have zeros in all entries. If an interpolation scheme has been
        !found to have appeared in previous basis, then record the basis number in the corresponding
        !xi direction. e.g. First basis: bi-quadratic, Second basis: quadratic-cubic, then sameBasis(3)
        !for the second basis will be [1,0,0]
        sameBasis=0
        DO xiIdx=1,meshDimension
          DO basisIdx=1,basisIndex-1
            IF(bases(basisIndex)%ptr%NUMBER_OF_NODES_XIC(xiIdx)== &
              & bases(basisIdx)%ptr%NUMBER_OF_NODES_XIC(xiIdx)) THEN
              sameBasis(xiIdx)=basisIndex
            ENDIF
          ENDDO !basisIdx
        ENDDO !xiIdx
        !Check if the interpolation scheme has appeared in previous basis
        basisAppeared=.FALSE.
        IF(sameBasis(1)/=0) THEN
          SELECT CASE(meshDimension)
          CASE(1)
            basisAppeared=.TRUE.
          CASE(2)
            IF(sameBasis(1)==sameBasis(2)) basisAppeared=.TRUE.
          CASE(3)
            IF(sameBasis(1)==sameBasis(2).AND.sameBasis(1)==sameBasis(3)) basisAppeared=.TRUE.
          END SELECT
        ENDIF
        IF(basisIndex==1) THEN
          !If this is the first basis, don't do anything
          DO nodeIdx=1,SIZE(nodeComponentNumbers)
            nodeUserNumbers(nodeIdx)=nodeComponentNumbers(nodeIdx)
          ENDDO !nodeIdx
        ELSE IF(basisAppeared) THEN
          !If the basis has appeared before, reuse node user numbers
          DO nodeIdx=1,SIZE(nodeComponentNumbers)
            nodeUserNumbers(nodeIdx)=generatedMesh%mesh%topology(sameBasis(1))% &
              & ptr%elements%elements(elementNumber)%USER_ELEMENT_NODES(nodeIdx)
          ENDDO !nodeIdx
        ELSE
          !If the basis has never appeared exactly in previous basis

          !Find corner node user number from the first basis
          basisFirstComponent=>bases(1)%ptr
          DO localNodeIdx3=1,2
            DO localNodeIdx2=1,2
              DO localNodeIdx1=1,2
                nodeIdxCurrent=localNodeIdx1
                nodeIdxFirst=localNodeIdx1
                IF(localNodeIdx1==2) THEN
                  nodeIdxCurrent=numberOfNodesXic(1)
                  nodeIdxFirst=basisFirstComponent%NUMBER_OF_NODES_XIC(1)
                ENDIF
                IF(meshDimension>1.AND.localNodeIdx2==2) THEN
                  nodeIdxCurrent=nodeIdxCurrent+(numberOfNodesXic(2)-1)*numberOfNodesXic(1)
                  nodeIdxFirst=nodeIdxFirst+(basisFirstComponent%NUMBER_OF_NODES_XIC(2)-1)* &
                    & basisFirstComponent%NUMBER_OF_NODES_XIC(1)
                ENDIF
                IF(meshDimension>2.AND.localNodeIdx3==2) THEN
                  nodeIdxCurrent=nodeIdxCurrent+numberOfNodesXic(1)*numberOfNodesXic(2)*(numberOfNodesXic(3)-1)
                  nodeIdxFirst=nodeIdxFirst+basisFirstComponent%NUMBER_OF_NODES_XIC(1)* &
                    & basisFirstComponent%NUMBER_OF_NODES_XIC(2)*(basisFirstComponent%NUMBER_OF_NODES_XIC(3)-1)
                ENDIF
                nodeUserNumbers(nodeIdxCurrent)=generatedMesh%mesh%topology(1)%ptr%elements%elements(elementNumber)% &
                  & GLOBAL_ELEMENT_NODES(nodeIdxFirst)
              ENDDO !localNodeIdx1
            ENDDO !localNodeIdx2
          ENDDO !localNodeIdx3

          !Find edge node user number from previous basis
          IF(sameBasis(1)/=0.AND.meshDimension>1) THEN !Do not consider 1D since it's a complete new basis
            basisPrevious=>bases(sameBasis(1))%ptr
            DO localNodeIdx3=1,2
              DO localNodeIdx2=1,2
                DO localNodeIdx1=2,numberOfNodesXic(1)-1
                  nodeIdxCurrent=localNodeIdx1
                  nodeIdxPrevious=localNodeIdx1
                  IF(localNodeIdx2==2) THEN
                    nodeIdxCurrent=nodeIdxCurrent+(numberOfNodesXic(2)-1)*numberOfNodesXic(1)
                    nodeIdxPrevious=nodeIdxPrevious+(basisPrevious%NUMBER_OF_NODES_XIC(2)-1)*basisPrevious%NUMBER_OF_NODES_XIC(1)
                  ENDIF
                  IF(meshDimension>2 .AND. localNodeIdx3==2) THEN
                    nodeIdxCurrent=nodeIdxCurrent+numberOfNodesXic(1)*numberOfNodesXic(2)* &
                      & (numberOfNodesXic(3)-1)
                    nodeIdxPrevious=nodeIdxPrevious+basisPrevious%NUMBER_OF_NODES_XIC(1)*basisPrevious% &
                      & NUMBER_OF_NODES_XIC(2)*(basisPrevious%NUMBER_OF_NODES_XIC(3)-1)
                  ENDIF
                  nodeUserNumbers(nodeIdxCurrent)=generatedMesh%mesh%topology(sameBasis(1))% &
                    & ptr%elements%elements(elementNumber)%GLOBAL_ELEMENT_NODES(nodeIdxPrevious)
                ENDDO !localNodeIdx1
              ENDDO !localNodeIdx2
            ENDDO !localNodeIdx3
          ENDIF
          IF(sameBasis(2)/=0) THEN
            basisPrevious=>bases(sameBasis(2))%ptr
            DO localNodeIdx3=1,2
              DO localNodeIdx2=2,numberOfNodesXic(2)-1
                DO localNodeIdx1=1,2
                  IF(localNodeIdx1==1) THEN
                    nodeIdxCurrent=localNodeIdx1+(localNodeIdx2-1)*numberOfNodesXic(1)
                    nodeIdxPrevious=localNodeIdx1+(localNodeIdx2-1)*basisPrevious%NUMBER_OF_NODES_XIC(1)
                  ELSE
                    nodeIdxCurrent=localNodeIdx2*numberOfNodesXic(1)
                    nodeIdxPrevious=localNodeIdx2*basisPrevious%NUMBER_OF_NODES_XIC(1)
                  ENDIF
                  IF(meshDimension>2 .AND. localNodeIdx3==2) THEN
                    nodeIdxCurrent=nodeIdxCurrent+numberOfNodesXic(1)*numberOfNodesXic(2)* &
                      & (numberOfNodesXic(3)-1)
                    nodeIdxPrevious=nodeIdxPrevious+basisPrevious%NUMBER_OF_NODES_XIC(1)*basisPrevious% &
                      & NUMBER_OF_NODES_XIC(2)*(basisPrevious%NUMBER_OF_NODES_XIC(3)-1)
                  ENDIF
                  nodeUserNumbers(nodeIdxCurrent)=generatedMesh%mesh%topology(sameBasis(2))% &
                    & ptr%elements%elements(elementNumber)%GLOBAL_ELEMENT_NODES(nodeIdxPrevious)
                ENDDO !localNodeIdx1
              ENDDO !localNodeIdx2
            ENDDO !localNodeIdx3
          ENDIF
          IF(sameBasis(3)/=0) THEN !Must be 3D
            basisPrevious=>bases(sameBasis(3))%ptr
            nodeIdxCurrent=0
            nodeIdxPrevious=0
            DO localNodeIdx3=2,numberOfNodesXic(3)-1
              DO localNodeIdx2=1,2
                IF(localNodeIdx2==2) THEN
                  nodeIdxCurrent=(numberOfNodesXic(2)-1)*numberOfNodesXic(1)+numberOfNodesXic(1)*numberOfNodesXic(2)* &
                    & (numberOfNodesXic(3)-1)
                  nodeIdxPrevious=(basisPrevious%NUMBER_OF_NODES_XIC(1)-1)*basisPrevious%NUMBER_OF_NODES_XIC(1)+ &
                    & basisPrevious%NUMBER_OF_NODES_XIC(1)*basisPrevious%NUMBER_OF_NODES_XIC(2)* &
                    & (basisPrevious%NUMBER_OF_NODES_XIC(3)-1)
                ENDIF
                DO localNodeIdx1=1,2
                  IF(localNodeIdx1==1) THEN
                    nodeIdxCurrent=1+nodeIdxCurrent
                    nodeIdxPrevious=1+nodeIdxPrevious
                  ELSE
                    nodeIdxCurrent=numberOfNodesXic(1)+nodeIdxCurrent
                    nodeIdxPrevious=basisPrevious%NUMBER_OF_NODES_XIC(1)+nodeIdxPrevious
                  ENDIF
                  nodeUserNumbers(nodeIdxCurrent)=generatedMesh%mesh%topology(sameBasis(3))% &
                    & ptr%elements%elements(elementNumber)%GLOBAL_ELEMENT_NODES(nodeIdxPrevious)
                ENDDO !localNodeIdx1
              ENDDO !localNodeIdx2
            ENDDO !localNodeIdx3
          ENDIF
          !The following code would only be executed if 3D (automatically satisfied, don't need to check,
          !since there must be at least 1 direction that has different interpolation scheme, if two direction
          ! has the same interpolation that has appeared before, then interpolation for the last direction
          ! must be different) and has same basis in 2 xi direction
          !i.e. find user node numbers for face nodes
          IF(sameBasis(1)==sameBasis(2).AND.sameBasis(1)/=0) THEN
            basisPrevious=>bases(sameBasis(1))%ptr
            DO localNodeIdx3=1,2
              DO localNodeIdx2=2,numberOfNodesXic(2)-1
                DO localNodeIdx1=2,numberOfNodesXic(1)-1
                  nodeIdxCurrent=localNodeIdx1+(localNodeIdx2-1)*numberOfNodesXic(1)
                  nodeIdxPrevious=localNodeIdx1+(localNodeIdx2-1)*basisPrevious%NUMBER_OF_NODES_XIC(1)
                  IF(localNodeIdx3==2) THEN
                    nodeIdxCurrent=nodeIdxCurrent+numberOfNodesXic(1)*numberOfNodesXic(2)* &
                      & (numberOfNodesXic(3)-1)
                    nodeIdxPrevious=nodeIdxPrevious+basisPrevious%NUMBER_OF_NODES_XIC(1)*basisPrevious% &
                      & NUMBER_OF_NODES_XIC(2)*(basisPrevious%NUMBER_OF_NODES_XIC(3)-1)
                  ENDIF
                  nodeUserNumbers(nodeIdxCurrent)=generatedMesh%mesh%topology(sameBasis(1))% &
                    & ptr%elements%elements(elementNumber)%GLOBAL_ELEMENT_NODES(nodeIdxPrevious)
                ENDDO !localNodeIdx1
              ENDDO !localNodeIdx2
            ENDDO !localNodeIdx3
          ELSE IF(sameBasis(1)==sameBasis(3).AND.sameBasis(1)/=0) THEN
            basisPrevious=>bases(sameBasis(1))%ptr
            nodeIdxCurrent=0
            nodeIdxPrevious=0
            DO localNodeIdx3=2,numberOfNodesXic(3)-1
              DO localNodeIdx2=1,2
                IF(localNodeIdx2==2) THEN
                  nodeIdxCurrent=(numberOfNodesXic(2)-1)*numberOfNodesXic(1)+numberOfNodesXic(1)*numberOfNodesXic(2)* &
                    & (localNodeIdx3-1)
                  nodeIdxPrevious=(basisPrevious%NUMBER_OF_NODES_XIC(2)-1)*basisPrevious%NUMBER_OF_NODES_XIC(1)+ &
                    & basisPrevious%NUMBER_OF_NODES_XIC(1)*basisPrevious%NUMBER_OF_NODES_XIC(2)*(localNodeIdx3-1)
                ENDIF
                DO localNodeIdx1=2,numberOfNodesXic(1)-1
                  nodeIdxCurrent=localNodeIdx1+nodeIdxCurrent
                  nodeIdxPrevious=localNodeIdx1+nodeIdxPrevious
                  nodeUserNumbers(nodeIdxCurrent)=generatedMesh%mesh%topology(sameBasis(1))% &
                    & ptr%elements%elements(elementNumber)%GLOBAL_ELEMENT_NODES(nodeIdxPrevious)
                ENDDO !localnodeIdx1
              ENDDO !localNodeIdx2
            ENDDO !localNodeIdx3
          ELSE IF(sameBasis(2)==sameBasis(3).AND.sameBasis(2)/=0) THEN
            basisPrevious=>bases(sameBasis(2))%ptr
            DO localNodeIdx3=2,numberOfNodesXic(3)-1
              DO localNodeIdx2=2,numberOfNodesXic(2)-1
                DO localNodeIdx1=1,2
                  IF(localNodeIdx1==1) THEN
                    nodeIdxCurrent=1+(localNodeIdx2-1)*numberOfNodesXic(1)+numberOfNodesXic(1)* &
                      & numberOfNodesXic(2)*(localNodeIdx3-1)
                    nodeIdxPrevious=1+(localNodeIdx2-1)*basisPrevious%NUMBER_OF_NODES_XIC(1)+basisPrevious%NUMBER_OF_NODES_XIC(1)* &
                      & basisPrevious%NUMBER_OF_NODES_XIC(2)*(localNodeIdx3-1)
                  ELSE
                    nodeIdxCurrent=localNodeIdx2*numberOfNodesXic(1)+numberOfNodesXic(1)* &
                      & numberOfNodesXic(2)*(localNodeIdx3-1)
                    nodeIdxPrevious=localNodeIdx2*basisPrevious%NUMBER_OF_NODES_XIC(1)+basisPrevious%NUMBER_OF_NODES_XIC(1)* &
                      & basisPrevious%NUMBER_OF_NODES_XIC(2)*(localNodeIdx3-1)
                  ENDIF
                  nodeUserNumbers(nodeIdxCurrent)=generatedMesh%mesh%topology(sameBasis(2))% &
                    & ptr%elements%elements(elementNumber)%GLOBAL_ELEMENT_NODES(nodeIdxPrevious)
                ENDDO
              ENDDO
            ENDDO
          ENDIF

          !Find the largest node user number in the previous basis
          nodeOffsetLastBasis=0
          lastElementNumber=generatedMesh%mesh%topology(1)%ptr%elements%NUMBER_OF_ELEMENTS !The mesh has the same topology regardless of mesh components
          DO basisIdx=1,basisIndex-1
            numberLocalNodes=SIZE(generatedMesh%MESH%TOPOLOGY(basisIdx)%ptr%elements%elements(lastElementNumber)% &
              & GLOBAL_ELEMENT_NODES,1)
            DO localNodeIdx=1,numberLocalNodes
              IF (generatedMesh%mesh%topology(basisIdx)%ptr%elements%elements(lastElementNumber)% &
                & GLOBAL_ELEMENT_NODES(localNodeIdx)>nodeOffsetLastBasis) THEN
                nodeOffsetLastBasis=generatedMesh%mesh%topology(basisIdx)%ptr%elements%elements(lastElementNumber)% &
                  & GLOBAL_ELEMENT_NODES(localNodeIdx)
              ENDIF
            ENDDO !localNodeIdx
          ENDDO !basisIdx

          !Calculate number of zeros nodes in different dimensions
          indexCount=1
          zeroCountXi1=0
          zeroCountXi12=0
          totalNumberZeroNodes=0
          edgeNode=0
          DO localNodeIdx3=1,numberOfNodesXic(3)
            DO localNodeIdx2=1,numberOfNodesXic(2)
              nodeCount=0
              DO localNodeIdx1=1,numberOfNodesXic(1)
                nodeIdx=(localNodeIdx3-1)*numberOfNodesXic(1)*numberOfNodesXic(2)+(localNodeIdx2-1)*numberOfNodesXic(1)+ &
                  & localNodeIdx1
                IF(nodeUserNumbers(nodeIdx)==0) THEN
                  nodeCount=nodeCount+1
                  totalNumberZeroNodes=totalNumberZeroNodes+1 !Total number of zeros in an element
                ENDIF
              ENDDO !localNodeIdx1
              zeroCountXi1(indexCount)=nodeCount !Total number of zero summed up across xi1 direction.
              IF(nodeCount==numberOfNodesXic(1)) edgeNode(indexCount)=1 !Shared edge node (with zero value) in xi1 direction (1 number for each node in xi2 direction)
              zeroCountXi12(localNodeIdx3)=zeroCountXi12(localNodeIdx3)+zeroCountXi1(indexCount) !Total number of zero summed on xi1-xi2 faces
              indexCount=indexCount+1
            ENDDO !localNodeIdx2
          ENDDO !localNodeIdx3
          
          !Calculate how many zero nodes has occurred in previous elements
          nodeOffsetElement=0
          IF(meshDimension==2.AND.elementIdx(2)/=1) THEN !Zero nodes occurred in the previous rows of elements
            offsetUnit=totalNumberZeroNodes-zeroCountXi1(1)-SUM(edgeNode(1:numberOfNodesXic(2)))+edgeNode(indexCount)
            !This is number of zero nodes in the elements before the current row of elements
            nodeOffsetElement=(elementIdx(2)-1)*numberOfElementsXi(1)*offsetUnit+(elementIdx(2)-1)* &
              & SUM(edgeNode(2:numberOfNodesXic(2)-1))
          ELSE IF(meshDimension==3.AND.elementIdx(3)/=1) THEN !Zero nodes occurred in the previous layer of elements
            nodeOffsetXi3Accumulated=0
            DO localNodeIdx3=1,numberOfNodesXic(3)-1
              offsetUnit=zeroCountXi12(localNodeIdx3)-zeroCountXi1((localNodeIdx3-1)*numberOfNodesXic(2)+1)- &
                & SUM(edgeNode((localNodeIdx3-1)*numberOfNodesXic(2)+1:localNodeIdx3*numberOfNodesXic(2)))+ &
                & edgeNode((localNodeIdx3-1)*numberOfNodesXic(2)+1)
              nodeOffsetXi3Accumulated=nodeOffsetXi3Accumulated+offsetUnit*numberOfElementsXi(1)*numberOfElementsXi(2)+ &
                & (numberOfElementsXi(1)-1)*(zeroCountXi1((localNodeIdx3-1)*numberOfNodesXic(2)+1)- &
                & edgeNode((localNodeIdx3-1)*numberOfNodesXic(2)+1))+zeroCountXi1((localNodeIdx3-1)*numberOfNodesXic(2)+1)+ &
                & SUM(edgeNode((localNodeIdx3-1)*numberOfNodesXic(2)+2:localNodeIdx3*numberOfNodesXic(2)))* &
                & numberOfElementsXi(2)
            ENDDO !localNodeIdx3
            nodeOffsetElement=(elementIdx(3)-1)*nodeOffsetXi3Accumulated
          ENDIF
          
          !Compute other nodes which haven't appeared in previous basis
          indexCount=1
          nodeOffsetElementXi12=0
          nodeOffsetXi2=0 !Number of zero nodes in the current row
          nodeOffsetXi3Accumulated=0 !Number of zero nodes in the layers in xi3 direction (localNodeIdx3)
          DO localNodeIdx3=1,numberOfNodesXic(3)
            nodeOffsetXi2Accumulated=0 !Number of zero nodes in the previous rows
            offsetUnit=zeroCountXi12(localNodeIdx3)-zeroCountXi1((localNodeIdx3-1)*numberOfNodesXic(2)+1)- &
              & SUM(edgeNode((localNodeIdx3-1)*numberOfNodesXic(2)+1:localNodeIdx3*numberOfNodesXic(2)))+ &
              & edgeNode((localNodeIdx3-1)*numberOfNodesXic(2)+1)
            IF(elementIdx(2)/=1.AND.meshDimension==3) THEN
              nodeOffsetElementXi12=offsetUnit*(elementIdx(2)-1)*numberOfElementsXi(1)+ &
                & (elementIdx(2)-1)*SUM(edgeNode((localNodeIdx3-1)*numberOfNodesXic(2)+2:localNodeIdx3*numberOfNodesXic(2)))
            ENDIF
            DO localNodeIdx2=1,numberOfNodesXic(2)
              nodeOffsetXi2=(zeroCountXi1(indexCount)-edgeNode(indexCount))*(elementIdx(1)-1)
              nodeOffset=nodeOffsetLastBasis+nodeOffsetElement+nodeOffsetXi3Accumulated+ &
                & nodeOffsetElementXi12+nodeOffsetXi2Accumulated+nodeOffsetXi2
              DO localNodeIdx1=1,numberOfNodesXic(1)
                !Local node index in the current element
                nodeIdx=(localNodeIdx3-1)*numberOfNodesXic(1)*numberOfNodesXic(2)+(localNodeIdx2-1)* &
                  & numberOfNodesXic(1)+localNodeIdx1
                IF(nodeUserNumbers(nodeIdx)==0) THEN
                  !This is for 2D case
                  nodeOffset=nodeOffset+1
                  nodeUserNumbers(nodeIdx)=nodeOffset
                ENDIF
              ENDDO !localNodeIdx1
              nodeOffsetXi2Accumulated=nodeOffsetXi2Accumulated+(zeroCountXi1(indexCount)-edgeNode(indexCount))* &
                & numberOfElementsXi(1)+edgeNode(indexCount)
              indexCount=indexCount+1
            ENDDO !localNodeIdx2
            IF(meshDimension==3) THEN
              nodeOffsetXi3Accumulated=nodeOffsetXi3Accumulated+offsetUnit*numberOfElementsXi(1)*numberOfElementsXi(2)+ &
                & (numberOfElementsXi(1)-1)*(zeroCountXi1((localNodeIdx3-1)*numberOfNodesXic(2)+1)- &
                & edgeNode((localNodeIdx3-1)*numberOfNodesXic(2)+1))+zeroCountXi1((localNodeIdx3-1)*numberOfNodesXic(2)+1)+ &
                & SUM(edgeNode((localNodeIdx3-1)*numberOfNodesXic(2)+2:localNodeIdx3*numberOfNodesXic(2)))* &
                & numberOfElementsXi(2)
            ENDIF
          ENDDO !localNodeIdx3
        ENDIF
      ELSE
<<<<<<< HEAD
        CALL FlagError("Generated mesh is not associated",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("nodeComponentNumbers and nodeUserNumbers arrays have different sizes.",err,error,*999)
    ENDIF
    CALL Exits("GeneratedMeshRegularComponentNodesToUserNumbers")
    RETURN
999 CALL Errors("GeneratedMeshRegularComponentNodesToUserNumbers",err,error)
    CALL Exits("GeneratedMeshRegularComponentNodesToUserNumbers")
    RETURN 1
  END SUBROUTINE GeneratedMeshRegularComponentNodesToUserNumbers
=======
        CALL FlagError("Generated mesh is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("NODE_COMPONENT_NUMBERS and NODE_USER_NUMBERS arrays have different sizes.",ERR,ERROR,*999)
    ENDIF
    EXITS("GeneratedMesh_RegularComponentNodesToUserNumbers")
    RETURN
999 ERRORS("GeneratedMesh_RegularComponentNodesToUserNumbers",ERR,ERROR)
    EXITS("GeneratedMesh_RegularComponentNodesToUserNumbers")
    RETURN 1
    
  END SUBROUTINE GeneratedMesh_RegularComponentNodesToUserNumbers
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b

  !
  !================================================================================================================================
  !

  !>Retrieve the user node number for a component number in a regular generated mesh
  !>This routine only works for Lagrange/Hermite elements
<<<<<<< HEAD
  SUBROUTINE GeneratedMeshRegularComponentNodeToUserNumber(generatedMesh,basisIndex,nodeComponentNumber,nodeUserNumber,err,error,*)
=======
  SUBROUTINE GeneratedMesh_RegularComponentNodeToUserNumber(GENERATED_MESH,BASIS_INDEX, &
      & NODE_COMPONENT_NUMBER,NODE_USER_NUMBER,ERR,ERROR,*)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh  !<A pointer to the generated mesh object
    INTEGER(INTG),INTENT(IN) :: basisIndex  !<The number of the basis being used
    INTEGER(INTG),INTENT(IN) :: nodeComponentNumber  !<The node numbers for this component basis
    INTEGER(INTG),INTENT(OUT) :: nodeUserNumber  !<On return, the corresponding user numbers
    INTEGER(INTG) :: ERR  !<The error code
    TYPE(VARYING_STRING) :: ERROR  !<The error string

    !Local variables
<<<<<<< HEAD
    TYPE(BASIS_PTR_TYPE), POINTER :: bases(:)
    INTEGER(INTG) :: meshDimension,elementNumber,localNodeNumber,numberofNodesLayer,xiIdx,xicIdx
    INTEGER(INTG) :: elementIdx(3),nodeIdx(3),numberOfNodesXic(3),numberOfElementsXi(3),reminderTemp

    CALL Enters("GeneratedMeshRegularComponentNodeToUserNumber",err,error,*999)

    IF(ASSOCIATED(generatedMesh)) THEN
      IF(ASSOCIATED(generatedMesh%regularMesh)) THEN
        meshDimension=generatedMesh%regularMesh%meshDimension
        bases=>generatedMesh%regularMesh%bases
        !Compute the number of elements in each xi direction
        numberOfElementsXi=1
        DO xiIdx=1,meshDimension
          numberOfElementsXi(xiIdx)=generatedMesh%regularMesh%numberOfElementsXi(xiIdx)
        ENDDO !xiIdx
        !Compute th number of nodes in each xic direction
        numberOfNodesXic=1
        DO xicIdx=1,meshDimension
          numberOfNodesXic(xicIdx)=bases(basisIndex)%ptr%NUMBER_OF_NODES_XIC(xicIdx)
        ENDDO !xicIdx
      ELSE
        CALL FlagError("The regular mesh for this generated mesh is not associated.",err,error,*999)
=======
    TYPE(BASIS_PTR_TYPE), POINTER :: BASES(:)
    INTEGER(INTG) :: NUM_BASES,NUM_DIMS,ELEMENT_NO,LOCAL_NODE_NO,NUMBER_OF_NODES_LAYER,xi_idx
    INTEGER(INTG) :: ELEM_IDX(3),NODE_IDX(3),NUMBER_OF_NODES_XIC(3),NUMBER_OF_ELEMENTS_XI(3),REMINDER_TEMP

    ENTERS("GeneratedMesh_RegularComponentNodeToUserNumber",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      IF(ASSOCIATED(GENERATED_MESH%REGULAR_MESH)) THEN
        NUM_BASES=SIZE(GENERATED_MESH%REGULAR_MESH%BASES)
        NUM_DIMS=GENERATED_MESH%REGULAR_MESH%MESH_DIMENSION
        BASES=>GENERATED_MESH%REGULAR_MESH%BASES
        NUMBER_OF_ELEMENTS_XI=1
        DO xi_idx=1,NUM_DIMS
          NUMBER_OF_ELEMENTS_XI(xi_idx)=GENERATED_MESH%REGULAR_MESH%NUMBER_OF_ELEMENTS_XI(xi_idx)
        ENDDO
        !Number of nodes in each xi direction
        NUMBER_OF_NODES_XIC=1
        DO xi_idx=1,NUM_DIMS
          NUMBER_OF_NODES_XIC(xi_idx)=BASES(BASIS_INDEX)%PTR%NUMBER_OF_NODES_XIC(xi_idx)
        ENDDO
      ELSE
        CALL FlagError("The regular mesh for this generated mesh is not associated.",ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
      ENDIF

      !Calculate current element/node indices/number
      reminderTemp=0;
      elementIdx=1;
      nodeIdx=1;
      SELECT CASE(meshDimension)
      CASE(1)
        !Calculate xi1 element index
        elementIdx(1)=(nodeComponentNumber-1)/(numberOfNodesXic(1)-1)+1
        nodeIdx(1)=MOD(nodeComponentNumber-1,numberOfNodesXic(1)-1)+1
        !If it's the last node in the line
        IF(elementIdx(1)>numberOfElementsXi(1)) THEN
          elementIdx(1)=elementIdx(1)-1
          nodeIdx(1)=numberOfNodesXic(1)
        ENDIF
        !Calculate element number
        elementNumber=elementIdx(1)
        localNodeNumber=nodeIdx(1)
      CASE(2)
        !Calculate xi2 element index
        numberOfNodesLayer=((numberOfNodesXic(1)-1)*numberOfElementsXi(1)+1)*(numberOfNodesXic(2)-1)
        elementIdx(2)=(nodeComponentNumber-1)/numberofNodesLayer+1
        reminderTemp=MOD(nodeComponentNumber-1,numberofNodesLayer)
        numberOfNodesLayer=((numberOfNodesXic(1)-1)*numberOfElementsXi(1)+1)
        nodeIdx(2)=reminderTemp/numberOfNodesLayer+1
        !If it's the last line of nodes in the line
        IF(elementIdx(2)>numberOfElementsXi(2)) THEN
          elementIdx(2)=elementIdx(2)-1
          nodeIdx(2)=numberOfNodesXic(2)
        ENDIF
        !Calculate xi1 element index
        reminderTemp=MOD(reminderTemp,numberOfNodesLayer)
        elementIdx(1)=reminderTemp/(numberOfNodesXic(1)-1)+1
        nodeIdx(1)=MOD(reminderTemp,numberOfNodesXic(1)-1)+1
        !If it's the last node in the line
        IF(elementIdx(1)>numberOfElementsXi(1)) THEN
          elementIdx(1)=elementIdx(1)-1
          nodeIdx(1)=numberOfNodesXic(1)
        ENDIF
        !Calculate element number
        elementNumber=(elementIdx(2)-1)*numberOfElementsXi(1)+elementIdx(1)
        localNodeNumber=(nodeIdx(2)-1)*numberOfNodesXic(1)+nodeIdx(1)
      CASE(3)
        !Calculate xi3 element index
        numberofNodesLayer=((numberOfNodesXic(1)-1)*numberOfElementsXi(1)+1)*((numberOfNodesXic(2)-1)* &
          & numberOfElementsXi(2)+1)*(numberOfNodesXic(3)-1) !Multiple planes of nodes
        elementIdx(3)=(nodeComponentNumber-1)/numberOfNodesLayer+1
        reminderTemp=MOD(nodeComponentNumber-1,numberOfNodesLayer) !Multiple planes of nodes
        numberOfNodesLayer=((numberOfNodesXic(1)-1)*numberOfElementsXi(1)+1)*((numberOfNodesXic(2)-1)* &
          & numberOfElementsXi(2)+1) !One plane of nodes
        nodeIdx(3)=reminderTemp/numberofNodesLayer+1
        !If it's the last node in the line
        IF(elementIdx(3)>numberOfElementsXi(3)) THEN
          elementIdx(3)=elementIdx(3)-1
          nodeIdx(3)=numberOfNodesXic(3)
        ENDIF
        reminderTemp=MOD(reminderTemp,numberOfNodesLayer) !One plane of nodes
        !Calculate xi2 element index
        numberOfNodesLayer=((numberOfNodesXic(1)-1)*numberOfElementsXi(1)+1)*(numberOfNodesXic(2)-1) !Multiple lines of nodes
        elementIdx(2)=reminderTemp/numberOfNodesLayer+1
        reminderTemp=MOD(reminderTemp,numberOfNodesLayer) !Multiple lines of nodes
        numberOfNodesLayer=(numberOfNodesXic(1)-1)*numberOfElementsXi(1)+1 !One line of nodes
        nodeIdx(2)=reminderTemp/numberOfNodesLayer+1
        reminderTemp=MOD(reminderTemp,numberOfNodesLayer) !One line of nodes
        !If it's the last node in the line
        IF(elementIdx(2)>numberOfElementsXi(2)) THEN
          elementIdx(2)=elementIdx(2)-1
          nodeIdx(2)=numberOfNodesXic(2)
        ENDIF
        !Calculate xi1 element index
        elementIdx(1)=reminderTemp/(numberOfNodesXic(1)-1)+1
        nodeIdx(1)=MOD(reminderTemp,numberOfNodesXic(1)-1)+1
        IF(elementIdx(1)>numberOfElementsXi(1)) THEN
          elementIdx(1)=elementIdx(1)-1
          nodeIdx(1)=numberOfNodesXic(1)
        ENDIF
        !Calculate element number
        elementNumber=(elementIdx(3)-1)*numberOfElementsXi(1)*numberOfElementsXi(2)+(elementIdx(2)-1)* &
          & numberOfElementsXi(1)+elementIdx(1)
        localNodeNumber=(nodeIdx(3)-1)*numberOfNodesXic(1)*numberOfNodesXic(2)+(nodeIdx(2)-1)*numberOfNodesXic(1)+ &
          & nodeIdx(1)
      END SELECT
      !Retrieve node user number
      IF(ASSOCIATED(generatedMesh%mesh)) THEN
        nodeUserNumber=generatedMesh%mesh%topology(basisIndex)%ptr%elements%elements(elementNumber)% &
          & USER_ELEMENT_NODES(localNodeNumber)
      ELSE
<<<<<<< HEAD
        CALL FlagError("The mesh for this generated mesh is not associated.",err,error,*999)
      ENDIF

    ELSE
      CALL FlagError("Generated mesh is not associated",err,error,*999)
    ENDIF

    CALL Exits("GeneratedMeshRegularComponentNodeToUserNumber")
    RETURN
999 CALL Errors("GeneratedMeshRegularComponentNodeToUserNumber",err,error)
    CALL Exits("GeneratedMeshRegularComponentNodeToUserNumber")
    RETURN 1
  END SUBROUTINE GeneratedMeshRegularComponentNodeToUserNumber
=======
        CALL FlagError("The mesh for this generated mesh is not associated.",ERR,ERROR,*999)
      ENDIF

    ELSE
      CALL FlagError("Generated mesh is not associated",ERR,ERROR,*999)
    ENDIF

    EXITS("GeneratedMesh_RegularComponentNodeToUserNumber")
    RETURN
999 ERRORS("GeneratedMesh_RegularComponentNodeToUserNumber",ERR,ERROR)
    EXITS("GeneratedMesh_RegularComponentNodeToUserNumber")
    RETURN 1
    
  END SUBROUTINE GeneratedMesh_RegularComponentNodeToUserNumber
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b

  !
  !================================================================================================================================
  !

  !>Calculates the user node number for a node numbered using one basis.
  !>This is currently only used for cylinder meshes, other mesh types don't require this.
<<<<<<< HEAD
  FUNCTION GeneratedMeshUserNumberToComponentNode(generatedMesh,basisIndex,nodeUserNumber,err,error)
    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh object
    INTEGER(INTG),INTENT(IN) :: basisIndex !<The number of the basis being used
    INTEGER(INTG),INTENT(IN) :: nodeUserNumber !<The corresponding user node number
    INTEGER(INTG) :: err !<The error code
    TYPE(VARYING_STRING) :: error !<The error string
    !Function variable
    INTEGER(INTG) :: GeneratedMeshUserNumberToComponentNode !<On return, the node number for this component basis
    !Local variables
    INTEGER(INTG) :: numberOfBases,meshDimension,basisIdx,xiIdx,remainder,temporaryTerm,numberCornerNodes,nodeOffset, &
      & basisNumberOfNodes,position(3),cornerNodeFactor(3),basisElementFactor(3),numberPreviousCorners,numberOfElementsXi(3)
    LOGICAL :: finishedCount,offEdge
    TYPE(BASIS_TYPE), POINTER :: basis
    TYPE(BASIS_PTR_TYPE), POINTER :: bases(:)
    TYPE(VARYING_STRING) :: localError

    NULLIFY(basis)
    NULLIFY(bases)
    
    CALL Enters("GeneratedMeshUserNumberToComponentNode",err,error,*999)

    numberCornerNodes=1
    remainder=nodeUserNumber-1 !use zero based numbering
    position=0

    IF(ASSOCIATED(generatedMesh)) THEN
=======
  FUNCTION USER_NUMBER_TO_COMPONENT_NODE(GENERATED_MESH,BASIS_INDEX,NODE_USER_NUMBER,ERR,ERROR)
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH        !<A pointer to the generated mesh object
    INTEGER(INTG),INTENT(IN) :: BASIS_INDEX                     !<The number of the basis being used
    INTEGER(INTG),INTENT(IN) :: NODE_USER_NUMBER                !<The corresponding user node number
    INTEGER(INTG) :: ERR          !<The error code
    TYPE(VARYING_STRING) :: ERROR !<The error string
    !function variable
    INTEGER(INTG) :: USER_NUMBER_TO_COMPONENT_NODE !<On return, the node number for this component basis
    !local variables
    INTEGER(INTG) :: NUM_BASES,NUM_DIMS,basis_idx,ni,REMAINDER,TEMP_TERM,NUM_CORNER_NODES,NODE_OFFSET,BASIS_NUM_NODES
    INTEGER(INTG) :: POS(3),CORNER_NODE_FACTOR(3),BASIS_ELEMENT_FACTOR(3),NUM_PREVIOUS_CORNERS
    INTEGER(INTG), POINTER :: NUMBER_OF_ELEMENTS_XI(:)
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(BASIS_PTR_TYPE), POINTER :: BASES(:)
    LOGICAL :: FINISHED_COUNT,OFF_EDGE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("USER_NUMBER_TO_COMPONENT_NODE",ERR,ERROR,*999)

    NULLIFY(BASIS)
    NULLIFY(BASES)
    NUM_CORNER_NODES=1
    REMAINDER=NODE_USER_NUMBER-1 !use zero based numbering
    POS=0

    IF(ASSOCIATED(GENERATED_MESH)) THEN
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
      !Only cylinder mesh type uses this now, although it was previously used by regular
      !meshes so some things relate to that.
      SELECT CASE(generatedMesh%generatedType)
      CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
<<<<<<< HEAD
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(GENERATED_MESH_POLAR_MESH_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
=======
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE(GENERATED_MESH_POLAR_MESH_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
      CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
        IF(ASSOCIATED(generatedMesh%cylinderMesh)) THEN
          numberOfBases=SIZE(generatedMesh%cylinderMesh%bases)
          meshDimension=generatedMesh%cylinderMesh%meshDimension
          bases=>generatedMesh%cylinderMesh%bases
          numberOfElementsXi(1:meshDimension)=generatedMesh%cylinderMesh%numberOfElementsXi(1:meshDimension)
        ELSE
<<<<<<< HEAD
          CALL FlagError("The cylinder mesh for this generated mesh is not associated.",err,error,*999)
        ENDIF
      CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="The generated mesh generated type of "//TRIM(NumberToVString(generatedMesh%generatedType,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
=======
          CALL FlagError("The cylinder mesh for this generated mesh is not associated.",ERR,ERROR,*999)
        ENDIF
      CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="The generated mesh generated type of "// &
            & TRIM(NUMBER_TO_VSTRING(GENERATED_MESH%GENERATED_TYPE,"*",ERR,ERROR))//" is invalid."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
      END SELECT
      IF(basisIndex<=numberOfBases) THEN
        IF(numberOfBases==1) THEN
          !If is the only basis, don't do anything
          GeneratedMeshUserNumberToComponentNode=nodeUserNumber
        ELSE
          temporaryTerm=1
          numberCornerNodes=1
          DO xiIdx=1,meshDimension
            numberCornerNodes=numberCornerNodes*(numberOfElementsXi(xiIdx)+1)
            cornerNodeFactor(xiIdx)=1
            IF(xiIdx>1) THEN
              temporaryTerm=temporaryTerm*(numberOfElementsXi(xiIdx-1)+1)
              cornerNodeFactor(xiIdx)=cornerNodeFactor(xiIdx)*temporaryTerm
            ENDIF
          ENDDO
          !Adjust for other mesh types
          IF(generatedMesh%generatedType==GENERATED_MESH_CYLINDER_MESH_TYPE) THEN
            cornerNodeFactor(3)=cornerNodeFactor(3)-numberOfElementsXi(1)-1
            numberCornerNodes=numberCornerNodes-(numberOfElementsXi(1)+1)*(numberOfElementsXi(3)+1)
          ENDIF
          nodeOffset=numberCornerNodes
          DO basisIdx=1,basisIndex-1
            basis=>bases(basisIdx)%ptr
            basisNumberOfNodes=1
            DO xiIdx=1,meshDimension
              basisNumberOfNodes=basisNumberOfNodes*(numberOfElementsXi(xiIdx)*(basis%NUMBER_OF_NODES_XIC(xiIdx)-1)+1)
            ENDDO
            !Adjust for other mesh types
            IF(generatedMesh%generatedType==GENERATED_MESH_CYLINDER_MESH_TYPE) THEN
              basisNumberOfNodes=basisNumberOfNodes-(numberOfElementsXi(1)+1)*(basis%NUMBER_OF_nodes_xic(1)-1)* &
                  & (numberOfElementsXi(3)+1)*(basis%NUMBER_OF_nodes_xic(3)-1)
            ENDIF
            nodeOffset=nodeOffset+basisNumberOfNodes-numberCornerNodes
          ENDDO
          basis=>bases(basisIndex)%ptr
          temporaryTerm=1
          DO xiIdx=1,meshDimension
            basisElementFactor(xiIdx)=basis%NUMBER_OF_NODES_XIC(xiIdx)-1
            IF(xiIdx>1) THEN
              temporaryTerm=temporaryTerm*((basis%NUMBER_OF_NODES_XIC(xiIdx-1)-1)*numberOfElementsXi(xiIdx-1)+1)
              basisElementFactor(xiIdx)=basisElementFactor(xiIdx)*temporaryTerm
            ENDIF
          ENDDO
          !Adjust for other mesh types
          IF(generatedMesh%generatedType==GENERATED_MESH_CYLINDER_MESH_TYPE) THEN
            !subtract nodes along line where y wraps around
            basisElementFactor(3)=basisElementFactor(3)-(numberOfElementsXi(1)*(basis%NUMBER_OF_NODES_XIC(1)-1)+1)* &
              & (basis%NUMBER_OF_NODES_XIC(3)-1)
          ENDIF
          IF(nodeUserNumber<=numberCornerNodes) THEN
            !we have a node on a corner
            IF(meshDimension>2) THEN
              position(3)=remainder/cornerNodeFactor(3)
              remainder=MOD(remainder,cornerNodeFactor(3))
            ENDIF
            IF(meshDimension>1) THEN
              position(2)=remainder/cornerNodeFactor(2)
              remainder=MOD(remainder,cornerNodeFactor(2))
            ENDIF
            position(1)=remainder/cornerNodeFactor(1)
            GeneratedMeshUserNumberToComponentNode=position(1)*basisElementFactor(1)+position(2)*basisElementFactor(2)+ &
              & position(3)*basisElementFactor(3)
            GeneratedMeshUserNumberToComponentNode=GeneratedMeshUserNumberToComponentNode+1
          ELSE IF(nodeUserNumber>nodeOffset) THEN
            remainder=remainder-nodeOffset
            DO xiIdx=1,meshDimension
              basisElementFactor(xiIdx)=basisElementFactor(xiIdx)-cornerNodeFactor(xiIdx)
            ENDDO
            numberPreviousCorners=0
            finishedCount=.FALSE.
            offEdge=.FALSE.
            IF(meshDimension>2) THEN
              IF(generatedMesh%generatedType==GENERATED_MESH_CYLINDER_MESH_TYPE.AND. &
                  & (MOD(remainder,basisElementFactor(3)) > basisElementFactor(2)*numberOfElementsXi(2)-1)) THEN
                offEdge=.TRUE.
              ELSE IF(generatedMesh%generatedType==GENERATED_MESH_REGULAR_MESH_TYPE.AND. &
                  & MOD(remainder,basisElementFactor(3)) > (basisElementFactor(2)*numberOfElementsXi(2)+ &
                  & basisElementFactor(1)*numberOfElementsXi(1)-1)) THEN
                offEdge=.TRUE.
              ENDIF
              IF(offEdge) THEN
                numberPreviousCorners=numberPreviousCorners+cornerNodeFactor(3)*(1+remainder/basisElementFactor(3))
                remainder=MOD(remainder,basisElementFactor(3))
                finishedCount=.TRUE.
              ELSE
                numberPreviousCorners=numberPreviousCorners+cornerNodeFactor(3)*(remainder/basisElementFactor(3))
                remainder=MOD(remainder,basisElementFactor(3))
              ENDIF
            ENDIF
            IF((meshDimension>1) .AND. (finishedCount.NEQV..TRUE.)) THEN
              IF(MOD(remainder,basisElementFactor(2)) > &
                  & basisElementFactor(1)*numberOfElementsXi(1)-1) THEN
                numberPreviousCorners=numberPreviousCorners+cornerNodeFactor(2)*(1+remainder/basisElementFactor(2))
                remainder=MOD(remainder,basisElementFactor(2))
                finishedCount=.TRUE.
              ELSE
                numberPreviousCorners=numberPreviousCorners+cornerNodeFactor(2)*(remainder/basisElementFactor(2))
                remainder=MOD(remainder,basisElementFactor(2))
              ENDIF
            ENDIF
            IF(finishedCount.NEQV..TRUE.) THEN
              numberPreviousCorners=numberPreviousCorners+cornerNodeFactor(1)*(remainder/basisElementFactor(1))+1
            ENDIF
            nodeOffset=nodeOffset-numberPreviousCorners
            GeneratedMeshUserNumberToComponentNode=nodeUserNumber-nodeOffset
          ELSE
<<<<<<< HEAD
            CALL FlagError("Invalid node number specified.",err,error,*999)
          ENDIF
        ENDIF
      ELSE
        localError="Mesh component must be less than or equal to "//(NumberToVString(numberOfBases,"*",err,error))// &
            & " but it is "//(NumberToVString(basisIndex,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Generated mesh is not associated",err,error,*999)
    ENDIF

    CALL Exits("GeneratedMeshUserNumberToComponentNode")
    RETURN
999 CALL Errors("GeneratedMeshUserNumberToComponentNode",err,error)
    CALL Exits("GeneratedMeshUserNumberToComponentNode")
=======
            CALL FlagError("Invalid node number specified.",ERR,ERROR,*999)
          ENDIF
        ENDIF
      ELSE
        LOCAL_ERROR="Mesh component must be less than or equal to "//(NUMBER_TO_VSTRING(NUM_BASES,"*",ERR,ERROR))// &
            & " but it is "//(NUMBER_TO_VSTRING(BASIS_INDEX,"*",ERR,ERROR))//"."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Generated mesh is not associated",ERR,ERROR,*999)
    ENDIF

    EXITS("USER_NUMBER_TO_COMPONENT_NODE")
    RETURN
999 ERRORSEXITS("USER_NUMBER_TO_COMPONENT_NODE",ERR,ERROR)
>>>>>>> 894d3af948ec146bfb9de99b726655046f5dd36b
    RETURN
    
  END FUNCTION GeneratedMeshUserNumberToComponentNode

  !
  !================================================================================================================================
  !

END MODULE GeneratedMeshRoutines



