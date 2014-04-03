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

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !> \addtogroup GeneratedMeshRoutines_GeneratedMeshTypes GeneratedMeshRoutines::GeneratedMeshTypes
  !> \brief Generated mesh types.
  !> \see GeneratedMeshRoutines,OPENCMISS_GeneratedMeshTypes
  !>@{
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_REGULAR_MESH_TYPE=1 !<A regular generated mesh. \see GeneratedMeshRoutines_GeneratedMeshTypes,GeneratedMeshRoutines
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_POLAR_MESH_TYPE=2 !<A polar generated mesh. \see GeneratedMeshRoutines_GeneratedMeshTypes,GeneratedMeshRoutines
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_FRACTAL_TREE_MESH_TYPE=3 !<A fractal tree generated mesh. \see GeneratedMeshRoutines_GeneratedMeshTypes,GeneratedMeshRoutines
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_CYLINDER_MESH_TYPE=4 !<A cylinder generated mesh. \see GeneratedMeshRoutines_GeneratedMeshTypes,GeneratedMeshRoutines
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_ELLIPSOID_MESH_TYPE=5 !<An ellipsoid generated mesh. \see GeneratedMeshRoutines_GeneratedMeshTypes,GeneratedMeshRoutines
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
  INTERFACE GeneratedMeshCreateStart
    MODULE PROCEDURE GeneratedMeshCreateStartInterface
    MODULE PROCEDURE GeneratedMeshCreateStartRegion
  END INTERFACE GeneratedMeshCreateStart

  !>Initialises the generated meshes for a region or interface.
  INTERFACE GeneratedMeshesInitialise
    MODULE PROCEDURE GeneratedMeshesInitialiseInterface
    MODULE PROCEDURE GeneratedMeshesInitialiseRegion 
  END INTERFACE GeneratedMeshesInitialise

  !>Finds a generated mesh in a list of generated meshes in a region or interface.
  INTERFACE GeneratedMeshUserNumberFind
    MODULE PROCEDURE GeneratedMeshUserNumberFindInterface
    MODULE PROCEDURE GeneratedMeshUserNumberFindRegion 
  END INTERFACE GeneratedMeshUserNumberFind

  PUBLIC GENERATED_MESH_REGULAR_MESH_TYPE,GENERATED_MESH_POLAR_MESH_TYPE
  PUBLIC GENERATED_MESH_FRACTAL_TREE_MESH_TYPE,GENERATED_MESH_CYLINDER_MESH_TYPE
  PUBLIC GENERATED_MESH_ELLIPSOID_MESH_TYPE
  
  PUBLIC GENERATED_MESH_CYLINDER_INNER_SURFACE,GENERATED_MESH_CYLINDER_OUTER_SURFACE
  PUBLIC GENERATED_MESH_CYLINDER_TOP_SURFACE,GENERATED_MESH_CYLINDER_BOTTOM_SURFACE
  PUBLIC GENERATED_MESH_ELLIPSOID_INNER_SURFACE,GENERATED_MESH_ELLIPSOID_OUTER_SURFACE
  PUBLIC GENERATED_MESH_ELLIPSOID_TOP_SURFACE
  PUBLIC GENERATED_MESH_REGULAR_LEFT_SURFACE,GENERATED_MESH_REGULAR_RIGHT_SURFACE
  PUBLIC GENERATED_MESH_REGULAR_TOP_SURFACE,GENERATED_MESH_REGULAR_BOTTOM_SURFACE
  PUBLIC GENERATED_MESH_REGULAR_FRONT_SURFACE,GENERATED_MESH_REGULAR_BACK_SURFACE
  PUBLIC GENERATED_MESHES_INITIALISE,GeneratedMeshesFinalise

  PUBLIC GeneratedMeshesInitialise

  PUBLIC GeneratedMeshBaseVectorsSet

  PUBLIC GeneratedMeshCoordinateSystemGet

  PUBLIC GeneratedMeshCreateStartGeneratedMeshCreateFinish

  PUBLIC GeneratedMeshDestroy

  PUBLIC GeneratedMeshBasisSet,GeneratedMeshBasisGet

  PUBLIC GeneratedMeshExtentSet,GeneratedMeshExtentGet

  PUBLIC GeneratedMeshNumberOfElementsSet,GeneratedMeshNumberOfElementsGet

  PUBLIC GeneratedMeshOriginSet,GeneratedMeshOriginGet

  PUBLIC GeneratedMeshTypeSet,GeneratedMeshTypeGet

  PUBLIC GeneratedMeshGeometricParametersCalculate

  PUBLIC GeneratedMeshRegionGet

  PUBLIC GeneratedMeshUserNumberFind

  PUBLIC GeneratedMeshUserNumberFind
  
  PUBLIC GENERATED_MESH_SURFACE_GET

CONTAINS

  !
  !================================================================================================================================
  !

  !>Gets the basis of a generated mesh. \see OPENCMISS::CMISSGeneratedMeshBasisGet
  SUBROUTINE GeneratedMeshBasisGet(generatedMesh,bases,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh to get the bases of
    TYPE(BASIS_PTR_TYPE) :: bases(:) !<BASES(basisIdx). On return, the bases of mesh to generate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: basisIdx,numberOfBases

    CALL Enters("GeneratedMeshBasisGet",err,error,*999)

    IF(ASSOCIATED(generatedMesh)) THEN
      IF(generatedMesh%generatedMeshFinished) THEN
        SELECT CASE(generatedMesh%generatedType)
        CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
          IF(ASSOCIATED(generatedMesh%regularMesh)) THEN
            IF(ALLOCATED(generatedMesh%regularMesh%bases)) THEN
              IF(SIZE(bases,1)>=SIZE(generatedMesh%regularMesh%bases,1)) THEN
                DO basisIdx=1,SIZE(generatedMesh%regularMesh%bases,1)
                  IF(ASSOCIATED(basis(basisIdx)%ptr)) THEN
                    localError="The pointer at location "//NumberToVString(TRIM(basisIdx,"*",err,error))// &
                      & " in the specified bases array is associated."
                    CALL FlagError(localError,err,error,*999)
                  ELSE                    
                    bases(basisIdx)%ptr=>generatedMesh%regularMesh%bases(basisIdx)%ptr
                  ENDIF
                ENDDO
              ELSE
                localError="The size of the specified bases array, "// &
                  & TRIM(NumberToVString(SIZE(bases,1),"*",err,error))// &
                  & ", is too small to hold the number of bases in the generated mesh of "// &
                  & TRIM(NumberToVString(SIZE(generatedMesh%regularMesh%bases,1),"*",err,error))//"."
                CALL FlagError(localError,err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Generated mesh bases are not allocated.",err,error,*999)
            ENDIF
          ELSE
            CALL FlagError("Generated mesh regular mesh is not associated.",err,error,*999)
          ENDIF
        CASE(GENERATED_MESH_POLAR_MESH_TYPE)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
          IF(ASSOCIATED(generatedMesh%cylinderMesh)) THEN
            IF(ALLOCATED(generatedMesh%cylinderMesh%bases)) THEN
              IF(SIZE(bases,1)>=SIZE(generatedMesh%cylinderMesh%bases,1)) THEN
                DO basisIdx=1,SIZE(generatedMesh%cylinderMesh%bases,1)
                  IF(ASSOCIATED(basis(basisIdx)%ptr)) THEN
                    localError="The pointer at location "//NumberToVString(TRIM(basisIdx,"*",err,error))// &
                      & " in the specified bases array is associated."
                    CALL FlagError(localError,err,error,*999)
                  ELSE                    
                    bases(basisIdx)%ptr=>generatedMesh%cylinderMesh%bases(basisIdx)%ptr
                  ENDIF
                ENDDO
              ELSE
                localError="The size of the specified bases array, "// &
                  & TRIM(NumberToVString(SIZE(bases,1),"*",err,error))// &
                  & ", is too small to hold the number of bases in the generated mesh of "// &
                  & TRIM(NumberToVString(SIZE(generatedMesh%cylinderMesh%bases,1),"*",err,error))//"."
                CALL FlagError(localError,err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Generated mesh bases are not allocated.",err,error,*999)
            ENDIF
          ELSE
            CALL FlagError("Generated mesh cylinder mesh is not associated.",err,error,*999)
          ENDIF
        CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
          IF(ASSOCIATED(generatedMesh%ellipsoidMesh)) THEN
            IF(ALLOCATED(generatedMesh%ellipsoidMesh%bases)) THEN
              IF(SIZE(bases,1)>=SIZE(generatedMesh%ellipsoidMesh%bases,1)) THEN
                DO basisIdx=1,SIZE(generatedMesh%ellipsoidMesh%bases,1)
                  IF(ASSOCIATED(basis(basisIdx)%ptr)) THEN
                    localError="The pointer at location "//NumberToVString(TRIM(basisIdx,"*",err,error))// &
                      & " in the specified bases array is associated."
                    CALL FlagError(localError,err,error,*999)
                  ELSE                    
                    bases(basisIdx)%ptr=>generatedMesh%ellipsoidMesh%bases(basisIdx)%ptr
                  ENDIF
                ENDDO
              ELSE
                localError="The size of the specified bases array, "// &
                  & TRIM(NumberToVString(SIZE(bases,1),"*",err,error))// &
                  & ", is too small to hold the number of bases in the generated mesh of "// &
                  & TRIM(NumberToVString(SIZE(generatedMesh%ellipsoidMesh%bases,1),"*",err,error))//"."
                CALL FlagError(localError,err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Generated mesh bases are not allocated.",err,error,*999)
            ENDIF
          ELSE
            CALL FlagError("Generated mesh ellipsoid mesh is not associated.",err,error,*999)
          ENDIF
        CASE DEFAULT
          localError="The generated mesh generated type of "// &
              & TRIM(NumberToVString(generatedMesh%generatedType,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        CALL FlagError("Generated mesh has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Generated mesh is already associated.",err,error,*999)
    ENDIF

    CALL Exits("GeneratedMeshBasisGet")
    RETURN
999 CALL Errors("GeneratedMeshBasisGet",err,error)
    CALL Exits("GeneratedMeshBasisGet")
    RETURN 1
  END SUBROUTINE GeneratedMeshBasisGet

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
    INTEGER(INTG) :: coordinateDimension,basisIdx,numberOfBases,numberOfXi,basisType
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: coordinateSystem
    TYPE(VARYING_STRING) :: localError

    NULLIFY(coordinateSystem)

    CALL Enters("GeneratedMeshBasisSet",err,error,*999)

    IF(ASSOCIATED(generatedMesh)) THEN
      IF(generatedMesh%generatedMeshFinished) THEN
        CALL FlagError("Generated mesh has been finished.",err,error,*999)
      ELSE
        CALL GeneratedMeshCoordinateSystemGet(generatedMesh,coordinateSystem,err,error,*999)
        CALL COORDINATE_SYSTEM_DIMENSION_GET(coordinateSystem,coordinateDimension,err,error,*999)
        numberOfBases=SIZE(bases,1)
        IF(numberOfBases >= 1) THEN
          !Check the supplied bases
          DO basisIdx=1,numberOfBases
            IF(ASSOCIATED(bases(basisIdx)%ptr)) THEN
              IF(bases(basisIdx)%ptr%NUMBER_OF_XI<=coordinateDimension) THEN
                IF(basisIdx == 1) THEN
                  numberOfXi=bases(1)%ptr%NUMBER_OF_XI
                  basisType=bases(1)%ptr%type
                ELSE                  
                  IF(bases(basisIdx)%ptr%NUMBER_OF_XI /= numberOfXi) THEN
                    CALL FlagError("All bases must have the same number of xi.",err,error,*999)
                  ENDIF
                  IF(bases(basesIdx)%ptr%type /= basisType) THEN
                    CALL FlagError("Using different basis types is not supported for generated meshes.",err,error,*999)
                  ENDIF
                ENDIF
                IF(.NOT.ALL(basis%COLLAPSED_XI==BASIS_NOT_COLLAPSED)) &
                  & CALL FlagError("Degenerate (collapsed) basis not implemented.",err,error,*999)
              ELSE
                localError="The basis number of xi dimensions of "// &
                  & TRIM(NumberToVString(bases(basis_idx)%ptr%NUMBER_OF_XI,"*",err,error))// &
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
          SELECT CASE(generatedMesh%generatedType)
          CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
            IF(ASSOCIATED(generatedMesh%regularMesh)) THEN
              !Store the bases
              IF(ALLOCATED(generatedMesh%regularMesh%bases)) DEALLOCATE(generatedMesh%regularMesh%bases)
              ALLOCATE(generatedMesh%regularMesh%bases(numberOfBases),STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate bases.",err,error,*999)
              DO basisIdx=1,numberOfBases
                generatedMesh%regularMesh%bases(basisIdx)%ptr=>bases(basisIdx)%ptr
              ENDDO !basisIdx
              !Reset the number of elements in each xi direction
              ALLOCATE(newNumberOfElementsXi(numberOfXi),STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate new number of elements xi.",err,error,*999)
              IF(generatedMesh%regularMesh%meshDimension==0) THEN
                !First time, default attributes
                newNumberOfElementsXi(1:numberOfXi)=1
              ELSE IF(generatedMesh%regularMesh%meshDimension>numberOfXi) THEN
                !New number of xi is less than the old number of xi
                newNumberOfElementsXi(1:numberOfXi)=generatedMesh%regularMesh%numberOfElementsXi(1:numberOfXi)
              ELSE
                !New number of xi is more than the old number of xi
                newNumberOfElementsXi(1:generatedMesh%regularMesh%meshDimension)= &
                  & generatedMesh%regularMesh%numberOfElementsXi(1:generatedMesh%regularMesh%meshDimension)
                newNumberOfElementsXi(generatedMesh%regularMesh%meshDimension+1:numberOfXi)= &
                  & generatedMesh%regularMesh%numberOfElementsXi(generatedMesh%regularMesh%meshDimension+1:numberOfXi)
              ENDIF
              CALL MOVE_ALLOC(newNumberOfElementsXi,generatedMesh%regularMesh%numberOfElementsXi)
              !Reset the number of elements in each xi direction
              ALLOCATE(newBaseVectors(coordinateDimension,numberOfXi),STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate new base vectors.",err,error,*999)
              IF(ALLOCATED(generatedMesh%regularMesh%baseVectors)) THEN
                IF(generatedMesh%regularMesh%meshDimension>numberOfXi) THEN
                  !New number of xi is less than the old number of xi
                  newBaseVectors(1:coordinateDimension,1:numberOfXi)= &
                    & generatedMesh%regularMesh%baseVectors(1:coordinateDimension,1:numberOfXi)
                ELSE
                  !New number of xi is more than the old number of xi
                  newBaseVectors(1:coordinateDimension,1:generatedMesh%regularMesh%meshDimension)= &
                    & generatedMesh%regularMesh%baseVectors(1:coordinateDimension,1:generatedMesh%regularMesh%meshDimension)
                  newBaseVectors(1:coordinateDimension,generatedMesh%regularMesh%meshDimension+1:numberOfXi)= &
                    & generatedMesh%regularMesh%baseVectors(1:coordinateDimension,generatedMesh%regularMesh% &
                    & meshDimension+1:numberOfXi)
                ENDIF
              ELSE
                !We don't have any base vectors defined so default the values to the extents
                IF(numberOfXi==1) THEN
                  !The base vector is just the extent vector
                  newBaseVectors(1:coordinateDimension,1)=generatedMesh%regularMesh%maximumExtent(1:coordinateDimension)
                ELSE
                  newBaseVectors=0.0_DP
                  IF(numberOfXi<coordinateDimension) THEN
                    !Find the first number of mesh dimensions for which the extent is non-zero.
                    count=0
                    coordinateIdx=1
                    DO xiIdx=1,numberOfXi
                      DO WHILE(ABS(generatedMesh%regularMesh%maximumExtent(coordinateIdx))<=ZERO_TOLERANCE)
                        coordinateIdx=coordinateIdxdx+1
                      ENDDO !While
                      newBaseVectors(coordinateIdx,xiIdx)=generatedMesh%regularMesh%maximumExtent(coordinateIdx)
                      coordinateIdx=coordinateIdx+1
                      count=count+1
                    ENDDO !xiIdx
                    IF(count/=numberOfXi)  &
                      & CALL FlagError("Invalid mesh extent. There number of non-zero components is < the mesh dimension.", &
                      & err,error,*999)
                  ELSE
                    !Number of xi is the same as the number of coordinates
                    !The default base vectors are aligned with the coordinate vectors
                    DO coordinateIdx=1,coordinateDimension
                      newBaseVectors(coordinateIdx,coordinateIdx)=generatedMesh%regularMesh%maximumExtent(coordinateIdx)
                    ENDDO !coordinate_idx
                  ENDIF
                ENDIF
              ENDIF
              CALL MOVE_ALLOC(newNumberOfElementsXi,generatedMesh%regularMesh%numberOfElementsXi)
              !Set the mesh dimension to be the new number of Xi
              generatedMesh%regularMesh%meshDimension=numberOfXi
            ELSE
              CALL FlagError("Regular generated mesh is not associated.",err,error,*999)
            ENDIF
          CASE(GENERATED_MESH_POLAR_MESH_TYPE)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
            IF(ASSOCIATED(generatedMesh%cylinderMesh)) THEN
              IF(ALLOCATED(generatedMesh%cylinderMesh%baseVectors)) THEN
                CALL FlagError("Can not reset the basis if base vectors have been specified.",err,error,*999)
              ELSE
                IF(ALLOCATED(generatedMesh%cylinderMesh%bases)) DEALLOCATE(generatedMesh%cylinderMesh%bases)
                ALLOCATE(generatedMesh%cylinderMesh%bases(numberOfBases),STAT=err)
                IF(err/=0) CALL FlagError("Could not allocate bases.",err,error,*999)
                DO basisIdx=1,numberOfBases
                  generatedMesh%cylinderMesh%bases(basisIdx)%ptr=>bases(basisIdx)%ptr
                ENDDO !basisIdx
                generatedMesh%cylinderMesh%meshDimension=numberOfXi
              ENDIF
            ELSE
              CALL FlagError("Cylinder generated mesh is not associated.",err,error,*999)
            ENDIF
          CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
            IF(ASSOCIATED(generatedMesh%ellipsoidMesh)) THEN
              IF(ALLOCATED(generatedMesh%ellipsoidMesh%bases)) DEALLOCATE(generatedMesh%ellipsoidMesh%bases)
              ALLOCATE(generatedMesh%ellipsoidMesh%bases(numberOfBases),STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate bases.",err,error,*999)
              DO basisIdx=1,numberOfBases
                generatedMesh%ellipsoidMesh%bases(basisIdx)%ptr=>bases(basisIdx)%ptr
              ENDDO !basisIdx
              generatedMesh%ellipsoidMesh%meshDimension=numberOfXi
            ELSE
              CALL FlagError("Ellpsoid generated mesh is not associated.",err,error,*999)
            ENDIF
          CASE DEFAULT
            localError="The generated mesh type of "//TRIM(NumberToVString(generatedMesh%generatedType,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ENDIF
      ELSEIF
        CALL FlagError("No bases where supplied.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Generated mesh is already associated.",err,error,*999)
    ENDIF

    CALL Exits("GeneratedMeshBasisSet")
    RETURN
999 CALL Errors("GeneratedMeshBasisSet",err,error)
    CALL Exits("GeneratedMeshBasisSet")
    RETURN 1
  END SUBROUTINE GeneratedMeshBasisSet

  !
  !================================================================================================================================
  !

  !>Sets/changes the base vectors of a generated mesh. \see OPENCMISS::CMISSGeneratedMeshBaseVectorsSet
  SUBROUTINE GeneratedMeshBaseVectorsSet(generatedMesh,baseVectors,err,error,*)

    !Argument variables
    TYPE(generatedMesh_TYPE), POINTER :: generatedMesh !<A pointer to the generated mesh to set the base vectors fo
    REAL(DP), INTENT(IN) :: baseVectors(:,:) !<baseVectors(coordinateIdx,xiIdx). The base vectors for the generated mesh to set.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: coordinateDimension
    TYPE(BASIS_TYPE), POINTER :: basis
    TYPE(BASIS_PTR_TYPE), POINTER :: bases(:)
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: coordinateSystem
    TYPE(VARYING_STRING) :: localError

    NULLIFY(coordinateSystem)
    NULLIFY(bases)
    NULLIFY(basis)
    
    CALL Enters("GeneratedMeshBaseVectorsSet",err,error,*999)

    IF(ASSOCIATED(generatedMesh)) THEN
      IF(generatedMesh%generatedMeshFinished) THEN
        CALL FlagError("Generated mesh has been finished.",err,error,*999)
      ELSE
        CALL GeneratedMeshCoordinateSystemGet(generatedMesh,coordianteSystem,err,error,*999)
        CALL COORDINATE_SYSTEM_DIMENSION_GET(coordinateSystem,coordinateDimension,err,error,*999)
        IF(SIZE(baseVectors,1)==coordinateDimension) THEN
          SELECT CASE(generatedMesh%generatedType)
          CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
            IF(ASSOCIATED(generatedMesh%regularMesh)) THEN
              bases=>generatedMesh%regularMesh%bases
              IF(ASSOCIATED(bases)) THEN
                basis=>bases(1)%ptr !Bases should all have same number of xi
                IF(ASSOCIATED(bases)) THEN
                  IF(SIZE(baseVectors,2)==basis%NUMBER_OF_XI) THEN
                    IF(ALLOCATED(generatedMesh%regularMesh%baseVectors)) DEALLOCATE(generatedMesh%regularMesh%baseVectors)
                    ALLOCATE(generatedMesh%regularMesh%baseVectors(SIZE(baseVectors,1),SIZE(baseVectors,2)),STAT=err)
                    IF(err/=0) CALL FlagError("Could not allocate base vectors.",err,error,*999)
                    generatedMesh%regularMesh%baseVectors=baseVectors
                  ELSE
                    localError="The size of the second dimension of base vectors of "// &
                      & TRIM(NumberToVString(SIZE(baseVectors,2),"*",err,error))// &
                      & " is invalid. The second dimension size must match the number of mesh dimensions of "// &
                      & TRIM(NumberToVString(basis%NUMBER_OF_XI,"*",err,error))//"."
                    CALL FlagError(localError,err,error,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Bases are not associated.",err,error,*999)
                ENDIF
              ELSE
                CALL FlagError("You must set the generated mesh basis before setting base vectors.",err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Regular generated mesh is not associated.",err,error,*999)
            ENDIF
          CASE(GENERATED_MESH_POLAR_MESH_TYPE)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
            IF(ASSOCIATED(generatedMesh%cylinderMesh)) THEN
              bases=>generatedMesh%cylinderMesh%bases
              IF(ASSOCIATED(bases)) THEN
                basis=>bases(1)%ptr !Bases should all have same number of xi
                IF(ASSOCIATED(bases)) THEN
                  IF(SIZE(baseVectors,2)==basis%NUMBER_OF_XI) THEN
                    IF(ALLOCATED(generatedMesh%cylinderMesh%baseVectors)) DEALLOCATE(generatedMesh%cylinderMesh%baseVectors)
                    ALLOCATE(generatedMesh%cylinderMesh%baseVectors(SIZE(baseVectors,1),SIZE(baseVectors,2)),STAT=err)
                    IF(err/=0) CALL FlagError("Could not allocate base vectors.",err,error,*999)
                    generatedMesh%cylinderMesh%baseVectors=baseVectors
                  ELSE
                    localError="The size of the second dimension of base vectors of "// &
                      & TRIM(NumberToVString(SIZE(baseVectors,2),"*",err,error))// &
                      & " is invalid. The second dimension size must match the number of mesh dimensions of "// &
                      & TRIM(NumberToVString(basis%NUMBER_OF_XI,"*",err,error))//"."
                    CALL FlagError(localError,err,error,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Bases are not associated.",err,error,*999)
                ENDIF
              ELSE
                CALL FlagError("You must set the generated mesh basis before setting base vectors.",err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Cylinder generated mesh is not associated.",err,error,*999)
            ENDIF
          CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The generated mesh mesh type of "//TRIM(NumberToVString(generatedMesh%generatedType,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ELSE
          localError="The size of the first dimension of base vectors of "// &
            & TRIM(NumberToVString(SIZE(baseVectors,1),"*",err,error))// &
            & " is invalid. The first dimension size must match the coordinate system dimension of "// &
            & TRIM(NumberToVString(coordinateDimension,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
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
          ENDIF
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Generated mesh is not associated.",err,error,*999)
    ENDIF

    CALL Exits("GeneratedMeshCoordinateSystemGet")
    RETURN
999 CALL Errors("GeneratedMeshCoordinateSystemGet",err,error)
    CALL Exits("GeneratedMeshCoordinateSystemGet")
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

    CALL Enters("GeneratedMeshCreateFinish",err,error,*999)

    IF(ASSOCIATED(generatedMesh)) THEN
      IF(generatedMesh%generatedMeshFinished) THEN
        CALL FlagError("Generated mesh has already been finished.",err,error,*999)
      ELSE
        IF(ASSOCIATED(mesh)) THEN
          CALL FlagError("Mesh is already associated.",err,error,*999)
        ELSE
          SELECT CASE(generatedMesh%generatedType)
          CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
            CALL GeneratedMeshRegularCreateFinish(generatedMesh,meshUserNumber,err,error,*999)
          CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
            CALL GENERATED_MESH_CYLINDER_CREATE_FINISH(generatedMesh,meshUserNumber,err,error,*999)
          CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
            CALL GeneratedMeshEllipsoidCreateFinish(generatedMesh,meshUserNumber,err,error,*999)
          CASE(GENERATED_MESH_POLAR_MESH_TYPE)
            CALL FlagError("Not implmented.",err,error,*999)
          CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The generated mesh mesh type of "// &
              & TRIM(NumberToVString(generatedMesh%generatedType,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
          !Return the pointers
          mesh=>generatedMesh%mesh
          mesh%GENERATED_MESH=>generatedMesh
          generatedMesh%generatedMeshFinished=.TRUE.
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Generated mesh is not associated.",err,error,*999)
    ENDIF

    CALL Exits("GeneratedMeshCreateFinish")
    RETURN
999 CALL Errors("GeneratedMeshCreateFinish",err,error)
    CALL Exits("GeneratedMeshCreateFinish")
    RETURN 1
    
  END SUBROUTINE GeneratedMeshCreateFinish

  !
  !================================================================================================================================
  !

  !>Starts the creation of a generic generated mesh.
  SUBROUTINE GeneratedMeshCreateStartGeneric(generatedMeshes,userNumber,generatedMesh,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshesType), POINTER :: generatedMeshes !<A pointer to the generated meshes
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the generated mesh to create
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<On exit, a pointer to the created generated mesh. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,generatedMeshIdx
    TYPE(GeneratedMeshType), POINTER :: newGeneratedMesh
    TYPE(GeneratedMeshPtrType), ALLOCATABLE :: newGeneratedMeshes(:)
    TYPE(VARYING_STRING) :: dummyError

    NULLIFY(newGeneratedMesh)

    CALL Enters("GeneratedMeshCreateStartGeneric",err,error,*997)

    IF(ASSOCIATED(generatedMeshes)) THEN
      IF(ASSOCIATED(generatedMesh)) THEN
        CALL FlagError("Generated mesh is already associated.",err,error,*997)
      ELSE
        !Initialise generated mesh
        CALL GeneratedMeshInitialise(newGeneratedMesh,err,error,*999)
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
    TYPE(VARYING_STRING) :: localError

    CALL Enters("GeneratedMeshCreateStartInterface",err,error,*999)

    IF(ASSOCIATED(interface)) THEN
      IF(ASSOCIATED(generatedMesh)) THEN
        CALL FlagError("Generated mesh is already associated.",err,error,*999)
      ELSE
        NULLIFY(generatedMesh)
        CALL GENERATED_MESH_USER_NUMBER_FIND(userNumber,interface,generatedMesh,err,error,*999)
        IF(ASSOCIATED(generatedMesh)) THEN
          localError="The specified user number of "//TRIM(NumberToVString(userNumber,"*",err,error))// &
            & " has already been used for a generated mesh."
          CALL FlagError(localError,err,error,*999)
        ELSE
          IF(ASSOCIATED(interface%GENERATED_MESHES)) THEN
            CALL GeneratedMeshCreateStartGeneric(interface%GENERATED_MESHES,userNumber,generatedMesh,err,error,*999)
            generatedMesh%interface=>interface
          ELSE
            CALL FlagError("Interface generated meshes is not associated.",err,error,*999)
          ENDIF
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Interface is not associated.",err,error,*999)
    ENDIF

    CALL Exits("GeneratedMeshCreateStartInterface")
    RETURN
999 CALL Errors("GeneratedMeshCreateStartInterface",err,error)
    CALL Exits("GeneratedMeshCreateStartInterface")
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
    TYPE(VARYING_STRING) :: localError

    CALL Enters("GeneratedMeshCreateStartRegion",err,error,*999)

    IF(ASSOCIATED(region)) THEN
      IF(ASSOCIATED(generatedMesh)) THEN
        CALL FlagError("Generated mesh is already associated.",err,error,*999)
      ELSE
        NULLIFY(generatedMesh)
        CALL GENERATED_MESH_USER_NUMBER_FIND(userNumber,region,generatedMesh,err,error,*999)
        IF(ASSOCIATED(generatedMesh)) THEN
          localError="The specified user number of "//TRIM(NumberToVString(userNumber,"*",err,error))// &
            & " has already been used for a generated mesh."
          CALL FlagError(localError,err,error,*999)
        ELSE
          IF(ASSOCIATED(region%GENERATED_MESHES)) THEN
            CALL GeneratedMeshCreateStartGeneric(region%GENERATED_MESHES,userNumber,generatedMesh,err,error,*999)
            generatedMesh%region=>region
          ELSE
            CALL FlagError("Region generated meshes is not associated.",err,error,*999)
          ENDIF
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Region is not associated.",err,error,*999)
    ENDIF

    CALL Exits("GeneratedMeshCreateStartRegion")
    RETURN
999 CALL Errors("GeneratedMeshCreateStartRegion",err,error)
    CALL Exits("GeneratedMeshCreateStartRegion")
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
    INTEGER(INTG) :: generatedMeshIdx,generatedMeshPosition
    TYPE(GeneratedMeshesType), POINTER :: generatedMeshes
    TYPE(GeneratedMeshPtrType), ALLOCATABLE :: newGeneratedMeshes(:)

    CALL Enters("GeneratedMeshDestroy",err,error,*998)

    IF(ASSOCIATED(generatedMesh)) THEN
      generatedMeshes=>generatedMesh%generatedMeshes
      IF(ASSOCIATED(generatedMeshes)) THEN
        IF(ASSOCIATED(generatedMeshes%generatedMeshes)) THEN
          generatedMeshPosition=generatedMesh%globalNumber
          CALL GeneratedMeshFinalise(generatedMesh,err,error,*999)
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
              ENDIF
            ENDDO !generatedMeshIdx
            CALL MOVE_ALLOC(newGeneratedMeshes,generatedMeshes%generatedMeshes)
            generatedMeshes%numberOfGeneratedMeshes=generatedMeshes%numberOfGeneratedMeshes-1
          ELSE
            DEALLOCATE(generatedMeshes%generatedMeshes)
            generatedMeshes%numberOfGeneratedMeshes=0
          ENDIF
        ELSE
          CALL FlagError("Generated meshes are not associated",err,error,*998)
        ENDIF
      ELSE
        CALL FlagError("Generated mesh generated meshes is not associated.",err,error,*998)
      ENDIF
    ELSE
      CALL FlagError("Generated mesh is not associated",err,error,*998)
    END IF

    CALL Exits("GeneratedMeshDestroy")
    RETURN
999 IF(ASSOCIATED(newGeneratedMeshes)) DEALLOCATE(newGeneratedMeshes)
998 CALL Errors("GeneratedMeshDestroy",err,error)
    CALL Exits("GeneratedMeshDestroy")
    RETURN 1
    
  END SUBROUTINE GeneratedMeshDestroy

  !
  !================================================================================================================================
  !

  !>Gets the extent of a generated mesh. \see OPENCMISS::CMISSGeneratedMeshExtentGet
  SUBROUTINE GeneratedMeshExtentGet(generatedMesh,extent,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh to get the type of
    REAL(DP), INTENT(OUT) :: extent(:) !<On return, the mesh extent. For regular meshes this is the maximum extent per axis. For cylindral meshes this is the inner & outer radii and length of cylinder
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    CALL Enters("GeneratedMeshExtentGet",err,error,*999)

    IF(ASSOCIATED(generatedMesh)) THEN
      SELECT CASE(generatedMesh%generatedType)
      CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
        IF(SIZE(extent,1)>=SIZE(generatedMesh%regularMesh%maximumExtent,1)) THEN
          extent(1:generatedMesh%regularMesh%coordinateDimension)= &
            & generatedMesh%regularMesh%maximumExtent(1:generatedMesh%regularMesh%coordinateDimension)
        ELSE
          localError="The size of extent is too small. The supplied size is "// &
            & TRIM(NumberToVString(SIZE(extent,1),"*",err,error))//" and it needs to be >= "// &
            & TRIM(NumberToVString(SIZE(generatedMesh%regularMesh%maximumExtent,1),"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      CASE(GENERATED_MESH_POLAR_MESH_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
        IF(SIZE(extent,1)>=SIZE(generatedMesh%cylinderMesh%cylinderExtent,1)) THEN
          extent(1:generatedMesh%cylinderMesh%coordinateDimension)= &
            & generatedMesh%cylinderMesh%cylinderExtent(1:generatedMesh%regularMesh%coordinateDimension)
        ELSE
          localError="The size of extent is too small. The supplied size is "// &
            & TRIM(NumberToVString(SIZE(extent,1),"*",err,error))//" and it needs to be >= "// &
            & TRIM(NumberToVString(SIZE(generatedMesh%cylinderMesh%cylinderExtent,1),"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
        IF(SIZE(extent,1)>=SIZE(generatedMesh%ellipsoidMesh%ellipsoidExtent,1)) THEN
          extent(1:generatedMesh%ellipsoidMesh%coordinateDimension)= &
            & generatedMesh%ellipsoidMesh%ellipsoidExtent(1:generatedMesh%ellipsoidMesh%coordinateDimension)
        ELSE
          localError="The size of extent is too small. The supplied size is "// &
            & TRIM(NumberToVString(SIZE(extent,1),"*",err,error))//" and it needs to be >= "// &
            & TRIM(NumberToVString(SIZE(generatedMesh%ellipsoidMesh%ellipsoidExtent,1),"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      CASE DEFAULT
        localError="The generated mesh mesh type of "//TRIM(NumberToVString(generatedMesh%generatedType,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Generated mesh is not associated",err,error,*999)
    ENDIF

    CALL Exits("GeneratedMeshExtentGet")
    RETURN
999 CALL Errors("GeneratedMeshExtentGet",err,error)
    CALL Exits("GeneratedMeshExtentGet")
    RETURN
  END SUBROUTINE GeneratedMeshExtentGet
  !
  !================================================================================================================================
  !

  !>Sets/changes the extent of a generated mesh. \see OPENCMISS::CMISSGeneratedMeshExtentSet
  SUBROUTINE GeneratedMeshExtentSet(generatedMesh,extent,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh to set the type of
    REAL(DP), INTENT(IN) :: extent(:) !<The extent of the generated mesh (MAXIMUM for regular type, CYLINDER for cylinder type)
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: coordinateDimension
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: coordinateSystem
    TYPE(VARYING_STRING) :: localError

    NULLIFY(coordinateSystem)
    
    CALL Enters("GeneratedMeshExtentSet",err,error,*999)

    IF(ASSOCIATED(generatedMesh)) THEN
      IF(generatedMesh%generatedMeshFinished) THEN
        CALL FlagError("Generated mesh has been finished.",err,error,*999)
      ELSE
        CALL GeneratedMeshCoordinateSystemGet(generatedMesh,coordinateSystem,err,error,*999)
        CALL COORDINATE_SYSTEM_DIMENSION_GET(coordinateSystem,coordinateDimension,err,error,*999)
        SELECT CASE(generatedMesh%generatedType)
        CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
          IF(SIZE(extent,1)>=coordinateDimension) THEN
            IF(ASSOCIATED(generatedMesh%regularMesh)) THEN
              IF(L2Norm(extent(1:coordinateDimension))>ZERO_TOLERANCE) THEN
                generatedMesh%regularMesh%maximumExtent(1:coordinateDimension)=extent(1:coordinateDimension)
              ELSE
                CALL FlagError("The norm of the mesh extent is zero.",err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Regular generated mesh is not associated.",err,error,*999)
            ENDIF
          ELSE
            localError="The extent size of "//TRIM(NumberToVString(SIZE(extent,1),"*",err,error))// &
                & " is invalid. The extent size must be >= the coordinate system dimension of "// &
                & TRIM(NumberToVString(coordinateDimension,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        CASE(GENERATED_MESH_POLAR_MESH_TYPE)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
          IF(ASSOCIATED(generatedMesh%cylinderMesh)) THEN
            meshDimension=generatedMesh%cylinderMesh%meshDimension
            IF(SIZE(extent,1)>=meshDimension) THEN
              generatedMesh%cylinderMesh%cylinderExtent(1:meshDimension)=extent(1:meshDimension)
            ELSE
              CALL FlagError("Cylinder generated mesh is not associated.",err,error,*999)
            ENDIF
          ELSE
            localError="The extent size of "//TRIM(NumberToVString(SIZE(extent,1),"*",err,error))// &
                & " is invalid. The extent size must >= the mesh dimension of "// &
                & TRIM(NumberToVString(meshDimension,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
          IF((SIZE(extent,1)-1)>=coordinateDimension) THEN
            IF(ASSOCIATED(generatedMesh%ellipsoidMesh)) THEN
              generatedMesh%ellipsoidMesh%ellipsoidExtent(1:coordinateDimension+1)=extent(1:coordinateDimension+1)
            ELSE
              CALL FlagError("Ellipsoid generated mesh is not associated.",err,error,*999)
            ENDIF
          ELSE
            localError="The extent size of "//TRIM(NumberToVString(SIZE(extent,1),"*",err,error))// &
                & " is invalid. The extent size must >= to one plus the coordinate system dimension of "// &
                & TRIM(NumberToVString(coordinateDimension,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        CASE DEFAULT
          localError="The generated mesh mesh type of "// &
              & TRIM(NumberToVString(generatedMesh%generatedType,"*",err,error))// &
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
    RETURN 1
  END SUBROUTINE GeneratedMeshExtentSet

  !
  !================================================================================================================================
  !

  !>Get one of the surfaces of a generated mesh. \see OPENCMISS::CMISSGeneratedMeshSurfaceGet
  SUBROUTINE GENERATED_MESH_SURFACE_GET(GENERATED_MESH,MESH_COMPONENT,SURFACE_TYPE,SURFACE_NODES,NORMAL_XI,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: GENERATED_MESH !<A pointer to the generated mesh to get the type of
    INTEGER(INTG), INTENT(IN) :: MESH_COMPONENT !<The mesh component to get the surface for.
    INTEGER(INTG), INTENT(IN) :: SURFACE_TYPE !<The surface you are interested in
    INTEGER(INTG), ALLOCATABLE, INTENT(OUT) :: SURFACE_NODES (:) !<The nodes on the specified surface
    INTEGER(INTG), INTENT(OUT) :: NORMAL_XI !<The normal outward pointing xi direction
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(GeneratedMeshEllipsoidType), POINTER :: ELLIPSOID_MESH
    TYPE(GeneratedMeshCylinderType), POINTER :: CYLINDER_MESH
    TYPE(GeneratedMeshRegularType), POINTER :: REGULAR_MESH
    TYPE(VARYING_STRING) :: localError
!     INTEGER(INTG), ALLOCATABLE :: NODES(:)


    CALL Enters("GENERATED_MESH_SURFACE_GET",err,error,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      SELECT CASE(GENERATED_MESH%generatedType)
      CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
        REGULAR_MESH=>GENERATED_MESH%regularMesh
        CALL GENERATED_MESH_REGULAR_SURFACE_GET(REGULAR_MESH,MESH_COMPONENT,SURFACE_TYPE,SURFACE_NODES,NORMAL_XI,err,error,*999)
      CASE(GENERATED_MESH_POLAR_MESH_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
        CYLINDER_MESH=>GENERATED_MESH%cylinderMesh
        CALL GENERATED_MESH_CYLINDER_SURFACE_GET(CYLINDER_MESH,MESH_COMPONENT,SURFACE_TYPE,SURFACE_NODES,NORMAL_XI,err,error,*999)
      CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
        ELLIPSOID_MESH=>GENERATED_MESH%ellipsoidMesh
        CALL GENERATED_MESH_ELLIPSOID_SURFACE_GET(ELLIPSOID_MESH,MESH_COMPONENT,SURFACE_TYPE,SURFACE_NODES,NORMAL_XI, &
            & err,error,*999)
      CASE DEFAULT
        localError="The generated mesh mesh type of "//TRIM(NumberToVString(GENERATED_MESH%generatedType,"*",err,error))// &
            & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Generated mesh is not associated",err,error,*999)
    ENDIF

    CALL Exits("GENERATED_MESH_SURFACE_GET")
    RETURN
999 CALL Errors("GENERATED_MESH_SURFACE_GET",err,error)
    CALL Exits("GENERATED_MESH_SURFACE_GET")
    RETURN
  END SUBROUTINE GENERATED_MESH_SURFACE_GET

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

    CALL Enters("GeneratedMeshFinalise",err,error,*999)

    IF(ASSOCIATED(generatedMesh)) THEN
      CALL GENERATED_MESH_REGULAR_FINALISE(generatedMesh%regularMesh,err,error,*999)
      CALL GENERATED_MESH_CYLINDER_FINALISE(generatedMesh%cylinderMesh,err,error,*999)
      CALL GENERATED_MESH_ELLIPSOID_FINALISE(generatedMesh%ellipsoidMesh,err,error,*999)
      DEALLOCATE(generatedMesh)
    ENDIF

    CALL Exits("GeneratedMeshFinalise")
    RETURN
999 CALL Errors("GeneratedMeshFinalise",err,error)
    CALL Exits("GeneratedMeshFinalise")
    RETURN 1
  END SUBROUTINE GeneratedMeshFinalise

  !
  !================================================================================================================================
  !

  !>Initialises a generated mesh.
  SUBROUTINE GeneratedMeshInitialise(generatedMesh,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    CALL Enters("GeneratedMeshInitialise",err,error,*999)

    IF(ASSOCIATED(generatedMesh)) THEN
      CALL FlagError("Generated mesh is already associated.",err,error,*999)
    ELSE
      ALLOCATE(generatedMesh,STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate generated mesh.",err,error,*999)
      generatedMesh%userNumber=0
      generatedMesh%globalNumber=0
      generatedMesh%generatedMeshFinished=.FALSE.
      NULLIFY(generatedMesh%region)
      NULLIFY(generatedMesh%interface)
      generatedMesh%generatedType=0
      NULLIFY(generatedMesh%regularMesh)
      NULLIFY(generatedMesh%cylinderMesh)
      NULLIFY(generatedMesh%ellipsoidMesh)
      NULLIFY(generatedMesh%mesh)
      !Default to a regular mesh.
      CALL GENERATED_MESH_REGULAR_INITIALISE(generatedMesh,err,error,*999)
    ENDIF

    CALL Exits("GeneratedMeshInitialise")
    RETURN
999 CALL Errors("GeneratedMeshInitialise",err,error)
    CALL Exits("GeneratedMeshInitialise")
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

    CALL Enters("GeneratedMeshNumberOfElementsGet",err,error,*999)

    IF(ASSOCIATED(generatedMesh)) THEN
      SELECT CASE(generatedMesh%generatedType)
      CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
        meshDimension=generatedMesh%regularMesh%meshDimension
        IF(SIZE(numberOfElements,1)>=meshDimension) THEN
          numberOfElements(1:meshDimension)=generatedMesh%regularMesh%numberOfElementsXi(1:meshDimension)
        ELSE
          localError="The size of number of elements is too small. The supplied size is "// &
            & TRIM(NumberToVString(SIZE(numberOfElements,1),"*",err,error))//" and it needs to be >= "// &
            & TRIM(NumberToVString(meshDimension,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      CASE(GENERATED_MESH_POLAR_MESH_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
        meshDimension=generatedMesh%cylinderMesh%meshDimension
        IF(SIZE(numberOfElements,1)>=meshDimension) THEN
          numberOfElements(1:meshDimension)=generatedMesh%cylinderMesh%numberOfElementsXi(1:meshDimension)
        ELSE
          localError="The size of number of elements is too small. The supplied size is "// &
            & TRIM(NumberToVString(SIZE(numberOfElements,1),"*",err,error))//" and it needs to be >= "// &
            & TRIM(NumberToVString(meshDimension,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
        meshDimension=generatedMesh%cylinderMesh%meshDimension
        IF(SIZE(numberOfElements,1)>=meshDimension) THEN
          numberOfElements(1:meshDimension)=generatedMesh%ellipsoidMesh%numberOfElementsXi(1:meshDimension)
        ELSE
          localError="The size of number of elements is too small. The supplied size is "// &
            & TRIM(NumberToVString(SIZE(numberOfElements,1),"*",err,error))//" and it needs to be >= "// &
            & TRIM(NumberToVString(meshDimension,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      CASE DEFAULT
        localError="The generated mesh mesh type of "//TRIM(NumberToVString(generatedMesh%generatedType,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Generated mesh is not associated.",err,error,*999)
    ENDIF

    CALL Exits("GeneratedMeshNumberOfElementsGet")
    RETURN
999 CALL Errors("GeneratedMeshNumberOfElementsGet",err,error)
    CALL Exits("GeneratedMeshNumberOfElementsGet")
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
    TYPE(GeneratedMeshEllipsoidType), POINTER :: cylinderMesh
    TYPE(GeneratedMeshRegularType), POINTER :: regularMesh
    TYPE(VARYING_STRING) :: localError

    CALL Enters("GeneratedMeshNumberOfElementsSet",err,error,*999)

    IF(ASSOCIATED(generatedMesh)) THEN
      IF(generatedMesh%generatedMeshFinished) THEN
        CALL FlagError("Generated mesh has been finished.",err,error,*999)
      ELSE
        SELECT CASE(generatedMesh%generatedType)
        CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
          regularMesh=>generatedMesh%regularMesh
          IF(ASSOCIATED(regularMesh)) THEN
            meshDimension=regularMesh%meshDimension
            IF(ASSOCIATED(regularMesh%bases)) THEN
              IF(SIZE(numberOfElements,1)>=meshDimension) THEN
                IF(ALL(numberOfElements>0)) THEN
                  regularMesh%numberOfElementsXi(1:meshDimension)=numberOfElements(1:meshDimension)
                ELSE
                  CALL FlagError("Must have 1 or more elements in all directions.",err,error,*999)
                ENDIF
              ELSE
                localError="The number of elements xi size of "// &
                  & TRIM(NumberToVString(SIZE(numberOfElements,1),"*",err,error))// &
                  & " is invalid. The number of elements size must be >= the mesh dimension of "// &
                  & TRIM(NumberToVString(meshDimension,"*",err,error))//"."
                CALL FlagError(localError,err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Must set the generated mesh basis before setting the number of elements.",err,error,*999)
            ENDIF
          ELSE
            CALL FlagError("Regular generated mesh is not associated.",err,error,*999)
          ENDIF
        CASE(GENERATED_MESH_POLAR_MESH_TYPE)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
          cylinderMesh=>generatedMesh%cylinderMesh
          IF(ASSOCIATED(cylinderMesh)) THEN
            meshDimension=cylinderMesh%meshDimension
            IF(ASSOCIATED(cylinderMesh%bases)) THEN
              IF(SIZE(numberOfElements,1)>=meshDimension) THEN
                IF(ALL(numberOfElements>0)) THEN
                  cylinderMesh%numberOfElementsXi(1:meshDimension)=numberOfElements(1:meshDimension)
                ELSE
                  CALL FlagError("Must have 1 or more elements in all directions.",err,error,*999)
                ENDIF
              ELSE
                localError="The number of elements xi size of "// &
                  & TRIM(NumberToVString(SIZE(numberOfElements,1),"*",err,error))// &
                  & " is invalid. The number of elements size must be >= the mesh dimension of "// &
                  & TRIM(NumberToVString(meshDimension,"*",err,error))//"."
                CALL FlagError(localError,err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Must set the generated mesh basis before setting the number of elements.",err,error,*999)
            ENDIF
          ELSE
            CALL FlagError("Cylinder generated mesh is not associated.",err,error,*999)
          ENDIF
        CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
          ellipsoidMesh=>generatedMesh%ellipsoidMesh
          IF(ASSOCIATED(ellipsoidMesh)) THEN
            meshDimension=ellipsoidMesh%meshDimension
            IF(ASSOCIATED(ellipsoidMesh%bases)) THEN
              IF(SIZE(numberOfElements,1)>=meshDimension) THEN
                IF(ALL(numberOfElements>0)) THEN
                  ellipsoidMesh%numberOfElementsXi(1:meshDimension)=numberOfElements(1:meshDimension)
                ELSE
                  CALL FlagError("Must have 1 or more elements in all directions.",err,error,*999)
                ENDIF
              ELSE
                localError="The number of elements xi size of "// &
                  & TRIM(NumberToVString(SIZE(numberOfElements,1),"*",err,error))// &
                  & " is invalid. The number of elements size must be >= the mesh dimension of "// &
                  & TRIM(NumberToVString(meshDimension,"*",err,error))//"."
                CALL FlagError(localError,err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Must set the generated mesh basis before setting the number of elements.",err,error,*999)
            ENDIF
          ELSE
            CALL FlagError("Ellipsoid generated mesh is not associated.",err,error,*999)
          ENDIF
        CASE DEFAULT
          localError="The generated mesh mesh type of "//TRIM(NumberToVString(generatedMesh%generatedType,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FlagError("Generated mesh is not associated.",err,error,*999)
    ENDIF

    CALL Exits("GeneratedMeshNumberOfElementsSet")
    RETURN
999 CALL Errors("GeneratedMeshNumberOfElementsSet",err,error)
    CALL Exits("GeneratedMeshNumberOfElementsSet")
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
    INTEGER(INTG) :: coordinateDimension
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: coordinateSystem
    TYPE(VARYING_STRING) :: localError

    NULLIFY(coordinateSystem)

    CALL Enters("GeneratedMeshOriginGet",err,error,*999)

    IF(ASSOCIATED(generatedMesh)) THEN
      CALL GeneratedMeshCoordinateSystemGet(generatedMesh,coordinateSystem,err,error,*999)
      CALL COORDINATE_SYSTEM_DIMENSION_GET(coordianteSystem,coordinateDimension,err,error,*999)
      IF(SIZE(origin,1)>=coordinateDimension) THEN
        SELECT CASE(generatedMesh%generatedType)
        CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
          origin(1:coordinateDimension)=generatedMesh%regularMesh%origin(1:coordinateDimension)
        CASE(GENERATED_MESH_POLAR_MESH_TYPE)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
          origin(1:coordinateDimension)=generatedMesh%cylinderMesh%origin(1:coordinateDimension)
        CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
          origin(1:coordinateDimension)=generatedMesh%ellipsoidMesh%origin(1:coordinateDimension)
        CASE DEFAULT
          localError="The generated mesh mesh type of "//TRIM(NumberToVString(generatedMesh%generatedType,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        localError="The size of origin is too small. The supplied size is "// &
          & TRIM(NumberToVString(SIZE(origin,1),"*",err,error))// &
          & " and it needs to be >= the coordinate dimension of "// &
          & TRIM(NumberToVString(coordinateDimension,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Generated mesh is not associated",err,error,*999)
    ENDIF

    CALL Exits("GeneratedMeshOriginGet")
    RETURN
999 CALL Errors("GeneratedMeshOriginGet",err,error)
    CALL Exits("GeneratedMeshOriginGet")
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
    INTEGER(INTG) :: coordianteDimension
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: coordinateSystem
    TYPE(VARYING_STRING) :: localError

    NULLIFY(coordinateSystem)
    
    CALL Enters("GeneratedMeshOriginSet",err,error,*999)

    IF(ASSOCIATED(generatedMesh)) THEN
      IF(generatedMesh%generatedMeshFinished) THEN
        CALL FlagError("Generated mesh has been finished.",err,error,*999)
      ELSE
        CALL GeneratedMeshCoordinateSystemGet(generatedMesh,coordianteSystem,err,error,*999)
        CALL COORDINATE_SYSTEM_DIMENSION_GET(coordinateSystem,coordinateDimension,err,error,*999)
        IF(SIZE(origin,1)>=coordinateDimension) THEN
          SELECT CASE(generatedMesh%generatedType)
          CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
            IF(ASSOCIATED(generatedMesh%regularMesh)) THEN
              generatedMesh%regularMesh%origin(1:coordinateDimension)=origin(1:coordinateDimension)
            ELSE
              CALL FlagError("Regular generated mesh is not associated.",err,error,*999)
            ENDIF
          CASE(GENERATED_MESH_POLAR_MESH_TYPE)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
            IF(ASSOCIATED(generatedMesh%cylinderMesh)) THEN
              generatedMesh%cylinderMesh%origin(1:coordinateDimension)=origin(1:coordinateDimension)
            ELSE
              CALL FlagError("Cylinder generated mesh is not associated.",err,error,*999)
            END IF
          CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
            IF(ASSOCIATED(generatedMesh%ellipsoidMesh)) THEN
             generatedMesh%ellipsoidMesh%origin(1:coordinateDimension)=origin(1:coordinateDimension)
            ELSE
              CALL FlagError("Ellipsoid generated mesh is not associated.",err,error,*999)
            END IF
          CASE DEFAULT
            localError="The generated mesh mesh type of "//TRIM(NumberToVString(generatedMesh%generatedType,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ELSE
          localError="The origin size of "//TRIM(NumberToVString(SIZE(origin,1),"*",err,error))// &
            & " is invalid. The origin size must be >= the coordinate system dimension of "// &
            & TRIM(NumberToVString(coordinateDimension,"*",err,error))//"."
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
    RETURN 1
  END SUBROUTINE GeneratedMeshOriginSet

  !
  !================================================================================================================================
  !

  !>Start to create the regular generated mesh type
  SUBROUTINE GeneratedMeshRegularCreateFinish(generatedMesh,meshUserNumber,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh
    INTEGER(INTG), INTENT(IN) :: meshUserNumber !<The user number for the mesh to generate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: coordinate_idx,COUNT,ELEMENT_FACTOR,grid_ne,GRID_NUMBER_OF_ELEMENTS,ni,ne,ne1,ne2,ne3,nn,nn1,nn2,nn3,np, &
      & NUMBER_OF_ELEMENTS_XI(3),TOTAL_NUMBER_OF_NODES_XI(3),TOTAL_NUMBER_OF_NODES,NUMBER_OF_CORNER_NODES, &
      & TOTAL_NUMBER_OF_ELEMENTS,xi_idx,NUM_BASES,basis_idx,BASIS_NUMBER_OF_NODES
    INTEGER(INTG), ALLOCATABLE :: elementNodes(:),elementNodesUserNumbers(:)
    TYPE(BASIS_TYPE), POINTER :: basis
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: coordinateSystem
    TYPE(GeneratedMeshRegularType), POINTER :: regularMesh
    TYPE(INTERFACE_TYPE), POINTER :: interface
    TYPE(MeshComponentElementsType), POINTER :: meshElements
    TYPE(NODES_TYPE), POINTER :: nodes
    TYPE(REGION_TYPE), POINTER :: region
    TYPE(VARYING_STRING) :: localError

    NULLIFY(coordinateSystem)
    
    CALL Enters("GeneratedMeshRegularCreateFinish",err,error,*999)

    IF(ASSOCIATED(generatedMesh)) THEN
      CALL GeneratedMeshCoordinateSystemGet(generatedMesh,coordinateSystem,err,error,*999)
      CALL COORDINATE_SYSTEM_DIMENSION_GET(coordinateSystem,coordinateDimension,err,error,*999)
      CALL COORDINATE_SYSTEM_TYPE_GET(coordinateSystem,coordinateType,err,error,*999)
      region=>generatedMesh%region
      interface=>generatedMesh%interface
      regularMesh=>generatedMesh%regularMesh
      IF(ASSOCIATED(regularMesh)) THEN
        SELECT CASE(coordinateType)
        CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
          IF(ALLOCATED(regularMesh%bases)) THEN
            !Use first basis to get number of xi
            basis=>REGULAR_MESH%bases(1)%PTR
            SELECT CASE(BASIS%TYPE)
            CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE,BASIS_SIMPLEX_TYPE)
              !Calculate the sizes of a regular grid of elements with the appropriate number of basis nodes in each dimension of
              !the grid element
              totalNumberOfNodes=1
              gridNumberOfElements=1
              numberOfElementsXi=1
              numberOfBases=SIZE(regularMesh%bases)
              DO xiIdx=1,regularMesh%meshDimension
                !Set total number of nodes to corner nodes only
                totalNumberOfNodes=totalNumberOfNodes*(regularMesh%numberOfElementsXi(xiIdx)+1)
                numberOfElementsXi(xiIdx)=regularMesh%numberOfElementsXi(xiIdx)
                gridNumberOfElements=gridNumberOfElements*regularMesh%numberOfElementsXi(xiIdx)
              ENDDO !xiIdx
              numberOfCornerNodes=totalNumberOfNodes
              !Add extra nodes for each basis
              !Will end up with some duplicate nodes if bases have the same interpolation in one direction
              DO basisIdx=1,numberOfBases
                basis=>regularMesh%basis(basisIdx)%ptr
                basisNumberOfNodes=1
                DO xiIdx=1,regularMesh%meshDimension
                  basisNumberofNodes=BasisNumberOfNodes*((basis%NUMBER_OF_NODES_XIC(xiIdx)-1)* &
                    & regularMesh%numberOfElementsXi(xiIdx)+1)
                ENDDO !xiIdx
                totalNumberOfNodes=totalNumberOfNodes+basisNumberOfNodes-numberOfCornerNodes
              ENDDO
              !Compute the element factor i.e., the number of sub elements each grid element will be split into.
              IF(basis%TYPE==BASIS_LAGRANGE_HERMITE_TP_TYPE) THEN
                elementFactor=1
              ELSE
                SELECT CASE(regularMesh%meshDimension)
                CASE(1)
                  elementFactor=1
                CASE(2)
                  elementFactor=2
                CASE(3)
                  elementFactor=6
                CASE DEFAULT
                  localError="The mesh dimension dimension of "// &
                    & TRIM(NumberToVString(regularMesh%meshDimension,"*",err,error))//" is invalid."
                  CALL FlagError(localError,err,error,*999)
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
                CALL MESH_CREATE_START(meshUserNumber,region,regularMesh%meshDimension,generatedMesh%mesh,err,error,*999)
              ELSE
                CALL MESH_CREATE_START(meshUserNumber,interface,regularMesh%meshDimension,generatedMesh%mesh,err,error,*999)
              ENDIF
              !Set the number of mesh components
              CALL MESH_NUMBER_OF_COMPONENTS_SET(generatedMesh%mesh,numberOfBases,err,error,*999)
              !Create the elements
              CALL MESH_NUMBER_OF_ELEMENTS_SET(generatedMesh%mesh,totalNumberOfElements,err,error,*999)
              DO basisIdx=1,numberOfBases
                basis=>regularMesh%bases(basisIdx)%ptr
                !Get number of nodes in each xi direction for this basis
                DO xiIdx=1,regularMesh%meshDimension
                  totalNumberOfNodesXi(xiIdx)=(basis%NUMBER_OF_NODES_XIC(xiIdx)-1)*regularMesh%numberOfElementsXi(xiIdx)+1
                ENDDO
                NULLIFY(meshElements)
                CALL MESH_TOPOLOGY_ELEMENTS_CREATE_START(generatedMesh%mesh,basisIdx,basis,meshElements,err,error,*999)
                !Set the elements for the regular mesh
                IF(ALLOCATED(elementNodes)) DEALLOCATE(elementNodes)
                ALLOCATE(elementNodes(basis%NUMBER_OF_NODES),STAT=err)
                IF (ALLOCATED(elementNodesUserNumbers)) DEALLOCATE(elementNodesUserNumbers)
                ALLOCATE(elementNodesUserNumbers(basis%NUMBER_OF_NODES),STAT=err)
                IF(err/=0) CALL FlagError("Could not allocate element nodes.",err,error,*999)
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
                                        locaNodeIdx=localNodeidx+1
                                        elementNodes(localNodeIdx)=nodeIdx+(localNodeIdx1-1)+(localnodeIdx2-1)* &
                                          & totalNumberOfNodes(1)+(localNodeIdx3-1)*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                                      ENDDO !localNodeIdx1
                                    ENDDO !localnodeIdx2
                                  ENDDO !localNodeIdx3
                                ENDIF
                              ENDIF
                              CALL GENERATED_MESH_REGULAR_COMPONENT_NODES_TO_USER_NUMBERS(regularMesh%generatedMesh, &
                                & basisIdx,elementNodes,elementNodesUserNumbers,err,error,*999)
                              CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(elementIdx,meshElements,elementNodesUserNumbers, &
                                & err,error,*999)
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
                                CALL COMPONENT_NODES_TO_USER_NUMBERS(regularMesh%generatedMesh,basisIdx,elementNodes, &
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
                                  CALL COMPONENT_NODES_TO_USER_NUMBERS(regularMesh%generatedMesh,basisIdx,elementNodes, &
                                    & elementNodesUserNumbers,err,error,*999)
                                  CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(elementIdx,meshElements,elementNodesUserNumbers, &
                                    & err,error,*999)
                                  !Second sub-element
                                  elementIdx=(gridElementIdx-1)*elementFactor+2
                                  elementNodes(1)=nodeIdx
                                  elementNodes(2)=nodeIdx+1+totalNumberOfNodesXi(1)
                                  elementNodes(3)=nodeIdx+totalNumberOfNodesXi(1)
                                  CALL COMPONENT_NODES_TO_USER_NUMBERS(regularMesh%generatedMesh,basisIdx,elementNodes, &
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
                                  CALL COMPONENT_NODES_TO_USER_NUMBERS(regularMesh%generatedMesh,basisIdx,elementNodes, &
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
                                  CALL COMPONENT_NODES_TO_USER_NUMBERS(regularMesh%generatedMesh,basisIdx,elementNodes, &
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
                                  CALL COMPONENT_NODES_TO_USER_NUMBERS(regularMesh%generatedMesh,basisIdx,elementNodes, &
                                    & elementNodesUserNumbers,err,error,*999)
                                  CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(elementIdx,meshElements,elementNodesUserNumbers, &
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
                                  CALL COMPONENT_NODES_TO_USER_NUMBERS(regularMesh%generatedMesh,basisIdx,elementNodes, &
                                    & elementNodesUserNumbers,err,error,*999)
                                  CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(elementIdx,meshElements,elementNodesUserNumbers, &
                                    & err,error,*999)
                                CASE DEFAULT
                                  localError="The simplex basis interpolation order of "// &
                                    & TRIM(NumberToVString(basis%INTERPOLATION_ORDER(1),"*",err,error))// &
                                    & " is invalid."
                                  CALL FlagError(localError,err,error,*999)
                                END SELECT
                              CASE(3)
                                !Tetrahedra element
                                !Break the grid cube element into 6 tetrahedra (so that we have a break down the main diagonal of the
                                !cube in order to allow for the middle node in quadratics to be included). The 6 tetrahedra are
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
                                  CALL COMPONENT_NODES_TO_USER_NUMBERS(regularMesh%generatedMesh,basisIdx,elementNodes, &
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
                                  CALL COMPONENT_NODES_TO_USER_NUMBERS(regularMesh%generatedMesh,basisIdx,elementNodes, &
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
                                  CALL COMPONENT_NODES_TO_USER_NUMBERS(regularMesh%generatedMesh,basisIdx,elementNodes, &
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
                                  CALL COMPONENT_NODES_TO_USER_NUMBERS(regularMesh%generatedMesh,basisIdx,elementNodes, &
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
                                  CALL COMPONENT_NODES_TO_USER_NUMBERS(regularMesh%generatedMesh,basisIdx,elementNodes, &
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
                                  CALL COMPONENT_NODES_TO_USER_NUMBERS(regularMesh%generatedMesh,basisIdx,elementNodes, &
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
                                  CALL COMPONENT_NODES_TO_USER_NUMBERS(regularMesh%generatedMesh,basisIdx,elementNodes, &
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
                                  CALL COMPONENT_NODES_TO_USER_NUMBERS(regularMesh%generatedMesh,basisIdx,elementNodes, &
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
                                  CALL COMPONENT_NODES_TO_USER_NUMBERS(regularMesh%generatedMesh,basisIdx,elementNodes, &
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
                                  CALL COMPONENT_NODES_TO_USER_NUMBERS(regularMesh%generatedMesh,basisIdx,elementNodes, &
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
                                  CALL COMPONENT_NODES_TO_USER_NUMBERS(regularMesh%generatedMesh,basisIdx,elementNodes, &
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
                                  CALL COMPONENT_NODES_TO_USER_NUMBERS(regularMesh%generatedMesh,basisIdx,elementNodes, &
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
                                  CALL COMPONENT_NODES_TO_USER_NUMBERS(regularMesh%generatedMesh,basisIdx,elementNodes, &
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
                                  CALL COMPONENT_NODES_TO_USER_NUMBERS(regularMesh%generatedMesh,basisIdx,elementNodes, &
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
                                  CALL COMPONENT_NODES_TO_USER_NUMBERS(regularMesh%generatedMesh,basisIdx,elementNodes, &
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
                                  CALL COMPONENT_NODES_TO_USER_NUMBERS(regularMesh%generatedMesh,basisIdx,elementNodes, &
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
                                  CALL COMPONENT_NODES_TO_USER_NUMBERS(regularMesh%generatedMesh,basisIdx,elementNodes, &
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
                                  CALL COMPONENT_NODES_TO_USER_NUMBERS(regularMesh%generatedMesh,basisIdx,elementNodes, &
                                    & elementNodesUserNumbers,err,error,*999)
                                  CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(elementIdx,meshElements,elementNodesUserNumbers, &
                                    & err,error,*999)
                                CASE DEFAULT
                                  localError="The simplex basis interpolation order of "// &
                                    & TRIM(NumberToVString(basis%INTERPOLATION_ORDER(1),"*",err,error))// &
                                    & " is invalid."
                                  CALL FlagError(localError,err,error,*999)
                                END SELECT
                              CASE DEFAULT
                                localError="The simplex number of xi directions of "// &
                                  & TRIM(NumberToVString(basis%NUMBER_OF_XI,"*",err,error))// &
                                  & " is invalid."
                                CALL FlagError(localError,err,error,*999)
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
        CASE DEFAULT
          localError="The coordinate system type of "//TRIM(NumberToVString(coordinateType,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        CALL FlagError("Regular mesh is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Generated Mesh is not associated.",err,error,*999)
    ENDIF

    IF(ALLOCATED(elementNodes)) DEALLOCATE(elementNodes)

    CALL Exits("GeneratedMeshRegularCreateFinish")
    RETURN
    ! TODO invalidate other associations
999 IF(ALLOCATED(elementNodes)) DEALLOCATE(elementNodes)
    CALL Errors("GeneratedMeshRegularCreateFinish",err,error)
    CALL Exits("GeneratedMeshRegularCreateFinish")
    RETURN 1
  END SUBROUTINE GeneratedMeshRegularCreateFinish

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
    INTEGER(INTG) :: totalNumberOfNodes,totalNumberOfElements,meshDimension
    INTEGER(INTG) :: basisNumberOfNodes,cornerNumberOfNodes
    INTEGER(INTG) :: elementIdx1,elementIdx2,elementIdx3,localNodeIdx1,localNodeIdx2,localNodeIdx3,from1,from2,from3, &
      & localNodeIdx,elementIdx,mc
    INTEGER(INTG), ALLOCATABLE :: APEX_ELEMENT_NODES(:), WALL_ELEMENT_NODES(:)
    INTEGER(INTG), ALLOCATABLE :: APEX_ELEMENT_NODES_USER_NUMBERS(:), WALL_ELEMENT_NODES_USER_NUMBERS(:)
    INTEGER(INTG), ALLOCATABLE :: NIDX(:,:,:),CORNER_NODES(:,:,:),EIDX(:,:,:)
    INTEGER(INTG), ALLOCATABLE :: numberOfElementsXi(:)
    REAL(DP) :: DELTA(3),DELTAi(3)
    TYPE(BASIS_TYPE), POINTER :: basis1,basis2
    TYPE(GeneratedMeshEllipsoidType), POINTER :: ellipsoidMesh
    TYPE(MeshComponentElementsType), POINTER :: meshElements
    TYPE(NODES_TYPE), POINTER :: nodes
    TYPE(REGION_TYPE), POINTER :: region

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
                              IF(ALLOCATED(NIDX)) DEALLOCATE(NIDX)
                              IF(ALLOCATED(EIDX)) DEALLOCATE(EIDX)
                              IF(ALLOCATED(CORNER_NODES)) DEALLOCATE(CORNER_NODES)
                              CALL GENERATED_MESH_ELLIPSOID_BUILD_NODE_INDICES(numberOfElementsXi,basis1% &
                                  NUMBER_OF_NODES_XIC, ellipsoidMesh%ellipsoidExtent, totalNumberOfNodes, &
                                  totalNumberOfElements, NIDX,CORNER_NODES,EIDX,DELTA,DELTAi,err,error,*999)
                              IF(mc==1) THEN
                                CALL MESH_NUMBER_OF_ELEMENTS_SET(generatedMesh%MESH,totalNumberOfElements, &
                                  & err,error,*999)
                              ENDIF

                              !Create the default node set
                              !TODO we finish create after the nodes are initialised?

                              NULLIFY(meshElements)
                              CALL MESH_TOPOLOGY_ELEMENTS_CREATE_START(generatedMesh%MESH,mc/2+1,basis1,meshElements, &
                                  ERR, ERROR,*999)
                              !Set the elements for the ellipsoid mesh
                              IF(ALLOCATED(WALL_ELEMENT_NODES)) DEALLOCATE(WALL_ELEMENT_NODES)
                              IF(ALLOCATED(APEX_ELEMENT_NODES)) DEALLOCATE(APEX_ELEMENT_NODES)
                              IF(ALLOCATED(WALL_ELEMENT_NODES_USER_NUMBERS)) DEALLOCATE(WALL_ELEMENT_NODES_USER_NUMBERS)
                              IF(ALLOCATED(APEX_ELEMENT_NODES_USER_NUMBERS)) DEALLOCATE(APEX_ELEMENT_NODES_USER_NUMBERS)
                              ALLOCATE(WALL_ELEMENT_NODES(basis1%NUMBER_OF_NODES),STAT=err)
                              IF(err/=0) CALL FlagError("Could not allocate wall element nodes.",err,error,*999)
                              ALLOCATE(APEX_ELEMENT_NODES(basis2%NUMBER_OF_NODES),STAT=err)
                              IF(err/=0) CALL FlagError("Could not allocate apex element nodes.",err,error,*999)
                              ALLOCATE(WALL_ELEMENT_NODES_USER_NUMBERS(basis1%NUMBER_OF_NODES),STAT=err)
                              IF(err/=0) CALL FlagError("Could not allocate wall element nodes.",err,error,*999)
                              ALLOCATE(APEX_ELEMENT_NODES_USER_NUMBERS(basis2%NUMBER_OF_NODES),STAT=err)
                              IF(err/=0) CALL FlagError("Could not allocate apex element nodes.",err,error,*999)
                              ! calculate element topology (nodes per each element)
                              ! the idea is to translate given (r,theta,z) to NIDX equivalents, which include interior nodes
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
                                    APEX_ELEMENT_NODES(localNodeIdx)=NIDX(localNodeIdx1,localNodeIdx2,localNodeIdx3)
                                    DO localNodeIdx2=from2+1,from2+basis2%NUMBER_OF_NODES_XIC(2)-1
                                      DO localNodeIdx1=from1,from1+basis2%NUMBER_OF_NODES_XIC(1)-1
                                        localNodeIdx=localNodeIdx+1
                                        ! circumferential loop-around
                                        IF(localNodeIdx1>SIZE(NIDX,1)) THEN
                                          APEX_ELEMENT_NODES(localNodeIdx)=NIDX(1,localNodeIdx2,localNodeIdx3)
                                        ELSE
                                          APEX_ELEMENT_NODES(localNodeIdx)=NIDX(localNodeIdx1,localNodeIdx2,localNodeIdx3)
                                        ENDIF
                                      ENDDO ! localNodeIdx1
                                    ENDDO ! localNodeIdx2
                                  ENDDO ! localNodeIdx3
                                  elementIdx=elementIdx+1
                                  CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_SET(ne,meshElements,basis2,err,error,*999)
                                  CALL COMPONENT_NODES_TO_USER_NUMBERS(ellipsoidMesh%generatedMesh,mc,APEX_ELEMENT_NODES, &
                                      & APEX_ELEMENT_NODES_USER_NUMBERS,err,error,*999)
                                  CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(ne,meshElements, &
                                    APEX_ELEMENT_NODES_USER_NUMBERS,err,error,*999)
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
                                          IF(localNodeIdx1>SIZE(NIDX,1)) THEN
                                            WALL_ELEMENT_NODES(localNodeIdx)=NIDX(1,localNodeIdx2,localNodeIdx3)
                                          ELSE
                                            WALL_ELEMENT_NODES(localNodeIdx)=NIDX(localNodeIdx1,localNodeIdx2,localNodeIdx3)
                                          ENDIF
                                        ENDDO ! localNodeIdx1
                                      ENDDO ! localNodeIdx2
                                    ENDDO ! localNodeIdx3
                                    elementIdx=elementIdx+1
                                    CALL COMPONENT_NODES_TO_USER_NUMBERS(ellipsoidMesh%generatedMesh,mc,WALL_ELEMENT_NODES, &
                                        & WALL_ELEMENT_NODES_USER_NUMBERS,err,error,*999)
                                    CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(ne,meshElements, &
                                        & WALL_ELEMENT_NODES_USER_NUMBERS,err,error,*999)
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
                  CALL FlagError("The number of dimensions of the given regular mesh does not match the size of &
                      &the origin.",err,error,*999)
                ENDIF
              ELSE
                CALL FlagError("Ellipsoid mesh requires a 3 dimensional coordinate system.",err,error,*999)
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
        CALL FlagError("Ellipsoid mesh is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Generated Mesh is not associated.",err,error,*999)
    ENDIF

    CALL Exits("GeneratedMeshEllipsoidCreateFinish")
    RETURN
    ! TODO invalidate other associations
999 IF(ALLOCATED(NIDX)) DEALLOCATE(NIDX)
    IF(ALLOCATED(EIDX)) DEALLOCATE(EIDX)
    IF(ALLOCATED(CORNER_NODES)) DEALLOCATE(CORNER_NODES)
    IF(ALLOCATED(numberOfElementsXi)) DEALLOCATE(numberOfElementsXi)
    IF(ALLOCATED(WALL_ELEMENT_NODES)) DEALLOCATE(WALL_ELEMENT_NODES)
    IF(ALLOCATED(APEX_ELEMENT_NODES)) DEALLOCATE(APEX_ELEMENT_NODES)
    CALL Errors("GeneratedMeshEllipsoidCreateFinish",err,error)
    CALL Exits("GeneratedMeshEllipsoidCreateFinish")
    RETURN 1
  END SUBROUTINE GeneratedMeshEllipsoidCreateFinish
  !
  !================================================================================================================================
  !

  !>Start to create the regular generated mesh type
  SUBROUTINE GENERATED_MESH_CYLINDER_CREATE_FINISH(GENERATED_MESH,MESH_USER_NUMBER,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: GENERATED_MESH !<A pointer to the generated mesh
    INTEGER(INTG), INTENT(IN) :: MESH_USER_NUMBER !<The user number for the mesh to generate.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(GeneratedMeshCylinderType), POINTER :: CYLINDER_MESH
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
    TYPE(MeshComponentElementsType), POINTER :: MESH_ELEMENTS

    CALL Enters("GENERATED_MESH_CYLINDER_CREATE_FINISH",err,error,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      CYLINDER_MESH=>GENERATED_MESH%cylinderMesh
      IF(ASSOCIATED(CYLINDER_MESH)) THEN
        REGION=>GENERATED_MESH%region
        IF(ASSOCIATED(REGION)) THEN
          IF(ASSOCIATED(REGION%COORDINATE_SYSTEM)) THEN
            !TODO is regular type only for COORDINATE_RECTANGULAR_CARTESIAN_TYPE?
            !If that, should we use IF rather than select?
            SELECT CASE(REGION%COORDINATE_SYSTEM%TYPE)
            CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
              !Determine the coordinate system and create the regular mesh for that system
              CYLINDER_MESH%MESH_DIMENSION=REGION%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS
              NUMBER_OF_DIMENSIONS=CYLINDER_MESH%MESH_DIMENSION
              IF(NUMBER_OF_DIMENSIONS==3) THEN ! hard-coded for 3D only
                IF(.NOT.ALLOCATED(CYLINDER_MESH%ORIGIN)) THEN
                  ALLOCATE(CYLINDER_MESH%ORIGIN(NUMBER_OF_DIMENSIONS),STAT=err)
                  IF(err/=0) CALL FlagError("Could not allocate origin.",err,error,*999)
                  CYLINDER_MESH%ORIGIN=0.0_DP
                ENDIF
                IF(SIZE(CYLINDER_MESH%ORIGIN)==CYLINDER_MESH%MESH_DIMENSION) THEN
                  IF(SIZE(CYLINDER_MESH%cylinderExtent)==CYLINDER_MESH%MESH_DIMENSION) THEN
                    IF(ALLOCATED(CYLINDER_MESH%bases)) THEN
                      ALLOCATE(NUMBER_ELEMENTS_XI(SIZE(CYLINDER_MESH%numberOfElementsXi)),STAT=err)
                      IF(err/=0) CALL FlagError("Could not allocate number of elements xi.",err,error,*999)
                      NUMBER_ELEMENTS_XI(1:SIZE(CYLINDER_MESH%numberOfElementsXi))= &
                        & CYLINDER_MESH%numberOfElementsXi(1:SIZE(CYLINDER_MESH%numberOfElementsXi))
                      CALL MESH_CREATE_START(MESH_USER_NUMBER,GENERATED_MESH%region,SIZE(NUMBER_ELEMENTS_XI,1), &
                        & GENERATED_MESH%MESH,err,error,*999)
                      CALL MESH_NUMBER_OF_COMPONENTS_SET(GENERATED_MESH%MESH,SIZE(CYLINDER_MESH%bases),err,error,*999)
                      !Calculate number of nodes
                      CORNER_NUMBER_OF_NODES=(NUMBER_ELEMENTS_XI(3)+1)*NUMBER_ELEMENTS_XI(2)*(NUMBER_ELEMENTS_XI(1)+1)
                      TOTAL_NUMBER_OF_NODES=CORNER_NUMBER_OF_NODES
                      DO basis_idx=1,SIZE(CYLINDER_MESH%bases)
                        BASIS=>CYLINDER_MESH%bases(basis_idx)%PTR
                        IF(ASSOCIATED(BASIS)) THEN
                          BASIS_NUMBER_OF_NODES=((BASIS%NUMBER_OF_NODES_XIC(3)-1)*NUMBER_ELEMENTS_XI(3)+1)* &
                              & ((BASIS%NUMBER_OF_NODES_XIC(2)-1)*NUMBER_ELEMENTS_XI(2))* &
                              & ((BASIS%NUMBER_OF_NODES_XIC(1)-1)*NUMBER_ELEMENTS_XI(1)+1)
                          TOTAL_NUMBER_OF_NODES=TOTAL_NUMBER_OF_NODES+BASIS_NUMBER_OF_NODES-CORNER_NUMBER_OF_NODES
                        ELSE
                          CALL FlagError("Basis is not associated.",err,error,*999)
                        ENDIF
                      ENDDO
                      NULLIFY(NODES)
                      CALL NODES_CREATE_START(REGION,TOTAL_NUMBER_OF_NODES,NODES,err,error,*999)
                      !Finish the nodes creation
                      CALL NODES_CREATE_FINISH(NODES,err,error,*999)
                      !Set the total number of elements
                      TOTAL_NUMBER_OF_ELEMENTS=NUMBER_ELEMENTS_XI(1)*NUMBER_ELEMENTS_XI(2)*NUMBER_ELEMENTS_XI(3)
                      CALL MESH_NUMBER_OF_ELEMENTS_SET(GENERATED_MESH%MESH,TOTAL_NUMBER_OF_ELEMENTS,err,error,*999)
                      DO basis_idx=1,SIZE(CYLINDER_MESH%bases)
                        BASIS=>CYLINDER_MESH%bases(basis_idx)%PTR
                        IF(ASSOCIATED(BASIS)) THEN
                          SELECT CASE(BASIS%TYPE)
                          CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
                            IF(BASIS%NUMBER_OF_XI==SIZE(NUMBER_ELEMENTS_XI,1)) THEN
                              IF(.NOT.ALL(NUMBER_ELEMENTS_XI>0)) &
                                & CALL FlagError("Must have 1 or more elements in all directions.",err,error,*999)
                              IF(NUMBER_ELEMENTS_XI(2)<3) &
                                CALL FlagError("Need >2 elements around the circumferential direction.",err,error,*999)
                              IF(.NOT.ALL(BASIS%COLLAPSED_XI==BASIS_NOT_COLLAPSED))  &
                                & CALL FlagError("Degenerate (collapsed) basis not implemented.",err,error,*999)
                              !Calculate nodes and element sizes
                              IF(ALLOCATED(NIDX)) DEALLOCATE(NIDX)
                              IF(ALLOCATED(EIDX)) DEALLOCATE(EIDX)
                              CALL GENERATED_MESH_CYLINDER_BUILD_NODE_INDICES(NUMBER_ELEMENTS_XI,BASIS%NUMBER_OF_NODES_XIC, &
                                & CYLINDER_MESH%cylinderExtent, TOTAL_NUMBER_OF_NODES,TOTAL_NUMBER_OF_ELEMENTS, &
                                & NIDX,EIDX,DELTA,DELTAi,err,error,*999)
                              !Set the elements for the cylinder mesh
                              IF(ALLOCATED(ELEMENT_NODES)) DEALLOCATE(ELEMENT_NODES)
                              IF(ALLOCATED(ELEMENT_NODES_USER_NUMBERS)) DEALLOCATE(ELEMENT_NODES_USER_NUMBERS)
                              ALLOCATE(ELEMENT_NODES_USER_NUMBERS(BASIS%NUMBER_OF_NODES),STAT=err)
                              ALLOCATE(ELEMENT_NODES(BASIS%NUMBER_OF_NODES),STAT=err)
                              IF(err/=0) CALL FlagError("Could not allocate element nodes.",err,error,*999)
                              !Create the elements
                              NULLIFY(MESH_ELEMENTS)
                              CALL MESH_TOPOLOGY_ELEMENTS_CREATE_START(GENERATED_MESH%MESH,basis_idx,BASIS,MESH_ELEMENTS, &
                                  & err,error,*999)
                              ! calculate element topology (nodes per each element)
                              ! the idea is to translate given (r,theta,z) to NIDX equivalents, which include interior nodes
                              elementIdx=0
                              DO ne3=1,NUMBER_ELEMENTS_XI(3)
                                from3=NINT(DELTA(3)*(ne3-1)/DELTAi(3)+1)
                                DO ne2=1,NUMBER_ELEMENTS_XI(2)
                                  from2=NINT(DELTA(2)*(ne2-1)/DELTAi(2)+1)
                                  DO ne1=1,NUMBER_ELEMENTS_XI(1)
                                    from1=NINT(DELTA(1)*(ne1-1)/DELTAi(1)+1)
                                    nn=0
                                    ! number of nodes in an element is dependent on basis used
                                    DO nn3=from3,from3+BASIS%NUMBER_OF_NODES_XIC(3)-1
                                      DO nn2=from2,from2+BASIS%NUMBER_OF_NODES_XIC(2)-1
                                        DO nn1=from1,from1+BASIS%NUMBER_OF_NODES_XIC(1)-1
                                          nn=nn+1
                                          ! compensate for circumferential loop-around
                                          IF(nn2>SIZE(NIDX,2)) THEN
                                            ! DEBUG: little check here
                                            IF(nn2>SIZE(NIDX,2)+1) CALL FlagError("NIDX needs debugging",err,error,*999)
                                            ELEMENT_NODES(nn)=NIDX(nn1,1,nn3)
                                          ELSE
                                            ELEMENT_NODES(nn)=NIDX(nn1,nn2,nn3)
                                          ENDIF
                                        ENDDO ! nn1
                                      ENDDO ! nn2
                                    ENDDO ! nn3
                                    elementIdx=elementIdx+1
                                    CALL COMPONENT_NODES_TO_USER_NUMBERS(CYLINDER_MESH%GENERATED_MESH,basis_idx,ELEMENT_NODES, &
                                        & ELEMENT_NODES_USER_NUMBERS,err,error,*999)
                                    CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(ne,MESH_ELEMENTS,ELEMENT_NODES_USER_NUMBERS, &
                                        & err,error,*999)
                                  ENDDO ! ne1
                                ENDDO ! ne2
                              ENDDO ! ne3
                              CALL MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH(MESH_ELEMENTS,err,error,*999)
                            ELSE
                              CALL FlagError("The number of xi directions of the given basis does not match the size of &
                                &the number of elements for the mesh.",err,error,*999)
                            ENDIF
                          CASE(BASIS_SIMPLEX_TYPE)
                            CALL FlagError("Cylinder meshes with simplex basis types is not implemented.",err,error,*999)
                          CASE DEFAULT
                            CALL FlagError("Basis type is either invalid or not implemented.",err,error,*999)
                          END SELECT
                        ELSE
                          CALL FlagError("Basis is not associated.",err,error,*999)
                        ENDIF
                      ENDDO
                      !Finish the mesh
                      CALL MESH_CREATE_FINISH(GENERATED_MESH%MESH,err,error,*999)
                    ELSE
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

    CALL Exits("GENERATED_MESH_CYLINDER_CREATE_FINISH")
    RETURN
    ! TODO invalidate other associations
999 IF(ALLOCATED(NIDX)) DEALLOCATE(NIDX)
    IF(ALLOCATED(EIDX)) DEALLOCATE(EIDX)
    IF(ALLOCATED(NUMBER_ELEMENTS_XI)) DEALLOCATE(NUMBER_ELEMENTS_XI)
    IF(ALLOCATED(ELEMENT_NODES)) DEALLOCATE(ELEMENT_NODES)
    CALL Errors("GENERATED_MESH_CYLINDER_CREATE_FINISH",err,error)
    CALL Exits("GENERATED_MESH_CYLINDER_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE GENERATED_MESH_CYLINDER_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Finalise the cylinder mesh type
  SUBROUTINE GENERATED_MESH_CYLINDER_FINALISE(CYLINDER_MESH,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshCylinderType), POINTER :: CYLINDER_MESH !<A pointer to the generated mesh
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    CALL Enters("GENERATED_MESH_CYLINDER_FINALISE",err,error,*999)

    IF(ASSOCIATED(CYLINDER_MESH)) THEN
      IF(ALLOCATED(CYLINDER_MESH%ORIGIN)) DEALLOCATE(CYLINDER_MESH%ORIGIN)
      IF(ALLOCATED(CYLINDER_MESH%cylinderExtent)) DEALLOCATE(CYLINDER_MESH%cylinderExtent)
      IF(ALLOCATED(CYLINDER_MESH%numberOfElementsXi)) DEALLOCATE(CYLINDER_MESH%numberOfElementsXi)
      IF(ALLOCATED(CYLINDER_MESH%bases)) DEALLOCATE(CYLINDER_MESH%bases)
      DEALLOCATE(CYLINDER_MESH)
    ENDIF

    CALL Exits("GENERATED_MESH_CYLINDER_FINALISE")
    RETURN
    ! TODO invalidate other associations
999 CALL Errors("GENERATED_MESH_CYLINDER_FINALISE",err,error)
    CALL Exits("GENERATED_MESH_CYLINDER_FINALISE")
    RETURN 1
  END SUBROUTINE GENERATED_MESH_CYLINDER_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialise the cylinder generated mesh type
  SUBROUTINE GENERATED_MESH_CYLINDER_INITIALISE(GENERATED_MESH,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: GENERATED_MESH !<A pointer to the generated mesh
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    CALL Enters("GENERATED_MESH_CYLINDER_INITIALISE",err,error,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      IF(ASSOCIATED(GENERATED_MESH%cylinderMesh)) THEN
        CALL FlagError("Cylinder mesh is already associated for this generated mesh.",err,error,*998)
      ELSE
        ALLOCATE(GENERATED_MESH%cylinderMesh,STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate cylinder generated mesh.",err,error,*999)
        GENERATED_MESH%cylinderMesh%GENERATED_MESH=>GENERATED_MESH
        GENERATED_MESH%generatedType=GENERATED_MESH_CYLINDER_MESH_TYPE
      ENDIF
    ELSE
      CALL FlagError("Generated mesh is not associated.",err,error,*998)
    ENDIF

    CALL Exits("GENERATED_MESH_CYLINDER_INITIALISE")
    RETURN
999 CALL GENERATED_MESH_CYLINDER_FINALISE(GENERATED_MESH%cylinderMesh,dummyErr,dummyError,*998)
998 CALL Errors("GENERATED_MESH_CYLINDER_INITIALISE",err,error)
    CALL Exits("GENERATED_MESH_CYLINDER_INITIALISE")
    RETURN 1
  END SUBROUTINE GENERATED_MESH_CYLINDER_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalise the regular mesh type
  SUBROUTINE GENERATED_MESH_REGULAR_FINALISE(REGULAR_MESH,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshRegularType), POINTER :: REGULAR_MESH !<A pointer to the generated mesh
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    CALL Enters("GENERATED_MESH_REGULAR_FINALISE",err,error,*999)

    IF(ASSOCIATED(REGULAR_MESH)) THEN
      IF(ALLOCATED(REGULAR_MESH%ORIGIN)) DEALLOCATE(REGULAR_MESH%ORIGIN)
      IF(ALLOCATED(REGULAR_MESH%maximumExtent)) DEALLOCATE(REGULAR_MESH%maximumExtent)
      IF(ALLOCATED(REGULAR_MESH%numberOfElementsXi)) DEALLOCATE(REGULAR_MESH%numberOfElementsXi)
      IF(ALLOCATED(REGULAR_MESH%baseVectors)) DEALLOCATE(REGULAR_MESH%baseVectors)
      IF(ALLOCATED(REGULAR_MESH%bases)) DEALLOCATE(REGULAR_MESH%bases)
      DEALLOCATE(REGULAR_MESH)
    ENDIF

    CALL Exits("GENERATED_MESH_REGULAR_FINALISE")
    RETURN
    ! TODO invalidate other associations
999 CALL Errors("GENERATED_MESH_REGULAR_FINALISE",err,error)
    CALL Exits("GENERATED_MESH_REGULAR_FINALISE")
    RETURN 1
  END SUBROUTINE GENERATED_MESH_REGULAR_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialise the regular generated mesh type
  SUBROUTINE GENERATED_MESH_REGULAR_INITIALISE(GENERATED_MESH,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: GENERATED_MESH !<A pointer to the generated mesh
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
    
    NULLIFY(COORDINATE_SYSTEM)

    CALL Enters("GENERATED_MESH_REGULAR_INITIALISE",err,error,*998)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      IF(ASSOCIATED(GENERATED_MESH%regularMesh)) THEN
        CALL FlagError("Regular mesh is already associated for this generated mesh.",err,error,*998)
      ELSE
      CALL GeneratedMeshCoordinateSystemGet(GENERATED_MESH,COORDINATE_SYSTEM,err,error,*999)
        ALLOCATE(GENERATED_MESH%regularMesh,STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate regular generated mesh.",err,error,*999)
        GENERATED_MESH%regularMesh%GENERATED_MESH=>GENERATED_MESH
        ALLOCATE(generatedMesh%regularMesh%origin(coordinateDimension),STAT=err)
        IF(err/=) CALL FlagError("Could not allocate origin.",err,error,*998)
        generatedMesh%regularMesh%origin=0.0_DP
        ALLOCATE(generatedMesh%regularMesh%maximumExtent(coordinateDimension),STAT=err)
        IF(err/=) CALL FlagError("Could not allocate maximum extent.",err,error,*998)
        generatedMesh%regularMesh%maximumExtent=1.0_DP
        GENERATED_MESH%generatedType=GENERATED_MESH_REGULAR_MESH_TYPE
      ENDIF
    ELSE
      CALL FlagError("Generated mesh is not associated.",err,error,*998)
    ENDIF

    CALL Exits("GENERATED_MESH_REGULAR_INITIALISE")
    RETURN
999 CALL GENERATED_MESH_REGULAR_FINALISE(GENERATED_MESH%regularMesh,dummyErr,dummyError,*998)
998 CALL Errors("GENERATED_MESH_REGULAR_INITIALISE",err,error)
    CALL Exits("GENERATED_MESH_REGULAR_INITIALISE")
    RETURN 1
  END SUBROUTINE GENERATED_MESH_REGULAR_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalise ellipsoid mesh type
  SUBROUTINE GENERATED_MESH_ELLIPSOID_FINALISE(ELLIPSOID_MESH,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshEllipsoidType), POINTER :: ELLIPSOID_MESH !<A pointer to the generated mesh
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    CALL Enters("GENERATED_MESH_ELLIPSOID_FINALISE",err,error,*999)

    IF(ASSOCIATED(ELLIPSOID_MESH)) THEN
      IF(ALLOCATED(ELLIPSOID_MESH%ORIGIN)) DEALLOCATE(ELLIPSOID_MESH%ORIGIN)
      IF(ALLOCATED(ELLIPSOID_MESH%ellipsoidExtent)) DEALLOCATE(ELLIPSOID_MESH%ellipsoidExtent)
      IF(ALLOCATED(ELLIPSOID_MESH%numberOfElementsXi)) DEALLOCATE(ELLIPSOID_MESH%numberOfElementsXi)
      IF(ALLOCATED(ELLIPSOID_MESH%bases)) DEALLOCATE(ELLIPSOID_MESH%bases)
      DEALLOCATE(ELLIPSOID_MESH)
    ENDIF

    CALL Exits("GENERATED_MESH_ELLIPSOID_FINALISE")
    RETURN
    ! TODO invalidate other associations
999 CALL Errors("GENERATED_MESH_ELLIPSOID_FINALISE",err,error)
    CALL Exits("GENERATED_MESH_ELLIPSOID_FINALISE")
    RETURN 1
  END SUBROUTINE GENERATED_MESH_ELLIPSOID_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialise the ellipsoid generated mesh type
  SUBROUTINE GENERATED_MESH_ELLIPSOID_INITIALISE(GENERATED_MESH,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: GENERATED_MESH !<A pointer to the generated mesh
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    CALL Enters("GENERATED_MESH_ELLIPSOID_INITIALISE",err,error,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      IF(ASSOCIATED(GENERATED_MESH%ellipsoidMesh)) THEN
        CALL FlagError("Ellipsoid mesh is already associated for this generated mesh.",err,error,*998)
      ELSE
        ALLOCATE(GENERATED_MESH%ellipsoidMesh,STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate ellipsoid generated mesh.",err,error,*999)
        GENERATED_MESH%ellipsoidMesh%GENERATED_MESH=>GENERATED_MESH
        GENERATED_MESH%generatedType=GENERATED_MESH_ELLIPSOID_MESH_TYPE
      ENDIF
    ELSE
      CALL FlagError("Generated mesh is not associated.",err,error,*998)
    ENDIF

    CALL Exits("GENERATED_MESH_ELLIPSOID_INITIALISE")
    RETURN
999 CALL GENERATED_MESH_ELLIPSOID_FINALISE(GENERATED_MESH%ellipsoidMesh,dummyErr,dummyError,*998)
998 CALL Errors("GENERATED_MESH_ELLIPSOID_INITIALISE",err,error)
    CALL Exits("GENERATED_MESH_ELLIPSOID_INITIALISE")
    RETURN 1
  END SUBROUTINE GENERATED_MESH_ELLIPSOID_INITIALISE

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

    CALL Enters("GeneratedMeshTypeGet",err,error,*999)

    IF(ASSOCIATED(generatedMesh)) THEN
      meshType=generatedMesh%generatedType
    ELSE
      CALL FlagError("Generated mesh is not associated.",err,error,*999)
    ENDIF

    CALL Exits("GeneratedMeshTypeGet")
    RETURN
999 CALL Errors("GeneratedMeshTypeGet",err,error)
    CALL Exits("GeneratedMeshTypeGet")
    RETURN 1
    
  END SUBROUTINE GeneratedMeshTypeGet

  !
  !================================================================================================================================
  !

  !>Sets/changes the type of a generated mesh. \see OPENCMISS::CMISSGeneratedMeshTypeSet
  SUBROUTINE GeneratedMeshTypeSet(generatedMesh,meshType,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: gneratedMesh !<A pointer to the generated mesh to set the type of
    INTEGER(INTG), INTENT(IN) :: meshType !<The type of mesh to generate \see GeneratedMeshRoutines_GeneratedMeshTypes,GeneratedMeshRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: oldMeshType
    TYPE(VARYING_STRING) :: localError

    CALL Enters("GeneratedMeshTypeSet",err,error,*999)

    IF(ASSOCIATED(generatedMesh)) THEN
      IF(generatedMesh%generatedMeshFinished) THEN
        CALL FlagError("Generated mesh has already been finished.",err,error,*999)
      ELSE
        oldMeshType=generatedMesh%generatedType
        IF(oldMeshType/=meshType) THEN
          !Initialise the new generated mesh type
          SELECT CASE(meshType)
          CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
            CALL GENERATED_MESH_REGULAR_INITIALISE(generatedMesh,err,error,*999)
          CASE(GENERATED_MESH_POLAR_MESH_TYPE)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
            CALL GENERATED_MESH_CYLINDER_INITIALISE(generatedMesh,err,error,*999)
          CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
            CALL GENERATED_MESH_ELLIPSOID_INITIALISE(generatedMesh,err,error,*999)
          CASE DEFAULT
            localError="The specified generated mesh mesh type of "//TRIM(NumberToVString(meshType,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
          !Finalise the new generated mesh type
          SELECT CASE(oldMeshType)
          CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
            CALL GENERATED_MESH_REGULAR_FINALISE(generatedMesh%regularMesh,err,error,*999)
          CASE(GENERATED_MESH_POLAR_MESH_TYPE)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
            CALL GENERATED_MESH_CYLINDER_FINALISE(generatedMesh%cylinderMesh,err,error,*999)
          CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
            CALL GENERATED_MESH_ELLIPSOID_FINALISE(generatedMesh%ellipsoidMesh,err,error,*999)
          CASE DEFAULT
            localError="The generated mesh mesh type of "//TRIM(NumberToVString(oldMeshType,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Generated mesh is not associated.",err,error,*999)
    ENDIF

    CALL Exits("GeneratedMeshTypeSet")
    RETURN
999 CALL Errors("GeneratedMeshTypeSet",err,error)
    CALL Exits("GeneratedMeshTypeSet")
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

    CALL Enters("GeneratedMeshUserNumberFindGeneric",err,error,*999)

    IF(ASSOCIATED(generatedMeshes)) THEN
      IF(ASSOCIATED(generatedMesh)) THEN
        CALL FlagError("Generated mesh is already associated.",err,error,*999)
      ELSE
        NULLIFY(generatedMesh)
        generatedMeshIdx=1
        DO WHILE(generatedMeshIdx<=generatedMeshes%numberOfGeneratedMeshes.AND..NOT.ASSOCIATED(generatedMesh))
          IF(generatedMeshes%generatedMeshes(generatedMeshIdx)%PTR%userNumber==userNumber) THEN
            generatedMesh=>generatedMeshes%generatedMeshes(generatedMeshIdx)%PTR
            EXIT
          ELSE
            generatedMeshIdx=generatedMeshIdx+1
          ENDIF
        ENDDO
      ENDIF
    ELSE
      CALL FlagError("Generated meshes is not associated.",err,error,*999)
    ENDIF

    CALL Exits("GeneratedMeshUserNumberFindGeneric")
    RETURN
999 CALL Errors("GeneratedMeshUserNumberFindGeneric",err,error)
    CALL Exits("GeneratedMeshUserNumberFindGeneric")
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

    CALL Enters("GeneratedMeshUserNumberFindInterface",err,error,*999)

    IF(ASSOCIATED(interface)) THEN
      CALL GeneratedMeshUserNumberFindGeneric(userNumber,interface%GENERATED_MESHES,generatedMesh,err,error,*999)
    ELSE
      CALL FlagError("Interface is not associated.",err,error,*999)
    ENDIF

    CALL Exits("GeneratedMeshUserNumberFindInterface")
    RETURN
999 CALL Errors("GeneratedMeshUserNumberFindInterface",err,error)
    CALL Exits("GeneratedMeshUserNumberFindInterface")
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

    CALL Enters("GeneratedMeshUserNumberFindRegion",err,error,*999)

    IF(ASSOCIATED(region)) THEN
      CALL GeneratedMeshUserNumberFindGeneric(userNumber,region%GENERATED_MESHES,generatedMesh,err,error,*999)
    ELSE
      CALL FlagError("Region is not associated.",err,error,*999)
    ENDIF

    CALL Exits("GeneratedMeshUserNumberFindRegion")
    RETURN
999 CALL Errors("GeneratedMeshUserNumberFindRegion",err,error)
    CALL Exits("GeneratedMeshUserNumberFindRegion")
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

    CALL Enters("GeneratedMeshesFinalise",err,error,*999)

    IF(ASSOCIATED(generatedMeshes)) THEN
      DO WHILE(generatedMeshes%numberOfGeneratedMeshes>0)
        generatedMesh=>generatedMeshes%generatedMeshes(1)%ptr
        CALL GeneratedMeshDestroy(generatedMesh,err,error,*999)
      ENDDO !generatedMeshIdx
      DEALLOCATE(generatedMeshes)
    ENDIF

    CALL Exits("GeneratedMeshesFinalise")
    RETURN
999 CALL Errors("GeneratedMeshesFinalise",err,error)
    CALL Exits("GeneratedMeshesFinalise")
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

    CALL Enters("GeneratedMeshesInitialiseGeneric",err,error,*998)

    IF(ASSOCIATED(generatedMeshes)) THEN
      CALL FlagError("Generated meshes is already associated.",err,error,*998)
    ELSE
      ALLOCATE(generatedMeshes,STAT=err)
      IF(err/=0) CALL FlagError("Generated meshes is not associated.",err,error,*999)
      generatedMeshes%numberOfGeneratedMeshes=0
      NULLIFY(generatedMeshes%generatedMeshes)
      NULLIFY(generatedMeshes%region)
      NULLIFY(generatedMeshes%interface)
    ENDIF

    CALL Exits("GeneratedMeshesInitialiseGeneric")
    RETURN
999 CALL GeneratedMeshesFinalise(generatedMeshes,dummyErr,dummyError,*998)
998 CALL Errors("GeneratedMeshesInitialiseGeneric",err,error)
    CALL Exits("GeneratedMeshesInitialiseGeneric")
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

    CALL Enters("GeneratedMeshesInitialiseInterface",err,error,*999)

    IF(ASSOCIATED(interface)) THEN
      IF(ASSOCIATED(interface%GENERATED_MESHES)) THEN
        CALL FlagError("Interface generated meshes is already associated.",err,error,*999)
      ELSE
        CALL GeneratedMeshesInitialiseGeneric(interface%GENERATED_MESHES,err,error,*999)
        interface%GENERATED_MESHES%interface=>interface
      ENDIF
    ELSE
      CALL FlagError("Interface is not associated.",err,error,*999)
    ENDIF

    CALL Exits("GeneratedMeshesInitialiseInterface")
    RETURN
999 CALL Errors("GeneratedMeshesInitialiseInterface",err,error)
    CALL Exits("GeneratedMeshesInitialiseInterface")
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

    CALL Enters("GeneratedMeshesInitialiseRegion",err,error,*999)

    IF(ASSOCIATED(region)) THEN
      IF(ASSOCIATED(region%GENERATED_MESHES)) THEN
        CALL FlagError("Region generated meshes is already associated.",err,error,*999)
      ELSE
        CALL GeneratedMeshesInitialiseGeneric(region%GENERATED_MESHES,err,error,*999)
        region%GENERATED_MESHES%region=>region
      ENDIF
    ELSE
      CALL FlagError("Region is not associated.",err,error,*999)
    ENDIF

    CALL Exits("GeneratedMeshesInitialiseRegion")
    RETURN
999 CALL Errors("GeneratedMeshesInitialiseRegion",err,error)
    CALL Exits("GeneratedMeshesInitialiseRegion")
    RETURN 1
  END SUBROUTINE GeneratedMeshesInitialiseRegion

  !
  !================================================================================================================================
  !

  !>Updates the geometric field parameters from the initial nodal positions of the mesh. \see OPENCMISS::CMISSGeneratedMeshGeometricParametersCalculate
  SUBROUTINE GeneratedMeshGeometricParametersCalculate(generatedMesh,geometricField,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<The mesh which is generated by the generated mesh \todo is this necessary???
    TYPE(FIELD_TYPE), POINTER :: geometricField !<A pointer to the field to update the geometric parameters for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(REGION_TYPE), POINTER :: meshRegion,fieldRegion
    TYPE(VARYING_STRING) :: localError

    NULLIFY(meshRegion)
    NULLIFY(fieldRegion)

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
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
                  CALL GeneratedMeshCylinderGeometricParametersCalculate(generatedMesh%cylinderMesh,geometricField,err,error,*999)
                CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
                  CALL GENERATED_MESH_ELLIPSOID_GEOMETRIC_PARAMETERS_CALCULATE(generatedMesh%ellipsoidMesh,geometricField,err,error,*999)
                CASE DEFAULT
                  localError="The generated mesh mesh type of "// &
                    & TRIM(NumberToVString(generatedMesh%generatedType,"*",err,error))// &
                    & " is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              ELSE
                localError="The generated mesh region user number of "// &
                  & TRIM(NumberToVString(meshRegion%USER_NUMBER,"*",err,error))// &
                  & " does not match the geometric field region user number of "//&
                  & TRIM(NumberToVString(fieldRegion%USER_NUMBER,"*",err,error))//"."
                CALL FlagError(localErr,err,error,*999)
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

    CALL Enters("GeneratedMeshRegionGet",err,error,*999)

    IF(ASSOCIATED(generatedMesh)) THEN
      IF(ASSOCIATED(region)) THEN
        CALL FlagError("Region is already associated.",err,error,*999)
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
              localError="The parent region not associated for generated mesh number "// &
                & TRIM(NumberToVString(generatedMesh%userNumber,"*",err,error))//" of interface number "// &
                & TRIM(NumberToVString(interface%USER_NUMBER,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ELSE
            localError="The region or interface is not associated for generated mesh number "// &
              & TRIM(NumberToVString(generatedMesh%userNumber,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Generated mesh is not associated.",err,error,*999)
    ENDIF

    CALL Exits("GeneratedMeshRegionGet")
    RETURN
999 CALL Errors("GeneratedMeshRegionGet",err,error)
    CALL Exits("GeneratedMeshRegionGet")
    RETURN 1
    
  END SUBROUTINE GeneratedMeshRegionGet

  !
  !
  !================================================================================================================================
  !

  !>Updates the geometric field parameters from the initial nodal positions of the regular mesh.
  SUBROUTINE GeneratedMeshRegularGeometricParametersCalculate(regularMesh,geometricField,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshRegularType), POINTER :: regularMesh !<A pointer to the regular mesh object
    TYPE(FIELD_TYPE), POINTER :: geometricField !<A pointer to the field to update the geometric parameters for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    !Local variables
    INTEGER(INTG) :: componentIdx,derivativeIdx,componentNode,meshComponent,nodeIdx,nodePositionIdx(3), &
      & totalNumberOfNodesXi(3),xiIdx,nodeUserNumber
    REAL(DP) :: deltaCoordinate(3,3),MY_ORIGIN(3),VALUE
    REAL(DP) :: derivativeValues(MAXIMUM_GLOBAL_DERIV_NUMBER)
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: coordinateSystem
    TYPE(DOMAIN_TYPE), POINTER :: domain
    TYPE(DOMAIN_NODES_TYPE), POINTER :: domainNodes
    TYPE(FIELD_VARIABLE_COMPONENT_TYPE), POINTER :: fieldVariableComponent
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable
    TYPE(VARYING_STRING) :: localError
    LOGICAL :: nodeExists,ghostNode

    NULLIFY(coordinateSystem)

    CALL Enters("GeneratedMeshRegularGeometricParametersCalculate",err,error,*999)

    IF(ASSOCIATED(regularMesh)) THEN
      IF(ASSOCIATED(geometricField)) THEN
        CALL GeneratedMeshCoordinateSystemGet(generatedMesh,coordinateSystem,err,error,*999)
        CALL COORDINATE_SYSTEM_DIMENSION_GET(coordinateSystem,coordinateDimension,err,error,*999)
        CALL COORDINATE_SYSTEM_TYPE_GET(coordinateSystem,coordinateType,err,error,*999)
        IF(coordinateType==COORDINATE_RECTANGULAR_CARTESIAN_TYPE) THEN
          deltaCoordinate=0.0_DP
          totalNumberOfNodesXi=1
          CALL FIELD_VARAIBLE_GET(geometricField,FIELD_U_VARIABLE_TYPE,fieldVariable,err,error,*999)
          DO componentIdx=1,fieldVariable%NUMBER_OF_COMPONENTS
            fieldVariableComponent=>fieldVariable%components(componentIdx)
            meshComponent=fieldVariableComponent%MESH_COMPONENT_NUMBER
            basis=>regularMesh%bases(meshComponent)%ptr
            !Calculate the total number of nodes in each xi direction
            DO xiIdx=1,regularMesh%meshDimension
              totalNumberOfNodesXi(xiIdx)=(basis%NUMBER_OF_NODES_XIC(xiIdx)-2)*regularMesh%numberOfElementsXi(xiIdx)+ &
                & regularMesh%numberOfElementsXi(xiIdx)+1
            ENDDO !xiIdx
            !Calculate delta
            DO xiIdx=1,regularMesh%meshDimension
              deltaCoordinate(1:coordinateDimension,xiIdx)=regularMesh%baseVectors(1:coordinateDimension,xiIdx)/ &
                & REAL(totalNumberOfNodesXi(xiIdx)-1,DP)
            ENDDO !xiIdx
            derivativeValues=0.0_DP
            CALL FIELD_SCALING_TYPE_GET(geometricField,scalingType,err,error,*999)
            SELECT CASE(scalingType)
            CASE(FIELD_NO_SCALING,FIELD_UNIT_SCALING)
              IF(regularMesh%numberOfElementsXi(1)>0) THEN
                derivativeValues(GLOBAL_DERIV_S1)=regularMesh%baseVectors(componentIdx,1)/regularMesh%numberOfElementsXi(1)
              END IF
              IF(regularMesh%meshDimension>1) THEN
                IF(regularMesh%numberOfElementsXi(2)>0) THEN
                  derivativeValues(GLOBAL_DERIV_S2)=regularMesh%baseVectors(componentIdx,2)/regularMesh%numberOfElementsXi(2)
                END IF
              ENDIF
              IF(regularMesh%meshDimension>2) THEN
                IF(regularMesh%numberOfElementsXi(3)>0) THEN
                  derivativeValues(GLOBAL_DERIV_S3)=regularMesh%baseVectors(componentIdx,3)/regularMesh%numberOfElementsXi(3)
                END IF
              ENDIF
            CASE(FIELD_ARC_LENGTH_SCALING,FIELD_ARITHMETIC_MEAN_SCALING, &
              & FIELD_GEOMETRIC_MEAN_SCALING,FIELD_HARMONIC_MEAN_SCALING)
              !Arc length scalings
              IF(regularMesh%numberOfElementsXi(1)>0) THEN
                derivativeValues(GLOBAL_DERIV_S1)=regularMesh%baseVectors(componentIdx,1)/L2Norm(regularMesh%baseVectors(:,1))
              END IF
              IF(regularMesh%meshDimension>1) THEN
                IF(regularMesh%numberOfElementsXi(2)>0) THEN
                  derivativeValues(GLOBAL_DERIV_S2)=regularMesh%baseVectors(componentIdx,2)/L2Norm(regularMesh%baseVectors(:,2))
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
                CALL GENERATED_MESH_REGULAR_COMPONENT_NODE_TO_USER_NUMBER(regularMesh%generatedMesh,meshComponent, &
                  & componentNode,nodeUserNumber,err,error,*999)
              ELSE
                nodeUserNumber=COMPONENT_NODE_TO_USER_NUMBER(regularMesh%GENERATED_MESH,MESH_COMPONENT, &
                  & componentNode,err,error)
              END IF
              CALL DOMAIN_TOPOLOGY_NODE_CHECK_EXISTS(fieldVariableComponent%domain%topology,nodeUserNumber, &
                nodeExists,nodeIdx,ghostNode,err,error,*999)
              IF(nodeExists.AND..NOT.ghostNode) THEN
                nodePositionIdx(3)=(componentNode-1)/(totalNumberOfNodesXi(2)*totalNumberOfNodesXi(1))+1
                nodePositionIdx(2)=MOD(componentNode-1,totalNumberOfNodesXi(2)*totalNumberOfNodesXi(1))/ &
                  & totalNumberOfNodesXi(1)+1
                nodePositionIdx(1)=MOD(MOD(componentNode-1,totalNumberOfNodesXi(2)*totalNumberOfNodesXi(1)), &
                  & totalNumberOfNodesXi(1))+1
                VALUE=0.0_DP
                DO xiIdx=1,regularMesh%meshDimension
                  VALUE=VALUE+REAL(nodePositionIdx(xiIdx)-1,DP)*deltaCoordinate(componentIdx,xiIdx)
                ENDDO !xiIdx
                VALUE=regularMesh%origin(componentIdx)+VALUE
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

  !>Updates the geometric field parameters from the initial nodal positions of the mesh. Derivatives are averaged via straight line approximation, except for circumferential component
  SUBROUTINE GeneratedMeshCylinderGeometricParametersCalculate(cylinderMesh,geometricField,err,error,*)

    ! Argument variables
    TYPE(GeneratedMeshCylinderType), POINTER :: cylinderMesh !<A pointer to the cylinder mesh object
    TYPE(FIELD_TYPE), POINTER :: geometricField !<A pointer to the field to update the geometric parameters for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    ! Local variables
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(DOMAIN_TYPE),POINTER :: DOMAIN
    TYPE(DOMAIN_NODES_TYPE), POINTER :: domainNodes
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable
    TYPE(FIELD_VARIABLE_COMPONENT_TYPE), POINTER :: fieldVariableComponent
    INTEGER(INTG) :: NUMBER_ELEMENTS_XI(3),NUMBER_OF_NODES_XIC(3)
    INTEGER(INTG) :: totalNumberOfNodesXi(3),INTERPOLATION_TYPES(3)
    INTEGER(INTG) :: componentIdx,xiIdx
    INTEGER(INTG) :: nodeIdx,globalNodeIdx,componentNodeIdx,ny,nk
    INTEGER(INTG) :: numberOfPlanarNodes,SCALING_TYPE,MESH_COMPONENT
    INTEGER(INTG), ALLOCATABLE :: NIDX(:,:,:),EIDX(:,:,:)
    INTEGER(INTG) :: nodePositionIdx(3) ! holds r,theta,z indices
    REAL(DP) :: deltaCoordinate(3),deltaCoordinateXi(3),polarCoordinates(3),rectangularCoordinates(3)
    REAL(DP) :: CYLINDER_EXTENT(3),DERIV
    TYPE(VARYING_STRING) :: localError

    NULLIFY(coordinateSystem)

    CALL Enters("GeneratedMeshCylinderGeometricParametersCalculate",err,error,*999)

    IF(ASSOCIATED(cylinderMesh)) THEN
      IF(ASSOCIATED(geometricField)) THEN
        CALL GeneratedMeshCoordinateSystemGet(generatedMesh,coordinateSystem,err,error,*999)
        CALL COORDINATE_SYSTEM_DIMENSION_GET(coordinateSystem,coordinateDimension,err,error,*999)
        CALL COORDINATE_SYSTEM_TYPE_GET(coordinateSystem,coordinateType,err,error,*999)
        IF(coordinateType==COORDINATE_RECTANGULAR_CARTESIAN_TYPE) THEN
          fieldVariable=>geometricField%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(fieldVariable)) THEN
            IF(fieldVariable%NUMBER_OF_COMPONENTS==3) THEN
              DO componentIdx=1,fieldVariable%NUMBER_OF_COMPONENTS
                fieldVariableComponent=>fieldVariable%components(componentIdx)
                meshComponent=fieldVariableComponent%MESH_COMPONENT_NUMBER
                basis=>cylinderMesh%bases(mechComponent)%ptr
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
                domain=>fieldVariable%COMPONENTS(MESH_COMPONENT)%domain
                domainNodes=>domain%TOPOLOGY%NODES
                DO nodeIdx=1,domainNodes%NUMBER_OF_NODES
                  globalNodeIdx=domainNodes%NODES(nodeIdx)%GLOBAL_NUMBER
                  componentNodeIdx=USER_NUMBER_TO_COMPONENT_NODE(cylinderMesh%generatedMeshes,meshComponent,globalNodeIdx,err,error)
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
                ny=fieldVariableComponent%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(nodeIdx)%DERIVATIVES(1)%VERSIONS(1)
                CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,ny, &
                  & rectangularCoordinates(componentIdx),err,error,*999)
                ! Do derivatives: if there are derivatives, we can assume it's cubic hermite
                !   given that quadratic hermites are only used for collapsed hex elements,
                !   but NB mixed bases have to be handled (e.g. CH-CH-linear combinations)
                IF(domainNodes%NODES(nodeIdx)%NUMBER_OF_DERIVATIVES>1) THEN
                  ! Since I decided how xi 1,2,3 line up with the cylinder polar coordinates,
                  ! we know a priori that only some of the derivatives are nonzero (analytically).
                  ! NOTE: if hermite type used, should assign FIELD_UNIT_SCALING type for this to work
                  DO nk=2,fieldVariableComponent%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(nodeIdx)%NUMBER_OF_DERIVATIVES
                CALL FIELD_SCALING_TYPE_GET(geometricField,scalingType,err,error,*999)
              SELECT CASE(scalingType)
              CASE(FIELD_NO_SCALING,FIELD_UNIT_SCALING)
                    SELECT CASE(domainNodes%NODES(nodeIdx)%DERIVATIVES(nk)%GLOBAL_DERIVATIVE_INDEX)
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
                    ! assign derivative
                    !Default to version 1 of each node derivative
                    ny=fieldVariableComponent%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(nodeIdx)%DERIVATIVES(nk)% &
                      & VERSIONS(1)
                    CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                         & ny,DERIV,err,error,*999)
                  ENDDO !nk
                ENDIF !derivatives
              ENDDO !nodeIdx
            ENDDO !componentIdx
          ELSE
            CALL FlagError("All field variable components must have node-based interpolation.",err,error,*999)
          ENDIF
          CALL FIELD_PARAMETER_SET_UPDATE_START(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
          CALL FIELD_PARAMETER_SET_UPDATE_FINISH(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
        ELSE
          CALL FlagError("Geometric field must be three dimensional.",err,error,*999)
        ENDIF
      ELSE
        localError="The standard field variable is not associated for field number "// &
          & TRIM(NumberToVString(geometricField%USER_NUMBER,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ELSE
      localError="Field number "//TRIM(NumberToVString(geometricField%USER_NUMBER,"*",err,error))//" is not a geometric field."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    ! all done
    IF(ALLOCATED(NIDX)) DEALLOCATE(NIDX)
    IF(ALLOCATED(EIDX)) DEALLOCATE(EIDX)

    CALL Exits("GeneratedMeshCylinderGeometricParametersCalculate")
    RETURN
999 IF(ALLOCATED(NIDX)) DEALLOCATE(NIDX)
    IF(ALLOCATED(EIDX)) DEALLOCATE(EIDX)
    CALL Errors("GeneratedMeshCylinderGeometricParametersCalculate",err,error)
    CALL Exits("GeneratedMeshCylinderGeometricParametersCalculate")
    RETURN 1

  END SUBROUTINE GeneratedMeshCylinderGeometricParametersCalculate

  !
  !================================================================================================================================
  !

  !>Updates the geometric field parameters from the initial nodal positions of the mesh.
  !>Derivatives are averaged via straight line approximation, except for circumferential component
  SUBROUTINE GENERATED_MESH_ELLIPSOID_GEOMETRIC_PARAMETERS_CALCULATE(ELLIPSOID_MESH,FIELD,err,error,*)

    ! Argument variables
    TYPE(GeneratedMeshEllipsoidType), POINTER :: ELLIPSOID_MESH !<A pointer to the ellipsoid mesh object
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to update the geometric parameters for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    ! Local variables
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(DOMAIN_TYPE),POINTER :: DOMAIN
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(FIELD_VARIABLE_COMPONENT_TYPE), POINTER :: FIELD_VARIABLE_COMPONENT
    INTEGER(INTG) :: MY_COMPUTATIONAL_NODE,DOMAIN_NUMBER,MESH_COMPONENT,basis_idx
    INTEGER(INTG) :: NUMBER_ELEMENTS_XIQ(3),NUMBER_OF_NODES_XICQ(3)
    INTEGER(INTG) :: TOTAL_NUMBER_NODES_XI(3),INTERPOLATION_TYPES(3)
    INTEGER(INTG) :: component_idx,xi_idx
    INTEGER(INTG) :: np,npg,i,j,k, local_node
    INTEGER(INTG) :: SCALING_TYPE!,NUMBER_OF_PLANAR_NODES
    INTEGER(INTG), ALLOCATABLE :: NIDX(:,:,:),EIDX(:,:,:)
    !INTEGER(INTG) :: node_idx(3) ! holds r,theta,z indices
    REAL(DP) :: DELTA(3),DELTAi(3),RECT_COORDS(3),t,phi,alpha,xi,nu,x,y,z
    REAL(DP) :: ELLIPSOID_EXTENT(4)
    TYPE(VARYING_STRING) :: localError

    NULLIFY(BASIS,DOMAIN,DECOMPOSITION,DOMAIN_NODES,FIELD_VARIABLE,FIELD_VARIABLE_COMPONENT)

    CALL Enters("GENERATED_MESH_ELLIPSOID_GEOMETRIC_PARAMETERS_CALCULATE",err,error,*999)

    MY_COMPUTATIONAL_NODE=COMPUTATIONAL_NODE_NUMBER_GET(err,error)

    ! assign to the field
    np=0
    IF(FIELD%TYPE==FIELD_GEOMETRIC_TYPE) THEN
       FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR
       IF(ASSOCIATED(FIELD_VARIABLE)) THEN
          IF(FIELD_VARIABLE%NUMBER_OF_COMPONENTS==3) THEN
             MESH_COMPONENT=FIELD_VARIABLE%COMPONENTS(1)%MESH_COMPONENT_NUMBER
             DO component_idx=2,3
                IF(FIELD_VARIABLE%COMPONENTS(component_idx)%MESH_COMPONENT_NUMBER/=MESH_COMPONENT) THEN
                   CALL FlagError("Multiple mesh components for geometric components is not implemented.",err,error,*999)
                ENDIF
             ENDDO
             basis_idx=MESH_COMPONENT*2-1
             !< Ellipsoid_extent= inner long axis, inner short axis, wall thickness, top angle (from 0)
             ! calculate the total number of nodes in each xi direction
             IF(ALLOCATED(ELLIPSOID_MESH%bases)) THEN
                !Check that the all geometric bases use the same mesh component
                BASIS=>ELLIPSOID_MESH%bases(basis_idx)%PTR
                NUMBER_ELEMENTS_XI=ELLIPSOID_MESH%numberOfElementsXi
                NUMBER_OF_NODES_XIC=BASIS%NUMBER_OF_NODES_XIC
                DO xi_idx=1,3
                   TOTAL_NUMBER_NODES_XI(xi_idx)=(NUMBER_OF_NODES_XIC(xi_idx)-1)*NUMBER_ELEMENTS_XI(xi_idx)+1
                ENDDO
                TOTAL_NUMBER_NODES_XI(1)=TOTAL_NUMBER_NODES_XI(1)-1 ! theta loops around so slightly different
                ! calculate DELTAi now
                ELLIPSOID_EXTENT=ELLIPSOID_MESH%ellipsoidExtent
                DELTA(1)=TWOPI/NUMBER_ELEMENTS_XI(1)
                DELTA(2)=(PI-ELLIPSOID_EXTENT(4))/NUMBER_ELEMENTS_XI(2)
                DELTA(3)=ELLIPSOID_EXTENT(3)/NUMBER_ELEMENTS_XI(3)
                DO xi_idx=1,3
                   DELTAi(xi_idx)=DELTA(xi_idx)/(NUMBER_OF_NODES_XIC(xi_idx)-1)
                ENDDO
             ELSE
                CALL FlagError("Ellipsoid mesh does not have bases allocated.",err,error,*999)
             ENDIF
             CALL FIELD_SCALING_TYPE_GET(FIELD,SCALING_TYPE,err,error,*999)
             IF(SCALING_TYPE/=FIELD_UNIT_SCALING) &
                  & CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"  Note: If the ellipsoid looks wonky, set field scaling to &
                  & unit scaling type.",err,error,*999)
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
                   npg=COMPONENT_NODE_TO_USER_NUMBER(ELLIPSOID_MESH%GENERATED_MESH,basis_idx,np,err,error)
                   CALL DECOMPOSITION_NODE_DOMAIN_GET(DECOMPOSITION,npg,MESH_COMPONENT,DOMAIN_NUMBER,err,error,*999)
                   IF(DOMAIN_NUMBER==MY_COMPUTATIONAL_NODE) THEN
                      RECT_COORDS(1)=0
                      RECT_COORDS(2)=0
                      RECT_COORDS(3)=-ELLIPSOID_EXTENT(1)
                      DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                         !Default to version 1 of each node derivative
                         CALL FIELD_PARAMETER_SET_UPDATE_NODE(FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,npg, &
                              & component_idx,RECT_COORDS(component_idx),err,error,*999)
                         local_node=DOMAIN%MAPPINGS%NODES%GLOBAL_TO_LOCAL_MAP(npg)%local_number(1)
                         IF(DOMAIN_NODES%NODES(local_node)%NUMBER_OF_DERIVATIVES>1) THEN
                            CALL FlagError("Not generalized to hermittean elements.",err,error,*999)
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
                         npg=COMPONENT_NODE_TO_USER_NUMBER(ELLIPSOID_MESH%GENERATED_MESH,basis_idx,np,err,error)
                         CALL DECOMPOSITION_NODE_DOMAIN_GET(DECOMPOSITION,npg,MESH_COMPONENT,DOMAIN_NUMBER,err,error,*999)
                         IF(DOMAIN_NUMBER==MY_COMPUTATIONAL_NODE) THEN
                            DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                               CALL FIELD_PARAMETER_SET_UPDATE_NODE(FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,npg, &
                                    & component_idx,RECT_COORDS(component_idx),err,error,*999)
                               local_node=DOMAIN%MAPPINGS%NODES%GLOBAL_TO_LOCAL_MAP(npg)%local_number(1)
                               IF(DOMAIN_NODES%NODES(local_node)%NUMBER_OF_DERIVATIVES>1) THEN
                                  CALL FlagError("Not generalized to hermittean elements.",err,error,*999)
                               ENDIF !derivatives
                            ENDDO
                         ENDIF
                      ENDDO
                   ENDDO

                   DO k=2,TOTAL_NUMBER_NODES_XI(3)
                      !transmural loop
                      j=1
                      !apex nodes
                      RECT_COORDS(1)=0
                      RECT_COORDS(2)=0
                      RECT_COORDS(3)=-ELLIPSOID_EXTENT(1)-(k-1)*(DELTAi(3))
                      np=np+1
                      npg=COMPONENT_NODE_TO_USER_NUMBER(ELLIPSOID_MESH%GENERATED_MESH,basis_idx,np,err,error)
                      CALL DECOMPOSITION_NODE_DOMAIN_GET(DECOMPOSITION,npg,MESH_COMPONENT,DOMAIN_NUMBER,err,error,*999)
                      IF(DOMAIN_NUMBER==MY_COMPUTATIONAL_NODE) THEN
                         DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                            CALL FIELD_PARAMETER_SET_UPDATE_NODE(FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,npg, &
                                 & component_idx,RECT_COORDS(component_idx),err,error,*999)
                            local_node=DOMAIN%MAPPINGS%NODES%GLOBAL_TO_LOCAL_MAP(npg)%local_number(1)
                            IF(DOMAIN_NODES%NODES(local_node)%NUMBER_OF_DERIVATIVES>1) THEN
                               CALL FlagError("Not generalized to hermittean elements.",err,error,*999)
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
                            npg=COMPONENT_NODE_TO_USER_NUMBER(ELLIPSOID_MESH%GENERATED_MESH,basis_idx,np,err,error)
                            CALL DECOMPOSITION_NODE_DOMAIN_GET(DECOMPOSITION,npg,MESH_COMPONENT,DOMAIN_NUMBER,err,error,*999)
                            IF(DOMAIN_NUMBER==MY_COMPUTATIONAL_NODE) THEN
                               DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                                  CALL FIELD_PARAMETER_SET_UPDATE_NODE(FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,npg, &
                                       & component_idx,RECT_COORDS(component_idx),err,error,*999)
                                  local_node=DOMAIN%MAPPINGS%NODES%GLOBAL_TO_LOCAL_MAP(npg)%local_number(1)
                                  IF(DOMAIN_NODES%NODES(local_node)%NUMBER_OF_DERIVATIVES>1) THEN
                                     CALL FlagError("Not generalized to hermittean elements.",err,error,*999)
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
                      !transmural loop
                      alpha=ELLIPSOID_EXTENT(1)+(k-1)*(DELTAi(3))
                      j=1
                      !apex nodes
                      RECT_COORDS(1)=0
                      RECT_COORDS(2)=0
                      RECT_COORDS(3)=-alpha
                      np=np+1
                      npg=COMPONENT_NODE_TO_USER_NUMBER(ELLIPSOID_MESH%GENERATED_MESH,basis_idx,np,err,error)
                      CALL DECOMPOSITION_NODE_DOMAIN_GET(DECOMPOSITION,npg,MESH_COMPONENT,DOMAIN_NUMBER,err,error,*999)
                      IF(DOMAIN_NUMBER==MY_COMPUTATIONAL_NODE) THEN
                         DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                            CALL FIELD_PARAMETER_SET_UPDATE_NODE(FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,npg, &
                                 & component_idx,RECT_COORDS(component_idx),err,error,*999)
                            local_node=DOMAIN%MAPPINGS%NODES%GLOBAL_TO_LOCAL_MAP(npg)%local_number(1)
                            IF(DOMAIN_NODES%NODES(local_node)%NUMBER_OF_DERIVATIVES>1) THEN
                               CALL FlagError("Not generalized to hermittean elements.",err,error,*999)
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
                            npg=COMPONENT_NODE_TO_USER_NUMBER(ELLIPSOID_MESH%GENERATED_MESH,basis_idx,np,err,error)
                            CALL DECOMPOSITION_NODE_DOMAIN_GET(DECOMPOSITION,npg,MESH_COMPONENT,DOMAIN_NUMBER,err,error,*999)
                            IF(DOMAIN_NUMBER==MY_COMPUTATIONAL_NODE) THEN
                               DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                                  CALL FIELD_PARAMETER_SET_UPDATE_NODE(FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,npg, &
                                       & component_idx,RECT_COORDS(component_idx),err,error,*999)
                                  local_node=DOMAIN%MAPPINGS%NODES%GLOBAL_TO_LOCAL_MAP(npg)%local_number(1)
                                  IF(DOMAIN_NODES%NODES(local_node)%NUMBER_OF_DERIVATIVES>1) THEN
                                     CALL FlagError("Not generalized to hermittean elements.",err,error,*999)
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
                   npg=COMPONENT_NODE_TO_USER_NUMBER(ELLIPSOID_MESH%GENERATED_MESH,basis_idx,np,err,error)
                   CALL DECOMPOSITION_NODE_DOMAIN_GET(DECOMPOSITION,npg,MESH_COMPONENT,DOMAIN_NUMBER,err,error,*999)
                   IF(DOMAIN_NUMBER==MY_COMPUTATIONAL_NODE) THEN
                      RECT_COORDS(1)=0
                      RECT_COORDS(2)=0
                      RECT_COORDS(3)=-ELLIPSOID_EXTENT(1)
                      DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                         CALL FIELD_PARAMETER_SET_UPDATE_NODE(FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,npg, &
                              & component_idx,RECT_COORDS(component_idx),err,error,*999)
                         local_node=DOMAIN%MAPPINGS%NODES%GLOBAL_TO_LOCAL_MAP(npg)%local_number(1)
                         IF(DOMAIN_NODES%NODES(local_node)%NUMBER_OF_DERIVATIVES>1) THEN
                            CALL FlagError("Not generalized to hermittean elements.",err,error,*999)
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
                         npg=COMPONENT_NODE_TO_USER_NUMBER(ELLIPSOID_MESH%GENERATED_MESH,basis_idx,np,err,error)
                         CALL DECOMPOSITION_NODE_DOMAIN_GET(DECOMPOSITION,npg,MESH_COMPONENT,DOMAIN_NUMBER,err,error,*999)
                         IF(DOMAIN_NUMBER==MY_COMPUTATIONAL_NODE) THEN
                            DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                               CALL FIELD_PARAMETER_SET_UPDATE_NODE(FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,npg, &
                                    & component_idx,RECT_COORDS(component_idx),err,error,*999)
                               local_node=DOMAIN%MAPPINGS%NODES%GLOBAL_TO_LOCAL_MAP(npg)%local_number(1)
                               IF(DOMAIN_NODES%NODES(local_node)%NUMBER_OF_DERIVATIVES>1) THEN
                                  CALL FlagError("Not generalized to hermittean elements.",err,error,*999)
                               ENDIF !derivatives
                            ENDDO
                         ENDIF
                      ENDDO
                   ENDDO

                   DO k=2,TOTAL_NUMBER_NODES_XI(3)
                      !transmural loop
                      j=1
                      !apex nodes
                      RECT_COORDS(1)=0
                      RECT_COORDS(2)=0
                      RECT_COORDS(3)=-ELLIPSOID_EXTENT(1)-(k-1)*(DELTAi(3))
                      np=np+1
                      npg=COMPONENT_NODE_TO_USER_NUMBER(ELLIPSOID_MESH%GENERATED_MESH,basis_idx,np,err,error)
                      CALL DECOMPOSITION_NODE_DOMAIN_GET(DECOMPOSITION,npg,MESH_COMPONENT,DOMAIN_NUMBER,err,error,*999)
                      IF(DOMAIN_NUMBER==MY_COMPUTATIONAL_NODE) THEN
                         DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                            CALL FIELD_PARAMETER_SET_UPDATE_NODE(FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,npg, &
                                 & component_idx,RECT_COORDS(component_idx),err,error,*999)
                            local_node=DOMAIN%MAPPINGS%NODES%GLOBAL_TO_LOCAL_MAP(npg)%local_number(1)
                            IF(DOMAIN_NODES%NODES(local_node)%NUMBER_OF_DERIVATIVES>1) THEN
                               CALL FlagError("Not generalized to hermittean elements.",err,error,*999)
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
                            npg=COMPONENT_NODE_TO_USER_NUMBER(ELLIPSOID_MESH%GENERATED_MESH,basis_idx,np,err,error)
                            CALL DECOMPOSITION_NODE_DOMAIN_GET(DECOMPOSITION,npg,MESH_COMPONENT,DOMAIN_NUMBER,err,error,*999)
                            IF(DOMAIN_NUMBER==MY_COMPUTATIONAL_NODE) THEN
                               DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                                  CALL FIELD_PARAMETER_SET_UPDATE_NODE(FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,npg, &
                                       & component_idx,RECT_COORDS(component_idx),err,error,*999)
                                  local_node=DOMAIN%MAPPINGS%NODES%GLOBAL_TO_LOCAL_MAP(npg)%local_number(1)
                                  IF(DOMAIN_NODES%NODES(local_node)%NUMBER_OF_DERIVATIVES>1) THEN
                                     CALL FlagError("Not generalized to hermittean elements.",err,error,*999)
                                  ENDIF !derivatives
                               ENDDO
                            ENDIF
                         ENDDO
                      ENDDO
                   ENDDO
                ELSE
                   CALL FlagError("Not valid long axis - short axis relation",err,error,*999)
                ENDIF
             ELSE
                CALL FlagError("All field variable components must have node-based interpolation.",err,error,*999)
             ENDIF
             CALL FIELD_PARAMETER_SET_UPDATE_START(FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
             CALL FIELD_PARAMETER_SET_UPDATE_FINISH(FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
          ELSE
             CALL FlagError("Geometric field must be three dimensional.",err,error,*999)
          ENDIF
       ELSE
          localError="The standard field variable is not associated for field number "// &
               & TRIM(NumberToVString(FIELD%USER_NUMBER,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
       ENDIF
    ELSE
       localError="Field number "//TRIM(NumberToVString(FIELD%USER_NUMBER,"*",err,error))//" is not a geometric field."
       CALL FlagError(localError,err,error,*999)
    ENDIF

    ! all done
    IF(ALLOCATED(NIDX)) DEALLOCATE(NIDX)
    IF(ALLOCATED(EIDX)) DEALLOCATE(EIDX)

    CALL Exits("GENERATED_MESH_ELLIPSOID_GEOMETRIC_PARAMETERS_CALCULATE")
    RETURN
999 IF(ALLOCATED(NIDX)) DEALLOCATE(NIDX)
    IF(ALLOCATED(EIDX)) DEALLOCATE(EIDX)
    CALL Errors("GENERATED_MESH_ELLIPSOID_GEOMETRIC_PARAMETERS_CALCULATE",err,error)
    CALL Exits("GENERATED_MESH_ELLIPSOID_GEOMETRIC_PARAMETERS_CALCULATE")
    RETURN 1

  END SUBROUTINE GENERATED_MESH_ELLIPSOID_GEOMETRIC_PARAMETERS_CALCULATE

  !
  !================================================================================================================================
  !

  !>Provides an easy way to grab surfaces for boundary condition assignment
  SUBROUTINE GENERATED_MESH_REGULAR_SURFACE_GET(REGULAR_MESH,MESH_COMPONENT,SURFACE_TYPE,SURFACE_NODES,NORMAL_XI,err,error,*)

    ! Argument variables
    TYPE(GeneratedMeshRegularType), POINTER :: REGULAR_MESH !<A pointer to the regular mesh object
    INTEGER(INTG), INTENT(IN) :: MESH_COMPONENT !<The mesh component to get the surface for.
    INTEGER(INTG), INTENT(IN) :: SURFACE_TYPE !<A constant identifying the type of surface to get \see GeneratedMeshRoutines_GeneratedMeshRegularSurfaces,GeneratedMeshRoutines
    INTEGER(INTG), ALLOCATABLE, INTENT(OUT) :: SURFACE_NODES(:) !<On exit, contains the list of nodes belonging to the surface
    INTEGER(INTG), INTENT(OUT) :: NORMAL_XI !<On exit, contains the xi direction of the outward pointing normal of the surface
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    ! Local variables
    TYPE(BASIS_TYPE), POINTER :: BASIS
    INTEGER(INTG),ALLOCATABLE :: NIDX(:,:,:),EIDX(:,:,:)
    INTEGER(INTG) :: NUMBER_OF_ELEMENTS_XI(3) !Specified number of elements in each xi direction
    INTEGER(INTG) :: NUMBER_OF_NODES_XI(3) ! Number of nodes per element in each xi direction (basis property)
    INTEGER(INTG) :: num_dims,TOTAL_NUMBER_OF_NODES,TOTAL_NUMBER_OF_ELEMENTS,NODE_USER_NUMBER
    REAL(DP) :: DELTA(3),DELTAI(3)
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: node_counter,i,j,k

    CALL Enters("GENERATED_MESH_REGULAR_SURFACE_GET",err,error,*999)

    IF(ALLOCATED(REGULAR_MESH%numberOfElementsXi)) THEN
      num_dims=SIZE(REGULAR_MESH%numberOfElementsXi)
      IF(num_dims==2) THEN
        NUMBER_OF_ELEMENTS_XI(1:2)=REGULAR_MESH%numberOfElementsXi(1:2)
        NUMBER_OF_ELEMENTS_XI(3)=1
      ELSE IF (num_dims==1) THEN
        NUMBER_OF_ELEMENTS_XI(1)=REGULAR_MESH%numberOfElementsXi(1)
        NUMBER_OF_ELEMENTS_XI(2)=1
        NUMBER_OF_ELEMENTS_XI(3)=1
      ELSE
        NUMBER_OF_ELEMENTS_XI=REGULAR_MESH%numberOfElementsXi
      ENDIF
      IF(ASSOCIATED(REGULAR_MESH%bases(MESH_COMPONENT)%PTR)) THEN
        BASIS=>REGULAR_MESH%bases(MESH_COMPONENT)%PTR
        IF(.NOT.ALLOCATED(SURFACE_NODES)) THEN
          !Node that simplex bases have an extra area coordinate so size of number_of_nodes_xic=num_dims+1
          NUMBER_OF_NODES_XI(1:num_dims)=BASIS%NUMBER_OF_NODES_XIC(1:num_dims)
          NUMBER_OF_NODES_XI(num_dims+1:3)=1

          ! build indices first (some of these are dummy arguments)
          CALL GENERATED_MESH_REGULAR_BUILD_NODE_INDICES(NUMBER_OF_ELEMENTS_XI,NUMBER_OF_NODES_XI, &
              & REGULAR_MESH%maximumExtent,TOTAL_NUMBER_OF_NODES,TOTAL_NUMBER_OF_ELEMENTS,NIDX,EIDX,DELTA,DELTAI,err,error,*999)
          node_counter=0
          SELECT CASE(SURFACE_TYPE)
          CASE(GENERATED_MESH_REGULAR_LEFT_SURFACE)
            ALLOCATE(SURFACE_NODES((SIZE(NIDX,2))*(SIZE(NIDX,3))),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate NODES array.",err,error,*999)
            DO k=1,SIZE(NIDX,3)
              DO j=1,SIZE(NIDX,2)
                node_counter=node_counter+1
                SURFACE_NODES(node_counter)=NIDX(1,j,k)
              ENDDO
            ENDDO
            NORMAL_XI=-1
          CASE(GENERATED_MESH_REGULAR_RIGHT_SURFACE)
            ALLOCATE(SURFACE_NODES((SIZE(NIDX,2))*(SIZE(NIDX,3))),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate NODES array.",err,error,*999)
            DO k=1,SIZE(NIDX,3)
              DO j=1,SIZE(NIDX,2)
                node_counter=node_counter+1
                SURFACE_NODES(node_counter)=NIDX(SIZE(NIDX,1),j,k)
              ENDDO
            ENDDO
            NORMAL_XI=1
          CASE(GENERATED_MESH_REGULAR_TOP_SURFACE)
            ALLOCATE(SURFACE_NODES((SIZE(NIDX,1))*(SIZE(NIDX,2))),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate NODES array.",err,error,*999)
            DO j=1,SIZE(NIDX,2)
              DO i=1,SIZE(NIDX,1)
                node_counter=node_counter+1
                SURFACE_NODES(node_counter)=NIDX(i,j,SIZE(NIDX,3))
              ENDDO
            ENDDO
            NORMAL_XI=3
          CASE(GENERATED_MESH_REGULAR_BOTTOM_SURFACE)
            ALLOCATE(SURFACE_NODES((SIZE(NIDX,1))*(SIZE(NIDX,2))),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate NODES array.",err,error,*999)
            DO j=1,SIZE(NIDX,2)
              DO i=1,SIZE(NIDX,1)
                node_counter=node_counter+1
                SURFACE_NODES(node_counter)=NIDX(i,j,1)
              ENDDO
            ENDDO
            NORMAL_XI=-3
          CASE(GENERATED_MESH_REGULAR_FRONT_SURFACE)
            ALLOCATE(SURFACE_NODES((SIZE(NIDX,1))*(SIZE(NIDX,3))),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate NODES array.",err,error,*999)
            DO j=1,SIZE(NIDX,3)
              DO i=1,SIZE(NIDX,1)
                node_counter=node_counter+1
                SURFACE_NODES(node_counter)=NIDX(i,1,j)
              ENDDO
            ENDDO
            NORMAL_XI=-2
          CASE(GENERATED_MESH_REGULAR_BACK_SURFACE)
            ALLOCATE(SURFACE_NODES((SIZE(NIDX,1))*(SIZE(NIDX,3))),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate NODES array.",err,error,*999)
            DO j=1,SIZE(NIDX,3)
              DO i=1,SIZE(NIDX,1)
                node_counter=node_counter+1
                SURFACE_NODES(node_counter)=NIDX(i,SIZE(NIDX,2),j)
              ENDDO
            ENDDO
            NORMAL_XI=2
          CASE DEFAULT
            localError="The specified surface type of "//TRIM(NumberToVString(SURFACE_TYPE,"*",err,error))// &
              & " is invalid for a regular mesh."
            CALL FlagError(localError,err,error,*999)
          END SELECT
          !Now convert the component node numbering to user numbers if a mesh has multiple components
          DO node_counter=1,SIZE(SURFACE_NODES,1)
            SELECT CASE(REGULAR_MESH%bases(MESH_COMPONENT)%PTR%TYPE)
            CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
              CALL GENERATED_MESH_REGULAR_COMPONENT_NODE_TO_USER_NUMBER(REGULAR_MESH%GENERATED_MESH,MESH_COMPONENT, &
                  & SURFACE_NODES(node_counter),NODE_USER_NUMBER,err,error,*999)
              SURFACE_NODES(node_counter)=NODE_USER_NUMBER
            CASE(BASIS_SIMPLEX_TYPE)
              SURFACE_NODES(node_counter)=COMPONENT_NODE_TO_USER_NUMBER(REGULAR_MESH%GENERATED_MESH,MESH_COMPONENT, &
                  & SURFACE_NODES(node_counter),err,error)
              IF(err/=0) GOTO 999
            CASE DEFAULT
              CALL FlagError("The basis type of "//TRIM(NumberToVString(REGULAR_MESH%bases(MESH_COMPONENT)%PTR%TYPE, &
                & "*",err,error))//" is not implemented when getting a regular mesh surface.",err,error,*999)
            END SELECT
          END DO
        ELSE
          CALL FlagError("Output SURFACE_NODES array is already allocated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Regular mesh object does not have a basis associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Regular mesh object does not have number of elements property specified.",err,error,*999)
    ENDIF

    CALL Exits("GENERATED_MESH_REGULAR_SURFACE_GET")
    RETURN
999 CALL Errors("GENERATED_MESH_REGULAR_SURFACE_GET",err,error)
    CALL Exits("GENERATED_MESH_REGULAR_SURFACE_GET")
    RETURN 1
  END SUBROUTINE GENERATED_MESH_REGULAR_SURFACE_GET

  !
  !================================================================================================================================
  !

  !>Provides an easy way to grab surfaces for boundary condition assignment
  SUBROUTINE GENERATED_MESH_CYLINDER_SURFACE_GET(CYLINDER_MESH,MESH_COMPONENT,SURFACE_TYPE,SURFACE_NODES,NORMAL_XI,err,error,*)

    ! Argument variables
    TYPE(GeneratedMeshCylinderType), POINTER :: CYLINDER_MESH !<A pointer to the cylinder mesh object
    INTEGER(INTG), INTENT(IN) :: MESH_COMPONENT !<The mesh component to get the surface for.
    INTEGER(INTG), INTENT(IN) :: SURFACE_TYPE !<A constant identifying the type of surface to get \see GeneratedMeshRoutines_GeneratedMeshCylinderSurfaces,GeneratedMeshRoutines
    INTEGER(INTG), ALLOCATABLE, INTENT(OUT) :: SURFACE_NODES(:) !<On exit, contains the list of nodes belonging to the surface
    INTEGER(INTG), INTENT(OUT) :: NORMAL_XI !<On exit, contains the xi direction of the outward pointing normal of the surface
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    ! Local variables
    TYPE(BASIS_TYPE), POINTER :: BASIS
    INTEGER(INTG),ALLOCATABLE :: NIDX(:,:,:),EIDX(:,:,:)
    INTEGER(INTG) :: NUMBER_OF_ELEMENTS_XI(3) !Specified number of elements in each xi direction
    INTEGER(INTG) :: NUMBER_OF_NODES_XI(3) ! Number of nodes per element in each xi direction (basis property)
    INTEGER(INTG) :: total_number_of_nodes,total_number_of_elements
    REAL(DP) :: delta(3),deltai(3)
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: node_counter,i, j, k

    CALL Enters("GENERATED_MESH_CYLINDER_SURFACE_GET",err,error,*999)

    ! let's go
    IF(ALLOCATED(CYLINDER_MESH%numberOfElementsXi)) THEN
      NUMBER_OF_ELEMENTS_XI=CYLINDER_MESH%numberOfElementsXi
      IF(ASSOCIATED(CYLINDER_MESH%bases(MESH_COMPONENT)%PTR)) THEN
        BASIS=>CYLINDER_MESH%bases(MESH_COMPONENT)%PTR
        IF(.NOT.ALLOCATED(SURFACE_NODES)) THEN
          NUMBER_OF_NODES_XI=BASIS%NUMBER_OF_NODES_XIC
          ! build indices first (some of these are dummy arguments)
          CALL GENERATED_MESH_CYLINDER_BUILD_NODE_INDICES(NUMBER_OF_ELEMENTS_XI,NUMBER_OF_NODES_XI, &
              & cylinder_mesh%cylinder_extent,total_number_of_nodes,total_number_of_elements,NIDX,EIDX, &
              & delta,deltai,err,error,*999)
          node_counter=0
          SELECT CASE(SURFACE_TYPE)
          CASE(GENERATED_MESH_CYLINDER_INNER_SURFACE)
            ALLOCATE(SURFACE_NODES((SIZE(NIDX,2))*(SIZE(NIDX,3))),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate NODES array.",err,error,*999)
            DO k=1,SIZE(NIDX,3)
              DO j=1,SIZE(NIDX,2)
                node_counter=node_counter+1
                SURFACE_NODES(node_counter)=COMPONENT_NODE_TO_USER_NUMBER(CYLINDER_MESH%GENERATED_MESH,MESH_COMPONENT, &
                    & NIDX(1,j,k),err,error)
              ENDDO
            ENDDO
            NORMAL_XI=-1
          CASE(GENERATED_MESH_CYLINDER_OUTER_SURFACE)
            ALLOCATE(SURFACE_NODES((SIZE(NIDX,2))*(SIZE(NIDX,3))),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate NODES array.",err,error,*999)
            DO k=1,SIZE(NIDX,3)
              DO j=1,SIZE(NIDX,2)
                node_counter=node_counter+1
                SURFACE_NODES(node_counter)=COMPONENT_NODE_TO_USER_NUMBER(CYLINDER_MESH%GENERATED_MESH,MESH_COMPONENT, &
                    & NIDX(SIZE(NIDX,1),j,k),err,error)
              ENDDO
            ENDDO
            NORMAL_XI=1
          CASE(GENERATED_MESH_CYLINDER_TOP_SURFACE)
            ALLOCATE(SURFACE_NODES((SIZE(NIDX,1))*(SIZE(NIDX,2))),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate NODES array.",err,error,*999)
            DO j=1,SIZE(NIDX,2)
              DO i=1,SIZE(NIDX,1)
                node_counter=node_counter+1
                SURFACE_NODES(node_counter)=COMPONENT_NODE_TO_USER_NUMBER(CYLINDER_MESH%GENERATED_MESH,MESH_COMPONENT, &
                    & NIDX(i,j,SIZE(NIDX,3)),err,error)
              ENDDO
            ENDDO
            NORMAL_XI=3
          CASE(GENERATED_MESH_CYLINDER_BOTTOM_SURFACE)
            ALLOCATE(SURFACE_NODES((SIZE(NIDX,1))*(SIZE(NIDX,2))),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate NODES array.",err,error,*999)
            DO j=1,SIZE(NIDX,2)
              DO i=1,SIZE(NIDX,1)
                node_counter=node_counter+1
                SURFACE_NODES(node_counter)=COMPONENT_NODE_TO_USER_NUMBER(CYLINDER_MESH%GENERATED_MESH,MESH_COMPONENT, &
                    & NIDX(i,j,1),err,error)
              ENDDO
            ENDDO
            NORMAL_XI=-3
          CASE DEFAULT
            localError="The specified surface type of "//TRIM(NumberToVString(SURFACE_TYPE,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ELSE
          CALL FlagError("Output SURFACE_NODES array is already allocated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Cylinder mesh object does not have a basis associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Cylinder mesh object does not have number of elements property specified.",err,error,*999)
    ENDIF

    CALL Exits("GENERATED_MESH_CYLINDER_SURFACE_GET")
    RETURN
999 CALL Errors("GENERATED_MESH_CYLINDER_SURFACE_GET",err,error)
    CALL Exits("GENERATED_MESH_CYLINDER_SURFACE_GET")
    RETURN 1
  END SUBROUTINE GENERATED_MESH_CYLINDER_SURFACE_GET
  !
  !================================================================================================================================
  !

  !>Provides an easy way to grab surfaces for boundary condition assignment
  SUBROUTINE GENERATED_MESH_ELLIPSOID_SURFACE_GET(ELLIPSOID_MESH,MESH_COMPONENT,SURFACE_TYPE,SURFACE_NODES,NORMAL_XI,err,error,*)

    ! Argument variables
    TYPE(GeneratedMeshEllipsoidType), POINTER :: ELLIPSOID_MESH !<A pointer to the ellipsoid mesh object
    INTEGER(INTG), INTENT(IN) :: MESH_COMPONENT !<The mesh component to get the surface for.
    INTEGER(INTG), INTENT(IN) :: SURFACE_TYPE !<A constant identifying the type of surface to get \see GeneratedMeshRoutines_GeneratedMeshEllipsoidSurfaces,GeneratedMeshRoutines
    INTEGER(INTG), ALLOCATABLE, INTENT(OUT) :: SURFACE_NODES(:) !<On exit, contains the list of nodes belonging to the surface
    INTEGER(INTG), INTENT(OUT) :: NORMAL_XI !<On exit, contains the xi direction of the outward pointing normal of the surface
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    ! Local variables
    TYPE(BASIS_TYPE), POINTER :: BASIS
    INTEGER(INTG),ALLOCATABLE :: NIDX(:,:,:),EIDX(:,:,:),CORNER_NODES(:,:,:)
    INTEGER(INTG) :: NUMBER_OF_ELEMENTS_XI(3) !Specified number of elements in each xi direction
    INTEGER(INTG) :: NUMBER_OF_NODES_XI(3) ! Number of nodes per element in each xi direction (basis property)
    INTEGER(INTG) :: total_number_of_nodes,total_number_of_elements
    REAL(DP) :: delta(3),deltai(3)
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: node_counter,i, j, k

    CALL Enters("GENERATED_MESH_ELLIPSOID_SURFACE_GET",err,error,*999)

    ! let's go
    IF(ALLOCATED(ELLIPSOID_MESH%numberOfElementsXi)) THEN
      NUMBER_OF_ELEMENTS_XI=ELLIPSOID_MESH%numberOfElementsXi
      IF(ALLOCATED(ELLIPSOID_MESH%bases)) THEN

!         !Below, there is an issue:
!         !  BASIS=>ELLIPSOID_MESH%bases(MESH_COMPONENT)%PTR does not account for the fact that:
!         !  in 'GeneratedMeshEllipsoidCreateFinish' the following is done:
!         !  CALL MESH_NUMBER_OF_COMPONENTS_SET(GENERATED_MESH%MESH,SIZE(ELLIPSOID_MESH%bases)/2,err,error,*999)
!         !A temporary work around is the following (although this bug may need to be fixed in several places):
!
!         IF(MESH_COMPONENT==2) THEN
!           BASIS_COMPONENT = MESH_COMPONENT + 1
!         ELSE
!           BASIS_COMPONENT = MESH_COMPONENT
!         ENDIF
!
!         IF(ASSOCIATED(ELLIPSOID_MESH%bases(BASIS_COMPONENT)%PTR)) THEN
!           BASIS=>ELLIPSOID_MESH%bases(BASIS_COMPONENT)%PTR

        IF(ASSOCIATED(ELLIPSOID_MESH%bases(MESH_COMPONENT)%PTR)) THEN
          BASIS=>ELLIPSOID_MESH%bases(MESH_COMPONENT)%PTR
          IF(.NOT.ALLOCATED(SURFACE_NODES)) THEN
            NUMBER_OF_NODES_XI=BASIS%NUMBER_OF_NODES_XIC
            ! build indices first (some of these are dummy arguments)
            CALL GENERATED_MESH_ELLIPSOID_BUILD_NODE_INDICES(NUMBER_OF_ELEMENTS_XI,NUMBER_OF_NODES_XI, &
                & ellipsoid_mesh%ellipsoid_extent,total_number_of_nodes,total_number_of_elements,NIDX, &
                & CORNER_NODES,EIDX,delta,deltai,err,error,*999)
            node_counter=0

            SELECT CASE(SURFACE_TYPE)

            CASE(GENERATED_MESH_ELLIPSOID_INNER_SURFACE)
              ALLOCATE(SURFACE_NODES((SIZE(NIDX,1))*(SIZE(NIDX,2)-1)+1),STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate NODES array.",err,error,*999)
              j=1
              i=1
              node_counter=node_counter+1
              SURFACE_NODES(node_counter)=COMPONENT_NODE_TO_USER_NUMBER(ELLIPSOID_MESH%GENERATED_MESH,MESH_COMPONENT, &
                  NIDX(i,j,1),err,error)
              DO j=2,SIZE(NIDX,2)
                DO i=1, SIZE(NIDX,1)
                  node_counter=node_counter+1
                  IF (NIDX(i,j,1)/=0) THEN
                    SURFACE_NODES(node_counter)=COMPONENT_NODE_TO_USER_NUMBER(ELLIPSOID_MESH%GENERATED_MESH,MESH_COMPONENT, &
                        & NIDX(i,j,1),err,error)
                  ELSE
                    node_counter=node_counter-1
                  ENDIF
                ENDDO
              ENDDO
              NORMAL_XI=-3

            CASE(GENERATED_MESH_ELLIPSOID_OUTER_SURFACE)
              ALLOCATE(SURFACE_NODES((SIZE(NIDX,1))*(SIZE(NIDX,2)-1)+1),STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate NODES array.",err,error,*999)
              j=1
              i=1
              node_counter=node_counter+1
              SURFACE_NODES(node_counter)=COMPONENT_NODE_TO_USER_NUMBER(ELLIPSOID_MESH%GENERATED_MESH,MESH_COMPONENT, &
                  & NIDX(i,j,SIZE(NIDX,3)),err,error)
              DO j=2,SIZE(NIDX,2)
                DO i=1, SIZE(NIDX,1)
                  node_counter=node_counter+1
                  IF (NIDX(i,j,SIZE(NIDX,3))/=0) THEN
                    SURFACE_NODES(node_counter)=COMPONENT_NODE_TO_USER_NUMBER(ELLIPSOID_MESH%GENERATED_MESH,MESH_COMPONENT, &
                        & NIDX(i,j,SIZE(NIDX,3)),err,error)
                  ELSE
                    node_counter=node_counter-1
                  ENDIF
                ENDDO
              ENDDO
              NORMAL_XI=3

            CASE(GENERATED_MESH_ELLIPSOID_TOP_SURFACE)
              ALLOCATE(SURFACE_NODES((SIZE(NIDX,1))*(SIZE(NIDX,3))),STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate NODES array.",err,error,*999)
              DO k=1,SIZE(NIDX,3)
                DO i=1, SIZE(NIDX,1)
                  node_counter=node_counter+1
                  IF (NIDX(i,SIZE(NIDX,2),k)/=0) THEN
                    SURFACE_NODES(node_counter)=COMPONENT_NODE_TO_USER_NUMBER(ELLIPSOID_MESH%GENERATED_MESH,MESH_COMPONENT, &
                        & NIDX(i,SIZE(NIDX,2),k),err,error)
                  ELSE
                    node_counter=node_counter-1
                  ENDIF
                ENDDO
              ENDDO
              NORMAL_XI=2
            CASE DEFAULT
              localError="The specified surface type of "//TRIM(NumberToVString(SURFACE_TYPE,"*",err,error))//" is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          ELSE
            CALL FlagError("Output SURFACE_NODES array is already allocated.",err,error,*999)
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

    CALL Exits("GENERATED_MESH_ELLIPSOID_SURFACE_GET")
    RETURN
999 CALL Errors("GENERATED_MESH_ELLIPSOID_SURFACE_GET",err,error)
    CALL Exits("GENERATED_MESH_ELLIPSOID_SURFACE_GET")
    RETURN 1
  END SUBROUTINE GENERATED_MESH_ELLIPSOID_SURFACE_GET

  !
  !================================================================================================================================
  !

  !>Calculates the mesh topology information for a given regular mesh (Not to be called by user)
  SUBROUTINE GENERATED_MESH_REGULAR_BUILD_NODE_INDICES(NUMBER_ELEMENTS_XI,NUMBER_OF_NODES_XIC,MAXIMUM_EXTENT, &
      & TOTAL_NUMBER_OF_NODES,TOTAL_NUMBER_OF_ELEMENTS,NIDX,EIDX,DELTA,DELTAi,err,error,*)

    ! Argument variables
    INTEGER(INTG),INTENT(IN) :: NUMBER_ELEMENTS_XI(3) !<Specified number of elements in each xi direction
    INTEGER(INTG),INTENT(IN) :: NUMBER_OF_NODES_XIC(3) !<Number of nodes per element in each xi direction for this basis
    REAL(DP),INTENT(IN) :: MAXIMUM_EXTENT(3)         !<width, length and height of regular mesh
    INTEGER(INTG),INTENT(OUT) :: TOTAL_NUMBER_OF_NODES    !<On exit, contains total number of nodes in regular mesh component
    INTEGER(INTG),INTENT(OUT) :: TOTAL_NUMBER_OF_ELEMENTS !<On exit, contains total number of elements in regular mesh
    INTEGER(INTG),ALLOCATABLE,INTENT(OUT) :: NIDX(:,:,:)  !<Mapping array to find a node number for a given (x,y,z)
    INTEGER(INTG),ALLOCATABLE,INTENT(OUT) :: EIDX(:,:,:)  !<Mapping array to find an element number for a given (x,y,z)
    REAL(DP),INTENT(OUT) :: DELTA(3)  !<Step sizes in each of (x,y,z) for elements
    REAL(DP),INTENT(OUT) :: DELTAi(3) !<Step sizes in each of (x,y,z) for node (identical to DELTA if 2 nodes per xi direction)
    INTEGER(INTG) :: ERR !<The error code
    TYPE(VARYING_STRING) :: ERROR !<The error string

    ! Local variables
    INTEGER(INTG) :: xi_idx,ne1,ne2,ne3,nn1,nn2,nn3,NN,NE
    INTEGER(INTG) :: TOTAL_NUMBER_NODES_XI(3) !<Total number of nodes per element in each xi direction for this basis

    CALL Enters("GENERATED_MESH_REGULAR_BUILD_NODE_INDICES",err,error,*999)

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
        ALLOCATE(NIDX(TOTAL_NUMBER_NODES_XI(1),TOTAL_NUMBER_NODES_XI(2),TOTAL_NUMBER_NODES_XI(3)),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate NIDX array.",err,error,*999)
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
        ALLOCATE(EIDX(NUMBER_ELEMENTS_XI(1),NUMBER_ELEMENTS_XI(2),NUMBER_ELEMENTS_XI(3)),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate EIDX array.",err,error,*999)
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
        CALL FlagError("NIDX array is already allocated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("EIDX array is already allocated.",ERR,error,*999)
    ENDIF

    CALL Exits("GENERATED_MESH_REGULAR_BUILD_NODE_INDICES")
    RETURN
999 CALL Errors("GENERATED_MESH_REGULAR_BUILD_NODE_INDICES",err,error)
    CALL Exits("GENERATED_MESH_REGULAR_BUILD_NODE_INDICES")
    RETURN 1
  END SUBROUTINE GENERATED_MESH_REGULAR_BUILD_NODE_INDICES

  !
  !================================================================================================================================
  !

  !>Calculates the mesh topology information for a given cylinder (Not to be called by user)
  SUBROUTINE GENERATED_MESH_CYLINDER_BUILD_NODE_INDICES(NUMBER_ELEMENTS_XI,NUMBER_OF_NODES_XIC,CYLINDER_EXTENT, &
      & TOTAL_NUMBER_OF_NODES,TOTAL_NUMBER_OF_ELEMENTS,NIDX,EIDX,DELTA,DELTAi,err,error,*)

    ! Argument variables
    INTEGER(INTG),INTENT(IN) :: NUMBER_ELEMENTS_XI(3) !<Specified number of elements in each xi direction
    INTEGER(INTG),INTENT(IN) :: NUMBER_OF_NODES_XIC(3) !<Number of nodes per element in each xi direction (basis property)
    REAL(DP),INTENT(IN) :: CYLINDER_EXTENT(3)         !<inner & outer radii and height of cylinder
    INTEGER(INTG),INTENT(OUT) :: TOTAL_NUMBER_OF_NODES    !<On exit, contains total number of nodes in cylinder mesh
    INTEGER(INTG),INTENT(OUT) :: TOTAL_NUMBER_OF_ELEMENTS !<On exit, contains total number of elements in cylinder mesh
    INTEGER(INTG),ALLOCATABLE,INTENT(OUT) :: NIDX(:,:,:)  !<Mapping array to find a node number for a given (r,theta,z)
    INTEGER(INTG),ALLOCATABLE,INTENT(OUT) :: EIDX(:,:,:)  !<Mapping array to find an element number for a given (r,theta,z)
    REAL(DP),INTENT(OUT) :: DELTA(3)  !<Step sizes in each of (r,theta,z) for elements
    REAL(DP),INTENT(OUT) :: DELTAi(3) !<Step sizes in each of (r,theta,z) for node (identical to DELTA if 2 nodes per xi direction)
    INTEGER(INTG) :: ERR !<The error code
    TYPE(VARYING_STRING) :: ERROR !<The error string

    ! Local variables
    INTEGER(INTG) :: xi_idx,ne1,ne2,ne3,nn1,nn2,nn3,NN,NE
    INTEGER(INTG) :: TOTAL_NUMBER_NODES_XI(3) ! total number of nodes in each xi direction

    CALL Enters("GENERATED_MESH_CYLINDER_BUILD_NODE_INDICES",err,error,*999)

    ! Can skip most of the testing as this subroutine is only to be called by
    ! GENERATED_MESH_CYLINDER_CREATE_FINISH, which tests the input params.
    IF(.NOT.ALLOCATED(NIDX)) THEN
      IF(.NOT.ALLOCATED(EIDX)) THEN
        ! calculate DELTA and DELTAi
        DELTA(1)=(CYLINDER_EXTENT(2)-CYLINDER_EXTENT(1))/NUMBER_ELEMENTS_XI(1)
        DELTA(2)=TWOPI/NUMBER_ELEMENTS_XI(2)
        DELTA(3)=CYLINDER_EXTENT(3)/NUMBER_ELEMENTS_XI(3)
        DO xi_idx=1,3
          DELTAi(xi_idx)=DELTA(xi_idx)/(NUMBER_OF_NODES_XIC(xi_idx)-1)
        ENDDO

        ! calculate total elements and nodes
        DO xi_idx=1,3
          TOTAL_NUMBER_NODES_XI(xi_idx)=(NUMBER_OF_NODES_XIC(xi_idx)-1)*NUMBER_ELEMENTS_XI(xi_idx)+1
        ENDDO
        TOTAL_NUMBER_NODES_XI(2)=TOTAL_NUMBER_NODES_XI(2)-1 ! theta loops around so slightly different
        !TOTAL_NUMBER_OF_ELEMENTS=PRODUCT(NUMBER_ELEMENTS_XI)

        ! calculate NIDX first
        ALLOCATE(NIDX(TOTAL_NUMBER_NODES_XI(1),TOTAL_NUMBER_NODES_XI(2),TOTAL_NUMBER_NODES_XI(3)),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate NIDX array.",err,error,*999)
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
        ALLOCATE(EIDX(NUMBER_ELEMENTS_XI(1),NUMBER_ELEMENTS_XI(2),NUMBER_ELEMENTS_XI(3)),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate EIDX array.",err,error,*999)
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
        CALL FlagError("NIDX array is already allocated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("EIDX array is already allocated.",ERR,error,*999)
    ENDIF

    CALL Exits("GENERATED_MESH_CYLINDER_BUILD_NODE_INDICES")
    RETURN
999 CALL Errors("GENERATED_MESH_CYLINDER_BUILD_NODE_INDICES",err,error)
    CALL Exits("GENERATED_MESH_CYLINDER_BUILD_NODE_INDICES")
    RETURN 1
  END SUBROUTINE GENERATED_MESH_CYLINDER_BUILD_NODE_INDICES

  !
  !================================================================================================================================
  !

  !>Calculate the mesh topology information for a given ellipsoid (Not to be called by user)
  SUBROUTINE GENERATED_MESH_ELLIPSOID_BUILD_NODE_INDICES(NUMBER_ELEMENTS_XI,NUMBER_OF_NODES_XI,ELLIPSOID_EXTENT, &
    & TOTAL_NUMBER_OF_NODES,TOTAL_NUMBER_OF_ELEMENTS,NIDX,CORNER_NODES,EIDX,DELTA,DELTAi,err,error,*)

    ! Argument variables
    INTEGER(INTG),INTENT(IN) :: NUMBER_ELEMENTS_XI(3) !<Specified number of elements in each xi direction
    INTEGER(INTG),INTENT(IN) :: NUMBER_OF_NODES_XI(3) !<Number of nodes per element in each xi direction (basis property)
    REAL(DP),INTENT(IN) :: ELLIPSOID_EXTENT(4)         !< long axis, short axis, wall thickness, top angle
    INTEGER(INTG),INTENT(OUT) :: TOTAL_NUMBER_OF_NODES    !<On exit, contains total number of nodes in ellipsoid mesh
    INTEGER(INTG),INTENT(OUT) :: TOTAL_NUMBER_OF_ELEMENTS !<On exit, contains total number of elements in ellipsoid mesh
    INTEGER(INTG),ALLOCATABLE,INTENT(OUT) :: NIDX(:,:,:)  !<Mapping array to find a node number for a given (r,theta,z)
    INTEGER(INTG),ALLOCATABLE,INTENT(OUT) :: CORNER_NODES(:,:,:) ! Returns the array of corner nodes numbered
    INTEGER(INTG),ALLOCATABLE,INTENT(OUT) :: EIDX(:,:,:)  !<Mapping array to find an element number for a given (r,theta,z)
    REAL(DP),INTENT(OUT) :: DELTA(3)  !<Step sizes in each of (r,theta,z) for elements
    REAL(DP),INTENT(OUT) :: DELTAi(3) !<Step sizes in each of (r,theta,z) for node (identical to DELTA if 2 nodes per xi direction)
    INTEGER(INTG) :: ERR !<The error code
    TYPE(VARYING_STRING) :: ERROR !<The error string

    ! Local variables
    INTEGER(INTG) :: xi_idx,ne1,ne2,ne3,nn1,nn2,nn3,tn1,tn2,tn3,NN,NE
    INTEGER(INTG) :: TOTAL_NUMBER_NODES_XI(3) ! total number of nodes in each xi direction

    CALL Enters("GENERATED_MESH_ELLIPSOID_BUILD_NODE_INDICES",err,error,*999)

    ! Can skip most of the testing as this subroutine is only to be called by
    ! GeneratedMeshEllipsoidCreateFinish, which tests the input params.
    IF(.NOT.ALLOCATED(NIDX)) THEN
      IF(.NOT.ALLOCATED(EIDX)) THEN
        ! calculate DELTA and DELTAi
        DELTA(1)=TWOPI/NUMBER_ELEMENTS_XI(1)
        DELTA(2)=(PI-ELLIPSOID_EXTENT(4))/NUMBER_ELEMENTS_XI(2)
        DELTA(3)=ELLIPSOID_EXTENT(3)/NUMBER_ELEMENTS_XI(3)
        DO xi_idx=1,3
          DELTAi(xi_idx)=DELTA(xi_idx)/(NUMBER_OF_NODES_XI(xi_idx)-1)
        ENDDO

        ! calculate total elements and nodes
        DO xi_idx=1,3
          TOTAL_NUMBER_NODES_XI(xi_idx)=(NUMBER_OF_NODES_XI(xi_idx)-1)*NUMBER_ELEMENTS_XI(xi_idx)+1
        ENDDO
        TOTAL_NUMBER_NODES_XI(1)=TOTAL_NUMBER_NODES_XI(1)-1 ! circumferential loops around so slightly different
        TOTAL_NUMBER_OF_ELEMENTS=PRODUCT(NUMBER_ELEMENTS_XI)

        ! calculate NIDX first
        ALLOCATE(NIDX(TOTAL_NUMBER_NODES_XI(1),TOTAL_NUMBER_NODES_XI(2),TOTAL_NUMBER_NODES_XI(3)),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate NIDX array.",err,error,*999)
        ALLOCATE(CORNER_NODES(NUMBER_ELEMENTS_XI(1),NUMBER_ELEMENTS_XI(2)+1,NUMBER_ELEMENTS_XI(3)+1),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate NIDX array.",err,error,*999)

        !nn: node number inside element in certain direction
        !ne: element number in certain direction
        !tn: global node number in certain direction
        !NN: Node counter
        !Due to one more corner node than elements in transmural direction, first shell is taken separatly
        NN=0
        ne3=1
        nn3=1
        !Due to one more corner node than elements in longitudinal direction, apex elements are taken separatly
        ne2=1
        nn2=1
        ne1=1
        nn1=1
        !apex nodes
        NN=NN+1
        tn3=1
        tn2=1
        tn1=1
        NIDX(tn1,tn2,tn3)=NN
        CORNER_NODES(ne1,ne2,ne3)=NN
        DO ne2=1,NUMBER_ELEMENTS_XI(2)
          DO nn2=2,(NUMBER_OF_NODES_XI(2))
            tn2=tn2+1
            tn1=0
            DO ne1=1,NUMBER_ELEMENTS_XI(1)
              DO nn1=1,(NUMBER_OF_NODES_XI(1)-1)
                tn1=tn1+1
                NN=NN+1
                NIDX(tn1,tn2,tn3)=NN
                IF ((nn1==1).AND.(nn2==NUMBER_OF_NODES_XI(2))) THEN
                  CORNER_NODES(ne1,ne2+1,ne3)=NN
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        DO ne3=1,NUMBER_ELEMENTS_XI(3)
          DO nn3=2,NUMBER_OF_NODES_XI(3)
            ne2=1
            nn2=1
            ne1=1
            nn1=1
            !apex nodes
            NN=NN+1
            tn3=tn3+1
            tn2=1
            tn1=1
            NIDX(tn1,tn2,tn3)=NN
            IF (nn3==NUMBER_OF_NODES_XI(3)) THEN
              CORNER_NODES(ne1,ne2,ne3+1)=NN
            ENDIF
            DO ne2=1,NUMBER_ELEMENTS_XI(2)
              DO nn2=2,(NUMBER_OF_NODES_XI(2))
                tn2=tn2+1
                tn1=0
                DO ne1=1,NUMBER_ELEMENTS_XI(1)
                  DO nn1=1,(NUMBER_OF_NODES_XI(1)-1)
                    tn1=tn1+1
                    NN=NN+1
                    NIDX(tn1,tn2,tn3)=NN
                    IF ((nn1==1).AND.(nn3==NUMBER_OF_NODES_XI(3)).AND.(nn2==NUMBER_OF_NODES_XI(2))) THEN
                      CORNER_NODES(ne1,ne2+1,ne3+1)=NN
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        TOTAL_NUMBER_OF_NODES=NN

        ! now do EIDX
        ALLOCATE(EIDX(NUMBER_ELEMENTS_XI(1),NUMBER_ELEMENTS_XI(2),NUMBER_ELEMENTS_XI(3)),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate EIDX array.",err,error,*999)
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
        CALL FlagError("NIDX array is already allocated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("EIDX array is already allocated.",ERR,error,*999)
    ENDIF

    CALL Exits("GENERATED_MESH_ELLIPSOID_BUILD_NODE_INDICES")
    RETURN
999 CALL Errors("GENERATED_MESH_ELLIPSOID_BUILD_NODE_INDICES",err,error)
    CALL Exits("GENERATED_MESH_ELLIPSOID_BUILD_NODE_INDICES")
    RETURN 1
  END SUBROUTINE GENERATED_MESH_ELLIPSOID_BUILD_NODE_INDICES

  !
  !================================================================================================================================
  !

  !>Calculates the user node numbers for an array of nodes numbered using one basis
  SUBROUTINE COMPONENT_NODES_TO_USER_NUMBERS(GENERATED_MESH,BASIS_INDEX,NODE_COMPONENT_NUMBERS, &
      & NODE_USER_NUMBERS,err,error,*)

    TYPE(GeneratedMeshType), POINTER :: GENERATED_MESH   !<A pointer to the generated mesh object
    INTEGER(INTG),INTENT(IN) :: BASIS_INDEX                !<The number of the basis being used
    INTEGER(INTG),INTENT(IN) :: NODE_COMPONENT_NUMBERS(:)  !<The node numbers for this component basis
    INTEGER(INTG),INTENT(INOUT) :: NODE_USER_NUMBERS(:)    !<On return, the corresponding user numbers
    INTEGER(INTG) :: ERR          !<The error code
    TYPE(VARYING_STRING) :: ERROR !<The error string
    !local variables
    INTEGER(INTG) :: node_idx

    CALL Enters("COMPONENT_NODES_TO_USER_NUMBERS",err,error,*999)

    IF(SIZE(NODE_USER_NUMBERS)==SIZE(NODE_COMPONENT_NUMBERS)) THEN
      DO node_idx=1,SIZE(NODE_COMPONENT_NUMBERS)
        NODE_USER_NUMBERS(node_idx)=COMPONENT_NODE_TO_USER_NUMBER(GENERATED_MESH,BASIS_INDEX, &
            NODE_COMPONENT_NUMBERS(node_idx),err,error)
      ENDDO
    ELSE
      CALL FlagError("NODE_COMPONENT_NUMBERS and NODE_USER_NUMBERS arrays have different sizes.",err,error,*999)
    ENDIF

    CALL Exits("COMPONENT_NODES_TO_USER_NUMBERS")
    RETURN
999 CALL Errors("COMPONENT_NODES_TO_USER_NUMBERS",err,error)
    CALL Exits("COMPONENT_NODES_TO_USER_NUMBERS")
    RETURN 1
  END SUBROUTINE COMPONENT_NODES_TO_USER_NUMBERS

  !
  !================================================================================================================================
  !

  !>Calculates the user node number for a node numbered using one basis
  FUNCTION COMPONENT_NODE_TO_USER_NUMBER(GENERATED_MESH,BASIS_INDEX,NODE_COMPONENT_NUMBER,err,error)
    TYPE(GeneratedMeshType), POINTER :: GENERATED_MESH  !<A pointer to the generated mesh object
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
    TYPE(VARYING_STRING) :: localError

    CALL Enters("COMPONENT_NODE_TO_USER_NUMBER",err,error,*999)

    NULLIFY(BASIS)
    NULLIFY(BASES)
    NUM_CORNER_NODES=1
    REMAINDER=NODE_COMPONENT_NUMBER-1 !use zero based numbering
    REMAINDER2=NODE_COMPONENT_NUMBER-1
    COMPONENT_NODE_TO_USER_NUMBER=0
    POS=0
    POS2=0

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      SELECT CASE(GENERATED_MESH%generatedType)
      CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
        IF(ASSOCIATED(GENERATED_MESH%regularMesh)) THEN
          NUM_BASES=SIZE(GENERATED_MESH%regularMesh%bases)
          NUM_DIMS=GENERATED_MESH%regularMesh%meshDimension
          BASES=>GENERATED_MESH%regularMesh%bases
          NUMBER_OF_ELEMENTS_XI=>GENERATED_MESH%regularMesh%numberOfElementsXi
        ELSE
          CALL FlagError("The regular mesh for this generated mesh is not associated.",err,error,*999)
        ENDIF
      CASE(GENERATED_MESH_POLAR_MESH_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
        IF(ASSOCIATED(GENERATED_MESH%cylinderMesh)) THEN
          NUM_BASES=SIZE(GENERATED_MESH%cylinderMesh%bases)
          NUM_DIMS=GENERATED_MESH%cylinderMesh%meshDimension
          BASES=>GENERATED_MESH%cylinderMesh%bases
          NUMBER_OF_ELEMENTS_XI=>GENERATED_MESH%cylinderMesh%numberOfElementsXi
        ELSE
          CALL FlagError("The cylinder mesh for this generated mesh is not associated.",err,error,*999)
        ENDIF
      CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
        IF(ASSOCIATED(GENERATED_MESH%ellipsoidMesh)) THEN
          NUM_BASES=SIZE(GENERATED_MESH%ellipsoidMesh%bases)
          NUM_DIMS=GENERATED_MESH%ellipsoidMesh%meshDimension
          BASES=>GENERATED_MESH%ellipsoidMesh%bases
          NUMBER_OF_ELEMENTS_XI=>GENERATED_MESH%ellipsoidMesh%numberOfElementsXi
        ELSE
        CALL FlagError("The ellipsoid mesh for this generated mesh is not associated.",err,error,*999)
        ENDIF
      CASE DEFAULT
        localError="The generated mesh generated type of "// &
            & TRIM(NumberToVString(GENERATED_MESH%generatedType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      IF(BASIS_INDEX<=NUM_BASES) THEN
        IF(NUM_BASES==1) THEN
          !If is the only basis, don't do anything
          COMPONENT_NODE_TO_USER_NUMBER=NODE_COMPONENT_NUMBER
        ELSE
          TEMP_TERM=1
          NUM_CORNER_NODES=1
          DO ni=1,NUM_DIMS
            NUM_CORNER_NODES=NUM_CORNER_NODES*(NUMBER_OF_ELEMENTS_XI(ni)+1)
            CORNER_NODE_FACTOR(ni)=1
            IF(ni>1) THEN
              TEMP_TERM=TEMP_TERM*(NUMBER_OF_ELEMENTS_XI(ni-1)+1)
              CORNER_NODE_FACTOR(ni)=CORNER_NODE_FACTOR(ni)*TEMP_TERM
            ENDIF
          ENDDO
          !Adjust for other mesh types
          IF(GENERATED_MESH%generatedType==GENERATED_MESH_CYLINDER_MESH_TYPE) THEN
            CORNER_NODE_FACTOR(3)=CORNER_NODE_FACTOR(3)-NUMBER_OF_ELEMENTS_XI(1)-1
            NUM_CORNER_NODES=NUM_CORNER_NODES-(NUMBER_OF_ELEMENTS_XI(1)+1)*(NUMBER_OF_ELEMENTS_XI(3)+1)
          ELSE IF(GENERATED_MESH%generatedType==GENERATED_MESH_ELLIPSOID_MESH_TYPE) THEN
            CORNER_NODE_FACTOR(3)=CORNER_NODE_FACTOR(3)-NUMBER_OF_ELEMENTS_XI(1)-NUMBER_OF_ELEMENTS_XI(2)
            CORNER_NODE_FACTOR(2)=CORNER_NODE_FACTOR(2)-1
            NUM_CORNER_NODES=NUM_CORNER_NODES-(NUMBER_OF_ELEMENTS_XI(2)+1)*(NUMBER_OF_ELEMENTS_XI(3)+1)- &
                & (NUMBER_OF_ELEMENTS_XI(1)-1)*(NUMBER_OF_ELEMENTS_XI(3)+1)
          ENDIF
          NODE_OFFSET=NUM_CORNER_NODES
          IF(GENERATED_MESH%generatedType==GENERATED_MESH_ELLIPSOID_MESH_TYPE) THEN
            !Every second mesh component is the collapsed node version
            STEP=2
          ELSE
            STEP=1
          ENDIF
          DO basis_idx=1,BASIS_INDEX-1,STEP
            BASIS=>BASES(basis_idx)%PTR
            BASIS_NUM_NODES=1
            DO ni=1,NUM_DIMS
              BASIS_NUM_NODES=BASIS_NUM_NODES*(NUMBER_OF_ELEMENTS_XI(ni)*(BASIS%NUMBER_OF_NODES_XIC(ni)-1)+1)
            ENDDO
            !Adjust for other mesh types
            IF(GENERATED_MESH%generatedType==GENERATED_MESH_CYLINDER_MESH_TYPE) THEN
              BASIS_NUM_NODES=BASIS_NUM_NODES-(NUMBER_OF_ELEMENTS_XI(1)+1)*(BASIS%NUMBER_OF_nodes_xic(1)-1)* &
                  & (NUMBER_OF_ELEMENTS_XI(3)+1)*(BASIS%NUMBER_OF_nodes_xic(3)-1)
            ELSE IF(GENERATED_MESH%generatedType==GENERATED_MESH_ELLIPSOID_MESH_TYPE) THEN
              BASIS_NUM_NODES=BASIS_NUM_NODES-(NUMBER_OF_ELEMENTS_XI(2)*(BASIS%NUMBER_OF_NODES_XIC(2)-1)+1)* &
                & (NUMBER_OF_ELEMENTS_XI(3)*(BASIS%NUMBER_OF_NODES_XIC(3)-1)+1)- &
                & (NUMBER_OF_ELEMENTS_XI(1)*(BASIS%NUMBER_OF_NODES_XIC(1)-1)-1)* &
                & (NUMBER_OF_ELEMENTS_XI(3)*(BASIS%NUMBER_OF_NODES_XIC(3)-1)+1)
            ENDIF
            NODE_OFFSET=NODE_OFFSET+BASIS_NUM_NODES-NUM_CORNER_NODES
          ENDDO
          BASIS=>BASES(BASIS_INDEX)%PTR
          TEMP_TERM=1
          DO ni=1,NUM_DIMS
            BASIS_NODE_FACTOR(ni)=1
            BASIS_ELEMENT_FACTOR(ni)=BASIS%NUMBER_OF_NODES_XIC(ni)-1
            IF(ni>1) THEN
              TEMP_TERM=TEMP_TERM*((BASIS%NUMBER_OF_NODES_XIC(ni-1)-1)*NUMBER_OF_ELEMENTS_XI(ni-1)+1)
              BASIS_NODE_FACTOR(ni)=BASIS_NODE_FACTOR(ni)*TEMP_TERM
              BASIS_ELEMENT_FACTOR(ni)=BASIS_ELEMENT_FACTOR(ni)*TEMP_TERM
            ENDIF
          ENDDO
          !Adjust for other mesh types
          IF(GENERATED_MESH%generatedType==GENERATED_MESH_CYLINDER_MESH_TYPE) THEN
            !subtract nodes along line where y wraps around
            BASIS_NODE_FACTOR(3)=BASIS_NODE_FACTOR(3)-NUMBER_OF_ELEMENTS_XI(1)*(BASIS%NUMBER_OF_NODES_XIC(1)-1)-1
            BASIS_ELEMENT_FACTOR(3)=BASIS_ELEMENT_FACTOR(3)-(NUMBER_OF_ELEMENTS_XI(1)* &
                & (BASIS%NUMBER_OF_NODES_XIC(1)-1)+1)*(BASIS%NUMBER_OF_NODES_XIC(3)-1)
          ELSE IF(GENERATED_MESH%generatedType==GENERATED_MESH_ELLIPSOID_MESH_TYPE) THEN
            !subtract missing nodes at apex
            BASIS_NODE_FACTOR(3)=BASIS_NODE_FACTOR(3)-NUMBER_OF_ELEMENTS_XI(1)*(BASIS%NUMBER_OF_NODES_XIC(1)-1)+1
            BASIS_ELEMENT_FACTOR(3)=BASIS_ELEMENT_FACTOR(3)-(NUMBER_OF_ELEMENTS_XI(1)* &
                & (BASIS%NUMBER_OF_NODES_XIC(1)-1)+1)*(BASIS%NUMBER_OF_NODES_XIC(3)-1)
            !subtract nodes along line where x wraps around
            BASIS_NODE_FACTOR(3)=BASIS_NODE_FACTOR(3)-NUMBER_OF_ELEMENTS_XI(2)*(BASIS%NUMBER_OF_NODES_XIC(2)-1)-1
            BASIS_ELEMENT_FACTOR(3)=BASIS_ELEMENT_FACTOR(3)-(NUMBER_OF_ELEMENTS_XI(2)*(BASIS%NUMBER_OF_NODES_XIC(2)-1)-1)* &
                & (BASIS%NUMBER_OF_NODES_XIC(3)-1)
            BASIS_NODE_FACTOR(2)=BASIS_NODE_FACTOR(2)-1
            BASIS_ELEMENT_FACTOR(2)=BASIS_ELEMENT_FACTOR(2)-(BASIS%NUMBER_OF_NODES_XIC(2)-1)
          ENDIF
          !Work out if we have a corner node, otherwise add node numbers used by corners and
          !previous basis interpolations and subtract number of corner nodes used before the
          !given component node number to get the user number
          CORNER_NODE=.TRUE.
          IF(NUM_DIMS>2) THEN
            POS(3)=REMAINDER/BASIS_NODE_FACTOR(3)
            POS2(3)=REMAINDER2/BASIS_ELEMENT_FACTOR(3)
            REMAINDER=MOD(REMAINDER,BASIS_NODE_FACTOR(3))
            REMAINDER2=MOD(REMAINDER2,BASIS_ELEMENT_FACTOR(3))
            IF(MOD(POS(3),BASIS%NUMBER_OF_NODES_XIC(3)-1)/=0) CORNER_NODE=.FALSE.
          ENDIF
          IF(NUM_DIMS>1) THEN
            IF(GENERATED_MESH%generatedType==GENERATED_MESH_ELLIPSOID_MESH_TYPE) THEN
              !Need to account for missing nodes at apex
              IF(REMAINDER>0) THEN
                REMAINDER=REMAINDER+NUMBER_OF_ELEMENTS_XI(1)*(BASIS%NUMBER_OF_NODES_XIC(1)-1)-1
                REMAINDER2=REMAINDER2+NUMBER_OF_ELEMENTS_XI(1)*(BASIS%NUMBER_OF_NODES_XIC(1)-1)-1
              ENDIF
            ENDIF
            POS(2)=REMAINDER/BASIS_NODE_FACTOR(2)
            POS2(2)=REMAINDER2/BASIS_ELEMENT_FACTOR(2)
            REMAINDER=MOD(REMAINDER,BASIS_NODE_FACTOR(2))
            REMAINDER2=MOD(REMAINDER2,BASIS_ELEMENT_FACTOR(2))
            IF(MOD(POS(2),BASIS%NUMBER_OF_NODES_XIC(2)-1)/=0) CORNER_NODE=.FALSE.
          ENDIF
          POS(1)=REMAINDER/BASIS_NODE_FACTOR(1)
          POS2(1)=REMAINDER2/BASIS_ELEMENT_FACTOR(1)
          IF(MOD(POS(1),BASIS%NUMBER_OF_NODES_XIC(1)-1)/=0) CORNER_NODE=.FALSE.
          IF(CORNER_NODE) THEN
            COMPONENT_NODE_TO_USER_NUMBER=POS2(1)*CORNER_NODE_FACTOR(1)+POS2(2)*CORNER_NODE_FACTOR(2)+ &
                & POS2(3)*CORNER_NODE_FACTOR(3)
            IF(GENERATED_MESH%generatedType==GENERATED_MESH_ELLIPSOID_MESH_TYPE.AND.POS2(2)/=0) THEN
              !Subtract off non-existent nodes at apex
              COMPONENT_NODE_TO_USER_NUMBER=COMPONENT_NODE_TO_USER_NUMBER-(NUMBER_OF_ELEMENTS_XI(1)-1)
            ENDIF
            COMPONENT_NODE_TO_USER_NUMBER=COMPONENT_NODE_TO_USER_NUMBER+1
          ELSE
            !subtract previous corner nodes from node offset
            NUM_PREVIOUS_CORNERS=0
            FINISHED_COUNT=.FALSE.
            IF(NUM_DIMS>2) THEN
              IF(MOD(POS(3),BASIS%NUMBER_OF_NODES_XIC(3)-1)/=0) THEN
                NUM_PREVIOUS_CORNERS=NUM_PREVIOUS_CORNERS+CORNER_NODE_FACTOR(3)*(POS2(3)+1)
                FINISHED_COUNT=.TRUE.
              ELSE
                NUM_PREVIOUS_CORNERS=NUM_PREVIOUS_CORNERS+CORNER_NODE_FACTOR(3)*POS2(3)
              ENDIF
            ENDIF
            IF((NUM_DIMS>1) .AND. (FINISHED_COUNT.NEQV..TRUE.)) THEN
              IF(MOD(POS(2),BASIS%NUMBER_OF_NODES_XIC(2)-1)/=0) THEN
                NUM_PREVIOUS_CORNERS=NUM_PREVIOUS_CORNERS+CORNER_NODE_FACTOR(2)*(POS2(2)+1)
                FINISHED_COUNT=.TRUE.
              ELSE
                NUM_PREVIOUS_CORNERS=NUM_PREVIOUS_CORNERS+CORNER_NODE_FACTOR(2)*POS2(2)
              ENDIF
              IF(GENERATED_MESH%generatedType==GENERATED_MESH_ELLIPSOID_MESH_TYPE) THEN
                NUM_PREVIOUS_CORNERS=NUM_PREVIOUS_CORNERS-(NUMBER_OF_ELEMENTS_XI(1)-1)
              ENDIF
            ENDIF
            IF(FINISHED_COUNT.NEQV..TRUE.) THEN
              IF(MOD(POS(1),BASIS%NUMBER_OF_NODES_XIC(1)-1)/=0) THEN
                NUM_PREVIOUS_CORNERS=NUM_PREVIOUS_CORNERS+CORNER_NODE_FACTOR(1)*(POS2(1)+1)
              ELSE
                NUM_PREVIOUS_CORNERS=NUM_PREVIOUS_CORNERS+CORNER_NODE_FACTOR(1)*POS2(1)
              ENDIF
            ENDIF
            NODE_OFFSET=NODE_OFFSET-NUM_PREVIOUS_CORNERS
            COMPONENT_NODE_TO_USER_NUMBER=NODE_OFFSET+NODE_COMPONENT_NUMBER
          ENDIF
        ENDIF
      ELSE
        localError="Mesh component must be less than or equal to "//(NumberToVString(NUM_BASES,"*",err,error))// &
            & " but it is "//(NumberToVString(BASIS_INDEX,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Generated mesh is not associated",err,error,*999)
    ENDIF

    CALL Exits("COMPONENT_NODE_TO_USER_NUMBER")
    RETURN
999 CALL Errors("COMPONENT_NODE_TO_USER_NUMBER",err,error)
    CALL Exits("COMPONENT_NODE_TO_USER_NUMBER")
    RETURN
  END FUNCTION COMPONENT_NODE_TO_USER_NUMBER

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

  SUBROUTINE GENERATED_MESH_REGULAR_COMPONENT_NODES_TO_USER_NUMBERS(GENERATED_MESH,BASIS_INDEX, &
      & NODE_COMPONENT_NUMBERS,NODE_USER_NUMBERS,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: GENERATED_MESH   !<A pointer to the generated mesh object
    INTEGER(INTG),INTENT(IN) :: BASIS_INDEX                !<The number of the basis being used
    INTEGER(INTG),INTENT(IN) :: NODE_COMPONENT_NUMBERS(:)  !<The node numbers for this component basis
    INTEGER(INTG),INTENT(INOUT) :: NODE_USER_NUMBERS(:)    !<On return, the corresponding user numbers
    INTEGER(INTG) :: ERR          !<The error code
    TYPE(VARYING_STRING) :: ERROR !<The error string
    !Local variables

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

    CALL Enters("GENERATED_MESH_REGULAR_COMPONENT_NODES_TO_USER_NUMBERS",err,error,*999)

    IF(SIZE(NODE_USER_NUMBERS)==SIZE(NODE_COMPONENT_NUMBERS)) THEN
      NODE_USER_NUMBERS=0
      IF(ASSOCIATED(GENERATED_MESH)) THEN
        IF(ASSOCIATED(GENERATED_MESH%regularMesh)) THEN
          NUM_BASES=SIZE(GENERATED_MESH%regularMesh%bases)
          NUM_DIMS=GENERATED_MESH%regularMesh%meshDimension
          BASES=>GENERATED_MESH%regularMesh%bases
          NUMBER_OF_ELEMENTS_XI=1
          DO xi_idx=1,NUM_DIMS
            NUMBER_OF_ELEMENTS_XI(xi_idx)=GENERATED_MESH%regularMesh%numberOfElementsXi(xi_idx)
          ENDDO
        ELSE
        CALL FlagError("The regular mesh for this generated mesh is not associated.",err,error,*999)
        ENDIF

        !Number of nodes in each xi direction
        NUMBER_OF_NODES_XIC=1
        DO xi_idx=1,NUM_DIMS
          NUMBER_OF_NODES_XIC(xi_idx)=BASES(BASIS_INDEX)%PTR%NUMBER_OF_NODES_XIC(xi_idx)
        ENDDO

        !Calculate current element indices and number
        REMINDER_TEMP=0;
        ELEM_IDX=1;
        SELECT CASE(NUM_DIMS)
        CASE(1)
          !Calculate xi1 element index
          ELEM_IDX(1)=(NODE_COMPONENT_NUMBERS(1)-1)/(NUMBER_OF_NODES_XIC(1)-1)+1
          !Calculate element number
          ELEMENT_NO=ELEM_IDX(1)
        CASE(2)
          !Calculate xi2 element index
          NUMBER_OF_NODES_LAYER=((NUMBER_OF_NODES_XIC(1)-1)*NUMBER_OF_ELEMENTS_XI(1)+1)*(NUMBER_OF_NODES_XIC(2)-1)
          ELEM_IDX(2)=NODE_COMPONENT_NUMBERS(1)/NUMBER_OF_NODES_LAYER+1
          REMINDER_TEMP=MOD(NODE_COMPONENT_NUMBERS(1),NUMBER_OF_NODES_LAYER)
          !Calculate xi1 element index
          ELEM_IDX(1)=(REMINDER_TEMP-1)/(NUMBER_OF_NODES_XIC(1)-1)+1
          !Calculate element number
          ELEMENT_NO=(ELEM_IDX(2)-1)*NUMBER_OF_ELEMENTS_XI(1)+ELEM_IDX(1)
        CASE(3)
          !Calculate xi3 element index
          NUMBER_OF_NODES_LAYER=((NUMBER_OF_NODES_XIC(1)-1)*NUMBER_OF_ELEMENTS_XI(1)+1)*((NUMBER_OF_NODES_XIC(2)-1)* &
            & NUMBER_OF_ELEMENTS_XI(2)+1)*(NUMBER_OF_NODES_XIC(3)-1)
          ELEM_IDX(3)=NODE_COMPONENT_NUMBERS(1)/NUMBER_OF_NODES_LAYER+1
         REMINDER_TEMP=MOD(NODE_COMPONENT_NUMBERS(1),NUMBER_OF_NODES_LAYER)
          !Calculate xi2 element index
          NUMBER_OF_NODES_LAYER=((NUMBER_OF_NODES_XIC(1)-1)*NUMBER_OF_ELEMENTS_XI(1)+1)*(NUMBER_OF_NODES_XIC(2)-1)
          ELEM_IDX(2)=REMINDER_TEMP/NUMBER_OF_NODES_LAYER+1
          REMINDER_TEMP=MOD(REMINDER_TEMP,NUMBER_OF_NODES_LAYER)
          !Calculate xi1 element index
          ELEM_IDX(1)=(REMINDER_TEMP-1)/(NUMBER_OF_NODES_XIC(1)-1)+1
          !Calculate element number
          ELEMENT_NO=(ELEM_IDX(3)-1)*NUMBER_OF_ELEMENTS_XI(1)*NUMBER_OF_ELEMENTS_XI(2)+ &
            & (ELEM_IDX(2)-1)*NUMBER_OF_ELEMENTS_XI(1)+ELEM_IDX(1)
        END SELECT


        !If not the first basis, check if previous basis have same interpolation order in each xi direction
        !SAME_BASIS(3) is initialised to have zeros in all entries. If an interpolation scheme has been
        !found to have appeared in previous basis, then record the basis number in the corresponding
        !xi direction. e.g. First basis: bi-quadratic, Second basis: quadratic-cubic, then SAME_BASIS(3)
        !for the second basis will be [1,0,0]
        SAME_BASIS=0
        DO xi_idx=1,NUM_DIMS
          DO basis_idx=1,BASIS_INDEX-1
            IF(BASES(BASIS_INDEX)%PTR%NUMBER_OF_NODES_XIC(xi_idx)== &
              & BASES(basis_idx)%PTR%NUMBER_OF_NODES_XIC(xi_idx)) THEN
              SAME_BASIS(xi_idx)=basis_idx
            ENDIF
          ENDDO
        ENDDO
        !Check if the interpolation scheme has appeared in previous basis
        BASIS_APPEARED=.FALSE.
        IF(SAME_BASIS(1)/=0) THEN
          SELECT CASE(NUM_DIMS)
          CASE(1)
            BASIS_APPEARED=.TRUE.
          CASE(2)
            IF(SAME_BASIS(1)==SAME_BASIS(2)) BASIS_APPEARED=.TRUE.
          CASE(3)
            IF(SAME_BASIS(1)==SAME_BASIS(2) .AND. SAME_BASIS(1)==SAME_BASIS(3)) THEN
             BASIS_APPEARED=.TRUE.
            ENDIF
          END SELECT
        ENDIF
        IF(BASIS_INDEX==1) THEN
          !If this is the first basis, don't do anything
          DO node_idx=1,SIZE(NODE_COMPONENT_NUMBERS)
            NODE_USER_NUMBERS(node_idx)=NODE_COMPONENT_NUMBERS(node_idx)
          ENDDO
        ELSEIF(BASIS_APPEARED) THEN
          !If the basis has appeared before, reuse node user numbers
          DO node_idx=1,SIZE(NODE_COMPONENT_NUMBERS)
            NODE_USER_NUMBERS(node_idx)=GENERATED_MESH%MESH%TOPOLOGY(SAME_BASIS(1))% &
              & PTR%ELEMENTS%ELEMENTS(ELEMENT_NO)%USER_ELEMENT_NODES(node_idx)
          ENDDO
        ELSE
          !If the basis has never appeared exactly in previous basis

          !Find corner node user number from the first basis
          BASIS_FIRST_COMP=>BASES(1)%PTR
          DO nn3=1,2
            DO nn2=1,2
              DO nn1=1,2
                NODE_IDX_CUR=nn1
                NODE_IDX_FIRST=nn1
                IF(nn1==2) THEN
                  NODE_IDX_CUR=NUMBER_OF_NODES_XIC(1)
                  NODE_IDX_FIRST=BASIS_FIRST_COMP%NUMBER_OF_NODES_XIC(1)
                ENDIF
                IF(NUM_DIMS>1 .AND. nn2==2) THEN
                  NODE_IDX_CUR=NODE_IDX_CUR+(NUMBER_OF_NODES_XIC(2)-1)*NUMBER_OF_NODES_XIC(1)
                  NODE_IDX_FIRST=NODE_IDX_FIRST+(BASIS_FIRST_COMP%NUMBER_OF_NODES_XIC(2)-1)* &
                    & BASIS_FIRST_COMP%NUMBER_OF_NODES_XIC(1)
                ENDIF
                IF(NUM_DIMS>2 .AND. nn3==2) THEN
                  NODE_IDX_CUR=NODE_IDX_CUR+NUMBER_OF_NODES_XIC(1)* &
                    & NUMBER_OF_NODES_XIC(2)*(NUMBER_OF_NODES_XIC(3)-1)
                  NODE_IDX_FIRST=NODE_IDX_FIRST+BASIS_FIRST_COMP%NUMBER_OF_NODES_XIC(1)* &
                    & BASIS_FIRST_COMP%NUMBER_OF_NODES_XIC(2)*(BASIS_FIRST_COMP%NUMBER_OF_NODES_XIC(3)-1)
                ENDIF
                NODE_USER_NUMBERS(NODE_IDX_CUR)=GENERATED_MESH%MESH%TOPOLOGY(1)%PTR%ELEMENTS% &
                & ELEMENTS(ELEMENT_NO)%GLOBAL_ELEMENT_NODES(NODE_IDX_FIRST)
              ENDDO
            ENDDO
          ENDDO

          !Find edge node user number from previous basis
          IF(SAME_BASIS(1)/=0 .AND. NUM_DIMS>1) THEN !Do not consider 1D since it's a complete new basis
            BASIS_PRE=>BASES(SAME_BASIS(1))%PTR
            DO nn3=1,2
              DO nn2=1,2
                DO nn1=2,NUMBER_OF_NODES_XIC(1)-1
                  NODE_IDX_CUR=nn1
                  NODE_IDX_PRE=nn1
                  IF(nn2==2) THEN
                    NODE_IDX_CUR=NODE_IDX_CUR+(NUMBER_OF_NODES_XIC(2)-1)*NUMBER_OF_NODES_XIC(1)
                    NODE_IDX_PRE=NODE_IDX_PRE+(BASIS_PRE%NUMBER_OF_NODES_XIC(2)-1)*BASIS_PRE%NUMBER_OF_NODES_XIC(1)
                  ENDIF
                  IF(NUM_DIMS>2 .AND. nn3==2) THEN
                    NODE_IDX_CUR=NODE_IDX_CUR+NUMBER_OF_NODES_XIC(1)*NUMBER_OF_NODES_XIC(2)* &
                      & (NUMBER_OF_NODES_XIC(3)-1)
                    NODE_IDX_PRE=NODE_IDX_PRE+BASIS_PRE%NUMBER_OF_NODES_XIC(1)*BASIS_PRE% &
                      & NUMBER_OF_NODES_XIC(2)*(BASIS_PRE%NUMBER_OF_NODES_XIC(3)-1)
                  ENDIF
                  NODE_USER_NUMBERS(NODE_IDX_CUR)=GENERATED_MESH%MESH%TOPOLOGY(SAME_BASIS(1))% &
                    & PTR%ELEMENTS%ELEMENTS(ELEMENT_NO)%GLOBAL_ELEMENT_NODES(NODE_IDX_PRE)
                ENDDO
              ENDDO
            ENDDO
          ENDIF
          IF(SAME_BASIS(2)/=0) THEN
            BASIS_PRE=>BASES(SAME_BASIS(2))%PTR
            DO nn3=1,2
              DO nn2=2,NUMBER_OF_NODES_XIC(2)-1
                DO nn1=1,2
                  IF(nn1==1) THEN
                    NODE_IDX_CUR=nn1+(nn2-1)*NUMBER_OF_NODES_XIC(1)
                    NODE_IDX_PRE=nn1+(nn2-1)*BASIS_PRE%NUMBER_OF_NODES_XIC(1)
                  ELSE
                    NODE_IDX_CUR=nn2*NUMBER_OF_NODES_XIC(1)
                    NODE_IDX_PRE=nn2*BASIS_PRE%NUMBER_OF_NODES_XIC(1)
                  ENDIF
                  IF(NUM_DIMS>2 .AND. nn3==2) THEN
                    NODE_IDX_CUR=NODE_IDX_CUR+NUMBER_OF_NODES_XIC(1)*NUMBER_OF_NODES_XIC(2)* &
                      & (NUMBER_OF_NODES_XIC(3)-1)
                    NODE_IDX_PRE=NODE_IDX_PRE+BASIS_PRE%NUMBER_OF_NODES_XIC(1)*BASIS_PRE% &
                      & NUMBER_OF_NODES_XIC(2)*(BASIS_PRE%NUMBER_OF_NODES_XIC(3)-1)
                  ENDIF
                  NODE_USER_NUMBERS(NODE_IDX_CUR)=GENERATED_MESH%MESH%TOPOLOGY(SAME_BASIS(2))% &
                    & PTR%ELEMENTS%ELEMENTS(ELEMENT_NO)%GLOBAL_ELEMENT_NODES(NODE_IDX_PRE)
                ENDDO
              ENDDO
            ENDDO
          ENDIF
          IF(SAME_BASIS(3)/=0) THEN !Must be 3D
            BASIS_PRE=>BASES(SAME_BASIS(3))%PTR
            NODE_IDX_CUR=0
            NODE_IDX_PRE=0
            DO nn3=2,NUMBER_OF_NODES_XIC(3)-1
              DO nn2=1,2
                IF(nn2==2) THEN
                  NODE_IDX_CUR=(NUMBER_OF_NODES_XIC(2)-1)*NUMBER_OF_NODES_XIC(1)+NUMBER_OF_NODES_XIC(1)* &
                    & NUMBER_OF_NODES_XIC(2)*(NUMBER_OF_NODES_XIC(3)-1)
                  NODE_IDX_PRE=(BASIS_PRE%NUMBER_OF_NODES_XIC(1)-1)*BASIS_PRE%NUMBER_OF_NODES_XIC(1)+ &
                    & BASIS_PRE%NUMBER_OF_NODES_XIC(1)*BASIS_PRE%NUMBER_OF_NODES_XIC(2)* &
                    & (BASIS_PRE%NUMBER_OF_NODES_XIC(3)-1)
                ENDIF
                DO nn1=1,2
                  IF(nn1==1) THEN
                    NODE_IDX_CUR=1+NODE_IDX_CUR
                    NODE_IDX_PRE=1+NODE_IDX_PRE
                  ELSE
                    NODE_IDX_CUR=NUMBER_OF_NODES_XIC(1)+NODE_IDX_CUR
                    NODE_IDX_PRE=BASIS_PRE%NUMBER_OF_NODES_XIC(1)+NODE_IDX_PRE
                  ENDIF
                  NODE_USER_NUMBERS(NODE_IDX_CUR)=GENERATED_MESH%MESH%TOPOLOGY(SAME_BASIS(3))% &
                    & PTR%ELEMENTS%ELEMENTS(ELEMENT_NO)%GLOBAL_ELEMENT_NODES(NODE_IDX_PRE)
                ENDDO
              ENDDO
            ENDDO
          ENDIF
          !The following code would only be executed if 3D (automatically satisfied, don't need to check,
          !since there must be at least 1 direction that has different interpolation scheme, if two direction
          ! has the same interpolation that has appeared before, then interpolation for the last direction
          ! must be different) and has same basis in 2 xi direction
          !i.e. find user node numbers for face nodes
          IF(SAME_BASIS(1)==SAME_BASIS(2) .AND. SAME_BASIS(1)/=0) THEN
            BASIS_PRE=>BASES(SAME_BASIS(1))%PTR
            DO nn3=1,2
              DO nn2=2,NUMBER_OF_NODES_XIC(2)-1
                DO nn1=2,NUMBER_OF_NODES_XIC(1)-1
                  NODE_IDX_CUR=nn1+(nn2-1)*NUMBER_OF_NODES_XIC(1)
                  NODE_IDX_PRE=nn1+(nn2-1)*BASIS_PRE%NUMBER_OF_NODES_XIC(1)
                  IF(nn3==2) THEN
                    NODE_IDX_CUR=NODE_IDX_CUR+NUMBER_OF_NODES_XIC(1)*NUMBER_OF_NODES_XIC(2)* &
                      & (NUMBER_OF_NODES_XIC(3)-1)
                    NODE_IDX_PRE=NODE_IDX_PRE+BASIS_PRE%NUMBER_OF_NODES_XIC(1)*BASIS_PRE% &
                      & NUMBER_OF_NODES_XIC(2)*(BASIS_PRE%NUMBER_OF_NODES_XIC(3)-1)
                  ENDIF
                  NODE_USER_NUMBERS(NODE_IDX_CUR)=GENERATED_MESH%MESH%TOPOLOGY(SAME_BASIS(1))% &
                    & PTR%ELEMENTS%ELEMENTS(ELEMENT_NO)%GLOBAL_ELEMENT_NODES(NODE_IDX_PRE)
                ENDDO
              ENDDO
            ENDDO
          ELSE IF(SAME_BASIS(1)==SAME_BASIS(3) .AND. SAME_BASIS(1)/=0) THEN
            BASIS_PRE=>BASES(SAME_BASIS(1))%PTR
            NODE_IDX_CUR=0
            NODE_IDX_PRE=0
            DO nn3=2,NUMBER_OF_NODES_XIC(3)-1
              DO nn2=1,2
                IF(nn2==2) THEN
                  NODE_IDX_CUR=(NUMBER_OF_NODES_XIC(2)-1)*NUMBER_OF_NODES_XIC(1)+NUMBER_OF_NODES_XIC(1)* &
                    & NUMBER_OF_NODES_XIC(2)*(nn3-1)
                  NODE_IDX_PRE=(BASIS_PRE%NUMBER_OF_NODES_XIC(2)-1)*BASIS_PRE%NUMBER_OF_NODES_XIC(1)+ &
                    & BASIS_PRE%NUMBER_OF_NODES_XIC(1)*BASIS_PRE%NUMBER_OF_NODES_XIC(2)*(nn3-1)
                ENDIF
                DO nn1=2,NUMBER_OF_NODES_XIC(1)-1
                  NODE_IDX_CUR=nn1+NODE_IDX_CUR
                  NODE_IDX_PRE=nn1+NODE_IDX_PRE
                  NODE_USER_NUMBERS(NODE_IDX_CUR)=GENERATED_MESH%MESH%TOPOLOGY(SAME_BASIS(1))% &
                    & PTR%ELEMENTS%ELEMENTS(ELEMENT_NO)%GLOBAL_ELEMENT_NODES(NODE_IDX_PRE)
                ENDDO
              ENDDO
            ENDDO
          ELSE IF(SAME_BASIS(2)==SAME_BASIS(3) .AND. SAME_BASIS(2)/=0) THEN
            BASIS_PRE=>BASES(SAME_BASIS(2))%PTR
            DO nn3=2,NUMBER_OF_NODES_XIC(3)-1
              DO nn2=2,NUMBER_OF_NODES_XIC(2)-1
                DO nn1=1,2
                  IF(nn1==1) THEN
                    NODE_IDX_CUR=1+(nn2-1)*NUMBER_OF_NODES_XIC(1)+NUMBER_OF_NODES_XIC(1)* &
                      & NUMBER_OF_NODES_XIC(2)*(nn3-1)
                    NODE_IDX_PRE=1+(nn2-1)*BASIS_PRE%NUMBER_OF_NODES_XIC(1)+BASIS_PRE%NUMBER_OF_NODES_XIC(1)* &
                      & BASIS_PRE%NUMBER_OF_NODES_XIC(2)*(nn3-1)
                  ELSE
                    NODE_IDX_CUR=nn2*NUMBER_OF_NODES_XIC(1)+NUMBER_OF_NODES_XIC(1)* &
                      & NUMBER_OF_NODES_XIC(2)*(nn3-1)
                    NODE_IDX_PRE=nn2*BASIS_PRE%NUMBER_OF_NODES_XIC(1)+BASIS_PRE%NUMBER_OF_NODES_XIC(1)* &
                      & BASIS_PRE%NUMBER_OF_NODES_XIC(2)*(nn3-1)
                  ENDIF
                  NODE_USER_NUMBERS(NODE_IDX_CUR)=GENERATED_MESH%MESH%TOPOLOGY(SAME_BASIS(2))% &
                    & PTR%ELEMENTS%ELEMENTS(ELEMENT_NO)%GLOBAL_ELEMENT_NODES(NODE_IDX_PRE)
                ENDDO
              ENDDO
            ENDDO
          ENDIF

          !Find the largest node user number in the previous basis
          NODE_OFFSET_LAST_BASIS=0
          LAST_ELEM_NO=GENERATED_MESH%MESH%TOPOLOGY(1)%PTR%ELEMENTS%NUMBER_OF_ELEMENTS !The mesh has the same topology regardless of mesh components
          DO basis_idx=1,BASIS_INDEX-1
          number_of_nodes_temp=SIZE(GENERATED_MESH%MESH%TOPOLOGY(basis_idx)%PTR%ELEMENTS% &
            & ELEMENTS(LAST_ELEM_NO)%GLOBAL_ELEMENT_NODES,1)
            DO node_index_temp=1,number_of_nodes_temp
              IF (GENERATED_MESH%MESH%TOPOLOGY(basis_idx)%PTR%ELEMENTS%ELEMENTS(LAST_ELEM_NO)% &
                & GLOBAL_ELEMENT_NODES(node_index_temp)>NODE_OFFSET_LAST_BASIS) THEN
                NODE_OFFSET_LAST_BASIS=GENERATED_MESH%MESH%TOPOLOGY(basis_idx)%PTR%ELEMENTS%ELEMENTS(LAST_ELEM_NO)% &
                  &GLOBAL_ELEMENT_NODES(node_index_temp)
              ENDIF
            ENDDO !node_index_temp
          ENDDO !basis_idx

          !Calculate number of zeros nodes in different dimensions
          INDEX_COUNT=1
          ZERO_COUNT_XI1=0
          ZERO_COUNT_XI12=0
          TOTAL_ZERO_NODE=0
          EDGE_NODE=0
          DO nn3=1,NUMBER_OF_NODES_XIC(3)
            DO nn2=1,NUMBER_OF_NODES_XIC(2)
              NODE_COUNT=0
              DO nn1=1,NUMBER_OF_NODES_XIC(1)
                NODE_IDX=(nn3-1)*NUMBER_OF_NODES_XIC(1)*NUMBER_OF_NODES_XIC(2)+(nn2-1)* &
                  & NUMBER_OF_NODES_XIC(1)+nn1
                IF(NODE_USER_NUMBERS(NODE_IDX)==0) THEN
                  NODE_COUNT=NODE_COUNT+1
                  TOTAL_ZERO_NODE=TOTAL_ZERO_NODE+1 !Total number of zeros in an element
                ENDIF
              ENDDO !nn1
              ZERO_COUNT_XI1(INDEX_COUNT)=NODE_COUNT !Total number of zero summed up across xi1 direction.
              IF(NODE_COUNT==NUMBER_OF_NODES_XIC(1)) EDGE_NODE(INDEX_COUNT)=1 !Shared edge node (with zero value) in xi1 direction (1 number for each node in xi2 direction)
              ZERO_COUNT_XI12(nn3)=ZERO_COUNT_XI12(nn3)+ZERO_COUNT_XI1(INDEX_COUNT) !Total number of zero summed on xi1-xi2 faces
              INDEX_COUNT=INDEX_COUNT+1
            ENDDO !nn2
          ENDDO !nn3

         !Calculate how many zero nodes has occurred in previous elements
         NODE_OFFSET_ELEM=0
          IF(NUM_DIMS==2 .AND. ELEM_IDX(2)/=1) THEN !Zero nodes occurred in the previous rows of elements
            OFFSET_UNIT=TOTAL_ZERO_NODE-ZERO_COUNT_XI1(1)-SUM(EDGE_NODE(1:NUMBER_OF_NODES_XIC(2)))+EDGE_NODE(INDEX_COUNT)
            !This is number of zero nodes in the elements before the current row of elements
            NODE_OFFSET_ELEM=(ELEM_IDX(2)-1)*NUMBER_OF_ELEMENTS_XI(1)*OFFSET_UNIT+(ELEM_IDX(2)-1)* &
              & SUM(EDGE_NODE(2:NUMBER_OF_NODES_XIC(2)-1))
          ELSEIF(NUM_DIMS==3 .AND. ELEM_IDX(3)/=1) THEN !Zero nodes occurred in the previous layer of elements
            NODE_OFFSET_XI3_ACCUM=0
            DO nn3=1,NUMBER_OF_NODES_XIC(3)-1
              OFFSET_UNIT=ZERO_COUNT_XI12(nn3)-ZERO_COUNT_XI1((nn3-1)*NUMBER_OF_NODES_XIC(2)+1)- &
                & SUM(EDGE_NODE((nn3-1)*NUMBER_OF_NODES_XIC(2)+1:nn3*NUMBER_OF_NODES_XIC(2)))+ &
                & EDGE_NODE((nn3-1)*NUMBER_OF_NODES_XIC(2)+1)
              NODE_OFFSET_XI3_ACCUM=NODE_OFFSET_XI3_ACCUM+OFFSET_UNIT*NUMBER_OF_ELEMENTS_XI(1)*NUMBER_OF_ELEMENTS_XI(2)+ &
                & (NUMBER_OF_ELEMENTS_XI(1)-1)*(ZERO_COUNT_XI1((nn3-1)*NUMBER_OF_NODES_XIC(2)+1)- &
                & EDGE_NODE((nn3-1)*NUMBER_OF_NODES_XIC(2)+1))+ZERO_COUNT_XI1((nn3-1)*NUMBER_OF_NODES_XIC(2)+1)+ &
                & SUM(EDGE_NODE((nn3-1)*NUMBER_OF_NODES_XIC(2)+2:nn3*NUMBER_OF_NODES_XIC(2)))* &
                & NUMBER_OF_ELEMENTS_XI(2)
            ENDDO
            NODE_OFFSET_ELEM=(ELEM_IDX(3)-1)*NODE_OFFSET_XI3_ACCUM
          ENDIF

          !Compute other nodes which haven't appeared in previous basis
          INDEX_COUNT=1
          NODE_OFFSET_ELEM_XI12=0
          NODE_OFFSET_XI2=0 !Number of zero nodes in the current row
          NODE_OFFSET_XI3_ACCUM=0 !Number of zero nodes in the layers in xi3 direction (nn3)
          DO nn3=1,NUMBER_OF_NODES_XIC(3)
            NODE_OFFSET_XI2_ACCUM=0 !Number of zero nodes in the previous rows
            OFFSET_UNIT=ZERO_COUNT_XI12(nn3)-ZERO_COUNT_XI1((nn3-1)*NUMBER_OF_NODES_XIC(2)+1)- &
                & SUM(EDGE_NODE((nn3-1)*NUMBER_OF_NODES_XIC(2)+1:nn3*NUMBER_OF_NODES_XIC(2)))+ &
                & EDGE_NODE((nn3-1)*NUMBER_OF_NODES_XIC(2)+1)
            IF(ELEM_IDX(2)/=1 .AND. NUM_DIMS==3) THEN
              NODE_OFFSET_ELEM_XI12=OFFSET_UNIT*(ELEM_IDX(2)-1)*NUMBER_OF_ELEMENTS_XI(1)+ &
                & (ELEM_IDX(2)-1)*SUM(EDGE_NODE((nn3-1)*NUMBER_OF_NODES_XIC(2)+2:nn3*NUMBER_OF_NODES_XIC(2)))
            ENDIF
            DO nn2=1,NUMBER_OF_NODES_XIC(2)
              NODE_OFFSET_XI2=(ZERO_COUNT_XI1(INDEX_COUNT)-EDGE_NODE(INDEX_COUNT))*(ELEM_IDX(1)-1)
              NODE_OFFSET=NODE_OFFSET_LAST_BASIS+NODE_OFFSET_ELEM+NODE_OFFSET_XI3_ACCUM+ &
                & NODE_OFFSET_ELEM_XI12+NODE_OFFSET_XI2_ACCUM+NODE_OFFSET_XI2
              DO nn1=1,NUMBER_OF_NODES_XIC(1)
                !Local node index in the current element
                NODE_IDX=(nn3-1)*NUMBER_OF_NODES_XIC(1)*NUMBER_OF_NODES_XIC(2)+(nn2-1)* &
                  & NUMBER_OF_NODES_XIC(1)+nn1
                IF(NODE_USER_NUMBERS(NODE_IDX)==0) THEN
                  !This is for 2D case
                  NODE_OFFSET=NODE_OFFSET+1
                  NODE_USER_NUMBERS(NODE_IDX)=NODE_OFFSET
                ENDIF
              ENDDO !nn1
              NODE_OFFSET_XI2_ACCUM=NODE_OFFSET_XI2_ACCUM+(ZERO_COUNT_XI1(INDEX_COUNT)-EDGE_NODE(INDEX_COUNT))* &
                & NUMBER_OF_ELEMENTS_XI(1)+EDGE_NODE(INDEX_COUNT)
              INDEX_COUNT=INDEX_COUNT+1
            ENDDO !nn2
            IF(NUM_DIMS==3) THEN
              NODE_OFFSET_XI3_ACCUM=NODE_OFFSET_XI3_ACCUM+OFFSET_UNIT*NUMBER_OF_ELEMENTS_XI(1)*NUMBER_OF_ELEMENTS_XI(2)+ &
                & (NUMBER_OF_ELEMENTS_XI(1)-1)*(ZERO_COUNT_XI1((nn3-1)*NUMBER_OF_NODES_XIC(2)+1)- &
                & EDGE_NODE((nn3-1)*NUMBER_OF_NODES_XIC(2)+1))+ZERO_COUNT_XI1((nn3-1)*NUMBER_OF_NODES_XIC(2)+1)+ &
                & SUM(EDGE_NODE((nn3-1)*NUMBER_OF_NODES_XIC(2)+2:nn3*NUMBER_OF_NODES_XIC(2)))* &
                & NUMBER_OF_ELEMENTS_XI(2)
            ENDIF
          ENDDO !nn3
        ENDIF
      ELSE
        CALL FlagError("Generated mesh is not associated",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("NODE_COMPONENT_NUMBERS and NODE_USER_NUMBERS arrays have different sizes.",err,error,*999)
    ENDIF
    CALL Exits("GENERATED_MESH_REGULAR_COMPONENT_NODES_TO_USER_NUMBERS")
    RETURN
999 CALL Errors("GENERATED_MESH_REGULAR_COMPONENT_NODES_TO_USER_NUMBERS",err,error)
    CALL Exits("GENERATED_MESH_REGULAR_COMPONENT_NODES_TO_USER_NUMBERS")
    RETURN 1
  END SUBROUTINE GENERATED_MESH_REGULAR_COMPONENT_NODES_TO_USER_NUMBERS

  !
  !================================================================================================================================
  !

  !>Retrieve the user node number for a component number in a regular generated mesh
  !>This routine only works for Lagrange/Hermite elements
  SUBROUTINE GENERATED_MESH_REGULAR_COMPONENT_NODE_TO_USER_NUMBER(GENERATED_MESH,BASIS_INDEX, &
      & NODE_COMPONENT_NUMBER,NODE_USER_NUMBER,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: GENERATED_MESH  !<A pointer to the generated mesh object
    INTEGER(INTG),INTENT(IN) :: BASIS_INDEX  !<The number of the basis being used
    INTEGER(INTG),INTENT(IN) :: NODE_COMPONENT_NUMBER  !<The node numbers for this component basis
    INTEGER(INTG),INTENT(OUT) :: NODE_USER_NUMBER  !<On return, the corresponding user numbers
    INTEGER(INTG) :: ERR  !<The error code
    TYPE(VARYING_STRING) :: ERROR  !<The error string

    !Local variables
    TYPE(BASIS_PTR_TYPE), POINTER :: BASES(:)
    INTEGER(INTG) :: NUM_BASES,NUM_DIMS,ELEMENT_NO,LOCAL_NODE_NO,NUMBER_OF_NODES_LAYER,xi_idx
    INTEGER(INTG) :: ELEM_IDX(3),NODE_IDX(3),NUMBER_OF_NODES_XIC(3),NUMBER_OF_ELEMENTS_XI(3),REMINDER_TEMP

    CALL Enters("GENERATED_MESH_REGULAR_COMPONENT_NODE_TO_USER_NUMBER",err,error,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      IF(ASSOCIATED(GENERATED_MESH%regularMesh)) THEN
        NUM_BASES=SIZE(GENERATED_MESH%regularMesh%bases)
        NUM_DIMS=GENERATED_MESH%regularMesh%meshDimension
        BASES=>GENERATED_MESH%regularMesh%bases
        NUMBER_OF_ELEMENTS_XI=1
        DO xi_idx=1,NUM_DIMS
          NUMBER_OF_ELEMENTS_XI(xi_idx)=GENERATED_MESH%regularMesh%numberOfElementsXi(xi_idx)
        ENDDO
        !Number of nodes in each xi direction
        NUMBER_OF_NODES_XIC=1
        DO xi_idx=1,NUM_DIMS
          NUMBER_OF_NODES_XIC(xi_idx)=BASES(BASIS_INDEX)%PTR%NUMBER_OF_NODES_XIC(xi_idx)
        ENDDO
      ELSE
        CALL FlagError("The regular mesh for this generated mesh is not associated.",err,error,*999)
      ENDIF

      !Calculate current element/node indices/number
      REMINDER_TEMP=0;
      ELEM_IDX=1;
      NODE_IDX=1;
      SELECT CASE(NUM_DIMS)
      CASE(1)
        !Calculate xi1 element index
        ELEM_IDX(1)=(NODE_COMPONENT_NUMBER-1)/(NUMBER_OF_NODES_XIC(1)-1)+1
        NODE_IDX(1)=MOD(NODE_COMPONENT_NUMBER-1,NUMBER_OF_NODES_XIC(1)-1)+1
        !If it's the last node in the line
        IF (ELEM_IDX(1)>NUMBER_OF_ELEMENTS_XI(1)) THEN
          ELEM_IDX(1)=ELEM_IDX(1)-1
          NODE_IDX(1)=NUMBER_OF_NODES_XIC(1)
        ENDIF
        !Calculate element number
        ELEMENT_NO=ELEM_IDX(1)
        LOCAL_NODE_NO=NODE_IDX(1)
      CASE(2)
        !Calculate xi2 element index
        NUMBER_OF_NODES_LAYER=((NUMBER_OF_NODES_XIC(1)-1)*NUMBER_OF_ELEMENTS_XI(1)+1)*(NUMBER_OF_NODES_XIC(2)-1)
        ELEM_IDX(2)=(NODE_COMPONENT_NUMBER-1)/NUMBER_OF_NODES_LAYER+1
        REMINDER_TEMP=MOD(NODE_COMPONENT_NUMBER-1,NUMBER_OF_NODES_LAYER)
        NUMBER_OF_NODES_LAYER=((NUMBER_OF_NODES_XIC(1)-1)*NUMBER_OF_ELEMENTS_XI(1)+1)
        NODE_IDX(2)=REMINDER_TEMP/NUMBER_OF_NODES_LAYER+1
        !If it's the last line of nodes in the line
        IF (ELEM_IDX(2)>NUMBER_OF_ELEMENTS_XI(2)) THEN
          ELEM_IDX(2)=ELEM_IDX(2)-1
          NODE_IDX(2)=NUMBER_OF_NODES_XIC(2)
        ENDIF
        !Calculate xi1 element index
        REMINDER_TEMP=MOD(REMINDER_TEMP,NUMBER_OF_NODES_LAYER)
        ELEM_IDX(1)=REMINDER_TEMP/(NUMBER_OF_NODES_XIC(1)-1)+1
        NODE_IDX(1)=MOD(REMINDER_TEMP,NUMBER_OF_NODES_XIC(1)-1)+1
        !If it's the last node in the line
        IF (ELEM_IDX(1)>NUMBER_OF_ELEMENTS_XI(1)) THEN
          ELEM_IDX(1)=ELEM_IDX(1)-1
          NODE_IDX(1)=NUMBER_OF_NODES_XIC(1)
        ENDIF
        !Calculate element number
        ELEMENT_NO=(ELEM_IDX(2)-1)*NUMBER_OF_ELEMENTS_XI(1)+ELEM_IDX(1)
        LOCAL_NODE_NO=(NODE_IDX(2)-1)*NUMBER_OF_NODES_XIC(1)+NODE_IDX(1)
      CASE(3)
        !Calculate xi3 element index
        NUMBER_OF_NODES_LAYER=((NUMBER_OF_NODES_XIC(1)-1)*NUMBER_OF_ELEMENTS_XI(1)+1)*((NUMBER_OF_NODES_XIC(2)-1)* &
          & NUMBER_OF_ELEMENTS_XI(2)+1)*(NUMBER_OF_NODES_XIC(3)-1) !Multiple planes of nodes
        ELEM_IDX(3)=(NODE_COMPONENT_NUMBER-1)/NUMBER_OF_NODES_LAYER+1
        REMINDER_TEMP=MOD(NODE_COMPONENT_NUMBER-1,NUMBER_OF_NODES_LAYER) !Multiple planes of nodes
        NUMBER_OF_NODES_LAYER=((NUMBER_OF_NODES_XIC(1)-1)*NUMBER_OF_ELEMENTS_XI(1)+1)*((NUMBER_OF_NODES_XIC(2)-1)* &
          & NUMBER_OF_ELEMENTS_XI(2)+1) !One plane of nodes
        NODE_IDX(3)=REMINDER_TEMP/NUMBER_OF_NODES_LAYER+1
        IF (ELEM_IDX(3)>NUMBER_OF_ELEMENTS_XI(3)) THEN
          ELEM_IDX(3)=ELEM_IDX(3)-1
          NODE_IDX(3)=NUMBER_OF_NODES_XIC(3)
        ENDIF
        REMINDER_TEMP=MOD(REMINDER_TEMP,NUMBER_OF_NODES_LAYER) !One plane of nodes
        !Calculate xi2 element index
        NUMBER_OF_NODES_LAYER=((NUMBER_OF_NODES_XIC(1)-1)*NUMBER_OF_ELEMENTS_XI(1)+1)*(NUMBER_OF_NODES_XIC(2)-1) !Multiple lines of nodes
        ELEM_IDX(2)=REMINDER_TEMP/NUMBER_OF_NODES_LAYER+1
        REMINDER_TEMP=MOD(REMINDER_TEMP,NUMBER_OF_NODES_LAYER) !Multiple lines of nodes
        NUMBER_OF_NODES_LAYER=(NUMBER_OF_NODES_XIC(1)-1)*NUMBER_OF_ELEMENTS_XI(1)+1 !One line of nodes
        NODE_IDX(2)=REMINDER_TEMP/NUMBER_OF_NODES_LAYER+1
        REMINDER_TEMP=MOD(REMINDER_TEMP,NUMBER_OF_NODES_LAYER) !One line of nodes
        IF (ELEM_IDX(2)>NUMBER_OF_ELEMENTS_XI(2)) THEN
          ELEM_IDX(2)=ELEM_IDX(2)-1
          NODE_IDX(2)=NUMBER_OF_NODES_XIC(2)
        ENDIF
        !Calculate xi1 element index
        ELEM_IDX(1)=REMINDER_TEMP/(NUMBER_OF_NODES_XIC(1)-1)+1
        NODE_IDX(1)=MOD(REMINDER_TEMP,NUMBER_OF_NODES_XIC(1)-1)+1
        IF (ELEM_IDX(1)>NUMBER_OF_ELEMENTS_XI(1)) THEN
          ELEM_IDX(1)=ELEM_IDX(1)-1
          NODE_IDX(1)=NUMBER_OF_NODES_XIC(1)
        ENDIF
        !Calculate element number
        ELEMENT_NO=(ELEM_IDX(3)-1)*NUMBER_OF_ELEMENTS_XI(1)*NUMBER_OF_ELEMENTS_XI(2)+ &
          & (ELEM_IDX(2)-1)*NUMBER_OF_ELEMENTS_XI(1)+ELEM_IDX(1)
        LOCAL_NODE_NO=(NODE_IDX(3)-1)*NUMBER_OF_NODES_XIC(1)*NUMBER_OF_NODES_XIC(2)+(NODE_IDX(2)-1)*NUMBER_OF_NODES_XIC(1)+ &
          & NODE_IDX(1)
      END SELECT
      !Retrieve node user number
      IF(ASSOCIATED(GENERATED_MESH%MESH)) THEN
        NODE_USER_NUMBER=GENERATED_MESH%MESH%TOPOLOGY(BASIS_INDEX)%PTR%ELEMENTS%ELEMENTS(ELEMENT_NO)% &
          & USER_ELEMENT_NODES(LOCAL_NODE_NO)
      ELSE
        CALL FlagError("The mesh for this generated mesh is not associated.",err,error,*999)
      ENDIF

    ELSE
        CALL FlagError("Generated mesh is not associated",err,error,*999)
    ENDIF

    CALL Exits("GENERATED_MESH_REGULAR_COMPONENT_NODE_TO_USER_NUMBER")
    RETURN
999 CALL Errors("GENERATED_MESH_REGULAR_COMPONENT_NODE_TO_USER_NUMBER",err,error)
    CALL Exits("GENERATED_MESH_REGULAR_COMPONENT_NODE_TO_USER_NUMBER")
    RETURN 1
  END SUBROUTINE GENERATED_MESH_REGULAR_COMPONENT_NODE_TO_USER_NUMBER

  !
  !================================================================================================================================
  !

  !>Calculates the user node number for a node numbered using one basis.
  !>This is currently only used for cylinder meshes, other mesh types don't require this.
  FUNCTION USER_NUMBER_TO_COMPONENT_NODE(GENERATED_MESH,BASIS_INDEX,NODE_USER_NUMBER,err,error)
    TYPE(GeneratedMeshType), POINTER :: GENERATED_MESH        !<A pointer to the generated mesh object
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
    TYPE(VARYING_STRING) :: localError

    CALL Enters("USER_NUMBER_TO_COMPONENT_NODE",err,error,*999)

    NULLIFY(BASIS)
    NULLIFY(BASES)
    NUM_CORNER_NODES=1
    REMAINDER=NODE_USER_NUMBER-1 !use zero based numbering
    POS=0

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      !Only cylinder mesh type uses this now, although it was previously used by regular
      !meshes so some things relate to that.
      SELECT CASE(GENERATED_MESH%generatedType)
      CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(GENERATED_MESH_POLAR_MESH_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
        IF(ASSOCIATED(GENERATED_MESH%cylinderMesh)) THEN
          NUM_BASES=SIZE(GENERATED_MESH%cylinderMesh%bases)
          NUM_DIMS=GENERATED_MESH%cylinderMesh%meshDimension
          BASES=>GENERATED_MESH%cylinderMesh%bases
          NUMBER_OF_ELEMENTS_XI=>GENERATED_MESH%cylinderMesh%numberOfElementsXi
        ELSE
          CALL FlagError("The cylinder mesh for this generated mesh is not associated.",err,error,*999)
        ENDIF
      CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="The generated mesh generated type of "// &
            & TRIM(NumberToVString(GENERATED_MESH%generatedType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      IF(BASIS_INDEX<=NUM_BASES) THEN
        IF(NUM_BASES==1) THEN
          !If is the only basis, don't do anything
          USER_NUMBER_TO_COMPONENT_NODE=NODE_USER_NUMBER
        ELSE
          TEMP_TERM=1
          NUM_CORNER_NODES=1
          DO ni=1,NUM_DIMS
            NUM_CORNER_NODES=NUM_CORNER_NODES*(NUMBER_OF_ELEMENTS_XI(ni)+1)
            CORNER_NODE_FACTOR(ni)=1
            IF(ni>1) THEN
              TEMP_TERM=TEMP_TERM*(NUMBER_OF_ELEMENTS_XI(ni-1)+1)
              CORNER_NODE_FACTOR(ni)=CORNER_NODE_FACTOR(ni)*TEMP_TERM
            ENDIF
          ENDDO
          !Adjust for other mesh types
          IF(GENERATED_MESH%generatedType==GENERATED_MESH_CYLINDER_MESH_TYPE) THEN
            CORNER_NODE_FACTOR(3)=CORNER_NODE_FACTOR(3)-NUMBER_OF_ELEMENTS_XI(1)-1
            NUM_CORNER_NODES=NUM_CORNER_NODES-(NUMBER_OF_ELEMENTS_XI(1)+1)*(NUMBER_OF_ELEMENTS_XI(3)+1)
          ENDIF
          NODE_OFFSET=NUM_CORNER_NODES
          DO basis_idx=1,BASIS_INDEX-1
            BASIS=>BASES(basis_idx)%PTR
            BASIS_NUM_NODES=1
            DO ni=1,NUM_DIMS
              BASIS_NUM_NODES=BASIS_NUM_NODES*(NUMBER_OF_ELEMENTS_XI(ni)*(BASIS%NUMBER_OF_NODES_XIC(ni)-1)+1)
            ENDDO
            !Adjust for other mesh types
            IF(GENERATED_MESH%generatedType==GENERATED_MESH_CYLINDER_MESH_TYPE) THEN
              BASIS_NUM_NODES=BASIS_NUM_NODES-(NUMBER_OF_ELEMENTS_XI(1)+1)*(BASIS%NUMBER_OF_nodes_xic(1)-1)* &
                  & (NUMBER_OF_ELEMENTS_XI(3)+1)*(BASIS%NUMBER_OF_nodes_xic(3)-1)
            ENDIF
            NODE_OFFSET=NODE_OFFSET+BASIS_NUM_NODES-NUM_CORNER_NODES
          ENDDO
          BASIS=>BASES(BASIS_INDEX)%PTR
          TEMP_TERM=1
          DO ni=1,NUM_DIMS
            BASIS_ELEMENT_FACTOR(ni)=BASIS%NUMBER_OF_NODES_XIC(ni)-1
            IF(ni>1) THEN
              TEMP_TERM=TEMP_TERM*((BASIS%NUMBER_OF_NODES_XIC(ni-1)-1)*NUMBER_OF_ELEMENTS_XI(ni-1)+1)
              BASIS_ELEMENT_FACTOR(ni)=BASIS_ELEMENT_FACTOR(ni)*TEMP_TERM
            ENDIF
          ENDDO
          !Adjust for other mesh types
          IF(GENERATED_MESH%generatedType==GENERATED_MESH_CYLINDER_MESH_TYPE) THEN
            !subtract nodes along line where y wraps around
            BASIS_ELEMENT_FACTOR(3)=BASIS_ELEMENT_FACTOR(3)-(NUMBER_OF_ELEMENTS_XI(1)* &
                & (BASIS%NUMBER_OF_NODES_XIC(1)-1)+1)*(BASIS%NUMBER_OF_NODES_XIC(3)-1)
          ENDIF
          IF(NODE_USER_NUMBER<=NUM_CORNER_NODES) THEN
            !we have a node on a corner
            IF(NUM_DIMS>2) THEN
              POS(3)=REMAINDER/CORNER_NODE_FACTOR(3)
              REMAINDER=MOD(REMAINDER,CORNER_NODE_FACTOR(3))
            ENDIF
            IF(NUM_DIMS>1) THEN
              POS(2)=REMAINDER/CORNER_NODE_FACTOR(2)
              REMAINDER=MOD(REMAINDER,CORNER_NODE_FACTOR(2))
            ENDIF
            POS(1)=REMAINDER/CORNER_NODE_FACTOR(1)
            USER_NUMBER_TO_COMPONENT_NODE=POS(1)*BASIS_ELEMENT_FACTOR(1)+POS(2)*BASIS_ELEMENT_FACTOR(2)+ &
                & POS(3)*BASIS_ELEMENT_FACTOR(3)
            USER_NUMBER_TO_COMPONENT_NODE=USER_NUMBER_TO_COMPONENT_NODE+1
          ELSE IF(NODE_USER_NUMBER>NODE_OFFSET) THEN
            REMAINDER=REMAINDER-NODE_OFFSET
            DO ni=1,NUM_DIMS
              BASIS_ELEMENT_FACTOR(ni)=BASIS_ELEMENT_FACTOR(ni)-CORNER_NODE_FACTOR(ni)
            ENDDO
            NUM_PREVIOUS_CORNERS=0
            FINISHED_COUNT=.FALSE.
            OFF_EDGE=.FALSE.
            IF(NUM_DIMS>2) THEN
              IF(GENERATED_MESH%generatedType==GENERATED_MESH_CYLINDER_MESH_TYPE.AND. &
                  & (MOD(REMAINDER,BASIS_ELEMENT_FACTOR(3)) > BASIS_ELEMENT_FACTOR(2)*NUMBER_OF_ELEMENTS_XI(2)-1)) THEN
                OFF_EDGE=.TRUE.
              ELSE IF(GENERATED_MESH%generatedType==GENERATED_MESH_REGULAR_MESH_TYPE.AND. &
                  & MOD(REMAINDER,BASIS_ELEMENT_FACTOR(3)) > (BASIS_ELEMENT_FACTOR(2)*NUMBER_OF_ELEMENTS_XI(2)+ &
                  & BASIS_ELEMENT_FACTOR(1)*NUMBER_OF_ELEMENTS_XI(1)-1)) THEN
                OFF_EDGE=.TRUE.
              ENDIF
              IF(OFF_EDGE) THEN
                NUM_PREVIOUS_CORNERS=NUM_PREVIOUS_CORNERS+CORNER_NODE_FACTOR(3)*(1+REMAINDER/BASIS_ELEMENT_FACTOR(3))
                REMAINDER=MOD(REMAINDER,BASIS_ELEMENT_FACTOR(3))
                FINISHED_COUNT=.TRUE.
              ELSE
                NUM_PREVIOUS_CORNERS=NUM_PREVIOUS_CORNERS+CORNER_NODE_FACTOR(3)*(REMAINDER/BASIS_ELEMENT_FACTOR(3))
                REMAINDER=MOD(REMAINDER,BASIS_ELEMENT_FACTOR(3))
              ENDIF
            ENDIF
            IF((NUM_DIMS>1) .AND. (FINISHED_COUNT.NEQV..TRUE.)) THEN
              IF(MOD(REMAINDER,BASIS_ELEMENT_FACTOR(2)) > &
                  & BASIS_ELEMENT_FACTOR(1)*NUMBER_OF_ELEMENTS_XI(1)-1) THEN
                NUM_PREVIOUS_CORNERS=NUM_PREVIOUS_CORNERS+CORNER_NODE_FACTOR(2)*(1+REMAINDER/BASIS_ELEMENT_FACTOR(2))
                REMAINDER=MOD(REMAINDER,BASIS_ELEMENT_FACTOR(2))
                FINISHED_COUNT=.TRUE.
              ELSE
                NUM_PREVIOUS_CORNERS=NUM_PREVIOUS_CORNERS+CORNER_NODE_FACTOR(2)*(REMAINDER/BASIS_ELEMENT_FACTOR(2))
                REMAINDER=MOD(REMAINDER,BASIS_ELEMENT_FACTOR(2))
              ENDIF
            ENDIF
            IF(FINISHED_COUNT.NEQV..TRUE.) THEN
              NUM_PREVIOUS_CORNERS=NUM_PREVIOUS_CORNERS+CORNER_NODE_FACTOR(1)*(REMAINDER/BASIS_ELEMENT_FACTOR(1))+1
            ENDIF
            NODE_OFFSET=NODE_OFFSET-NUM_PREVIOUS_CORNERS
            USER_NUMBER_TO_COMPONENT_NODE=NODE_USER_NUMBER-NODE_OFFSET
          ELSE
            CALL FlagError("Invalid node number specified.",err,error,*999)
          ENDIF
        ENDIF
      ELSE
        localError="Mesh component must be less than or equal to "//(NumberToVString(NUM_BASES,"*",err,error))// &
            & " but it is "//(NumberToVString(BASIS_INDEX,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Generated mesh is not associated",err,error,*999)
    ENDIF

    CALL Exits("USER_NUMBER_TO_COMPONENT_NODE")
    RETURN
999 CALL Errors("USER_NUMBER_TO_COMPONENT_NODE",err,error)
    CALL Exits("USER_NUMBER_TO_COMPONENT_NODE")
    RETURN
  END FUNCTION USER_NUMBER_TO_COMPONENT_NODE

  !
  !================================================================================================================================
  !

END MODULE GeneratedMeshRoutines

