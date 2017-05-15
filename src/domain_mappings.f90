!> \file
!> \author Chris Bradley
!> \brief This module handles all domain mappings routines.
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

!> This module handles all domain mappings routines.
MODULE DOMAIN_MAPPINGS

#ifndef NOMPIMOD
  use MPI
#endif
  USE BASE_ROUTINES
  USE COMP_ENVIRONMENT
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE KINDS
  USE LISTS
  USE STRINGS
  USE TYPES

#include "macros.h"


  IMPLICIT NONE
#ifdef NOMPIMOD
#include "mpif.h"
#endif

  PRIVATE

  !Module parameters

  !> \addtogroup DOMAIN_MAPPINGS_DomainType DOMAIN_MAPPINGS::DomainType
  !> \see DOMAIN_MAPPINGS
  !>@{
  INTEGER(INTG), PARAMETER :: DOMAIN_LOCAL_INTERNAL=1 !<The domain item is internal to the domain \see DOMAIN_MAPPINGS_DomainType,DOMAIN_MAPPINGS
  INTEGER(INTG), PARAMETER :: DOMAIN_LOCAL_BOUNDARY=2 !<The domain item is on the boundary of the domain \see DOMAIN_MAPPINGS_DomainType,DOMAIN_MAPPINGS
  INTEGER(INTG), PARAMETER :: DOMAIN_LOCAL_GHOST=3 !<The domain item is ghosted from another domain \see DOMAIN_MAPPINGS_DomainType,DOMAIN_MAPPINGS
  !>@}

  !Module types

  !Module variables

  !Interfaces

  PUBLIC DOMAIN_LOCAL_INTERNAL,DOMAIN_LOCAL_BOUNDARY,DOMAIN_LOCAL_GHOST

  PUBLIC DOMAIN_MAPPINGS_MAPPING_FINALISE,DOMAIN_MAPPINGS_MAPPING_INITIALISE,DOMAIN_MAPPINGS_MAPPING_GLOBAL_INITIALISE, &
       & DOMAIN_MAPPINGS_GLOBAL_TO_LOCAL_GET,DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE, &
       & DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE2, fill_domain_mapping

CONTAINS

!> Initialise and fill the components of a new DOMAIN MAPPING instance.
  subroutine fill_domain_mapping( subdomain, map, refMap, numDomains, internalVAL, boundaryVAL, ghostVAL, * )

    ! Argument variables
      type(DOMAIN_MAPPING_TYPE),     pointer :: map !<pointer to the new DOMAIN_MAPPING instance to be filled
      type(DOMAIN_MAPPING_TYPE),     pointer :: refMap !<pointer to a reference DOMAIN MAPPING instance.  This mapping is used to retrieve the # of ADJACENT DOMAINS
      integer(INTG),              intent(in) :: numDomains !<# of active sub-domains the new DOMAIN MAPPING covers
      integer(INTG),              intent(in) :: subdomain
      integer(INTG), allocatable, intent(in) :: internalVAL(:) !<array of INTERNAL mesh parameters (elements,nodes or DOFs) that the new DOMAIN MAPPING contains
      integer(INTG), allocatable, intent(in) :: boundaryVAL(:) !<array of BOUNDARY mesh parameters (elements,nodes or DOFs) that the new DOMAIN MAPPING contains
      integer(INTG), allocatable, intent(in) :: ghostVAL(:) !<array of GHOST mesh parameters (elements,nodes or DOFs) that the new DOMAIN MAPPING contains

    ! Local variables
      integer(INTG) :: m, n, nn, np, err, status(MPI_STATUS_SIZE), domain_no
      integer(INTG), allocatable :: displ(:), recv_cnt(:), temp(:), tmp(:)
      type(VARYING_STRING) :: error

      EXITS( "fill_domain_mapping" )

    ! Input array check
      if ( .not.allocated(internalVAL) ) call FlagError( "INTERNAL values array not allocated", err, error, *999 )
      if ( .not.allocated(boundaryVAL) ) call FlagError( "BOUNDARY values array not allocated", err, error, *999 )
      if ( .not.allocated(ghostVAL) ) call FlagError( "GHOST values array not allocated", err, error, *999 )

    ! set counts 
      map%NUMBER_OF_LOCAL = map%NUMBER_OF_INTERNAL + map%NUMBER_OF_BOUNDARY
      map%TOTAL_NUMBER_OF_LOCAL = map%NUMBER_OF_LOCAL + map%NUMBER_OF_GHOST
      map%NUMBER_OF_DOMAINS = numDomains
 
    ! Set start/end of each type (INTERNAL, BOUNDARY, GHOST)
      map%INTERNAL_START = 1
      map%INTERNAL_FINISH = map%NUMBER_OF_INTERNAL
      map%BOUNDARY_START = map%INTERNAL_FINISH
      if ( map%NUMBER_OF_BOUNDARY>0 ) map%BOUNDARY_START = map%BOUNDARY_START + 1
      map%BOUNDARY_FINISH = map%INTERNAL_FINISH + map%NUMBER_OF_BOUNDARY
      map%GHOST_START = map%BOUNDARY_FINISH + 1
      map%GHOST_FINISH = map%BOUNDARY_FINISH + map%NUMBER_OF_GHOST

    ! Allocate and fill LOCAL_TO_GLOBAL_MAP array
      allocate( map%LOCAL_TO_GLOBAL_MAP(map%TOTAL_NUMBER_OF_LOCAL), STAT=err )
      if ( err/=0 ) call FlagError( "could not allocate LOCAL_TO_GLOBAL_MAP", err, error, *999 )

      map%LOCAL_TO_GLOBAL_MAP( 1:map%NUMBER_OF_INTERNAL ) = internalVAL( 1:map%NUMBER_OF_INTERNAL )
      if ( map%NUMBER_OF_BOUNDARY>0 ) map%LOCAL_TO_GLOBAL_MAP( map%BOUNDARY_START:map%BOUNDARY_FINISH ) &
                                  & = boundaryVAL( 1:map%NUMBER_OF_BOUNDARY )
      map%LOCAL_TO_GLOBAL_MAP( map%GHOST_START:map%GHOST_FINISH ) = ghostVAL( 1:map%NUMBER_OF_GHOST )

    ! Allocate and fill the DOMAIN_LIST array
      allocate( map%DOMAIN_LIST( map%TOTAL_NUMBER_OF_LOCAL), STAT=err )
      if ( err/=0 ) call FlagError( "could not allocate DOMAIN_LIST", err, error, *999 )

      do n = 1,map%TOTAL_NUMBER_OF_LOCAL
         map%DOMAIN_LIST( n ) = n
      enddo

    ! Allocate and fill the LOCAL_TYPE array
      allocate( map%LOCAL_TYPE(map%TOTAL_NUMBER_OF_LOCAL), STAT=err )
      if ( err/=0 ) call FlagError( "could not allocate LOCAL_TYPE array", err, error, *999 )

      map%LOCAL_TYPE( 1:map%INTERNAL_FINISH ) = DOMAIN_LOCAL_INTERNAL
      if ( map%NUMBER_OF_BOUNDARY>0 ) then
         map%LOCAL_TYPE( map%BOUNDARY_START:map%BOUNDARY_FINISH ) = DOMAIN_LOCAL_BOUNDARY
      endif
      map%LOCAL_TYPE( map%GHOST_START:map%GHOST_FINISH ) = DOMAIN_LOCAL_GHOST

!
! PART TWO - DOMAIN ADJACENCY INFO
!            Determine what other sub-domains the local sub-domain must communicate with
!--------------------------------------------------------------------------------------------------------------------------------
    ! construct the ADJACENT_DOMAINS_PTR array
      allocate( map%ADJACENT_DOMAINS_PTR(0:numDomains-1), STAT=err )
      if (err/=0) call FlagError( "could not allocate ADJACENT_DOMAINS_PTR array", err, error, *999 )
      map%ADJACENT_DOMAINS_PTR(:) = refMap%ADJACENT_DOMAINS_PTR(:)

    ! contruct the ADJACENT_DOMAINS_LIST array
      allocate( map%ADJACENT_DOMAINS_LIST( size(refMap%ADJACENT_DOMAINS_LIST) ), STAT=err )
      if (err/=0) call FlagError( "could not allocate ADJACENT_DOMAINS_LIST array", err, error, *999 )
      map%ADJACENT_DOMAINS_LIST(:) = refMap%ADJACENT_DOMAINS_LIST(:)

    ! gather the # of local nodes/elements/DOFs on every sub-domain into a global array
      allocate( recv_cnt(numDomains) )
      allocate( displ(numDomains) )
      allocate( map%NUMBER_OF_DOMAIN_LOCAL(0:numDomains-1), STAT=err )
      if (err/=0) call FlagError( "could not allocate NUMBER_OF_DOMAIN_LOCAL", err, error, *999 )

      recv_cnt(:) = 1
      do n = 1,numDomains
         displ(n) = n - 1
      enddo

      call MPI_Allgatherv( map%NUMBER_OF_LOCAL, 1, MPI_INTEGER, map%NUMBER_OF_DOMAIN_LOCAL, recv_cnt, displ, &
                         & MPI_INTEGER, COMPUTATIONAL_ENVIRONMENT%MPI_COMM, err )

    ! gather the # of ghost nodes/elements/DOFs on every sub-domain into a global array
      allocate( map%NUMBER_OF_DOMAIN_GHOST( 0:numDomains-1 ), STAT=err )
      if (err/=0) call FlagError( "could not allocate NUMBER_OF_DOMAIN_GHOST", err, error, *999 )

      call MPI_Allgatherv( map%NUMBER_OF_GHOST, 1, MPI_INTEGER, map%NUMBER_OF_DOMAIN_GHOST, recv_cnt, displ, &
                         & MPI_INTEGER, COMPUTATIONAL_ENVIRONMENT%MPI_COMM, err )

    ! gather the # of boundary nodes/elements/DOFs on every sub-domain into a global array
      allocate( map%NUMBER_OF_DOMAIN_BOUNDARY(0:numDomains-1), STAT=err )
      if ( err/=0 ) call FlagError( "could not allocate NUMBER_OF_DOMAIN_BOUNDARY", err, error, *999 )
      call MPI_Allgatherv( map%NUMBER_OF_BOUNDARY, 1, MPI_INTEGER, &
                         & map%NUMBER_OF_DOMAIN_BOUNDARY, recv_cnt, displ, MPI_INTEGER, &
                         & COMPUTATIONAL_ENVIRONMENT%MPI_COMM, err )
      deallocate( recv_cnt,displ )

    ! allocate an array of ADJACENT_DOMAINS structures for the element mapping
      map%NUMBER_OF_ADJACENT_DOMAINS = refMap%NUMBER_OF_ADJACENT_DOMAINS
      allocate( map%ADJACENT_DOMAINS( map%NUMBER_OF_ADJACENT_DOMAINS ), STAT=err )
      if (err/=0) call FlagError( "could not allocate ADJACENT_DOMAINS structure", err, error, *999 )

    ! set the domain IDs and initialize the send/receive counts in the array of structures
      do n = 1,map%NUMBER_OF_ADJACENT_DOMAINS
         map%ADJACENT_DOMAINS( n )%DOMAIN_NUMBER = refMap%ADJACENT_DOMAINS(n)%DOMAIN_NUMBER
         map%ADJACENT_DOMAINS( n )%NUMBER_OF_SEND_GHOSTS = 0
         map%ADJACENT_DOMAINS( n )%NUMBER_OF_RECEIVE_GHOSTS = 0
      enddo

  !
  ! PART THREE - HALO EXCHANGE INFO (send)
  !              Determine the local IDs of the global ID values that need to be sent to adjacent sub-domains
  !--------------------------------------------------------------------------------------------------------------------------------
      call MPI_Barrier( COMPUTATIONAL_ENVIRONMENT%MPI_COMM, err )
      do domain_no = 0,numDomains-1
         allocate( temp(map%NUMBER_OF_DOMAIN_GHOST(domain_no) ), STAT=err )
         if (err/=0) call FlagError( "could not allocate temporary array", err, error, *999 )

       ! domain_no will send the global IDs of its GHOST IDs to all sub-domains adjacent to it
         if ( subdomain==domain_no ) then
            temp = ghostVAL( 1:map%NUMBER_OF_GHOST )
            do n = 1,map%NUMBER_OF_ADJACENT_DOMAINS
               call MPI_Send( temp, map%NUMBER_OF_GHOST, MPI_INTEGER, map%ADJACENT_DOMAINS(n)%DOMAIN_NUMBER, 0, &
                            & COMPUTATIONAL_ENVIRONMENT%MPI_COMM, err )
            enddo

       ! the other sub-domains check if they are adjacent to domain_no.
         else
            do n = 1,map%NUMBER_OF_ADJACENT_DOMAINS
               if ( map%ADJACENT_DOMAINS(n)%DOMAIN_NUMBER==domain_no ) then

                ! If they are, they receive the sent GHOST node IDs
                  call MPI_Recv( temp, map%NUMBER_OF_DOMAIN_GHOST(domain_no), MPI_INTEGER, domain_no, MPI_ANY_TAG, &
                               & COMPUTATIONAL_ENVIRONMENT%MPI_COMM, status, err )

                ! Update the NUMBER_OF_SEND_GHOSTS counter
                  allocate( tmp(map%NUMBER_OF_DOMAIN_GHOST(domain_no)), STAT=err )
                  if (err/=0) call FlagError( "could not allocate temporary array", err, error, *999 )

                  do np = 1,map%NUMBER_OF_LOCAL
                  do m = 1,map%NUMBER_OF_DOMAIN_GHOST(domain_no)
                     if ( map%LOCAL_TO_GLOBAL_MAP(np)==temp(m) ) then
                        map%ADJACENT_DOMAINS(n)%NUMBER_OF_SEND_GHOSTS = map%ADJACENT_DOMAINS(n)%NUMBER_OF_SEND_GHOSTS + 1
                        tmp( map%ADJACENT_DOMAINS(n)%NUMBER_OF_SEND_GHOSTS ) = np
                     endif
                  enddo
                  enddo

                ! Allocate memory and fill in the local IDs of the ghost values to send domain_no
                  m = map%ADJACENT_DOMAINS(n)%NUMBER_OF_SEND_GHOSTS
                  allocate( map%ADJACENT_DOMAINS(n)%LOCAL_GHOST_SEND_INDICES( m ), STAT=err )
                  if (err/=0) call FlagError( "could not allocate LOCAL_GHOST_SEND_INDICES", err, error, *999 )

                  map%ADJACENT_DOMAINS(n)%LOCAL_GHOST_SEND_INDICES = tmp( 1:map%ADJACENT_DOMAINS(n)%NUMBER_OF_SEND_GHOSTS )
                  deallocate( tmp )
                  exit

               endif
            enddo
         endif ! subdomain

         deallocate( temp )
      enddo

  !
  ! PART FOUR - HALO EXCHANGE INFO (receive)
  !             Determine the local IDs of the ID values that the local sub-domain receives from adjacent sub-domains
  !--------------------------------------------------------------------------------------------------------------------------------
      do domain_no = 0,numDomains-1
         allocate( temp(map%NUMBER_OF_DOMAIN_LOCAL(domain_no)) )

       ! domain_no will send the global IDs of its LOCAL nodes to all sub-domains adjacent to it
         if ( subdomain==domain_no ) then
            temp( 1:map%NUMBER_OF_INTERNAL ) = internalVAL( 1:map%NUMBER_OF_INTERNAL )
            temp( map%BOUNDARY_START:map%BOUNDARY_FINISH ) = boundaryVAL( 1:map%NUMBER_OF_BOUNDARY )
            do n = 1,map%NUMBER_OF_ADJACENT_DOMAINS
               call MPI_Send( temp, map%NUMBER_OF_DOMAIN_LOCAL(domain_no), MPI_INTEGER, &
                           &  map%ADJACENT_DOMAINS(n)%DOMAIN_NUMBER, 0, COMPUTATIONAL_ENVIRONMENT%MPI_COMM, err )
            enddo
         else

          ! the other sub-domains check if they are adjacent to domain_no.
            do n = 1,map%NUMBER_OF_ADJACENT_DOMAINS
               if ( map%ADJACENT_DOMAINS(n)%DOMAIN_NUMBER==domain_no ) then

                ! If they are, they receive the sent LOCAL node IDs
                  call MPI_Recv( temp, map%NUMBER_OF_DOMAIN_LOCAL(domain_no), MPI_INTEGER, domain_no, MPI_ANY_TAG, &
                               & COMPUTATIONAL_ENVIRONMENT%MPI_COMM, status, err )

                ! Update the NUMBER_OF_RECEIVE_GHOSTS counter
                  allocate( tmp(map%NUMBER_OF_DOMAIN_LOCAL(domain_no)), STAT=err )

                  do np = map%GHOST_START,map%GHOST_FINISH
                  do m = 1,map%NUMBER_OF_DOMAIN_LOCAL(domain_no)
                     if ( map%LOCAL_TO_GLOBAL_MAP(map%DOMAIN_LIST(np))==temp(m) ) then
                        map%ADJACENT_DOMAINS(n)%NUMBER_OF_RECEIVE_GHOSTS = map%ADJACENT_DOMAINS(n)%NUMBER_OF_RECEIVE_GHOSTS + 1
                        tmp(map%ADJACENT_DOMAINS(n)%NUMBER_OF_RECEIVE_GHOSTS) = map%DOMAIN_LIST(np)
                        exit
                     endif
                  enddo
                  enddo

                ! Allocate memory and fill in the local IDs of the ghost values to send domain_no
                  m = map%ADJACENT_DOMAINS(n)%NUMBER_OF_RECEIVE_GHOSTS
                  allocate( map%ADJACENT_DOMAINS(n)%LOCAL_GHOST_RECEIVE_INDICES( m ), STAT=err )
                  if (err/=0) call FlagError( "could not allocate nodal LOCAL_GHOST_RECEIVE_INDICES", err, error, *999 )

                  map%ADJACENT_DOMAINS(n)%LOCAL_GHOST_RECEIVE_INDICES = tmp( 1:m )
                  deallocate( tmp )
                  exit

               endif
            enddo
         endif

         deallocate( temp )
      enddo

   !*******************************************************************************************************************
   !                                      DEBUGGING  /  DIAGNOSTICS
   !*******************************************************************************************************************

      if ( DIAGNOSTICS2 ) then
         if ( DIAGNOSTICS1 ) then
            call WRITE_STRING( DIAGNOSTIC_OUTPUT_TYPE, "Fill Domain:", err, error, *999 )
            call WRITE_STRING_VALUE( DIAGNOSTIC_OUTPUT_TYPE,"  # of INTERNAL parameters = ", &
                                   & map%NUMBER_OF_INTERNAL, err, error, *999 )
            call WRITE_STRING_VALUE( DIAGNOSTIC_OUTPUT_TYPE,"  # of BOUNDARY parameters = ", &
                                   & map%NUMBER_OF_BOUNDARY, err, error, *999 )
            call WRITE_STRING_VALUE( DIAGNOSTIC_OUTPUT_TYPE,"  # of GHOST parameters = ", &
                                   & map%NUMBER_OF_GHOST, err, error, *999 )
         endif
         call WRITE_STRING_VALUE( DIAGNOSTIC_OUTPUT_TYPE,"  # of adjacent domains = ", &
                                & map%NUMBER_OF_ADJACENT_DOMAINS, err, error, *999 )
         do n = 1,map%NUMBER_OF_ADJACENT_DOMAINS
            call WRITE_STRING_VALUE( DIAGNOSTIC_OUTPUT_TYPE,"  Adjacent domain ID = ", &
                                   & map%ADJACENT_DOMAINS(n)%DOMAIN_NUMBER, err, error, *999 )
            call WRITE_STRING_VECTOR( DIAGNOSTIC_OUTPUT_TYPE, 1, 1, &
                                    & map%ADJACENT_DOMAINS(n)%NUMBER_OF_SEND_GHOSTS, 8, 8, &
                                    & map%ADJACENT_DOMAINS(n)%LOCAL_GHOST_SEND_INDICES, &
                                    & '("  Local Ghost Send Indices List :",8(X,I5))', '(33X,8(X,I5))', &
                                    & err, error, *999 )
            call WRITE_STRING_VECTOR( DIAGNOSTIC_OUTPUT_TYPE, 1, 1, &
                                    & map%ADJACENT_DOMAINS(n)%NUMBER_OF_RECEIVE_GHOSTS, 8, 8, &
                                    & map%ADJACENT_DOMAINS(n)%LOCAL_GHOST_RECEIVE_INDICES, &
                                    & '("  Local Ghost Receive Indices List :",8(X,I5))', '(36X,8(X,I5))', &
                                    & err, error, *999 )
         enddo
      endif

      return
 999  ERRORSEXITS( "fill_domain_mapping", err, error )
      return 1
  end subroutine fill_domain_mapping

  !
  !================================================================================================================================
  !

  !>Finalises the adjacent domain and deallocates all memory for a domain mapping.
  SUBROUTINE DOMAIN_MAPPINGS_ADJACENT_DOMAIN_FINALISE(ADJACENT_DOMAIN,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_ADJACENT_DOMAIN_TYPE) :: ADJACENT_DOMAIN !<The adjacent domain to finalise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_MAPPINGS_ADJACENT_DOMAIN_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(ADJACENT_DOMAIN%LOCAL_GHOST_SEND_INDICES)) DEALLOCATE(ADJACENT_DOMAIN%LOCAL_GHOST_SEND_INDICES)
    IF(ALLOCATED(ADJACENT_DOMAIN%LOCAL_GHOST_RECEIVE_INDICES)) DEALLOCATE(ADJACENT_DOMAIN%LOCAL_GHOST_RECEIVE_INDICES)
    ADJACENT_DOMAIN%NUMBER_OF_SEND_GHOSTS=0
    ADJACENT_DOMAIN%NUMBER_OF_RECEIVE_GHOSTS=0
    ADJACENT_DOMAIN%DOMAIN_NUMBER=0

    EXITS("DOMAIN_MAPPINGS_ADJACENT_DOMAIN_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_MAPPINGS_ADJACENT_DOMAIN_FINALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_MAPPINGS_ADJACENT_DOMAIN_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialise the adjacent domain for a domain mapping.
  SUBROUTINE DOMAIN_MAPPINGS_ADJACENT_DOMAIN_INITIALISE(ADJACENT_DOMAIN,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_ADJACENT_DOMAIN_TYPE) :: ADJACENT_DOMAIN !<The adjacent domain to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_MAPPINGS_ADJACENT_DOMAIN_INITIALISE",ERR,ERROR,*999)

    ADJACENT_DOMAIN%NUMBER_OF_SEND_GHOSTS=0
    ADJACENT_DOMAIN%NUMBER_OF_RECEIVE_GHOSTS=0
    ADJACENT_DOMAIN%DOMAIN_NUMBER=0

    EXITS("DOMAIN_MAPPINGS_ADJACENT_DOMAIN_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_MAPPINGS_ADJACENT_DOMAIN_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_MAPPINGS_ADJACENT_DOMAIN_INITIALISE

  !
  !================================================================================================================================
  !

  !>Returns the local number, if it exists on the rank, for the specifed global number
  SUBROUTINE DOMAIN_MAPPINGS_GLOBAL_TO_LOCAL_GET(DOMAIN_MAPPING,GLOBAL_NUMBER,LOCAL_EXISTS,LOCAL_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOMAIN_MAPPING !<A pointer to the domain mapping to get the local number from
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER !<The global number to get the local number for
    LOGICAL, INTENT(OUT) :: LOCAL_EXISTS !<On exit, is .TRUE. if the specifed global number exists on the local rank, .FALSE. if not
    INTEGER(INTG), INTENT(OUT) :: LOCAL_NUMBER !<On exit, the local number corresponding to the global number if it exists. If it doesn't exist then 0.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("DOMAIN_MAPPINGS_GLOBAL_TO_LOCAL_GET",ERR,ERROR,*999)

    LOCAL_EXISTS=.FALSE.
    LOCAL_NUMBER=0
    IF(ASSOCIATED(DOMAIN_MAPPING)) THEN
      IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=DOMAIN_MAPPING%NUMBER_OF_GLOBAL) THEN
        IF(DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(GLOBAL_NUMBER)%DOMAIN_NUMBER(1)== &
          & COMPUTATIONAL_ENVIRONMENT%MY_COMPUTATIONAL_NODE_NUMBER) THEN
          LOCAL_NUMBER=DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(GLOBAL_NUMBER)%LOCAL_NUMBER(1)
          LOCAL_EXISTS=.TRUE.
        ENDIF
      ELSE
        LOCAL_ERROR="The specified global number of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
          & " is invalid. The number must be between 1 and "// &
          & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING%NUMBER_OF_GLOBAL,"*",ERR,ERROR))//"."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Domain mapping is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("DOMAIN_MAPPINGS_GLOBAL_TO_LOCAL_GET")
    RETURN
999 ERRORSEXITS("DOMAIN_MAPPINGS_GLOBAL_TO_LOCAL_GET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_MAPPINGS_GLOBAL_TO_LOCAL_GET


  !
  !================================================================================================================================
  !

  !>Calculates the domain mappings local map from a domain mappings global map.
  SUBROUTINE DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE( DOMAIN_MAPPING, ERR, ERROR, * )

    !Argument variables
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOMAIN_MAPPING !<The domain mapping to calculate the local mappings
    INTEGER(INTG),         INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING),  INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    INTEGER(INTG) :: domain_idx,domain_idx2,domain_no,domain_no2,global_number,idx,local_number,local_number2,NUMBER_INTERNAL, &
                   & NUMBER_BOUNDARY,NUMBER_GHOST,my_computational_node_number,MY_DOMAIN_INDEX,TEMP,NUMBER_OF_ADJACENT_DOMAINS, &
                   & RECEIVE_FROM_DOMAIN,DUMMY_ERR,NUMBER_OF_GHOST_RECEIVE,NUMBER_OF_GHOST_SEND,local_type,COUNT, &
                   & TOTAL_NUMBER_OF_ADJACENT_DOMAINS, ierr, errorcode
    INTEGER(INTG),       ALLOCATABLE :: ADJACENT_DOMAIN_MAP(:),ADJACENT_DOMAINS(:,:),SEND_LIST(:),RECEIVE_LIST(:)
    TYPE(LIST_PTR_TYPE), ALLOCATABLE :: GHOST_SEND_LISTS(:),GHOST_RECEIVE_LISTS(:)

    LOGICAL :: OWNED_BY_ALL,SEND_GLOBAL
    TYPE(VARYING_STRING) :: LOCAL_ERROR,DUMMY_ERROR

    ENTERS("DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE",ERR,ERROR,*999)

    if (.not.ASSOCIATED(DOMAIN_MAPPING))  call FlagError( "Domain mapping is not associated.", ERR, ERROR, *999 )

    my_computational_node_number=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)
    IF(ERR/=0) GOTO 999

    DO global_number=1,DOMAIN_MAPPING%NUMBER_OF_GLOBAL
      !If necessary, reset global domain index so that my computational node is in the first index position
       MY_DOMAIN_INDEX=1
       DO domain_idx=2,DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%NUMBER_OF_DOMAINS
          domain_no = DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%DOMAIN_NUMBER(domain_idx)
          IF (domain_no==my_computational_node_number) THEN
             MY_DOMAIN_INDEX = domain_idx
             EXIT
          ENDIF
       ENDDO !domain_idx
       IF (MY_DOMAIN_INDEX/=1) THEN !Swap domain index in the global to local map
          TEMP=DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%LOCAL_NUMBER(1)
          DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%LOCAL_NUMBER(1) = &
            & DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%LOCAL_NUMBER(MY_DOMAIN_INDEX)
          DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%LOCAL_NUMBER(MY_DOMAIN_INDEX) = TEMP
          TEMP=DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%DOMAIN_NUMBER(1)
          DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%DOMAIN_NUMBER(1) = &
            & DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%DOMAIN_NUMBER(MY_DOMAIN_INDEX)
          DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%DOMAIN_NUMBER(MY_DOMAIN_INDEX) = TEMP
          TEMP=DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%LOCAL_TYPE(1)
          DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%LOCAL_TYPE(1) = &
            & DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%LOCAL_TYPE(MY_DOMAIN_INDEX)
          DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%LOCAL_TYPE(MY_DOMAIN_INDEX) = TEMP
       ENDIF
    ENDDO !global_number

    EXITS("DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE")
    RETURN
999 RETURN 1

  END SUBROUTINE DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE

  SUBROUTINE DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE2( DOMAIN_MAPPING, ERR, ERROR, * )

    !Argument variables
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOMAIN_MAPPING !<The domain mapping to calculate the local mappings
    INTEGER(INTG),         INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING),  INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    INTEGER(INTG) :: domain_idx,domain_idx2,domain_no,domain_no2,global_number,idx,local_number,local_number2,NUMBER_INTERNAL, &
                   & NUMBER_BOUNDARY,NUMBER_GHOST,my_computational_node_number,MY_DOMAIN_INDEX,TEMP,NUMBER_OF_ADJACENT_DOMAINS, &
                   & RECEIVE_FROM_DOMAIN,DUMMY_ERR,NUMBER_OF_GHOST_RECEIVE,NUMBER_OF_GHOST_SEND,local_type,COUNT, &
                   & TOTAL_NUMBER_OF_ADJACENT_DOMAINS, ierr, errorcode
    INTEGER(INTG),       ALLOCATABLE :: ADJACENT_DOMAIN_MAP(:),ADJACENT_DOMAINS(:,:),SEND_LIST(:),RECEIVE_LIST(:)
    TYPE(LIST_PTR_TYPE), ALLOCATABLE :: GHOST_SEND_LISTS(:),GHOST_RECEIVE_LISTS(:)

    LOGICAL :: OWNED_BY_ALL,SEND_GLOBAL
    TYPE(VARYING_STRING) :: LOCAL_ERROR,DUMMY_ERROR

    ENTERS("DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE2",ERR,ERROR,*999)

    if (.not.ASSOCIATED(DOMAIN_MAPPING))  call FlagError( "Domain mapping is not associated.", ERR, ERROR, *999 )

    my_computational_node_number=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)
    IF(ERR/=0) GOTO 999

    !Calculate local to global maps from global to local map
    if (.not.allocated(DOMAIN_MAPPING%NUMBER_OF_DOMAIN_LOCAL)) &
       & ALLOCATE(DOMAIN_MAPPING%NUMBER_OF_DOMAIN_LOCAL(0:DOMAIN_MAPPING%NUMBER_OF_DOMAINS-1),STAT=ERR)
    IF(ERR/=0) CALL FlagError("Could not allocate number of domain local.",ERR,ERROR,*999)
    DOMAIN_MAPPING%NUMBER_OF_DOMAIN_LOCAL=0
    if (.not.allocated(DOMAIN_MAPPING%NUMBER_OF_DOMAIN_GHOST)) &
       & ALLOCATE(DOMAIN_MAPPING%NUMBER_OF_DOMAIN_GHOST(0:DOMAIN_MAPPING%NUMBER_OF_DOMAINS-1),STAT=ERR)
    IF(ERR/=0) CALL FlagError("Could not allocate number of domain ghost.",ERR,ERROR,*999)
    DOMAIN_MAPPING%NUMBER_OF_DOMAIN_GHOST=0
    NUMBER_INTERNAL=0
    NUMBER_BOUNDARY=0
    NUMBER_GHOST=0
    ALLOCATE(ADJACENT_DOMAINS(0:DOMAIN_MAPPING%NUMBER_OF_DOMAINS-1,0:DOMAIN_MAPPING%NUMBER_OF_DOMAINS-1),STAT=ERR)
    IF(ERR/=0) CALL FlagError("Could not allocate adjacent domains.",ERR,ERROR,*999)
    ADJACENT_DOMAINS=0
    DO global_number=1,DOMAIN_MAPPING%NUMBER_OF_GLOBAL
      !If necessary, reset global domain index so that my computational node is in the first index position
       MY_DOMAIN_INDEX=1
       DO domain_idx=2,DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%NUMBER_OF_DOMAINS
          domain_no = DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%DOMAIN_NUMBER(domain_idx)
          IF (domain_no==my_computational_node_number) THEN
             MY_DOMAIN_INDEX = domain_idx
             EXIT
          ENDIF
       ENDDO !domain_idx
       IF (MY_DOMAIN_INDEX/=1) THEN !Swap domain index in the global to local map
          TEMP=DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%LOCAL_NUMBER(1)
          DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%LOCAL_NUMBER(1) = &
            & DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%LOCAL_NUMBER(MY_DOMAIN_INDEX)
          DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%LOCAL_NUMBER(MY_DOMAIN_INDEX) = TEMP
          TEMP=DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%DOMAIN_NUMBER(1)
          DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%DOMAIN_NUMBER(1) = &
            & DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%DOMAIN_NUMBER(MY_DOMAIN_INDEX)
          DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%DOMAIN_NUMBER(MY_DOMAIN_INDEX) = TEMP
          TEMP=DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%LOCAL_TYPE(1)
          DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%LOCAL_TYPE(1) = &
            & DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%LOCAL_TYPE(MY_DOMAIN_INDEX)
          DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%LOCAL_TYPE(MY_DOMAIN_INDEX) = TEMP
       ENDIF
       DO domain_idx = 1,DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%NUMBER_OF_DOMAINS
          domain_no=DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%DOMAIN_NUMBER(domain_idx)
          DO domain_idx2=1,DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%NUMBER_OF_DOMAINS
             domain_no2=DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%DOMAIN_NUMBER(domain_idx2)
             ADJACENT_DOMAINS(domain_no,domain_no2)=1
          ENDDO !domain_idx2
        ENDDO !domain_idx
        DO domain_idx=1,DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%NUMBER_OF_DOMAINS
          domain_no=DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%DOMAIN_NUMBER(domain_idx)
          local_type=DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%LOCAL_TYPE(domain_idx)
          IF(local_type==DOMAIN_LOCAL_GHOST) THEN
            DOMAIN_MAPPING%NUMBER_OF_DOMAIN_GHOST(domain_no)=DOMAIN_MAPPING%NUMBER_OF_DOMAIN_GHOST(domain_no)+1
          ELSE
            DOMAIN_MAPPING%NUMBER_OF_DOMAIN_LOCAL(domain_no)=DOMAIN_MAPPING%NUMBER_OF_DOMAIN_LOCAL(domain_no)+1
          ENDIF
          IF(domain_no==my_computational_node_number) THEN
            SELECT CASE(local_type)
            CASE(DOMAIN_LOCAL_INTERNAL)
              NUMBER_INTERNAL=NUMBER_INTERNAL+1
            CASE(DOMAIN_LOCAL_BOUNDARY)
              NUMBER_BOUNDARY=NUMBER_BOUNDARY+1
            CASE(DOMAIN_LOCAL_GHOST)
              NUMBER_GHOST=NUMBER_GHOST+1
            CASE DEFAULT
              LOCAL_ERROR="The domain local type of "//TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP( &
                & global_number)%LOCAL_TYPE(domain_idx),"*",ERR,ERROR))//" is invalid."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          ENDIF
        ENDDO !domain_idx
      ENDDO !global_number

      !!TODO: move adjacent domains calculation back to where the global to local array is set up????
      NUMBER_OF_ADJACENT_DOMAINS=0
      TOTAL_NUMBER_OF_ADJACENT_DOMAINS=0
      DO domain_no=0,DOMAIN_MAPPING%NUMBER_OF_DOMAINS-1
        DO domain_no2=0,DOMAIN_MAPPING%NUMBER_OF_DOMAINS-1
          IF(domain_no/=domain_no2) THEN
            IF(ADJACENT_DOMAINS(domain_no,domain_no2)>0) THEN
              TOTAL_NUMBER_OF_ADJACENT_DOMAINS=TOTAL_NUMBER_OF_ADJACENT_DOMAINS+1
              IF(domain_no==my_computational_node_number) NUMBER_OF_ADJACENT_DOMAINS=NUMBER_OF_ADJACENT_DOMAINS+1
            ENDIF
          ENDIF
        ENDDO !domain_no2
      ENDDO !domain_no
      if ( .not.allocated(DOMAIN_MAPPING%ADJACENT_DOMAINS_PTR) ) &
        & ALLOCATE(DOMAIN_MAPPING%ADJACENT_DOMAINS_PTR(0:DOMAIN_MAPPING%NUMBER_OF_DOMAINS),STAT=ERR)
      IF(ERR/=0) CALL FlagError("Could not allocate adjacent domains ptr.",ERR,ERROR,*999)
      if ( .not.allocated(DOMAIN_MAPPING%ADJACENT_DOMAINS_LIST) ) &
        &  ALLOCATE(DOMAIN_MAPPING%ADJACENT_DOMAINS_LIST(TOTAL_NUMBER_OF_ADJACENT_DOMAINS),STAT=ERR)
      IF(ERR/=0) CALL FlagError("Could not allocate adjacent domains list.",ERR,ERROR,*999)
      COUNT=1
      DO domain_no=0,DOMAIN_MAPPING%NUMBER_OF_DOMAINS-1
        DOMAIN_MAPPING%ADJACENT_DOMAINS_PTR(domain_no)=COUNT
        DO domain_no2=0,DOMAIN_MAPPING%NUMBER_OF_DOMAINS-1
          IF(domain_no/=domain_no2) THEN
            IF(ADJACENT_DOMAINS(domain_no,domain_no2)>0) THEN
              DOMAIN_MAPPING%ADJACENT_DOMAINS_LIST(COUNT)=domain_no2
              COUNT=COUNT+1
            ENDIF
          ENDIF
        ENDDO !domain_no2
      ENDDO !domain_no
      DOMAIN_MAPPING%ADJACENT_DOMAINS_PTR(DOMAIN_MAPPING%NUMBER_OF_DOMAINS)=COUNT
      DEALLOCATE(ADJACENT_DOMAINS)

      if ( .not.allocated(DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP) ) then
         ALLOCATE(DOMAIN_MAPPING%DOMAIN_LIST(NUMBER_INTERNAL+NUMBER_BOUNDARY+NUMBER_GHOST),STAT=ERR)
         IF(ERR/=0) CALL FlagError("Could not allocate domain map domain list.",ERR,ERROR,*999)
         ALLOCATE(DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(NUMBER_INTERNAL+NUMBER_BOUNDARY+NUMBER_GHOST),STAT=ERR)
         IF(ERR/=0) CALL FlagError("Could not allocate domain map local to global list.",ERR,ERROR,*999)
         DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL=NUMBER_INTERNAL+NUMBER_BOUNDARY+NUMBER_GHOST
         DOMAIN_MAPPING%NUMBER_OF_LOCAL=NUMBER_INTERNAL+NUMBER_BOUNDARY
         DOMAIN_MAPPING%NUMBER_OF_INTERNAL=NUMBER_INTERNAL
         DOMAIN_MAPPING%NUMBER_OF_BOUNDARY=NUMBER_BOUNDARY
         DOMAIN_MAPPING%NUMBER_OF_GHOST=NUMBER_GHOST
         DOMAIN_MAPPING%INTERNAL_START=1
         DOMAIN_MAPPING%INTERNAL_FINISH=NUMBER_INTERNAL
         DOMAIN_MAPPING%BOUNDARY_START=NUMBER_INTERNAL+1
         DOMAIN_MAPPING%BOUNDARY_FINISH=NUMBER_INTERNAL+NUMBER_BOUNDARY
         DOMAIN_MAPPING%GHOST_START=NUMBER_INTERNAL+NUMBER_BOUNDARY+1
         DOMAIN_MAPPING%GHOST_FINISH=NUMBER_INTERNAL+NUMBER_BOUNDARY+NUMBER_GHOST
      endif

      NUMBER_INTERNAL=0
      NUMBER_BOUNDARY=0
      NUMBER_GHOST=0
      if ( .not.allocated(DOMAIN_MAPPING%ADJACENT_DOMAINS) ) then
      ALLOCATE(DOMAIN_MAPPING%ADJACENT_DOMAINS(NUMBER_OF_ADJACENT_DOMAINS),STAT=ERR)
      IF(ERR/=0) CALL FlagError("Could not allocate adjacent domains.",ERR,ERROR,*999)
      DOMAIN_MAPPING%NUMBER_OF_ADJACENT_DOMAINS=NUMBER_OF_ADJACENT_DOMAINS
      ALLOCATE(ADJACENT_DOMAIN_MAP(0:DOMAIN_MAPPING%NUMBER_OF_DOMAINS-1),STAT=ERR)
      IF(ERR/=0) CALL FlagError("Could not allocate adjacent domain map.",ERR,ERROR,*999)
      ALLOCATE(GHOST_SEND_LISTS(DOMAIN_MAPPING%NUMBER_OF_ADJACENT_DOMAINS),STAT=ERR)
      IF(ERR/=0) CALL FlagError("Could not allocate ghost send list.",ERR,ERROR,*999)
      ALLOCATE(GHOST_RECEIVE_LISTS(DOMAIN_MAPPING%NUMBER_OF_ADJACENT_DOMAINS),STAT=ERR)
      IF(ERR/=0) CALL FlagError("Could not allocate ghost recieve list.",ERR,ERROR,*999)
      endif
      DO domain_idx=1,DOMAIN_MAPPING%NUMBER_OF_ADJACENT_DOMAINS
        CALL DOMAIN_MAPPINGS_ADJACENT_DOMAIN_INITIALISE(DOMAIN_MAPPING%ADJACENT_DOMAINS(domain_idx),ERR,ERROR,*999)
        domain_no= &
          & DOMAIN_MAPPING%ADJACENT_DOMAINS_LIST(DOMAIN_MAPPING%ADJACENT_DOMAINS_PTR(my_computational_node_number)+domain_idx-1)
        DOMAIN_MAPPING%ADJACENT_DOMAINS(domain_idx)%DOMAIN_NUMBER=domain_no
        ADJACENT_DOMAIN_MAP(domain_no)=domain_idx
        NULLIFY(GHOST_SEND_LISTS(domain_idx)%PTR)
        CALL LIST_CREATE_START(GHOST_SEND_LISTS(domain_idx)%PTR,ERR,ERROR,*999)
        CALL LIST_DATA_TYPE_SET(GHOST_SEND_LISTS(domain_idx)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
        CALL LIST_INITIAL_SIZE_SET(GHOST_SEND_LISTS(domain_idx)%PTR,MAX(DOMAIN_MAPPING%NUMBER_OF_GHOST,1),ERR,ERROR,*999)
        CALL LIST_CREATE_FINISH(GHOST_SEND_LISTS(domain_idx)%PTR,ERR,ERROR,*999)
        NULLIFY(GHOST_RECEIVE_LISTS(domain_idx)%PTR)
        CALL LIST_CREATE_START(GHOST_RECEIVE_LISTS(domain_idx)%PTR,ERR,ERROR,*999)
        CALL LIST_DATA_TYPE_SET(GHOST_RECEIVE_LISTS(domain_idx)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
        CALL LIST_INITIAL_SIZE_SET(GHOST_RECEIVE_LISTS(domain_idx)%PTR,MAX(DOMAIN_MAPPING%NUMBER_OF_GHOST,1),ERR,ERROR,*999)
        CALL LIST_CREATE_FINISH(GHOST_RECEIVE_LISTS(domain_idx)%PTR,ERR,ERROR,*999)
      ENDDO !domain_idx

      DO global_number=1,DOMAIN_MAPPING%NUMBER_OF_GLOBAL
        SEND_GLOBAL=.FALSE.
        IF(DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%NUMBER_OF_DOMAINS>1) THEN
          !Check if we have a special case where the global number is owned by all domains e.g., as in a constant field
          IF(DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%NUMBER_OF_DOMAINS==DOMAIN_MAPPING%NUMBER_OF_DOMAINS) THEN
            OWNED_BY_ALL=.TRUE.
            DO domain_idx=1,DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%NUMBER_OF_DOMAINS
              local_type=DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%LOCAL_TYPE(domain_idx)
              OWNED_BY_ALL=OWNED_BY_ALL.AND.local_type==DOMAIN_LOCAL_INTERNAL
            ENDDO !domain_idx
          ELSE
            OWNED_BY_ALL=.FALSE.
          ENDIF
          IF(.NOT.OWNED_BY_ALL) THEN
            RECEIVE_FROM_DOMAIN=-1
            DO domain_idx=1,DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%NUMBER_OF_DOMAINS
              domain_no=DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%DOMAIN_NUMBER(domain_idx)
              local_type=DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%LOCAL_TYPE(domain_idx)
              IF(local_type/=DOMAIN_LOCAL_GHOST) THEN
                IF(domain_no==my_computational_node_number) SEND_GLOBAL=.TRUE.
                IF(RECEIVE_FROM_DOMAIN==-1) THEN
                  RECEIVE_FROM_DOMAIN=domain_no
                ELSE
                  LOCAL_ERROR="Invalid domain mapping. Global number "//TRIM(NUMBER_TO_VSTRING(global_number,"*",ERR,ERROR))// &
                    & " is owned by domain number "//TRIM(NUMBER_TO_VSTRING(RECEIVE_FROM_DOMAIN,"*",ERR,ERROR))// &
                    & " as well as domain number "//TRIM(NUMBER_TO_VSTRING(domain_no,"*",ERR,ERROR))//"."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ENDIF
            ENDDO !domain_idx
            IF(RECEIVE_FROM_DOMAIN==-1) THEN
              LOCAL_ERROR="Invalid domain mapping. Global number "//TRIM(NUMBER_TO_VSTRING(global_number,"*",ERR,ERROR))// &
                & " is not owned by any domain."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ENDIF
        ENDIF
        DO domain_idx=1,DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%NUMBER_OF_DOMAINS
          domain_no=DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%DOMAIN_NUMBER(domain_idx)
          local_number=DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%LOCAL_NUMBER(domain_idx)
          local_type=DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%LOCAL_TYPE(domain_idx)
          IF(domain_no==my_computational_node_number) THEN
            DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_number)=global_number
            SELECT CASE(DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%LOCAL_TYPE(domain_idx))
            CASE(DOMAIN_LOCAL_INTERNAL)
              NUMBER_INTERNAL=NUMBER_INTERNAL+1
              DOMAIN_MAPPING%DOMAIN_LIST(NUMBER_INTERNAL)=local_number
            CASE(DOMAIN_LOCAL_BOUNDARY)
              NUMBER_BOUNDARY=NUMBER_BOUNDARY+1
              DOMAIN_MAPPING%DOMAIN_LIST(DOMAIN_MAPPING%INTERNAL_FINISH+NUMBER_BOUNDARY)=local_number
            CASE(DOMAIN_LOCAL_GHOST)
              NUMBER_GHOST=NUMBER_GHOST+1
              DOMAIN_MAPPING%DOMAIN_LIST(DOMAIN_MAPPING%BOUNDARY_FINISH+NUMBER_GHOST)=local_number
              CALL LIST_ITEM_ADD(GHOST_RECEIVE_LISTS(ADJACENT_DOMAIN_MAP(RECEIVE_FROM_DOMAIN))%PTR,local_number,ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="The domain local type of "//TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP( &
                & global_number)%LOCAL_TYPE(domain_idx),"*",ERR,ERROR))//" is invalid."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          ELSE IF(SEND_GLOBAL.AND.local_type==DOMAIN_LOCAL_GHOST) THEN
            local_number2=DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%LOCAL_NUMBER(1) !The local number for this node
            CALL LIST_ITEM_ADD(GHOST_SEND_LISTS(ADJACENT_DOMAIN_MAP(domain_no))%PTR,local_number2,ERR,ERROR,*999)
          ENDIF
        ENDDO !domain_idx
      ENDDO !global_number

      DO domain_idx=1,DOMAIN_MAPPING%NUMBER_OF_ADJACENT_DOMAINS
        CALL LIST_REMOVE_DUPLICATES(GHOST_SEND_LISTS(domain_idx)%PTR,ERR,ERROR,*999)
        CALL LIST_DETACH_AND_DESTROY(GHOST_SEND_LISTS(domain_idx)%PTR,NUMBER_OF_GHOST_SEND,SEND_LIST,ERR,ERROR,*999)
        ALLOCATE(DOMAIN_MAPPING%ADJACENT_DOMAINS(domain_idx)%LOCAL_GHOST_SEND_INDICES(NUMBER_OF_GHOST_SEND),STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate local ghost send inidices.",ERR,ERROR,*999)
        DOMAIN_MAPPING%ADJACENT_DOMAINS(domain_idx)%LOCAL_GHOST_SEND_INDICES(1:NUMBER_OF_GHOST_SEND)= &
          & SEND_LIST(1:NUMBER_OF_GHOST_SEND)
        DOMAIN_MAPPING%ADJACENT_DOMAINS(domain_idx)%NUMBER_OF_SEND_GHOSTS=NUMBER_OF_GHOST_SEND
        DEALLOCATE(SEND_LIST)
        CALL LIST_REMOVE_DUPLICATES(GHOST_RECEIVE_LISTS(domain_idx)%PTR,ERR,ERROR,*999)
        CALL LIST_DETACH_AND_DESTROY(GHOST_RECEIVE_LISTS(domain_idx)%PTR,NUMBER_OF_GHOST_RECEIVE,RECEIVE_LIST,ERR,ERROR,*999)
        ALLOCATE(DOMAIN_MAPPING%ADJACENT_DOMAINS(domain_idx)%LOCAL_GHOST_RECEIVE_INDICES(NUMBER_OF_GHOST_RECEIVE),STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate local ghost receive inidices.",ERR,ERROR,*999)
        DOMAIN_MAPPING%ADJACENT_DOMAINS(domain_idx)%LOCAL_GHOST_RECEIVE_INDICES(1:NUMBER_OF_GHOST_RECEIVE)= &
          & RECEIVE_LIST(1:NUMBER_OF_GHOST_RECEIVE)
        DOMAIN_MAPPING%ADJACENT_DOMAINS(domain_idx)%NUMBER_OF_RECEIVE_GHOSTS=NUMBER_OF_GHOST_RECEIVE
        DEALLOCATE(RECEIVE_LIST)
      ENDDO !domain_idx

      DEALLOCATE(ADJACENT_DOMAIN_MAP)
      DEALLOCATE(GHOST_SEND_LISTS)
      DEALLOCATE(GHOST_RECEIVE_LISTS)

    EXITS("DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE2")
    RETURN
999 IF(ALLOCATED(SEND_LIST)) DEALLOCATE(SEND_LIST)
    IF(ALLOCATED(RECEIVE_LIST)) DEALLOCATE(RECEIVE_LIST)
    IF(ALLOCATED(ADJACENT_DOMAIN_MAP)) DEALLOCATE(ADJACENT_DOMAIN_MAP)
    IF(ALLOCATED(ADJACENT_DOMAINS)) DEALLOCATE(ADJACENT_DOMAINS)
    IF(ALLOCATED(GHOST_SEND_LISTS)) THEN
      DO domain_idx=1,SIZE(GHOST_SEND_LISTS)
        IF(ASSOCIATED(GHOST_SEND_LISTS(domain_idx)%PTR)) &
          & CALL LIST_DESTROY(GHOST_SEND_LISTS(domain_idx)%PTR,DUMMY_ERR,DUMMY_ERROR,*998)
      ENDDO ! domain_idx
998   DEALLOCATE(GHOST_SEND_LISTS)
    ENDIF
    IF(ALLOCATED(GHOST_RECEIVE_LISTS)) THEN
      DO domain_idx=1,SIZE(GHOST_RECEIVE_LISTS)
        IF(ASSOCIATED(GHOST_RECEIVE_LISTS(domain_idx)%PTR)) &
          & CALL LIST_DESTROY(GHOST_RECEIVE_LISTS(domain_idx)%PTR,DUMMY_ERR,DUMMY_ERROR,*997)
      ENDDO ! domain_idx
997   DEALLOCATE(GHOST_RECEIVE_LISTS)
    ENDIF
    ERRORSEXITS("DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE2",ERR,ERROR)
    RETURN 1

  END SUBROUTINE DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE2


  !
  !================================================================================================================================
  !

  !>Finalises the mapping for a domain mappings mapping and deallocates all memory.
  SUBROUTINE DOMAIN_MAPPINGS_MAPPING_FINALISE(DOMAIN_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOMAIN_MAPPING !<A pointer to the domain mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: idx

    ENTERS("DOMAIN_MAPPINGS_MAPPING_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN_MAPPING)) THEN
      IF(ALLOCATED(DOMAIN_MAPPING%NUMBER_OF_DOMAIN_LOCAL)) DEALLOCATE(DOMAIN_MAPPING%NUMBER_OF_DOMAIN_LOCAL)
      IF(ALLOCATED(DOMAIN_MAPPING%NUMBER_OF_DOMAIN_GHOST)) DEALLOCATE(DOMAIN_MAPPING%NUMBER_OF_DOMAIN_GHOST)
      IF(ALLOCATED(DOMAIN_MAPPING%DOMAIN_LIST)) DEALLOCATE(DOMAIN_MAPPING%DOMAIN_LIST)
      IF(ALLOCATED(DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP)) DEALLOCATE(DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP)
      IF(ALLOCATED(DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP)) THEN
        DO idx=1,SIZE(DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP,1)
          CALL DOMAIN_MAPPINGS_MAPPING_GLOBAL_FINALISE(DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(idx),ERR,ERROR,*999)
        ENDDO !idx
        DEALLOCATE(DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP)
      ENDIF
      IF(ALLOCATED(DOMAIN_MAPPING%ADJACENT_DOMAINS_PTR)) DEALLOCATE(DOMAIN_MAPPING%ADJACENT_DOMAINS_PTR)
      IF(ALLOCATED(DOMAIN_MAPPING%ADJACENT_DOMAINS_LIST)) DEALLOCATE(DOMAIN_MAPPING%ADJACENT_DOMAINS_LIST)
      IF(ALLOCATED(DOMAIN_MAPPING%ADJACENT_DOMAINS)) THEN
        DO idx=1,SIZE(DOMAIN_MAPPING%ADJACENT_DOMAINS,1)
          CALL DOMAIN_MAPPINGS_ADJACENT_DOMAIN_FINALISE(DOMAIN_MAPPING%ADJACENT_DOMAINS(idx),ERR,ERROR,*999)
        ENDDO !idx
        DEALLOCATE(DOMAIN_MAPPING%ADJACENT_DOMAINS)
      ENDIF
      DEALLOCATE(DOMAIN_MAPPING)
    ENDIF

    EXITS("DOMAIN_MAPPINGS_MAPPING_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_MAPPINGS_MAPPING_FINALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_MAPPINGS_MAPPING_FINALISE

  !
  !================================================================================================================================
  !

  !> Finalises the global mapping in the given domain mappings.
  SUBROUTINE DOMAIN_MAPPINGS_MAPPING_GLOBAL_FINALISE(MAPPING_GLOBAL_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_GLOBAL_MAPPING_TYPE) :: MAPPING_GLOBAL_MAP !<The domain global mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<Th error string
    !Local Variables

    ENTERS("DOMAIN_MAPPINGS_MAPPING_GLOBAL_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(MAPPING_GLOBAL_MAP%LOCAL_NUMBER)) DEALLOCATE(MAPPING_GLOBAL_MAP%LOCAL_NUMBER)
    IF(ALLOCATED(MAPPING_GLOBAL_MAP%DOMAIN_NUMBER)) DEALLOCATE(MAPPING_GLOBAL_MAP%DOMAIN_NUMBER)
    IF(ALLOCATED(MAPPING_GLOBAL_MAP%LOCAL_TYPE)) DEALLOCATE(MAPPING_GLOBAL_MAP%LOCAL_TYPE)

    EXITS("DOMAIN_MAPPINGS_MAPPING_GLOBAL_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_MAPPINGS_MAPPING_GLOBAL_FINALISE",ERR,ERROR)
    RETURN 1

  END SUBROUTINE DOMAIN_MAPPINGS_MAPPING_GLOBAL_FINALISE

  !
  !================================================================================================================================
  !

  !>Finalises the global mapping in the given domain mappings.
  SUBROUTINE DOMAIN_MAPPINGS_MAPPING_GLOBAL_INITIALISE(MAPPING_GLOBAL_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_GLOBAL_MAPPING_TYPE) :: MAPPING_GLOBAL_MAP !<The domain global mapping to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_MAPPINGS_MAPPING_GLOBAL_INITIALISE",ERR,ERROR,*999)

    MAPPING_GLOBAL_MAP%NUMBER_OF_DOMAINS=0

    EXITS("DOMAIN_MAPPINGS_MAPPING_GLOBAL_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_MAPPINGS_MAPPING_GLOBAL_INITIALISE",ERR,ERROR)
    RETURN 1

  END SUBROUTINE DOMAIN_MAPPINGS_MAPPING_GLOBAL_INITIALISE

  !
  !================================================================================================================================
  !

  !>Initialises the mapping for a domain mappings mapping.
  SUBROUTINE DOMAIN_MAPPINGS_MAPPING_INITIALISE(DOMAIN_MAPPING,NUMBER_OF_DOMAINS,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOMAIN_MAPPING !<A pointer to the domain mapping to initialise the mappings for
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_DOMAINS !<The number of domains
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("DOMAIN_MAPPINGS_MAPPING_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN_MAPPING)) THEN
      IF(NUMBER_OF_DOMAINS>0) THEN
        DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL=0
        DOMAIN_MAPPING%NUMBER_OF_LOCAL=0
        DOMAIN_MAPPING%NUMBER_OF_GLOBAL=0
        DOMAIN_MAPPING%NUMBER_OF_DOMAINS=NUMBER_OF_DOMAINS
        DOMAIN_MAPPING%NUMBER_OF_INTERNAL=0
        DOMAIN_MAPPING%NUMBER_OF_BOUNDARY=0
        DOMAIN_MAPPING%NUMBER_OF_GHOST=0
        DOMAIN_MAPPING%INTERNAL_START=0
        DOMAIN_MAPPING%INTERNAL_FINISH=0
        DOMAIN_MAPPING%BOUNDARY_START=0
        DOMAIN_MAPPING%BOUNDARY_FINISH=0
        DOMAIN_MAPPING%GHOST_START=0
        DOMAIN_MAPPING%GHOST_FINISH=0
        DOMAIN_MAPPING%NUMBER_OF_ADJACENT_DOMAINS=0
      ELSE
        LOCAL_ERROR="The specified number of domains of "//TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DOMAINS,"*",ERR,ERROR))// &
          & " is invalid. The number of domains must be > 0."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ENDIF

    EXITS("DOMAIN_MAPPINGS_MAPPING_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_MAPPINGS_MAPPING_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_MAPPINGS_MAPPING_INITIALISE

  !
  !================================================================================================================================
  !

END MODULE DOMAIN_MAPPINGS
