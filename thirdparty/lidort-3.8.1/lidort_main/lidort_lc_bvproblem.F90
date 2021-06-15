! ###########################################################
! #                                                         #
! #                  LIDORT_3p8p1                           #
! #                                                         #
! #    (LInearized Discrete Ordinate Radiative Transfer)    #
! #     --         -        -        -         -            #
! #                                                         #
! ###########################################################

! ###########################################################
! #                                                         #
! #  Authors :     Robert  J. D. Spurr (1)                  #
! #                Matthew J. Christi                       #
! #                                                         #
! #  Address (1) : RT Solutions, Inc.                       #
! #                9 Channing Street                        #
! #                Cambridge, MA 02138, USA                 #
! #                                                         #
! #  Tel:          (617) 492 1183                           #
! #  Email :       rtsolutions@verizon.net                  #
! #                                                         #
! #  This Version :   LIDORT_3p8p1                          #
! #  Release Date :   30 June 2019                          #
! #                                                         #
! #  Previous LIDORT Versions under Standard GPL 3.0:       #
! #  ------------------------------------------------       #
! #                                                         #
! #      3.7   F90, released June  2014                     #
! #      3.8   F90, released March 2017                     #
! #                                                         #
! #  Features Summary of Recent LIDORT Versions             #
! #  ------------------------------------------             #
! #                                                         #
! #       NEW: THERMAL SUPPLEMENT INCLUDED    (3.2)         #
! #       NEW: OUTGOING SPHERICITY CORRECTION (3.2)         #
! #       NEW: TOTAL COLUMN JACOBIANS         (3.3)         #
! #       VLIDORT COMPATIBILITY               (3.4)         #
! #       THREADED/OPTIMIZED F90 code         (3.5)         #
! #       EXTERNAL SS / NEW I/O STRUCTURES    (3.6)         #
! #                                                         #
! #       Surface-leaving, BRDF Albedo-scaling     (3.7)    # 
! #       Taylor series, BBF Jacobians, ThreadSafe (3.7)    #
! #       New Water-Leaving Treatment              (3.8)    #
! #       BRDF-Telescoping, enabled                (3.8)    #
! #       Several Performance Enhancements         (3.8)    #
! #       Water-leaving coupled code               (3.8.1)  #
! #       Planetary problem, media properties      (3.8.1)  #
! #                                                         #
! ###########################################################

!    ###################################################################
!    #                                                                 #
!    # This is Version 3.8.1 of the LIDORT software library.           #
!    # This library comes with the Standard GNU General Public License,#
!    # Version 3.0, 29 June 2007. Please read this license carefully.  #
!    #                                                                 #
!    #      LIDORT Copyright (c) 1999-2019.                            #
!    #          Robert Spurr, RT Solutions, Inc.                       #
!    #          9 Channing Street, Cambridge, MA 02138, USA.           #
!    #                                                                 #
!    #                                                                 #
!    # This file is part of LIDORT_3p8p1 ( Version 3.8.1. )            #
!    #                                                                 #
!    # LIDORT_3p8p1 is free software: you can redistribute it          #
!    # and/or modify it under the terms of the Standard GNU GPL        #
!    # (General Public License) as published by the Free Software      #
!    # Foundation, either version 3.0 of this License, or any          #
!    # later version.                                                  #
!    #                                                                 #
!    # LIDORT_3p8p1 is distributed in the hope that it will be         #
!    # useful, but WITHOUT ANY WARRANTY; without even the implied      #
!    # warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR         #
!    # PURPOSE. See the Standard GNU General Public License (GPL)      #
!    # for more details.                                               #
!    #                                                                 #
!    # You should have received a copy of the Standard GNU General     #
!    # Public License (GPL) Version 3.0, along with the LIDORT_3p8p1   #
!    # code package. If not, see <http://www.gnu.org/licenses/>.       #
!    #                                                                 #
!    ###################################################################

! ###############################################################
! #                                                             #
! # Subroutines in this Module                                  #
! #                                                             #
! #  For linearizations (atmospheric) with regular BVP          #
! #            LC_BVP_SOLUTION_MASTER                           #
! #            LC_BVP_COLUMN_SETUP                              #
! #            LC_BVP_SURFACE_SETUP                             #
! #            LC_BEAMSOLUTION_NEQK                             #
! #                                                             #
! #  For linearizations (atmospheric) with telescoped BVP       #
! #            LC_BVPTEL_SOLUTION_MASTER                        #
! #            LC_BVPTEL_COLUMN_SETUP                           #
! #            LC_BVPTEL_SURFACE_SETUP                          #
! #                                                             #
! ###############################################################

!  Upgrade, Version 3.8.  Implementation of BRDFs in the Linearized Telescoped problem
!  Upgrade, Version 3.8.1 Some surface-leaving linearization changes

module lidort_lc_bvproblem_m

!  Parameter types

   USE LIDORT_PARS_m  , only : fpk

!  LIDORT Use dependencies

   USE lidort_aux_m,    only : DGBTRS, DGETRS
   USE lidort_Taylor_m, only : TAYLOR_SERIES_L_1

private :: LC_BVP_COLUMN_SETUP, LC_BVP_SURFACE_SETUP, LC_BEAMSOLUTION_NEQK, &
           LC_BVPTEL_COLUMN_SETUP, LC_BVPTEL_SURFACE_SETUP
public  :: LC_BVP_SOLUTION_MASTER, LC_BVPTEL_SOLUTION_MASTER

contains

SUBROUTINE LC_BVP_SOLUTION_MASTER &
           ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, DO_INCLUDE_DIRECTBEAM,   & ! Input, Flags
             DO_WATER_LEAVING, DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS,    & ! Input, FLags
             DO_LAYER_SCATTERING, TAYLOR_ORDER, FOURIER, IBEAM, NLAYERS,   & ! Input, Flags and order
             NTOTAL, N_SUB, N_SUP, NSTREAMS, NSTREAMS_2, N_WEIGHTFUNCS,    & ! Input, Numbers
             DELTAU_VERT, L_DELTAU_VERT, SOLARBEAM_CUTOFF, QUAD_STRMWTS,   & ! Input, optical and control
             SURFACE_FACTOR, ALBEDO, BRDF_F, RF_DIRECT_BEAM, SLTERM,       & ! Input, Surface Stuff
             DELTAU_SLANT, INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,    & ! Input, Beam Quantities
             BANDMAT2, SMAT2, IPIVOT, SIPIVOT, LCONMASK, MCONMASK,         & ! Input, BVP Bandmat
             T_DELT_EIGEN, XPOS, XNEG, LCON, MCON,                         & ! Input, Homogeneous
             GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE, & ! Input, Greens Function
             LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR,         & ! Input, Linearized Beam Quantities
             LC_TRANS_ATMOS, L_KEIGEN, L_T_DELT_EIGEN, L_XPOS, L_XNEG,     & ! Input, Linearized Homog solution 
             L_ATERM_SAVE, L_BTERM_SAVE, L_T_WUPPER, L_T_WLOWER,           & ! Input, Linearized Greens + Thermal
             NCON, PCON, NCON_XVEC, PCON_XVEC, L_WUPPER, L_WLOWER,         & ! Output - Linearized Constants + PI
             STATUS, MESSAGE, TRACE )                                        ! Output - Exception handling

!  Linearization of the Boundary Problem Solution

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAXSTREAMS, MAXSTREAMS_2, MAXMOMENTS, MAXLAYERS, &
                                MAXBEAMS, MAX_ATMOSWFS, MAXBANDTOTAL, MAXTOTAL,  &
                                LIDORT_SUCCESS, LIDORT_SERIOUS

      IMPLICIT NONE

!  Subroutine arguments
!  ====================

!  Control and Optical
!  -------------------

!  Surface BRDF and inclusion flags

      LOGICAL  , intent(in)  ::  DO_INCLUDE_DIRECTBEAM
      LOGICAL  , intent(in)  ::  DO_INCLUDE_SURFACE
      LOGICAL  , intent(in)  ::  DO_BRDF_SURFACE
      LOGICAL  , intent(in)  ::  DO_WATER_LEAVING

!  Emission and solar source flags

      LOGICAL  , intent(in)  ::  DO_INCLUDE_THERMEMISS
      LOGICAL  , intent(in)  ::  DO_SOLAR_SOURCES

!  Local flags for the solution saving option

      LOGICAL  , intent(in)  ::  DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER  , intent(in)   :: TAYLOR_ORDER

!  Fourier component, beam number

      INTEGER  , intent(in)  ::  FOURIER, IBEAM

!  Number of layers, BVP control

      INTEGER  , intent(in)  ::  NLAYERS, NTOTAL, N_SUB, N_SUP

!  Number of streams

      INTEGER  , intent(in)  ::  NSTREAMS, NSTREAMS_2

!  Linearization control

      INTEGER  , intent(in)  ::  N_WEIGHTFUNCS

!  Input optical properties after delta-M scaling. 
!     These are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: DELTAU_VERT ( MAXLAYERS )
      REAL(fpk), intent(in)  :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: SOLARBEAM_CUTOFF(MAXBEAMS)

!  Quadrature input

      REAL(fpk), intent(in)  :: QUAD_STRMWTS ( MAXSTREAMS )

!  surface stuff
!  -------------

!  Factor = 1+delta(m,0), albedo

      REAL(fpk), intent(in)  :: SURFACE_FACTOR, ALBEDO

!  Fourier components of BRDF, in the following order (same all threads)
!    incident quadrature streams, reflected quadrature streams

      REAL(fpk), intent(in)  :: BRDF_F ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTREAMS )

!  4/9/19  RF_DIRECT_BEAM is the reflected beam (excludes SLTERM)

      REAL(fpk), intent(in)  :: RF_DIRECT_BEAM ( MAXSTREAMS, MAXBEAMS )
      REAL(fpk), intent(in)  :: SLTERM         ( MAXSTREAMS )

!  Derived optical thickness inputs

      REAL(fpk), intent(in)  :: DELTAU_SLANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS )

!  Beam quantities
!  ---------------

!  Average-secants, Initial and average-secant transmittance factors.

      REAL(fpk), intent(in)  :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS )

!  BVProblem inputs
!  ----------------

!  Matrix, Band-matrix

      REAL(fpk), intent(in)  :: SMAT2   (MAXSTREAMS_2,MAXSTREAMS_2)
      REAL(fpk), intent(in)  :: BANDMAT2(MAXBANDTOTAL,MAXTOTAL)

!  Pivot matrices

      INTEGER  , intent(in)  ::  IPIVOT  (MAXTOTAL)
      INTEGER  , intent(in)  ::  SIPIVOT (MAXSTREAMS_2)

!  Masking

      INTEGER  , intent(in)  ::  LCONMASK(MAXSTREAMS,MAXLAYERS)
      INTEGER  , intent(in)  ::  MCONMASK(MAXSTREAMS,MAXLAYERS)

!  Homogeneous and thermal solution variables
!  ------------------------------------------

!  Eigensolutions, eigenstream transmittances

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Solution constants of integration

      REAL(fpk), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)

!  Green functions
!  ---------------

!  Green function Multipliers for solution

      REAL(fpk), intent(in)  :: GFUNC_UP(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GFUNC_DN(MAXSTREAMS,MAXLAYERS)

!  Holding arrays for Multiplier coefficients

      REAL(fpk), intent(in)  :: GAMMA_M(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GAMMA_P(MAXSTREAMS,MAXLAYERS)

!  Green's function particular integral arrays

      REAL(fpk), intent(in)  :: ATERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: BTERM_SAVE(MAXSTREAMS,MAXLAYERS)

!  Linearized
!  ----------

!  Linearized Average-secants, Initial and average-secant transmittance factors.
!     LC_INITIAL_TRANS are Logarithmic derivatives (DOUBLE NORMALIZED)
!  4/9/19 Add linearization of TRANS_ATMOS for the waterleaing contribution

      REAL(fpk), intent(in)  :: LC_INITIAL_TRANS     ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LC_AVERAGE_SECANT    ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LC_T_DELT_MUBAR      ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LC_TRANS_ATMOS       ( MAXBEAMS, MAX_ATMOSWFS )

!  Linearized Eigenvalues, Eigensolutions, eigenstream transmittances

      REAL(fpk), intent(in)  :: L_KEIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized General Thermal solutions at the Lower boundary

      REAL(fpk), intent(in)  :: L_T_WUPPER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Saved quantities for the Green function solution
!     These are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: L_ATERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_BTERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  output arguments
!  ----------------

!  Linearized Beam solutions at the Lower and Upper layer boundaries

      REAL(fpk), intent(inout) :: L_WUPPER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(inout) :: L_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Solution constants of integration, and related quantities

      REAL(fpk), intent(inout) :: NCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(inout) :: PCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(inout) :: NCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(inout) :: PCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Exception handling. Updated 18 May 2010.

      INTEGER      , intent(out) :: STATUS
      CHARACTER*(*), intent(out) :: MESSAGE, TRACE

!  Local variables
!  ---------------

!  Column vectors for solving linearized BCs

      REAL(fpk)  :: COL2_WF    (MAXTOTAL,    MAX_ATMOSWFS)
      REAL(fpk)  :: SCOL2_WF   (MAXSTREAMS_2,MAX_ATMOSWFS)

!  Error tracing variables

      INTEGER     :: INFO
      CHARACTER*3 :: CI

!  Other local help variables 

      INTEGER     :: I, I1, Q, N, AA
      INTEGER     :: NS1, NS2         !mick eff 3/22/2017

!  Initialise Exception handling

      STATUS = LIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!mick eff 3/22/2017
!  Define some proxies

      NS1 = NSTREAMS + 1 ; NS2 = NSTREAMS_2

!  Linearization of the regular BVP case
!  =====================================

!  Set up the column vectors for Column linearizations
!  ---------------------------------------------------

!  Bulk: Compute the main column B' where AX = B'
! 4/9/19. Add water-leaving control, reflected direct beam, surface leaving linearization contribution

      CALL LC_BVP_COLUMN_SETUP &
           ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, DO_INCLUDE_DIRECTBEAM,   & ! Input, Flags
             DO_WATER_LEAVING, DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS,    & ! Input, FLags
             DO_LAYER_SCATTERING, TAYLOR_ORDER, FOURIER, IBEAM,            & ! Input, Flags and order
             NLAYERS, NTOTAL, NSTREAMS, NSTREAMS_2, N_WEIGHTFUNCS,         & ! Input, Numbers
             DELTAU_VERT, L_DELTAU_VERT, SOLARBEAM_CUTOFF, QUAD_STRMWTS,   & ! Input, optical and control
             SURFACE_FACTOR, ALBEDO, BRDF_F, RF_DIRECT_BEAM, SLTERM,       & ! Input, Surface Stuff
             DELTAU_SLANT, INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,    & ! Input, Beam Quantities
             T_DELT_EIGEN, XPOS, XNEG, LCON, MCON,                         & ! Input, Homogeneous
             GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE, & ! Input, Greens Function
             LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR,         & ! Input, Linearized Beam Quantities
             LC_TRANS_ATMOS, L_KEIGEN, L_T_DELT_EIGEN, L_XPOS, L_XNEG,     & ! Input, Linearized Homogeneous
             L_ATERM_SAVE, L_BTERM_SAVE, L_T_WUPPER, L_T_WLOWER,           & ! Input, Linearized Greens + thermal
             L_WUPPER, L_WLOWER, COL2_WF, SCOL2_WF )                         ! Output

!  BVP back-substitution: With compression (multilayers)
!  -----------------------------------------------------

      IF ( NLAYERS .GT. 1 ) THEN

!  LAPACK substitution (DGBTRS) using RHS column vector COL2_WF
!  BV solution for perturbed integration constants
!    ( call to LAPACK solver routine for back substitution )

        CALL DGBTRS  ( 'n', NTOTAL, N_SUB, N_SUP, N_WEIGHTFUNCS, &
              BANDMAT2, MAXBANDTOTAL, IPIVOT,  COL2_WF, MAXTOTAL, INFO )

!  Exception handling

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'Atmos_WEIGHTFUNCS for DGBTRS call in LC_BVP_SOLUTION_MASTER'
          STATUS  = LIDORT_SERIOUS
          RETURN
        ENDIF

!  Set integration constants NCON and PCON for +/- eigensolutions

        DO Q = 1, N_WEIGHTFUNCS
          DO N = 1, NLAYERS
            DO I = 1, NSTREAMS
              NCON(I,N,Q) = COL2_WF(LCONMASK(I,N),Q)
              PCON(I,N,Q) = COL2_WF(MCONMASK(I,N),Q)
            ENDDO
          ENDDO
        ENDDO

!  Solve the boundary problem: No compression, Single Layer only
!  -------------------------------------------------------------

      ELSE IF ( NLAYERS .EQ. 1 ) THEN

!  LAPACK substitution (DGETRS) using RHS column vector SCOL2_WF

        CALL DGETRS ( 'N', NTOTAL, N_WEIGHTFUNCS, SMAT2, MAXSTREAMS_2, SIPIVOT, &
                       SCOL2_WF, MAXSTREAMS_2, INFO )

!  Exception handling

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'Atmos_WEIGHTFUNCS for 1-layer: DGETRS call in LC_BVP_SOLUTION_MASTER'
          STATUS  = LIDORT_SERIOUS
          RETURN
        ENDIF

!  Set integration constants NCON and PCON for +/- eigensolutions

        N = 1
        DO Q = 1, N_WEIGHTFUNCS
          NCON(1:NSTREAMS,N,Q) = SCOL2_WF(1:NSTREAMS,Q)
          PCON(1:NSTREAMS,N,Q) = SCOL2_WF(NS1:NS2,Q)
        ENDDO

      ENDIF

!  linearized BVP results
!  ======================

!  Associated quantities

      DO Q = 1, N_WEIGHTFUNCS
        DO N = 1, NLAYERS
          DO AA = 1, NSTREAMS
            NCON_XVEC(1:NSTREAMS_2,AA,N,Q) = XPOS(1:NSTREAMS_2,AA,N) * NCON(AA,N,Q)
            PCON_XVEC(1:NSTREAMS_2,AA,N,Q) = XNEG(1:NSTREAMS_2,AA,N) * PCON(AA,N,Q)
          ENDDO
        ENDDO
      ENDDO

!  Debug.
!      if ( fourier_component.eq.3.and.ibeam.eq.2 ) then
!         do n = 20, nlayers
!            write(*,*)'regular',n, NCON(2,N,1), PCON(2,N,1),LCON(3,N),MCON(2,n)
!         enddo
!      endif

!  finish

      RETURN
END SUBROUTINE LC_BVP_SOLUTION_MASTER

!

SUBROUTINE LC_BVP_COLUMN_SETUP &
           ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, DO_INCLUDE_DIRECTBEAM,   & ! Input, Flags
             DO_WATER_LEAVING, DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS,    & ! Input, FLags
             DO_LAYER_SCATTERING, TAYLOR_ORDER, FOURIER, IBEAM,            & ! Input, Flags and order
             NLAYERS, NTOTAL, NSTREAMS, NSTREAMS_2, N_WEIGHTFUNCS,         & ! Input, Numbers
             DELTAU_VERT, L_DELTAU_VERT, SOLARBEAM_CUTOFF, QUAD_STRMWTS,   & ! Input, optical and control
             SURFACE_FACTOR, ALBEDO, BRDF_F, RF_DIRECT_BEAM, SLTERM,       & ! Input, Surface Stuff
             DELTAU_SLANT, INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,    & ! Input, Beam Quantities
             T_DELT_EIGEN, XPOS, XNEG, LCON, MCON,                         & ! Input, Homogeneous
             GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE, & ! Input, Greens Function
             LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR,         & ! Input, Linearized Beam Quantities
             LC_TRANS_ATMOS, L_KEIGEN, L_T_DELT_EIGEN, L_XPOS, L_XNEG,     & ! Input, Linearized Homogeneous
             L_ATERM_SAVE, L_BTERM_SAVE, L_T_WUPPER, L_T_WLOWER,           & ! Input, Linearized Greens + thermal
             L_WUPPER, L_WLOWER, COL2_WF, SCOL2_WF )                         ! Output

!  Linearized column vector setup (bulk property linearization)
! 4/9/19. Add water-leaving control, reflected direct beam, surface leaving linearization contribution
  
!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAXSTREAMS, MAXSTREAMS_2, MAXMOMENTS, MAXLAYERS, &
                                MAXBEAMS, MAX_ATMOSWFS, MAXTOTAL, ZERO

      IMPLICIT NONE

!  Subroutine arguments
!  ====================

!  Control and Optical
!  -------------------

!  Surface BRDF and inclusion flags

      LOGICAL  , intent(in)  ::  DO_INCLUDE_DIRECTBEAM
      LOGICAL  , intent(in)  ::  DO_INCLUDE_SURFACE
      LOGICAL  , intent(in)  ::  DO_BRDF_SURFACE
      LOGICAL  , intent(in)  ::  DO_WATER_LEAVING

!  Emission and solar source flags

      LOGICAL  , intent(in)  ::  DO_INCLUDE_THERMEMISS
      LOGICAL  , intent(in)  ::  DO_SOLAR_SOURCES

!  Local flags for the solution saving option

      LOGICAL  , intent(in)  ::  DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER  , intent(in)  :: TAYLOR_ORDER

!  Number of layers

      INTEGER  , intent(in)  ::  NLAYERS, NTOTAL

!  Number of streams

      INTEGER  , intent(in)  ::  NSTREAMS, NSTREAMS_2

!  Linearization control

      INTEGER  , intent(in)  ::  N_WEIGHTFUNCS

!  Fourier component, beam number

      INTEGER  , intent(in)  ::  FOURIER, IBEAM

!  Input optical properties after delta-M scaling. 
!     These are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: DELTAU_VERT ( MAXLAYERS )
      REAL(fpk), intent(in)  :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: SOLARBEAM_CUTOFF(MAXBEAMS)

!  Quadrature input

      REAL(fpk), intent(in)  :: QUAD_STRMWTS ( MAXSTREAMS )

!  surface stuff
!  -------------

!  Factor, albedo

      REAL(fpk), intent(in)  :: SURFACE_FACTOR, ALBEDO

!  Fourier components of BRDF, in the following order (same all threads)
!    incident quadrature streams, reflected quadrature streams

      REAL(fpk), intent(in)  :: BRDF_F ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTREAMS )

!  4/9/19 This is the reflected beam (excludes SLTERM)

      REAL(fpk), intent(in)  :: RF_DIRECT_BEAM ( MAXSTREAMS, MAXBEAMS )
      REAL(fpk), intent(in)  :: SLTERM         ( MAXSTREAMS )

!  Derived optical thickness inputs

      REAL(fpk), intent(in)  :: DELTAU_SLANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS )

!  Beam quantities
!  ---------------

!  Average-secants, Initial and average-secant transmittance factors.

      REAL(fpk), intent(in)  :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS )

!  Homogeneous
!  -----------

!  Eigensolutions, eigenstream transmittances

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Solution constants of integration

      REAL(fpk), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)

!  Green functions
!  ---------------

!  Green function Multipliers for solution

      REAL(fpk), intent(in)  :: GFUNC_UP(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GFUNC_DN(MAXSTREAMS,MAXLAYERS)

!  Holding arrays for Multiplier coefficients

      REAL(fpk), intent(in)  :: GAMMA_M(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GAMMA_P(MAXSTREAMS,MAXLAYERS)

!  Green's function particular integral arrays

      REAL(fpk), intent(in)  :: ATERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: BTERM_SAVE(MAXSTREAMS,MAXLAYERS)

!  Linearized
!  ----------

!  Linearized Average-secants, Initial and average-secant transmittance factors.
!     LC_INITIAL_TRANS are Logarithmic derivatives (DOUBLE NORMALIZED)
!  4/9/19 Add linearization of TRANS_ATMOS for the waterleaving contribution
      
      REAL(fpk), intent(in)  :: LC_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LC_T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LC_TRANS_ATMOS    ( MAXBEAMS,  MAX_ATMOSWFS )

!  Linearized Eigenvalues, Eigensolutions, eigenstream transmittances

      REAL(fpk), intent(in)  :: L_KEIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Thermal solutions at the Lower boundary

      REAL(fpk), intent(in)  :: L_T_WUPPER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Saved quantities for the Green function solution
!     These are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: L_ATERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_BTERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Output arguments
!  ----------------

!  Linearized Beam solutions at the Lower and Upper layer boundaries

      REAL(fpk), intent(out) :: L_WUPPER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: L_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Column vectors for solving linearized BBV Problem

      REAL(fpk), intent(out) :: COL2_WF    (MAXTOTAL,    MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: SCOL2_WF   (MAXSTREAMS_2,MAX_ATMOSWFS)

!  Local variables
!  ---------------

!  Help variables

      INTEGER    :: Q, AA, N, N1, I, I1, CM, C0, K, NV, M, IB
      REAL(fpk)  :: CPOS, CNEG, L_HOM, L_PARTIC, LCTERM
      REAL(fpk)  :: L_BEAM, L_HOMD, L_HOMU, FAC, FAC3
      LOGICAL    :: MODIFIED_BOUNDARY

!  Local linearized reflectance arrays

      REAL(fpk)  :: R2_L_BEAM(MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk)  :: R2_L_HOMP(MAXSTREAMS,MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk)  :: R2_L_HOMM(MAXSTREAMS,MAXSTREAMS,MAX_ATMOSWFS)
!mick eff 3/22/2017
      REAL(fpk)  :: HELP1 (NSTREAMS), HELP2 (NSTREAMS)

!  Initialise
!  ----------

      M  = FOURIER
      IB = IBEAM

!  Zero the results vectors

      COL2_WF(1:NTOTAL,1:MAX_ATMOSWFS) = ZERO

!  Copy already existing thermal linearizations
!    This is a very important zeroing.................!!!!!

      IF ( DO_INCLUDE_THERMEMISS ) THEN
        L_WUPPER(1:NSTREAMS_2,1:NLAYERS,1:N_WEIGHTFUNCS) = L_T_WUPPER(1:NSTREAMS_2,1:NLAYERS,1:N_WEIGHTFUNCS)
        L_WLOWER(1:NSTREAMS_2,1:NLAYERS,1:N_WEIGHTFUNCS) = L_T_WLOWER(1:NSTREAMS_2,1:NLAYERS,1:N_WEIGHTFUNCS)
      ELSE
        L_WUPPER(1:NSTREAMS_2,1:NLAYERS,1:N_WEIGHTFUNCS) = ZERO
        L_WLOWER(1:NSTREAMS_2,1:NLAYERS,1:N_WEIGHTFUNCS) = ZERO
      ENDIF

!  Top of first layer (TOA), UPPER boundary condition
!  --------------------------------------------------

      N = 1

!  Get the linearized beam solution for the first layer
!   Adds to the thermal solution if already present

      IF ( DO_SOLAR_SOURCES ) THEN
        CALL LC_BEAMSOLUTION_NEQK &
           ( DO_LAYER_SCATTERING, TAYLOR_ORDER,                    & ! Input, Flags and order
             NSTREAMS, NSTREAMS_2, N_WEIGHTFUNCS, M, IB, N,        & ! Input, Numbers
             DELTAU_VERT, L_DELTAU_VERT, SOLARBEAM_CUTOFF,         & ! Input, optical and control
             INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,          & ! Input, Beam Quantities
             LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR, & ! Input, Beam Quantities (Linearized)
             T_DELT_EIGEN, XPOS, L_KEIGEN, L_T_DELT_EIGEN, L_XPOS, & ! Input, Homogeneous solution stuff
             GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P,                 & ! Input, Greens Function
             ATERM_SAVE, BTERM_SAVE, L_ATERM_SAVE, L_BTERM_SAVE,   & ! Input, Greens Function
             L_WUPPER, L_WLOWER )                                    ! Output
      ENDIF

!  .. contribution WVAR from beam solution variations
!  .. contribution HVAR homogeneous (eigenvalue) solution variations

      DO Q = 1, N_WEIGHTFUNCS
        DO I = 1, NSTREAMS
          L_PARTIC = - L_WUPPER(I,N,Q)
          HELP1    =   T_DELT_EIGEN(1:NSTREAMS,N)   * L_XNEG(I,1:NSTREAMS,N,Q) &
                   + L_T_DELT_EIGEN(1:NSTREAMS,N,Q) *   XNEG(I,1:NSTREAMS,N)
          L_HOM    = DOT_PRODUCT ( LCON(1:NSTREAMS,N), L_XPOS(I,1:NSTREAMS,N,Q) ) &
                   + DOT_PRODUCT ( MCON(1:NSTREAMS,N), HELP1 )
          COL2_WF(I,Q) = L_PARTIC - L_HOM
        ENDDO
      ENDDO

!  Intermediate boundary conditions
!  --------------------------------

      DO N = 1, NLAYERS - 1

!  N1 is the layer below, C0 is the offset

        N1 = N + 1
        C0 = N*NSTREAMS_2 - NSTREAMS

!  Get the linearized beam solution for the next layer
!   Adds to the thermal solution if already present

        IF ( DO_SOLAR_SOURCES ) THEN
          CALL LC_BEAMSOLUTION_NEQK &
           ( DO_LAYER_SCATTERING, TAYLOR_ORDER,                    & ! Input, Flags and order
             NSTREAMS, NSTREAMS_2, N_WEIGHTFUNCS, M, IB, N1,       & ! Input, Numbers
             DELTAU_VERT, L_DELTAU_VERT, SOLARBEAM_CUTOFF,         & ! Input, optical and control
             INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,          & ! Input, Beam Quantities
             LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR, & ! Input, Beam Quantities (Linearized)
             T_DELT_EIGEN, XPOS, L_KEIGEN, L_T_DELT_EIGEN, L_XPOS, & ! Input, Homogeneous solution stuff
             GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P,                 & ! Input, Greens Function
             ATERM_SAVE, BTERM_SAVE, L_ATERM_SAVE, L_BTERM_SAVE,   & ! Input, Greens Function
             L_WUPPER, L_WLOWER )                                    ! Output
        ENDIF

!  .. 2 contributions to L_BEAM, from variations L_WUPPER L_WLOWER 
!  .. 2 contributions to L_HOM,  from variations above and below

        DO Q = 1, N_WEIGHTFUNCS
          DO I = 1, NSTREAMS_2
            CM = C0 + I
            L_PARTIC = L_WUPPER(I,N1,Q) - L_WLOWER(I,N,Q)
            HELP1    =   T_DELT_EIGEN(1:NSTREAMS,N1)   * L_XNEG(I,1:NSTREAMS,N1,Q) &
                     + L_T_DELT_EIGEN(1:NSTREAMS,N1,Q) *   XNEG(I,1:NSTREAMS,N1)
            L_HOMU   = DOT_PRODUCT ( LCON(1:NSTREAMS,N1), L_XPOS(I,1:NSTREAMS,N1,Q) ) &
                     + DOT_PRODUCT ( MCON(1:NSTREAMS,N1), HELP1 )
            HELP2    =   T_DELT_EIGEN(1:NSTREAMS,N)   * L_XPOS(I,1:NSTREAMS,N,Q) &
                     + L_T_DELT_EIGEN(1:NSTREAMS,N,Q) *   XPOS(I,1:NSTREAMS,N)
            L_HOMD   = DOT_PRODUCT ( LCON(1:NSTREAMS,N), HELP2 ) &
                     + DOT_PRODUCT ( MCON(1:NSTREAMS,N), L_XNEG(I,1:NSTREAMS,N,Q) )
            L_HOM    = L_HOMU - L_HOMD
            COL2_WF(CM,Q) = L_PARTIC + L_HOM
          ENDDO
        ENDDO

!  End layer

      ENDDO

!  LOWER layer 
!  -----------

      N = NLAYERS
      MODIFIED_BOUNDARY = .TRUE.

!  get the linearized downward-reflected term

      CALL LC_BVP_SURFACE_SETUP &
          ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, NSTREAMS, & ! Input
            NLAYERS, MODIFIED_BOUNDARY, FOURIER,           & ! Input
            SURFACE_FACTOR, ALBEDO, BRDF_F, N_WEIGHTFUNCS, & ! Input
            QUAD_STRMWTS, T_DELT_EIGEN, L_T_DELT_EIGEN,    & ! Input
            XPOS, L_XPOS, L_XNEG, L_WLOWER,                & ! Input
            R2_L_BEAM, R2_L_HOMP, R2_L_HOMM )                ! Output 

!  Compute the solution

      C0 = (N-1)*NSTREAMS_2 + NSTREAMS
      DO Q = 1, N_WEIGHTFUNCS
        DO I = 1, NSTREAMS
          CM = C0 + I
          I1 = I + NSTREAMS
          L_PARTIC = L_WLOWER(I1,N,Q) - R2_L_BEAM(I,Q)
          HELP1 =   T_DELT_EIGEN(1:NSTREAMS,N)   * L_XPOS(I1,1:NSTREAMS,N,Q) &
                + L_T_DELT_EIGEN(1:NSTREAMS,N,Q) *   XPOS(I1,1:NSTREAMS,N) &
                - R2_L_HOMP(I,1:NSTREAMS,Q)
          HELP2 = L_XNEG(I1,1:NSTREAMS,N,Q) - R2_L_HOMM(I,1:NSTREAMS,Q)
          L_HOM = DOT_PRODUCT ( LCON(1:NSTREAMS,N), HELP1 ) &
                + DOT_PRODUCT ( MCON(1:NSTREAMS,N), HELP2 )
          COL2_WF(CM,Q) = - L_PARTIC - L_HOM
        ENDDO
      ENDDO

!  Add direct beam variation to Final boundary
!  -------------------------------------------

! 4/9/19. There was a bug here. Formula was only correct for the reflected beam
      
      IF ( DO_INCLUDE_DIRECTBEAM ) THEN
        DO Q = 1, N_WEIGHTFUNCS
          DO I = 1, NSTREAMS
            CM = C0 + I
            COL2_WF(CM,Q) = COL2_WF(CM,Q) &
              - RF_DIRECT_BEAM(I,IBEAM) * DOT_PRODUCT(L_DELTAU_VERT(Q,1:NLAYERS),DELTAU_SLANT(N,1:NLAYERS,IBEAM))
          ENDDO
        ENDDO
      ENDIF

!  4/9/19 Add the linearization due to Adjusted waterleaving term

      IF ( DO_WATER_LEAVING .and. FOURIER .eq. 0 ) THEN
         DO Q = 1, N_WEIGHTFUNCS
            LCTERM = LC_TRANS_ATMOS(IBEAM,Q)
            DO I = 1, NSTREAMS
               CM = C0 + I
               COL2_WF(CM,Q) = COL2_WF(CM,Q) + LCTERM * SLTERM(I)
            ENDDO
         ENDDO
      ENDIF
        
!  Single layer - just copy

      IF ( NLAYERS .eq. 1 ) THEN
        SCOL2_WF(1:NTOTAL,1:N_WEIGHTFUNCS) = COL2_WF(1:NTOTAL,1:N_WEIGHTFUNCS)
      ENDIF

!  debug
!      do i = 1, ntotal
!        write(58,*)i,COL2_WF(i,1)
!      enddo
!      stop'LC_BVP'
 
!  Finish

      RETURN
END SUBROUTINE LC_BVP_COLUMN_SETUP

!

SUBROUTINE LC_BVP_SURFACE_SETUP &
          ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, NSTREAMS, & ! Input
            NLAYERS, MODIFIED_BCL4, FOURIER,               & ! Input
            SURFACE_FACTOR, ALBEDO, BRDF_F, N_WEIGHTFUNCS, & ! Input
            QUAD_STRMWTS, T_DELT_EIGEN, L_T_DELT_EIGEN,    & ! Input
            XPOS, L_XPOS, L_XNEG, L_WLOWER,                & ! Input
            R2_L_BEAM, R2_L_HOMP, R2_L_HOMM )                ! Output 

!  Linearized surface reflectance terms

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAXSTREAMS, MAXSTREAMS_2, MAXMOMENTS, MAXLAYERS, &
                                MAX_ATMOSWFS, ZERO

      IMPLICIT NONE

!  input arguments
!  ---------------

!  BRDF flag

      LOGICAL  , intent(in)  ::  DO_BRDF_SURFACE

!  Number of streams and layers

      INTEGER  , intent(in)  ::  NSTREAMS, NLAYERS

!  Number of weighting functions

      INTEGER  , intent(in)  ::  N_WEIGHTFUNCS

!  Fourier component

      INTEGER  , intent(in)  ::  FOURIER

!  Flag for type of boundary condition

      LOGICAL  , intent(in)  ::  MODIFIED_BCL4

!  overall surface flag and surface factor, albedo

      LOGICAL  , intent(in)  :: DO_INCLUDE_SURFACE
      REAL(fpk), intent(in)  :: SURFACE_FACTOR, ALBEDO

!  Fourier components of BRDF, in the following order (same all threads)
!    ( New code, 23 March 2010 )
!    incident quadrature streams, reflected quadrature streams

      REAL(fpk), intent(in)  :: BRDF_F ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTREAMS )

!  Quadrature input

      REAL(fpk), intent(in)  :: QUAD_STRMWTS ( MAXSTREAMS )

!  Eigenvector solutions

      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  transmittance factors for +/- eigenvalues

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)

!  Linearized Eigensolutions

      REAL(fpk), intent(in)  :: L_XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized General beam solutions at the Lower boundary

      REAL(fpk), intent(in)  :: L_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized transmittances, homogeneous solutions

      REAL(fpk), intent(in)  :: L_T_DELT_EIGEN ( MAXSTREAMS, MAXLAYERS, MAX_ATMOSWFS )

!  Output arguments
!  ----------------

      REAL(fpk), intent(out) :: R2_L_BEAM(MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: R2_L_HOMP(MAXSTREAMS,MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: R2_L_HOMM(MAXSTREAMS,MAXSTREAMS,MAX_ATMOSWFS)

!  Local variables
!  ---------------

      REAL(fpk)  :: PV_W(MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk)  :: HV_P(MAXSTREAMS,MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk)  :: HV_M(MAXSTREAMS,MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk)  :: HSP_U, HSM_U, REFL_P, REFL_B, REFL_M
      REAL(fpk)  :: FACTOR
      INTEGER    :: AA, J, Q, N, I, M

!  Initial section
!  ---------------

!  Always zero the result to start

      R2_L_BEAM(1:NSTREAMS,1:N_WEIGHTFUNCS) = ZERO
      R2_L_HOMP(1:NSTREAMS,1:NSTREAMS,1:N_WEIGHTFUNCS) = ZERO
      R2_L_HOMM(1:NSTREAMS,1:NSTREAMS,1:N_WEIGHTFUNCS) = ZERO

!  Fourier component

      M = FOURIER

!  Return if no albedo

      IF ( .NOT. DO_INCLUDE_SURFACE ) RETURN

!  Set up Auxiliary arrays
!  -----------------------

      N = NLAYERS

!  Particular integral parts

      DO Q = 1, N_WEIGHTFUNCS
        PV_W(1:NSTREAMS,Q) = L_WLOWER(1:NSTREAMS,N,Q) * QUAD_STRMWTS(1:NSTREAMS)
      ENDDO

!    Modified boundary condition: homogeneous parts

      IF ( MODIFIED_BCL4 ) THEN
        DO Q = 1, N_WEIGHTFUNCS
          DO AA = 1, NSTREAMS
            HV_P(1:NSTREAMS,AA,Q) = QUAD_STRMWTS(1:NSTREAMS) &
                                  * ( L_XPOS(1:NSTREAMS,AA,N,Q) *   T_DELT_EIGEN(AA,N) &
                                      + XPOS(1:NSTREAMS,AA,N)   * L_T_DELT_EIGEN(AA,N,Q) )
            HV_M(1:NSTREAMS,AA,Q) = QUAD_STRMWTS(1:NSTREAMS)*L_XNEG(1:NSTREAMS,AA,N,Q)
          ENDDO
        ENDDO
      ENDIF

!  Integrated Downward reflection (Calculation, Lambertian case)
!  -------------------------------------------------------------

     IF ( .not. DO_BRDF_SURFACE ) THEN

!  Amplitude

        FACTOR = SURFACE_FACTOR * ALBEDO

!  Only a solution for the isotropic part

        IF ( FOURIER .EQ. 0 ) THEN
          DO Q = 1, N_WEIGHTFUNCS

!  Particular solution

            R2_L_BEAM(1:NSTREAMS,Q) = FACTOR * SUM( PV_W(1:NSTREAMS,Q) )

!  Homogeneous solutions (only for modified BC)

            IF ( MODIFIED_BCL4 ) THEN
              DO AA = 1, NSTREAMS
                R2_L_HOMP(1:NSTREAMS,AA,Q) = FACTOR * SUM ( HV_P(1:NSTREAMS,AA,Q) )
                R2_L_HOMM(1:NSTREAMS,AA,Q) = FACTOR * SUM ( HV_M(1:NSTREAMS,AA,Q) )
              ENDDO
            ENDIF

!  end parameter loop

          ENDDO
        ENDIF

!  Integrated Downward reflection (Calculation, Bidirectional case)
!  ----------------------------------------------------------------

      ELSE

!mick eff 3/22/2017 - move I loop
        DO Q = 1, N_WEIGHTFUNCS

!  Particular solutions
!     @@@ Rob Fix 2/3/11,  Reverse J,I ---> I,J (J is incident)

          DO I = 1, NSTREAMS
            R2_L_BEAM(I,Q) = SURFACE_FACTOR * DOT_PRODUCT ( PV_W(1:NSTREAMS,Q), BRDF_F(M,I,1:NSTREAMS) )
          ENDDO

!  Homogeneous solutions
!     @@@ Rob Fix 2/3/11,  Reverse J,I ---> I,J (J is incident)

          IF ( MODIFIED_BCL4 ) THEN
            DO AA = 1, NSTREAMS
              DO I = 1, NSTREAMS
                R2_L_HOMP(I,AA,Q) = SURFACE_FACTOR * DOT_PRODUCT ( HV_P(1:NSTREAMS,AA,Q), BRDF_F(M,I,1:NSTREAMS) )
                R2_L_HOMM(I,AA,Q) = SURFACE_FACTOR * DOT_PRODUCT ( HV_M(1:NSTREAMS,AA,Q), BRDF_F(M,I,1:NSTREAMS) )
              ENDDO
            ENDDO
          ENDIF

!  end parameter and stream loops

        ENDDO

!  End BRDF clause

      ENDIF

!  Finish

      RETURN
END SUBROUTINE LC_BVP_SURFACE_SETUP

!

SUBROUTINE LC_BEAMSOLUTION_NEQK &
       ( DO_LAYER_SCATTERING, TAYLOR_ORDER,                    & ! Input, Flag and order
         NSTREAMS, NSTREAMS_2, N_WEIGHTFUNCS, M, IB, N,        & ! Input, Numbers
         DELTAU_VERT, L_DELTAU_VERT, SOLARBEAM_CUTOFF,         & ! Input, optical and control
         INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,          & ! Input, Beam Quantities
         LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR, & ! Input, Beam Quantities (Linearized)
         T_DELT_EIGEN, XPOS, L_KEIGEN, L_T_DELT_EIGEN, L_XPOS, & ! Input, Homogeneous solution stuff
         GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P,                 & ! Input, Greens Function
         ATERM_SAVE, BTERM_SAVE, L_ATERM_SAVE, L_BTERM_SAVE,   & ! Input, Greens Function
         L_WUPPER, L_WLOWER )                                    ! Output

!  Linearization of beam particular integral in layer N, due to column.
!   This is the bulk property linearization

!  Completely rewritten, 11 October 2013. Version 3.7
!    * Final version with Taylor series expansions
!    * Removal of LC_GAMMA

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAXSTREAMS, MAXSTREAMS_2, MAXMOMENTS, MAXLAYERS, &
                                MAXBEAMS, MAX_ATMOSWFS, ZERO, ONE, TAYLOR_SMALL

      IMPLICIT NONE

!  subroutine arguments
!  ====================

!  Control and Optical
!  -------------------

!  Local flags for the solution saving option

      LOGICAL  , intent(in)  ::  DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER  , intent(in)  :: TAYLOR_ORDER

!  Number of streams

      INTEGER  , intent(in)  ::  NSTREAMS, NSTREAMS_2

!  number of varying parameters (input)

      INTEGER  , intent(in)  ::  N_WEIGHTFUNCS

!  Fourier number, beam index, layer index

      INTEGER  , intent(in)  ::  M, IB, N

!  Input optical properties after delta-M scaling. 
!     These are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: DELTAU_VERT ( MAXLAYERS )
      REAL(fpk), intent(in)  :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: SOLARBEAM_CUTOFF(MAXBEAMS)

!  Beam quantities
!  ---------------

!  Average-secants, Initial and average-secant transmittance factors.

      REAL(fpk), intent(in)  :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS )

!  Linearized Average-secants, Initial and average-secant transmittance factors.
!     LC_INITIAL_TRANS are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: LC_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LC_T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Homogeneous solution variables
!  ------------------------------

!  Eigensolutions XPOS, eigenstream transmittances

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Linearized Eigenvalues, Eigensolutions XPOS, eigenstream transmittances

      REAL(fpk), intent(in)  :: L_KEIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Green functions
!  ---------------

!  Green function Multipliers for solution

      REAL(fpk), intent(in)  :: GFUNC_UP(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GFUNC_DN(MAXSTREAMS,MAXLAYERS)

!  Holding arrays for Multiplier coefficients

      REAL(fpk), intent(in)  :: GAMMA_M(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GAMMA_P(MAXSTREAMS,MAXLAYERS)

!  Green's function particular integral arrays

      REAL(fpk), intent(in)  :: ATERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: BTERM_SAVE(MAXSTREAMS,MAXLAYERS)

!  Linearized Saved quantities for the Green function solution
!     These are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: L_ATERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_BTERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  output arguments
!  ----------------

!  Linearized beam solutions at the Lower and Upper layer boundaries

      REAL(fpk), intent(inout) :: L_WUPPER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(inout) :: L_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Local variables
!  ===============

!  Local linearized Green's function multipliers

      REAL(fpk)  :: L_GFUNC_UP(MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk)  :: L_GFUNC_DN(MAXSTREAMS,MAX_ATMOSWFS)

!  Help variables

      INTEGER    :: AA, I, I1, Q
      REAL(fpk)  :: S_P_U, S_P_L, S_M_U, S_M_L
      REAL(fpk)  :: CONST, WDEL, ZDEL, ZWDEL, AST, BST, EPS, DELTA, LAM, MULT
      REAL(fpk)  :: L_WDEL, L_ZDEL, L_LAM, L_KEG, L_DELTA, L_AST, L_BST, L_MULT(MAX_ATMOSWFS)
      REAL(fpk)  :: LBSOL(MAXSTREAMS_2,MAX_ATMOSWFS,2)

!  Green's function solution
!  =========================

!  No linearized particular solution beyond the cutoff layer. ALSO -
!  Nothing if layer is inactive (Does not depend on solution saving)
!    Exit Immediately (Solutions have already been zeroed).

!      IF ((DO_SOLUTION_SAVING.AND..NOT.DO_LAYER_SCATTERING(M,N)) .OR. (N .GT.SOLARBEAM_CUTOFF(IB))) RETURN
      IF ( .NOT.DO_LAYER_SCATTERING(M,N) .OR. N.GT.SOLARBEAM_CUTOFF(IB) ) RETURN

!  Linearizations of optical depth integrations (Linearized Green function multipliers)
!  ------------------------------------------------------------------------------------

      CONST   = INITIAL_TRANS(N,IB)
      WDEL    = T_DELT_MUBAR(N,IB)

!  Start discrete ordinate loop

      DO AA = 1, NSTREAMS

         ZDEL  = T_DELT_EIGEN(AA,N)
         ZWDEL = ZDEL * WDEL
         AST   = CONST * ATERM_SAVE(AA,N) 
         BST   = CONST * BTERM_SAVE(AA,N) 

!  Downwelling, Make allowances for Taylor series

        IF ( ABS(GAMMA_M(AA,N)) .LT. TAYLOR_SMALL ) THEN
           EPS   = GAMMA_M(AA,N)
           DELTA = DELTAU_VERT(N)
           LAM   = AVERAGE_SECANT(N,IB)
           DO Q = 1, N_WEIGHTFUNCS
              L_LAM  = LC_AVERAGE_SECANT(N,IB,Q)
              L_KEG  = L_KEIGEN(AA,N,Q)
              L_DELTA = L_DELTAU_VERT(Q,N) * DELTA ! Input is single normalized
              CALL TAYLOR_SERIES_L_1 ( TAYLOR_ORDER, EPS, DELTA, L_DELTA, L_KEG, L_LAM, WDEL, LAM, L_MULT(Q) )
           ENDDO
        ELSE
           MULT = ( ZDEL - WDEL ) / GAMMA_M(AA,N)
           DO Q = 1, N_WEIGHTFUNCS
              L_ZDEL =  L_T_DELT_EIGEN  (AA,N,Q) ; L_KEG  = L_KEIGEN(AA,N,Q)
              L_WDEL =  LC_T_DELT_MUBAR (N,IB,Q) ; L_LAM  = LC_AVERAGE_SECANT(N,IB,Q)
              L_MULT(Q) = ( ( L_ZDEL - L_WDEL ) - MULT * (L_LAM - L_KEG) ) / GAMMA_M(AA,N)
           ENDDO
        ENDIF
        DO Q = 1, N_WEIGHTFUNCS
           L_AST =  LC_INITIAL_TRANS(N,IB,Q)  + L_ATERM_SAVE(AA,N,Q)
           L_GFUNC_DN(AA,Q) = GFUNC_DN(AA,N) * L_AST + L_MULT(Q) * AST
        ENDDO

!  Upwelling

        MULT = ( ONE - ZWDEL ) / GAMMA_P(AA,N)
        DO Q = 1, N_WEIGHTFUNCS
           L_ZDEL =  L_T_DELT_EIGEN  (AA,N,Q) ; L_KEG  = L_KEIGEN(AA,N,Q)
           L_WDEL =  LC_T_DELT_MUBAR (N,IB,Q) ; L_LAM  = LC_AVERAGE_SECANT(N,IB,Q)
           L_MULT(Q) = - ( L_ZDEL*WDEL + L_WDEL*ZDEL + MULT*(L_LAM+L_KEG) ) / GAMMA_P(AA,N)
        ENDDO
        DO Q = 1, N_WEIGHTFUNCS
           L_BST =  LC_INITIAL_TRANS(N,IB,Q)  + L_BTERM_SAVE(AA,N,Q)
           L_GFUNC_UP(AA,Q) = GFUNC_UP(AA,N) * L_BST + L_MULT(Q) * BST
        ENDDO

!  End discrete ordinate loop

      ENDDO

!  Set linearized form of particular integral at boundaries
!  --------------------------------------------------------

!mick eff 3/22/2017
      DO Q = 1, N_WEIGHTFUNCS
         DO I = 1, NSTREAMS
           I1 = I + NSTREAMS
           !S_P_U = ZERO ; S_M_U = ZERO
           !S_M_L = ZERO ; S_P_L = ZERO
           !DO AA = 1, NSTREAMS
           !   S_P_U = S_P_U + L_GFUNC_UP(AA,Q) * XPOS(I1,AA,N) + GFUNC_UP(AA,N) * L_XPOS(I1,AA,N,Q)
           !   S_M_U = S_M_U + L_GFUNC_UP(AA,Q) * XPOS(I,AA,N)  + GFUNC_UP(AA,N) * L_XPOS(I,AA,N,Q)
           !   S_M_L = S_M_L + L_GFUNC_DN(AA,Q) * XPOS(I1,AA,N) + GFUNC_DN(AA,N) * L_XPOS(I1,AA,N,Q)
           !   S_P_L = S_P_L + L_GFUNC_DN(AA,Q) * XPOS(I,AA,N)  + GFUNC_DN(AA,N) * L_XPOS(I,AA,N,Q)
           !ENDDO
           !LBSOL(I,Q,1)  = S_P_U ; LBSOL(I1,Q,1) = S_M_U
           !LBSOL(I1,Q,2) = S_M_L ; LBSOL(I,Q,2)  = S_P_L
           LBSOL(I,Q,1)  = DOT_PRODUCT ( L_GFUNC_UP(1:NSTREAMS,Q),  XPOS(I1,1:NSTREAMS,N)   ) &
                         + DOT_PRODUCT (   GFUNC_UP(1:NSTREAMS,N),L_XPOS(I1,1:NSTREAMS,N,Q) )
           LBSOL(I1,Q,1) = DOT_PRODUCT ( L_GFUNC_UP(1:NSTREAMS,Q),  XPOS(I, 1:NSTREAMS,N)   ) &
                         + DOT_PRODUCT (   GFUNC_UP(1:NSTREAMS,N),L_XPOS(I, 1:NSTREAMS,N,Q) )
           LBSOL(I1,Q,2) = DOT_PRODUCT ( L_GFUNC_DN(1:NSTREAMS,Q),  XPOS(I1,1:NSTREAMS,N)   ) &
                         + DOT_PRODUCT (   GFUNC_DN(1:NSTREAMS,N),L_XPOS(I1,1:NSTREAMS,N,Q) )
           LBSOL(I,Q,2)  = DOT_PRODUCT ( L_GFUNC_DN(1:NSTREAMS,Q),  XPOS(I, 1:NSTREAMS,N)   ) &
                         + DOT_PRODUCT (   GFUNC_DN(1:NSTREAMS,N),L_XPOS(I, 1:NSTREAMS,N,Q) )
        ENDDO
      ENDDO

!  Add to existing solution

!mick eff 3/22/2017
      !DO Q = 1, N_WEIGHTFUNCS
      !   DO I = 1, NSTREAMS_2
      !      L_WUPPER(I,N,Q) = L_WUPPER(I,N,Q) + LBSOL(I,Q,1)
      !      L_WLOWER(I,N,Q) = L_WLOWER(I,N,Q) + LBSOL(I,Q,2)
      !   ENDDO
      !ENDDO
      DO Q = 1, N_WEIGHTFUNCS
        L_WUPPER(1:NSTREAMS_2,N,Q) = L_WUPPER(1:NSTREAMS_2,N,Q) + LBSOL(1:NSTREAMS_2,Q,1)
        L_WLOWER(1:NSTREAMS_2,N,Q) = L_WLOWER(1:NSTREAMS_2,N,Q) + LBSOL(1:NSTREAMS_2,Q,2)
      ENDDO

!  Finish

      RETURN
END SUBROUTINE LC_BEAMSOLUTION_NEQK

!

SUBROUTINE LC_BVPTEL_SOLUTION_MASTER &
       ( DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM, DO_LAYER_SCATTERING,           & ! Input, Flags
         TAYLOR_ORDER, FOURIER, IBEAM, NLAYERS, NSTREAMS, NSTREAMS_2, N_COLUMNWFS, & ! Input, Numbers
         N_SUB, N_SUP, ACTIVE_LAYERS, NLAYERS_TEL, N_BVTELMATRIX_SIZE,             & ! Input, Numbers
         DELTAU_VERT, L_DELTAU_VERT, SOLARBEAM_CUTOFF,                             & ! Input, optical and control
         SURFACE_FACTOR, QUAD_STRMWTS, BRDF_F, DIRECT_BEAM, DELTAU_SLANT,          & ! Input, Surface +DB Stuff
         INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR, T_DELT_DISORDS,              & ! Input, Beam Quantities
         T_DELT_EIGEN, XPOS, XNEG, WLOWER, LCON, MCON,                             & ! Input, Homogeneous + PI
         GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,             & ! Input, Greens Function
         BANDTELMAT2, SMAT2, IPIVOTTEL, SIPIVOT, LCON_XVEC, MCON_XVEC,             & ! Input, BVP Bandmat
         LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR, L_T_DELT_DISORDS,   & ! Input , Linearized Beam Quantities
         L_KEIGEN, L_T_DELT_EIGEN, L_XPOS, L_XNEG, L_ATERM_SAVE, L_BTERM_SAVE,     & ! Input , Linearized solutions
         NCON, PCON, NCON_XVEC, PCON_XVEC, L_WUPPER, L_WLOWER,                     & ! Output, Linearized Constants + PI
         STATUS, MESSAGE, TRACE )                                                    ! Output - Exception handling

!  Linearization of the Telescoped Boundary Problem Solution
!   Version 3.8. Major extension to BRDFs with TELESCOPING. Implemented May 2016.

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAXSTREAMS, MAXSTREAMS_2, MAXMOMENTS, MAXLAYERS, &
                                MAXBEAMS, MAX_ATMOSWFS, MAXBANDTOTAL, MAXTOTAL,  &
                                ZERO, LIDORT_SUCCESS, LIDORT_SERIOUS

      IMPLICIT NONE

!  input arguments
!  ===============

!  Control and Optical
!  -------------------

!  surface control. Always going to be BRDF here

      LOGICAL  , intent(in)  :: DO_INCLUDE_SURFACE
      LOGICAL  , intent(in)  :: DO_INCLUDE_DIRECTBEAM

!  Local flags for the solution saving option

      LOGICAL  , intent(in)  ::  DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER  , intent(in)  :: TAYLOR_ORDER

!  Fourier number, beam index

      INTEGER  , intent(in)  ::  FOURIER, IBEAM

!  Number of layers

      INTEGER  , intent(in)  ::  NLAYERS

!  Number of streams

      INTEGER  , intent(in)  ::  NSTREAMS, NSTREAMS_2

!  Linearization control: number of varying parameters (input)

      INTEGER  , intent(in)  ::  N_COLUMNWFS

!  BVProblem Band matrix control

      INTEGER  , intent(in)  ::  N_SUB, N_SUP

!  Number of telescoped layers

      INTEGER  , intent(in)  ::  NLAYERS_TEL

!  Active layers for telescoping

      INTEGER  , intent(in)  ::  ACTIVE_LAYERS ( MAXLAYERS )

!  Size of BVP matrix for telescoped 

      INTEGER  , intent(in)  ::  N_BVTELMATRIX_SIZE

!  Input optical properties after delta-M scaling. 
!     These are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: DELTAU_VERT ( MAXLAYERS )
      REAL(fpk), intent(in)  :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: SOLARBEAM_CUTOFF(MAXBEAMS)

!  surface control

      REAL(fpk), intent(in)  :: SURFACE_FACTOR

!  Quadrature

      REAL(fpk), intent(in)  :: QUAD_STRMWTS(MAXSTREAMS)

!  Fourier components of BRDF
!    incident quadrature streams, reflected quadrature streams

      REAL(fpk), intent(in)  :: BRDF_F ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTREAMS )

!  Direct beam

      REAL(fpk), intent(in)  :: DIRECT_BEAM ( MAXSTREAMS, MAXBEAMS )

!  Derived optical thickness inputs

      REAL(fpk), intent(in)  :: DELTAU_SLANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS )

!  Beam quantities
!  ---------------

!  Average-secants, Initial and average-secant transmittance factors.

      REAL(fpk), intent(in)  :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS )

!  BVProblem inputs
!  ----------------

!  Matrix, Band-matrix

      REAL(fpk), intent(in)  :: SMAT2       (MAXSTREAMS_2,MAXSTREAMS_2)
      REAL(fpk), intent(in)  :: BANDTELMAT2 (MAXBANDTOTAL,MAXTOTAL)

!  Pivot matrices

      INTEGER  , intent(in)  ::  IPIVOTTEL  (MAXTOTAL)
      INTEGER  , intent(in)  ::  SIPIVOT    (MAXSTREAMS_2)

!  discrete ordinate factors (BVP telescoping, solutions saving)

      REAL(fpk), intent(in)  :: T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: L_T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Homogeneous solution variables
!  ------------------------------

!  Eigensolutions, eigenstream transmittances

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Particular integral at lower boundary

      REAL(fpk), intent(in)  :: WLOWER ( MAXSTREAMS_2, MAXLAYERS )

!  Solution constants of integration

      REAL(fpk), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)

!  Solution constants of integration multiplied by homogeneous solutions

      REAL(fpk), intent(in)  :: LCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Green functions
!  ---------------

!  Green function Multipliers for solution

      REAL(fpk), intent(in)  :: GFUNC_UP(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GFUNC_DN(MAXSTREAMS,MAXLAYERS)

!  Holding arrays for Multiplier coefficients

      REAL(fpk), intent(in)  :: GAMMA_M(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GAMMA_P(MAXSTREAMS,MAXLAYERS)

!  Green's function particular integral arrays

      REAL(fpk), intent(in)  :: ATERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: BTERM_SAVE(MAXSTREAMS,MAXLAYERS)

!  Linearized
!  ----------

!  Linearized Average-secants, Initial and average-secant transmittance factors.
!     LC_INITIAL_TRANS are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: LC_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LC_T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized Eigenvalues, Eigensolutions, eigenstream transmittances

      REAL(fpk), intent(in)  :: L_KEIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Saved quantities for the Green function solution
!     These are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: L_ATERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_BTERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Output arguments
!  ----------------

!  Linearized beam solutions at the Lower and Upper layer boundaries

      REAL(fpk), intent(out) :: L_WUPPER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: L_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Solution constants of integration, and related quantities

      REAL(fpk), intent(out) :: NCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: PCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: NCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: PCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Exception handling. Updated 18 May 2010.

      INTEGER      , intent(out) :: STATUS
      CHARACTER*(*), intent(out) :: MESSAGE, TRACE

!  Local variables
!  ---------------

!  Column vectors for solving linearized BCs

      REAL(fpk)  :: COLTEL2_WF (MAXTOTAL,    MAX_ATMOSWFS)
      REAL(fpk)  :: SCOL2_WF   (MAXSTREAMS_2,MAX_ATMOSWFS)

!  Error tracing variables

      INTEGER     :: INFO
      CHARACTER*3 :: CI

!  Other local help variables 

      INTEGER    :: I, I1, K, K1, Q, M
      INTEGER    :: NS, N, N1, NAF, NAL, NAL1, AA, C0, CM, CP
      REAL(fpk)  :: SHOM, L_HOM1, L_HOM2
!mick eff 3/22/2017
      INTEGER    :: NS1, NS2, CM1, CM2, CP1, CP2
      REAL(fpk)  :: HELP1 (NSTREAMS)

!  initialise Exception handling

      STATUS = LIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  Wf count

      M = FOURIER
      Q = 0
!mick eff 3/22/2017
      NS1 = NSTREAMS + 1 ; NS2 = NSTREAMS_2

!  Set up linearized BVP, column vector  B' where AX = B'
!  ======================================================

!  Bulk: Compute the main column B' where AX = B'

      CALL LC_BVPTEL_COLUMN_SETUP &
           ( DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM, DO_LAYER_SCATTERING,  & ! Input, Flags and order
             TAYLOR_ORDER, NLAYERS, NSTREAMS, NSTREAMS_2, N_COLUMNWFS,        & ! Input, Numbers
             FOURIER, IBEAM, ACTIVE_LAYERS, NLAYERS_TEL, N_BVTELMATRIX_SIZE,  & ! Input, Numbers for Telescoping
             DELTAU_VERT, L_DELTAU_VERT, SOLARBEAM_CUTOFF,                    & ! Input, optical and control
             INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,                     & ! Input, Beam Quantities
             T_DELT_EIGEN, T_DELT_DISORDS, XPOS, XNEG, WLOWER, LCON, MCON,    & ! Input, Homogeneous + WLOWER
             SURFACE_FACTOR, QUAD_STRMWTS, BRDF_F, DIRECT_BEAM, DELTAU_SLANT, & ! Input, Surface + Direct-beam
             GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,    & ! Input, Greens Function
             LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR,            & ! Input, Linearized Beam Quantities
             L_KEIGEN, L_T_DELT_EIGEN, L_XPOS, L_XNEG,                        & ! Input, Linearized Homogeneous
             L_T_DELT_DISORDS, L_ATERM_SAVE, L_BTERM_SAVE,                    & ! Input, Linearized Greens Function
             L_WUPPER, L_WLOWER, COLTEL2_WF, SCOL2_WF )                         ! Output

!  debug
!       if ( fourier .eq. 3 .and. ibeam.eq.4 ) then
!          do n = 13, N_BVTELMATRIX_SIZE
!             write(*,*)n,coltel2_wf(n,1),coltel2_wf(n,2)
!          enddo
!       endif

!  Solve linearized BVP, several active layers
!  ===========================================

      IF ( NLAYERS_TEL .GT. 1 ) THEN

!  BV solution for linearized integration constants
!    ( call to LAPACK solver routine for back substitution )

        CALL DGBTRS  ( 'n', N_BVTELMATRIX_SIZE, N_SUB, N_SUP, N_COLUMNWFS, &
                       BANDTELMAT2, MAXBANDTOTAL, IPIVOTTEL, COLTEL2_WF, MAXTOTAL, INFO )

!  Exception handling

        IF ( INFO .LT. 0 ) THEN
         WRITE(CI, '(I3)' ) INFO
         MESSAGE = 'argument i illegal value, for i = '//CI
         TRACE   = 'DGBTRS call in L_BVPTEL_BEAMSOLUTION_MASTER'
         STATUS  = LIDORT_SERIOUS
         RETURN
        ENDIF

!  Set linearized integration constants, active layers

!mick eff 3/22/2017
        C0 = - NSTREAMS_2
        DO NS = 1, NLAYERS_TEL
          N = ACTIVE_LAYERS(NS)
          C0 = C0 + NSTREAMS_2
          !DO I = 1, NSTREAMS
          !  CM = C0 + I
          !  CP = CM + NSTREAMS
          !  DO Q = 1, N_COLUMNWFS
          !    NCON(I,N,Q) = COLTEL2_WF(CM,Q)
          !    PCON(I,N,Q) = COLTEL2_WF(CP,Q)
          !  ENDDO
          !ENDDO
          CM1 = C0 + 1         ; CM2 = C0 + NSTREAMS
          CP1 = CM1 + NSTREAMS ; CP2 = CM2 + NSTREAMS
          DO Q = 1, N_COLUMNWFS
            NCON(1:NSTREAMS,N,Q) = COLTEL2_WF(CM1:CM2,Q)
            PCON(1:NSTREAMS,N,Q) = COLTEL2_WF(CP1:CP2,Q)
          ENDDO
        ENDDO

!  Solve linearized BVP: Single Layer only
!  =======================================

      ELSE IF ( NLAYERS_TEL .EQ. 1 ) THEN

!  LAPACK substitution (DGETRS) using RHS column vector SCOL2_WF

        CALL DGETRS ( 'N', NSTREAMS_2, N_COLUMNWFS, SMAT2, MAXSTREAMS_2, SIPIVOT, &
                       SCOL2_WF, MAXSTREAMS_2, INFO )

!  Exception handling

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGETRS call in LC_BVPTEL_BEAMSOLUTION_MASTER'
          STATUS  = LIDORT_SERIOUS
          RETURN
        ENDIF

!  Set linearized integration constants for active layer

!mick eff 3/22/2017
        N = ACTIVE_LAYERS(1)
        !DO K = 1, NSTREAMS
        !  K1 = K + NSTREAMS 
        !  DO Q = 1, N_COLUMNWFS
        !    NCON(K,N,Q) = SCOL2_WF(K,Q)
        !    PCON(K,N,Q) = SCOL2_WF(K1,Q)
        !  ENDDO
        !ENDDO
        DO Q = 1, N_COLUMNWFS
          NCON(1:NSTREAMS,N,Q) = SCOL2_WF(1:NSTREAMS,Q)
          PCON(1:NSTREAMS,N,Q) = SCOL2_WF(NS1:NS2,Q)
        ENDDO

!  end clause for backsubstitution

      ENDIF

!  Associated quantities for active layers
!  ---------------------------------------

!mick eff 3/22/2017
      DO NS = 1, NLAYERS_TEL
        N = ACTIVE_LAYERS(NS)
        !DO I = 1, NSTREAMS_2
        !  DO K = 1, NSTREAMS
        !    DO Q = 1, N_COLUMNWFS
        !      NCON_XVEC(I,K,N,Q) = NCON(K,N,Q) * XPOS(I,K,N)
        !      PCON_XVEC(I,K,N,Q) = PCON(K,N,Q) * XNEG(I,K,N)
        !    ENDDO
        !  ENDDO
        !ENDDO
        DO Q = 1, N_COLUMNWFS
          DO K = 1, NSTREAMS
            NCON_XVEC(1:NSTREAMS_2,K,N,Q) = NCON(K,N,Q) * XPOS(1:NSTREAMS_2,K,N)
            PCON_XVEC(1:NSTREAMS_2,K,N,Q) = PCON(K,N,Q) * XNEG(1:NSTREAMS_2,K,N)
          ENDDO
        ENDDO
      ENDDO

!  Set linearized integration constants for non-active layers
!  ==========================================================

!  Now we propagate the results upwards and downwards through the
!  appropriate non-active layers where there is no scattering.

!  Transmittance layers ABOVE active layer(s)
!  -----------------------------------------

!   --NCON values are zero (no downwelling radiation)
!   --PCON values propagated upwards from top of first active layer

!  layer immediately above first active layer
!   --- Require linearized solutions at top of first active layer
!   --- Additional linearizations required, first layer is always active

!mick eff 3/22/2017
      NAF = ACTIVE_LAYERS(1)
      IF ( NAF .GT. 1 ) THEN
        N1 = NAF - 1
        !DO I = 1, NSTREAMS
        !  I1 = I + NSTREAMS
        !  DO Q = 1, N_COLUMNWFS
        !    SHOM = ZERO
        !    DO AA = 1, NSTREAMS
        !      L_HOM1 = NCON_XVEC(I1,AA,NAF,Q) + LCON(AA,NAF)*L_XPOS(I1,AA,NAF,Q)
        !      L_HOM2 =   T_DELT_EIGEN(AA,NAF) * ( PCON_XVEC(I1,AA,NAF,Q) + MCON(AA,NAF)*L_XNEG(I1,AA,NAF,Q) ) &
        !             + L_T_DELT_EIGEN(AA,NAF,Q) * MCON_XVEC(I1,AA,NAF)
        !      SHOM = SHOM + L_HOM1 + L_HOM2
        !    ENDDO
        !    PCON(I,N1,Q) = L_WUPPER(I1,NAF,Q) + SHOM
        !    NCON(I,N1,Q) = ZERO
        !  ENDDO
        !ENDDO
        DO Q = 1, N_COLUMNWFS
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            HELP1 = PCON_XVEC(I1,1:NSTREAMS,NAF,Q) + MCON(1:NSTREAMS,NAF)*L_XNEG(I1,1:NSTREAMS,NAF,Q)
            SHOM  = SUM ( NCON_XVEC(I1,1:NSTREAMS,NAF,Q) ) &
                  + DOT_PRODUCT ( LCON(1:NSTREAMS,NAF), L_XPOS(I1,1:NSTREAMS,NAF,Q) ) &
                  + DOT_PRODUCT (   T_DELT_EIGEN(1:NSTREAMS,NAF), HELP1 ) &
                  + DOT_PRODUCT ( L_T_DELT_EIGEN(1:NSTREAMS,NAF,Q), MCON_XVEC(I1,1:NSTREAMS,NAF) )
            PCON(I,N1,Q) = L_WUPPER(I1,NAF,Q) + SHOM
            NCON(I,N1,Q) = ZERO
          ENDDO
        ENDDO
      ENDIF

!  For remaining non-active atmospheric layers to TOA, propagate upwards.
!   Additional linearizations if you are passing through the varying layer.

!mick eff 3/22/2017
      DO N = NAF-2, 1, -1
        N1 = N + 1
        !DO I = 1, NSTREAMS
        !  DO Q = 1, N_COLUMNWFS
        !    NCON(I,N,Q) = ZERO
        !    PCON(I,N,Q) = T_DELT_DISORDS(I,N1)   * PCON(I,N1,Q) &
        !              + L_T_DELT_DISORDS(I,N1,Q) * MCON(I,N1)
        !  ENDDO
        !ENDDO
        DO Q = 1, N_COLUMNWFS
          NCON(1:NSTREAMS,N,Q) = ZERO
          PCON(1:NSTREAMS,N,Q) =   T_DELT_DISORDS(1:NSTREAMS,N1)   * PCON(1:NSTREAMS,N1,Q) &
                               + L_T_DELT_DISORDS(1:NSTREAMS,N1,Q) * MCON(1:NSTREAMS,N1)
        ENDDO
      ENDDO     

!  Transmittance layers below active layer(s)
!  -----------------------------------------

!       ** Only do this if active scattering is above (not adjacent to) the surface layer

!   -- NCON values are  propagated downwards from bottom of last active layer
!   -- PCON values also propagated downwards, BUT only present if surface condition
!  1.   Require linearized solutions at bottom of last active layer
!  2.   Set values for layer immediately below last active layer
!  3.   Remaining layers to bottom, just propagate using discrete-ordinate transmittances

      NAL = ACTIVE_LAYERS(NLAYERS_TEL) ; NAL1 = NAL + 1
      IF ( NAL .LT. NLAYERS ) THEN

!  N-constants, always required

!mick eff 3/22/2017
         !DO I = 1, NSTREAMS
         !  DO Q = 1, N_COLUMNWFS
         !    SHOM = ZERO
         !    DO AA = 1, NSTREAMS
         !      L_HOM2 = PCON_XVEC(I,AA,NAL,Q) + MCON(AA,NAL)*L_XNEG(I,AA,NAL,Q)
         !      L_HOM1 =   T_DELT_EIGEN(AA,NAL) * ( NCON_XVEC(I,AA,NAL,Q) + LCON(AA,NAL)*L_XPOS(I,AA,NAL,Q) ) &
         !             + L_T_DELT_EIGEN(AA,NAL,Q) * LCON_XVEC(I,AA,NAL)
         !      SHOM = SHOM + L_HOM1 + L_HOM2
         !    ENDDO
         !    NCON(I,NAL1,Q) = L_WLOWER(I,NAL,Q) + SHOM
         !  ENDDO
         !ENDDO
         DO Q = 1, N_COLUMNWFS
           DO I = 1, NSTREAMS
             HELP1 = NCON_XVEC(I,1:NSTREAMS,NAL,Q) + LCON(1:NSTREAMS,NAL)*L_XPOS(I,1:NSTREAMS,NAL,Q)
             SHOM  = SUM ( PCON_XVEC(I,1:NSTREAMS,NAL,Q) ) &
                   + DOT_PRODUCT ( MCON(1:NSTREAMS,NAL), L_XNEG(I,1:NSTREAMS,NAL,Q) ) &
                   + DOT_PRODUCT (   T_DELT_EIGEN(1:NSTREAMS,NAL), HELP1 ) &
                   + DOT_PRODUCT ( L_T_DELT_EIGEN(1:NSTREAMS,NAL,Q), LCON_XVEC(I,1:NSTREAMS,NAL) )
             NCON(I,NAL1,Q) = L_WLOWER(I,NAL,Q) + SHOM
           ENDDO
         ENDDO

!mick eff 3/22/2017
         DO N = NAL + 2, NLAYERS
           N1 = N - 1
           !DO I = 1, NSTREAMS
           !  DO Q = 1, N_COLUMNWFS
           !    NCON(I,N,Q) =   T_DELT_DISORDS(I,N1)   * NCON(I,N1,Q) &
           !                + L_T_DELT_DISORDS(I,N1,Q) * LCON(I,N1)
           !  ENDDO
           !ENDDO
           DO Q = 1, N_COLUMNWFS
             NCON(1:NSTREAMS,N,Q) =   T_DELT_DISORDS(1:NSTREAMS,N1)   * NCON(1:NSTREAMS,N1,Q) &
                                  + L_T_DELT_DISORDS(1:NSTREAMS,N1,Q) * LCON(1:NSTREAMS,N1)
           ENDDO
         ENDDO

!  P-Constants need to be determined if there is a surface condition. Otherwise zero.

         IF ( DO_INCLUDE_SURFACE ) THEN
!mick eff 3/22/2017
           !DO I = 1, NSTREAMS
           !  I1 = I + NSTREAMS
           !  DO Q = 1, N_COLUMNWFS
           !    SHOM = ZERO
           !    DO AA = 1, NSTREAMS
           !      L_HOM2 = PCON_XVEC(I1,AA,NAL,Q) + MCON(AA,NAL)*L_XNEG(I1,AA,NAL,Q)
           !      L_HOM1 = T_DELT_EIGEN(AA,NAL) * &
           !          ( NCON_XVEC(I1,AA,NAL,Q) + LCON(AA,NAL)*L_XPOS(I1,AA,NAL,Q) ) + &
           !           L_T_DELT_EIGEN(AA,NAL,Q) * LCON_XVEC(I1,AA,NAL)
           !      SHOM = SHOM + L_HOM1 + L_HOM2
           !    ENDDO
           !    PCON(I,NAL1,Q) = L_WLOWER(I1,NAL,Q) + SHOM
           !    PCON(I,NAL1,Q) = ( PCON(I,NAL1,Q) - L_T_DELT_DISORDS(I,NAL1,Q) * MCON(I,NAL1) ) / T_DELT_DISORDS(I,NAL1)
           !  ENDDO
           !ENDDO
           DO Q = 1, N_COLUMNWFS
             DO I = 1, NSTREAMS
               I1 = I + NSTREAMS
               HELP1 = NCON_XVEC(I1,1:NSTREAMS,NAL,Q) + LCON(1:NSTREAMS,NAL)*L_XPOS(I1,1:NSTREAMS,NAL,Q)
               SHOM  = SUM ( PCON_XVEC(I1,1:NSTREAMS,NAL,Q) ) &
                     + DOT_PRODUCT ( MCON(1:NSTREAMS,NAL), L_XNEG(I1,1:NSTREAMS,NAL,Q) ) &
                     + DOT_PRODUCT (   T_DELT_EIGEN(1:NSTREAMS,NAL), HELP1 ) &
                     + DOT_PRODUCT ( L_T_DELT_EIGEN(1:NSTREAMS,NAL,Q), LCON_XVEC(I1,1:NSTREAMS,NAL) )
               PCON(I,NAL1,Q) = L_WLOWER(I1,NAL,Q) + SHOM
               PCON(I,NAL1,Q) = ( PCON(I,NAL1,Q) - L_T_DELT_DISORDS(I,NAL1,Q) * MCON(I,NAL1) ) / T_DELT_DISORDS(I,NAL1)
             ENDDO
           ENDDO

!mick eff 3/22/2017
           !DO N = NAL + 2, NLAYERS
           !  N1 = N - 1
           !  DO I = 1, NSTREAMS
           !    DO Q = 1, N_COLUMNWFS
           !      PCON(I,N,Q) = ( PCON(I,N1,Q) - L_T_DELT_DISORDS(I,N,Q) * MCON(I,N) ) / T_DELT_DISORDS(I,N) 
           !    ENDDO
           !  ENDDO
           !ENDDO
           DO Q = 1, N_COLUMNWFS
             DO N = NAL + 2, NLAYERS
               N1 = N - 1
               PCON(1:NSTREAMS,N,Q) = ( PCON(1:NSTREAMS,N1,Q) - L_T_DELT_DISORDS(1:NSTREAMS,N,Q) * MCON(1:NSTREAMS,N) ) &
                                      / T_DELT_DISORDS(1:NSTREAMS,N) 
             ENDDO
           ENDDO
         ELSE
!mick eff 3/22/2017
           !DO N = NAL + 1, NLAYERS
           !  DO Q = 1, N_COLUMNWFS
           !    PCON(1:NSTREAMS,N,Q) = ZERO
           !  ENDDO
           !ENDDO
           PCON(1:NSTREAMS,NAL+1:NLAYERS,1:N_COLUMNWFS) = ZERO
         ENDIF

!  End clause for non-active layers below telescoped problem

      ENDIF

!  Associated quantities for inactive layers
!  -----------------------------------------

!  Atmosphere layers with no scattering

!mick eff 3/22/2017
      !DO N = 1, NLAYERS
      !  IF ( N.LT.NAF .OR. N.GT.NAL ) THEN
      !    DO I = 1, NSTREAMS_2
      !      DO AA = 1, NSTREAMS
      !       DO Q = 1, N_COLUMNWFS
      !        NCON_XVEC(I,AA,N,Q) = NCON(AA,N,Q) * XPOS(I,AA,N)
      !        PCON_XVEC(I,AA,N,Q) = PCON(AA,N,Q) * XNEG(I,AA,N)
      !       ENDDO
      !      ENDDO
      !    ENDDO
      !  ENDIF
      !ENDDO
      DO Q = 1, N_COLUMNWFS
        DO N = 1, NLAYERS
          IF ( N.LT.NAF .OR. N.GT.NAL ) THEN
            DO AA = 1, NSTREAMS
              NCON_XVEC(1:NSTREAMS_2,AA,N,Q) = NCON(AA,N,Q) * XPOS(1:NSTREAMS_2,AA,N)
              PCON_XVEC(1:NSTREAMS_2,AA,N,Q) = PCON(AA,N,Q) * XNEG(1:NSTREAMS_2,AA,N)
            ENDDO
          ENDIF
        ENDDO
      ENDDO

!  Debug

!      if ( fourier.eq.3.and.ibeam.eq.2 ) then
!         do n = 20, nlayers
!            write(*,*)'telescp',n, NCON(2,N,1), PCON(2,N,1),LCON(3,N),MCON(2,n)
!         enddo
!      endif

!  finish

      RETURN
END SUBROUTINE LC_BVPTEL_SOLUTION_MASTER

!

SUBROUTINE LC_BVPTEL_SURFACE_SETUP                                             &
          ( DO_INCLUDE_SURFACE, NSTREAMS, NLAYERS, N, MBOUNDARY,               & ! Input
            FOURIER, N_WEIGHTFUNCS, SURFACE_FACTOR, QUAD_STRMWTS, BRDF_F,      & ! Input
            T_DELT_EIGEN, L_T_DELT_EIGEN, T_DELT_DISORDS, L_T_DELT_DISORDS,    & ! Input
            XPOS, XNEG, WLOWER, L_XPOS, L_XNEG, L_WLOWER,                      & ! Input
            R2_L_BEAM, R2_L_HOMP, R2_L_HOMM, CUMTRANS, L_CUMTRANS )              ! Output 

!  Linearized surface reflectance terms, Telescoping
!    Version 3.8. Major extension to BRDFs with TELESCOPING. Implemented May 2016.

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAXSTREAMS, MAXSTREAMS_2, MAXMOMENTS, MAXLAYERS, &
                                MAX_ATMOSWFS, ZERO, ONE

      IMPLICIT NONE

!  input arguments
!  ---------------

!  Surface (BRDF) flag

      LOGICAL  , intent(in)  ::  DO_INCLUDE_SURFACE

!  Number of streams and layer numbers

      INTEGER  , intent(in)  ::  NSTREAMS, NLAYERS, N

!  Number of weighting functions

      INTEGER  , intent(in)  ::  N_WEIGHTFUNCS

!  Fourier component

      INTEGER  , intent(in)  ::  FOURIER

!  Flag for type of boundary condition

      LOGICAL  , intent(in)  ::  MBOUNDARY

!  surface factor

      REAL(fpk), intent(in)  :: SURFACE_FACTOR

!  Fourier components of BRDF, in the following order
!    incident quadrature streams, reflected quadrature streams

      REAL(fpk), intent(in)  :: BRDF_F ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTREAMS )

!  Quadrature input

      REAL(fpk), intent(in)  :: QUAD_STRMWTS ( MAXSTREAMS )

!  transmittance factors discrete ordinate streams, and Linearizations

      REAL(fpk), intent(in)  :: T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: L_T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Eigenvector solutions

      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  General beam solutions at the Lower boundary

      REAL(fpk), intent(in)  :: WLOWER(MAXSTREAMS_2,MAXLAYERS)

!  transmittance factors for +/- eigenvalues

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)

!  Linearized Eigensolutions

      REAL(fpk), intent(in)  :: L_XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized General beam solutions at the Lower boundary

      REAL(fpk), intent(in)  :: L_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized transmittances, homogeneous solutions

      REAL(fpk), intent(in)  :: L_T_DELT_EIGEN ( MAXSTREAMS, MAXLAYERS, MAX_ATMOSWFS )

!  Output arguments
!  ----------------

      REAL(fpk), intent(out) :: R2_L_BEAM(MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: R2_L_HOMP(MAXSTREAMS,MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: R2_L_HOMM(MAXSTREAMS,MAXSTREAMS,MAX_ATMOSWFS)

!  cumulative tranmsittance (and linearization) from bottom of lowest active layer to surface
!     Calculated here, but  could be done earlier and passed in

      REAL(fpk), intent(out) :: CUMTRANS(MAXSTREAMS)
      REAL(fpk), intent(out) :: L_CUMTRANS(MAXSTREAMS,MAX_ATMOSWFS)

!  Local variables
!  ---------------

      REAL(fpk)  :: QCUMTRANS(MAXSTREAMS)
      REAL(fpk)  :: L_QCUMTRANS(MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk)  :: PS_W(MAXSTREAMS), PV_W(MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk)  :: HS_P(MAXSTREAMS,MAXSTREAMS), HV_P(MAXSTREAMS,MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk)  :: HS_M(MAXSTREAMS,MAXSTREAMS), HV_M(MAXSTREAMS,MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk)  :: BEAM, L_BEAM, HOMP, HOMM, L_HOMP, L_HOMM
      INTEGER    :: AA, J, Q, I, M, N1
!mick eff 3/22/2017
      !REAL(fpk)  :: HSP, HSM, L_HSP, L_HSM
      REAL(fpk)  :: HSP(NSTREAMS), HSM(NSTREAMS), L_HSP(NSTREAMS), L_HSM(NSTREAMS)

!  Initial section
!  ---------------

!  Always zero the result to start

!mick eff 3/22/2017
      !R2_L_BEAM = ZERO
      !R2_L_HOMP = ZERO
      !R2_L_HOMM = ZERO

      !CUMTRANS   = ONE
      !L_CUMTRANS = ZERO

      R2_L_BEAM(1:NSTREAMS,1:N_WEIGHTFUNCS) = ZERO
      R2_L_HOMP(1:NSTREAMS,1:NSTREAMS,1:N_WEIGHTFUNCS) = ZERO
      R2_L_HOMM(1:NSTREAMS,1:NSTREAMS,1:N_WEIGHTFUNCS) = ZERO

      CUMTRANS(1:NSTREAMS) = ONE
      L_CUMTRANS(1:NSTREAMS,1:N_WEIGHTFUNCS) = ZERO

!  Proxy Fourier component

      M = FOURIER

!  Return if no albedo

      IF ( .NOT. DO_INCLUDE_SURFACE ) RETURN

!  Cumulative transmittance and its linearization

!mick eff 3/22/2017 - initializing turned off (already done above) 
      !CUMTRANS(1:NSTREAMS) = ONE ; L_CUMTRANS(1:NSTREAMS,:) = ZERO
      DO N1 = NLAYERS, N+1, -1
         CUMTRANS(1:NSTREAMS) = CUMTRANS(1:NSTREAMS) * T_DELT_DISORDS(1:NSTREAMS,N1)
         DO Q = 1, N_WEIGHTFUNCS
           L_CUMTRANS(1:NSTREAMS,Q) = L_CUMTRANS(1:NSTREAMS,Q) * T_DELT_DISORDS  (1:NSTREAMS,N1) &
                                      + CUMTRANS(1:NSTREAMS)   * L_T_DELT_DISORDS(1:NSTREAMS,N1,Q)
         ENDDO
      ENDDO

!  Stored

      QCUMTRANS(1:NSTREAMS) = CUMTRANS(1:NSTREAMS) * QUAD_STRMWTS(1:NSTREAMS)
      DO Q = 1, N_WEIGHTFUNCS
         L_QCUMTRANS(1:NSTREAMS,Q) = L_CUMTRANS(1:NSTREAMS,Q) * QUAD_STRMWTS(1:NSTREAMS)
      ENDDO

!  Particular integral

!mick eff 3/22/2017
      !DO J = 1, NSTREAMS
      !  PS_W(J) = WLOWER(J,N) * QCUMTRANS(J)
      !  DO Q = 1, N_WEIGHTFUNCS
      !    PV_W(J,Q) = L_WLOWER(J,N,Q) * QCUMTRANS(J) + WLOWER(J,N) * L_QCUMTRANS(J,Q)
      !  ENDDO
      !ENDDO
      PS_W(1:NSTREAMS) = WLOWER(1:NSTREAMS,N) * QCUMTRANS(1:NSTREAMS)
      DO Q = 1, N_WEIGHTFUNCS
        PV_W(1:NSTREAMS,Q) = L_WLOWER(1:NSTREAMS,N,Q) *   QCUMTRANS(1:NSTREAMS) &
                           +   WLOWER(1:NSTREAMS,N)   * L_QCUMTRANS(1:NSTREAMS,Q)
      ENDDO

!  Modified boundary condition: homogeneous parts

!mick eff 3/22/2017
      IF ( MBOUNDARY ) THEN
        !DO J = 1, NSTREAMS
        !  DO AA = 1, NSTREAMS
        !    HSP = XPOS(J,AA,N) * T_DELT_EIGEN(AA,N)
        !    HSM = XNEG(J,AA,N)
        !    HS_P(J,AA) = HSP * QCUMTRANS(J)
        !    HS_M(J,AA) = HSM * QCUMTRANS(J)
        !    DO Q = 1, N_WEIGHTFUNCS
        !      L_HSP = L_XPOS(J,AA,N,Q) * T_DELT_EIGEN(AA,N) + XPOS(J,AA,N) * L_T_DELT_EIGEN(AA,N,Q)
        !      L_HSM = L_XNEG(J,AA,N,Q)
        !      HV_P(J,AA,Q) = L_HSP * QCUMTRANS(J) + HSP * L_QCUMTRANS(J,Q)
        !      HV_M(J,AA,Q) = L_HSM * QCUMTRANS(J) + HSM * L_QCUMTRANS(J,Q)
        !    ENDDO
        !  ENDDO
        !ENDDO
        DO AA = 1, NSTREAMS
          HSP = XPOS(1:NSTREAMS,AA,N) * T_DELT_EIGEN(AA,N)
          HSM = XNEG(1:NSTREAMS,AA,N)
          HS_P(1:NSTREAMS,AA) = HSP * QCUMTRANS(1:NSTREAMS)
          HS_M(1:NSTREAMS,AA) = HSM * QCUMTRANS(1:NSTREAMS)
          DO Q = 1, N_WEIGHTFUNCS
            L_HSP = L_XPOS(1:NSTREAMS,AA,N,Q) * T_DELT_EIGEN(AA,N) + XPOS(1:NSTREAMS,AA,N) * L_T_DELT_EIGEN(AA,N,Q)
            L_HSM = L_XNEG(1:NSTREAMS,AA,N,Q)
            HV_P(1:NSTREAMS,AA,Q) = L_HSP * QCUMTRANS(1:NSTREAMS) + HSP * L_QCUMTRANS(1:NSTREAMS,Q)
            HV_M(1:NSTREAMS,AA,Q) = L_HSM * QCUMTRANS(1:NSTREAMS) + HSM * L_QCUMTRANS(1:NSTREAMS,Q)
          ENDDO
        ENDDO
      ENDIF

!  Integrated Downward reflection (Calculation, Bidirectional case)
!     homogeneous and particular solutions.
!     @@@ Rob Fix 2/3/11,  Reverse J,I ---> I,J (J is incident)

      DO I = 1, NSTREAMS
        BEAM = DOT_PRODUCT ( PS_W(1:NSTREAMS),BRDF_F(M,I,1:NSTREAMS) )
        DO Q = 1, N_WEIGHTFUNCS
          L_BEAM = DOT_PRODUCT ( PV_W(1:NSTREAMS,Q),BRDF_F(M,I,1:NSTREAMS) )
          R2_L_BEAM(I,Q) = SURFACE_FACTOR * ( L_BEAM * CUMTRANS(I) + BEAM * L_CUMTRANS(I,Q) )
        ENDDO
        IF ( MBOUNDARY ) THEN
          DO AA = 1, NSTREAMS
            HOMP = DOT_PRODUCT ( HS_P(1:NSTREAMS,AA),BRDF_F(M,I,1:NSTREAMS) )
            HOMM = DOT_PRODUCT ( HS_M(1:NSTREAMS,AA),BRDF_F(M,I,1:NSTREAMS) )
            DO Q = 1, N_WEIGHTFUNCS
              L_HOMP = DOT_PRODUCT ( HV_P(1:NSTREAMS,AA,Q),BRDF_F(M,I,1:NSTREAMS) )
              L_HOMM = DOT_PRODUCT ( HV_M(1:NSTREAMS,AA,Q),BRDF_F(M,I,1:NSTREAMS) )
              R2_L_HOMP(I,AA,Q) = SURFACE_FACTOR * ( L_HOMP * CUMTRANS(I) + HOMP * L_CUMTRANS(I,Q) )
              R2_L_HOMM(I,AA,Q) = SURFACE_FACTOR * ( L_HOMM * CUMTRANS(I) + HOMM * L_CUMTRANS(I,Q) )
            ENDDO
          ENDDO
        ENDIF
      ENDDO

!  Finish

      RETURN
END SUBROUTINE LC_BVPTEL_SURFACE_SETUP

!

SUBROUTINE LC_BVPTEL_COLUMN_SETUP &
           ( DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM, DO_LAYER_SCATTERING,  & ! Input, Flags and order
             TAYLOR_ORDER, NLAYERS, NSTREAMS, NSTREAMS_2, N_COLUMNWFS,        & ! Input, Numbers
             FOURIER, IBEAM, ACTIVE_LAYERS, NLAYERS_TEL, N_BVTELMATRIX_SIZE,  & ! Input, Numbers for Telescoping
             DELTAU_VERT, L_DELTAU_VERT, SOLARBEAM_CUTOFF,                    & ! Input, optical and control
             INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,                     & ! Input, Beam Quantities
             T_DELT_EIGEN, T_DELT_DISORDS, XPOS, XNEG, WLOWER, LCON, MCON,    & ! Input, Homogeneous + WLOWER
             SURFACE_FACTOR, QUAD_STRMWTS, BRDF_F, DIRECT_BEAM, DELTAU_SLANT, & ! Input, Surface + Direct-beam
             GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,    & ! Input, Greens Function
             LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR,            & ! Input, Linearized Beam Quantities
             L_KEIGEN, L_T_DELT_EIGEN, L_XPOS, L_XNEG,                        & ! Input, Linearized Homogeneous
             L_T_DELT_DISORDS, L_ATERM_SAVE, L_BTERM_SAVE,                    & ! Input, Linearized Greens Function
             L_WUPPER, L_WLOWER, COLTEL2_WF, SCOL2_WF )                         ! Output

!  Column setup for the linearized telescoped BVP

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAXSTREAMS, MAXSTREAMS_2, MAXMOMENTS, MAXLAYERS, &
                                MAXBEAMS, MAX_ATMOSWFS, MAXTOTAL, ZERO

      IMPLICIT NONE

!  input arguments
!  ===============

!  Control and Optical
!  -------------------

!  Surface (BRDF) and directbeam flag

      LOGICAL  , intent(in)  ::  DO_INCLUDE_SURFACE
      LOGICAL  , intent(in)  ::  DO_INCLUDE_DIRECTBEAM

!  Local flags for the solution saving option

      LOGICAL  , intent(in)  ::  DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER  , intent(in)  :: TAYLOR_ORDER

!  Number of layers

      INTEGER  , intent(in)  ::  NLAYERS

!  Number of streams

      INTEGER  , intent(in)  ::  NSTREAMS, NSTREAMS_2

!  Linearization control: number of varying parameters (input)

      INTEGER  , intent(in)  ::  N_COLUMNWFS

!  Fourier number, beam index

      INTEGER  , intent(in)  ::  FOURIER, IBEAM

!  Number of telescoped layers

      INTEGER  , intent(in)  ::  NLAYERS_TEL

!  Active layers for telescoping

      INTEGER  , intent(in)  ::  ACTIVE_LAYERS ( MAXLAYERS )

!  Size of BVP matrix for telescoped 

      INTEGER  , intent(in)  ::  N_BVTELMATRIX_SIZE

!  Input optical properties after delta-M scaling. 
!     These are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: DELTAU_VERT ( MAXLAYERS )
      REAL(fpk), intent(in)  :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: SOLARBEAM_CUTOFF(MAXBEAMS)

!  Surface Quantities
!  ------------------

!  surface factor

      REAL(fpk), intent(in)  :: SURFACE_FACTOR

!  Fourier components of BRDF, in the following order (same all threads)
!    incident quadrature streams, reflected quadrature streams

      REAL(fpk), intent(in)  :: BRDF_F ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTREAMS )

!  Quadrature input

      REAL(fpk), intent(in)  :: QUAD_STRMWTS ( MAXSTREAMS )

!  Direct beam solutions

      REAL(fpk), intent(in)  :: DIRECT_BEAM ( MAXSTREAMS, MAXBEAMS )

!  Derived optical thickness inputs

      REAL(fpk), intent(in)  :: DELTAU_SLANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS )

!  Beam quantities
!  ---------------

!  Average-secants, Initial and average-secant transmittance factors.

      REAL(fpk), intent(in)  :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS )

!  Solution variables
!  ------------------

!  transmittance factors discrete ordinate streams, and Linearizations

      REAL(fpk), intent(in)  :: T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: L_T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Eigensolutions, eigenstream transmittances

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Solution constants of integration

      REAL(fpk), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)

!  particular solutions at Lower boundaries

      REAL(fpk), intent(in)  :: WLOWER ( MAXSTREAMS_2, MAXLAYERS )

!  Green functions
!  ---------------

!  Green function Multipliers for solution

      REAL(fpk), intent(in)  :: GFUNC_UP(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GFUNC_DN(MAXSTREAMS,MAXLAYERS)

!  Holding arrays for Multiplier coefficients

      REAL(fpk), intent(in)  :: GAMMA_M(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GAMMA_P(MAXSTREAMS,MAXLAYERS)

!  Green's function particular integral arrays

      REAL(fpk), intent(in)  :: ATERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: BTERM_SAVE(MAXSTREAMS,MAXLAYERS)

!  Linearized
!  ----------

!  Linearized Average-secants, Initial and average-secant transmittance factors.
!     LC_INITIAL_TRANS are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: LC_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LC_T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized Eigenvalues, Eigensolutions, eigenstream transmittances

      REAL(fpk), intent(in)  :: L_KEIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Saved quantities for the Green function solution
!     These are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: L_ATERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_BTERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  output arguments
!  ----------------

!  Linearized beam solutions at the Lower and Upper layer boundaries

      REAL(fpk), intent(inout) :: L_WUPPER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(inout) :: L_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Column vectors for solving linearized BCs

      REAL(fpk), intent(out) :: COLTEL2_WF (MAXTOTAL,    MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: SCOL2_WF   (MAXSTREAMS_2,MAX_ATMOSWFS)

!  local variables
!  ---------------

!  Output arguments from the Surface setup (reflectances and cumulative transmittances)
!  cumulative tranmsittance (and linearization) from bottom of lowest active layer to surface

      REAL(fpk) :: R2_L_BEAM(MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk) :: R2_L_HOMP(MAXSTREAMS,MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk) :: R2_L_HOMM(MAXSTREAMS,MAXSTREAMS,MAX_ATMOSWFS)

      REAL(fpk) :: CUMTRANS(MAXSTREAMS)
      REAL(fpk) :: L_CUMTRANS(MAXSTREAMS,MAX_ATMOSWFS)

!  help variables

      LOGICAL     :: MBOUNDARY
      INTEGER     :: Q,AA,N,N1,NS,I,I1,CM,C0,NAF,K
      REAL(fpk)   :: BEAM, FAC3, CPOS, CNEG, L_HOM, L_BEAM, L_HOMD, L_HOMU
!mick eff 3/22/2017
      REAL(fpk)   :: HELP1 (NSTREAMS), HELP2 (NSTREAMS)

!  Try this safety-first zeroing

!mick eff 3/22/2017
      !DO I = 1, NSTREAMS_2
      !  DO Q = 1, N_COLUMNWFS
      !    DO N = 1, NLAYERS
      !      L_WUPPER(I,N,Q) = ZERO
      !      L_WLOWER(I,N,Q) = ZERO
      !    ENDDO
      !  ENDDO
      !ENDDO
      L_WUPPER(1:NSTREAMS_2,1:NLAYERS,1:N_COLUMNWFS) = ZERO
      L_WLOWER(1:NSTREAMS_2,1:NLAYERS,1:N_COLUMNWFS) = ZERO

!  Get the linearized solutions for all active layers
!    Always need this, regardless of number of active layers

      DO NS = 1, NLAYERS_TEL
        N = ACTIVE_LAYERS(NS)
        CALL LC_BEAMSOLUTION_NEQK &
           ( DO_LAYER_SCATTERING, TAYLOR_ORDER,                      & ! Input, Flags and order
             NSTREAMS, NSTREAMS_2, N_COLUMNWFS, FOURIER, IBEAM, N,   & ! Input, Numbers
             DELTAU_VERT, L_DELTAU_VERT, SOLARBEAM_CUTOFF,           & ! Input, optical and control
             INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,            & ! Input, Beam Quantities
             LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR,   & ! Input, Beam Quantities (Linearized)
             T_DELT_EIGEN, XPOS, L_KEIGEN, L_T_DELT_EIGEN, L_XPOS,   & ! Input, Homogeneous solution stuff
             GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P,                   & ! Input, Greens Function
             ATERM_SAVE, BTERM_SAVE, L_ATERM_SAVE, L_BTERM_SAVE,     & ! Input, Greens Function
             L_WUPPER, L_WLOWER )                                      ! Output
      ENDDO

!  Go to special case for only 1 active layer

      IF ( NLAYERS_TEL .GT. 1 ) THEN

!  zero column vector

!mick eff 3/22/2017
        !COLTEL2_WF(1:N_BVTELMATRIX_SIZE,:) = ZERO
        COLTEL2_WF(1:N_BVTELMATRIX_SIZE,1:N_COLUMNWFS) = ZERO

!  top of first active layer, first boundary condition
!  ---------------------------------------------------

        NS = 1
        N = ACTIVE_LAYERS(NS)
        NAF = N

!  Require homogeneous and beam solution linearizations

!mick eff 3/22/2017
        !DO I = 1, NSTREAMS
        !  DO Q = 1, N_COLUMNWFS
        !    L_BEAM = - L_WUPPER(I,N,Q)
        !    L_HOM  = ZERO
        !    DO AA = 1, NSTREAMS
        !      CPOS = L_XPOS(I,AA,N,Q)
        !      CNEG = T_DELT_EIGEN(AA,N)   * L_XNEG(I,AA,N,Q) +  &
        !           L_T_DELT_EIGEN(AA,N,Q) *   XNEG(I,AA,N)
        !      L_HOM = L_HOM + LCON(AA,N) * CPOS + MCON(AA,N) * CNEG
        !    ENDDO
        !    COLTEL2_WF(I,Q) = L_BEAM - L_HOM
        !  ENDDO
        !ENDDO
        DO Q = 1, N_COLUMNWFS
          DO I = 1, NSTREAMS
            L_BEAM = - L_WUPPER(I,N,Q)
            HELP1  =   T_DELT_EIGEN(1:NSTREAMS,N)   * L_XNEG(I,1:NSTREAMS,N,Q) &
                   + L_T_DELT_EIGEN(1:NSTREAMS,N,Q) *   XNEG(I,1:NSTREAMS,N)
            L_HOM  = DOT_PRODUCT ( LCON(1:NSTREAMS,N),L_XPOS(I,1:NSTREAMS,N,Q) ) &
                   + DOT_PRODUCT ( MCON(1:NSTREAMS,N),HELP1 )
            COLTEL2_WF(I,Q) = L_BEAM - L_HOM
          ENDDO
        ENDDO

!  Intermediate boundaries between active layers
!  ---------------------------------------------

        DO NS = 1, NLAYERS_TEL - 1

!  Offsets

          N  = ACTIVE_LAYERS(NS)
          N1 = N + 1
          C0 = NS*NSTREAMS_2 - NSTREAMS

!  Get the linearized beam solution for the next layer N1

!mick eff 3/22/2017
          !DO I = 1, NSTREAMS_2
          !  CM = C0 + I
          !  DO Q = 1, N_COLUMNWFS
          !    L_BEAM = L_WUPPER(I,N1,Q) - L_WLOWER(I,N,Q)      
          !    L_HOMU = ZERO
          !    L_HOMD = ZERO
          !    DO AA = 1, NSTREAMS
          !      CPOS = L_XPOS(I,AA,N1,Q)
          !      CNEG =   T_DELT_EIGEN(AA,N1)   * L_XNEG(I,AA,N1,Q) + &
          !             L_T_DELT_EIGEN(AA,N1,Q) *   XNEG(I,AA,N1)
          !      L_HOMU = L_HOMU + LCON(AA,N1) * CPOS + MCON(AA,N1) * CNEG
          !      CNEG = L_XNEG(I,AA,N,Q)
          !      CPOS =   T_DELT_EIGEN(AA,N)    * L_XPOS(I,AA,N,Q)  +  &
          !             L_T_DELT_EIGEN(AA,N,Q)  *   XPOS(I,AA,N)
          !      L_HOMD = L_HOMD + LCON(AA,N)  * CPOS + MCON(AA,N)  * CNEG
          !    ENDDO
          !    L_HOM            = L_HOMU - L_HOMD
          !    COLTEL2_WF(CM,Q) = L_BEAM + L_HOM
          !   ENDDO
          !ENDDO
          DO Q = 1, N_COLUMNWFS
            DO I = 1, NSTREAMS_2
              CM = C0 + I
              L_BEAM = L_WUPPER(I,N1,Q) - L_WLOWER(I,N,Q)
              HELP1    =   T_DELT_EIGEN(1:NSTREAMS,N1)   * L_XNEG(I,1:NSTREAMS,N1,Q) &
                       + L_T_DELT_EIGEN(1:NSTREAMS,N1,Q) *   XNEG(I,1:NSTREAMS,N1)
              L_HOMU   = DOT_PRODUCT ( LCON(1:NSTREAMS,N1), L_XPOS(I,1:NSTREAMS,N1,Q) ) &
                       + DOT_PRODUCT ( MCON(1:NSTREAMS,N1), HELP1 )
              HELP2    =   T_DELT_EIGEN(1:NSTREAMS,N)   * L_XPOS(I,1:NSTREAMS,N,Q) &
                       + L_T_DELT_EIGEN(1:NSTREAMS,N,Q) *   XPOS(I,1:NSTREAMS,N)
              L_HOMD   = DOT_PRODUCT ( LCON(1:NSTREAMS,N), HELP2 ) &
                       + DOT_PRODUCT ( MCON(1:NSTREAMS,N), L_XNEG(I,1:NSTREAMS,N,Q) )
              L_HOM    = L_HOMU - L_HOMD
              COLTEL2_WF(CM,Q) = L_BEAM + L_HOM
            ENDDO
          ENDDO

!  End loop over intermediate active layer boundaries

        ENDDO

!  Final boundary, bottom of lowest active layer
!  ---------------------------------------------

        NS = NLAYERS_TEL
        N  = ACTIVE_LAYERS(NS)      
        C0 = (NS-1)*NSTREAMS_2 + NSTREAMS
        MBOUNDARY = .TRUE.

!  get the linearized downward-reflected term. New 5/24/16 Rob Fix

        CALL LC_BVPTEL_SURFACE_SETUP                                           &
          ( DO_INCLUDE_SURFACE, NSTREAMS, NLAYERS, N, MBOUNDARY,               & ! Input
            FOURIER, N_COLUMNWFS, SURFACE_FACTOR, QUAD_STRMWTS, BRDF_F,        & ! Input
            T_DELT_EIGEN, L_T_DELT_EIGEN, T_DELT_DISORDS, L_T_DELT_DISORDS,    & ! Input
            XPOS, XNEG, WLOWER, L_XPOS, L_XNEG, L_WLOWER,                      & ! Input
            R2_L_BEAM, R2_L_HOMP, R2_L_HOMM, CUMTRANS, L_CUMTRANS )              ! Output 

!  last active layer is varying, always

        IF ( DO_INCLUDE_SURFACE ) THEN
!mick eff 3/22/2017
          !DO I = 1, NSTREAMS
          !  CM = C0 + I
          !  I1 = I + NSTREAMS
          !  DO Q = 1, N_COLUMNWFS
          !    L_BEAM = L_WLOWER(I1,N,Q) - R2_L_BEAM(I,Q)
          !    L_HOM  = ZERO
          !    DO AA = 1, NSTREAMS
          !      CPOS =   T_DELT_EIGEN(AA,N)   * L_XPOS(I1,AA,N,Q) + &
          !             L_T_DELT_EIGEN(AA,N,Q) *   XPOS(I1,AA,N)
          !      CPOS =        CPOS       - R2_L_HOMP(I,AA,Q)
          !      CNEG = L_XNEG(I1,AA,N,Q) - R2_L_HOMM(I,AA,Q)
          !      L_HOM = L_HOM + LCON(AA,N)*CPOS + MCON(AA,N)*CNEG
          !    ENDDO
          !    COLTEL2_WF(CM,Q) = - L_BEAM - L_HOM
          !  ENDDO
          !ENDDO
          DO Q = 1, N_COLUMNWFS
            DO I = 1, NSTREAMS
              CM = C0 + I
              I1 = I + NSTREAMS
              L_BEAM = L_WLOWER(I1,N,Q) - R2_L_BEAM(I,Q)
              HELP1  =   T_DELT_EIGEN(1:NSTREAMS,N)   * L_XPOS(I1,1:NSTREAMS,N,Q) &
                     + L_T_DELT_EIGEN(1:NSTREAMS,N,Q) *   XPOS(I1,1:NSTREAMS,N) &
                     - R2_L_HOMP(I,1:NSTREAMS,Q)
              HELP2  = L_XNEG(I1,1:NSTREAMS,N,Q) - R2_L_HOMM(I,1:NSTREAMS,Q)
              L_HOM  = DOT_PRODUCT ( LCON(1:NSTREAMS,N), HELP1 ) &
                     + DOT_PRODUCT ( MCON(1:NSTREAMS,N), HELP2 )
              COLTEL2_WF(CM,Q) = - L_BEAM - L_HOM
            ENDDO
          ENDDO
        ELSE
!mick eff 3/22/2017
          !DO I = 1, NSTREAMS
          !  CM = C0 + I
          !  I1 = I + NSTREAMS
          !  DO Q = 1, N_COLUMNWFS
          !    L_BEAM = L_WLOWER(I1,N,Q)
          !    L_HOM  = ZERO
          !    DO AA = 1, NSTREAMS
          !      CPOS = T_DELT_EIGEN(AA,N)   * L_XPOS(I1,AA,N,Q) + &
          !           L_T_DELT_EIGEN(AA,N,Q) *   XPOS(I1,AA,N)
          !      CNEG = L_XNEG(I1,AA,N,Q)
          !      L_HOM = L_HOM + LCON(AA,N)*CPOS + MCON(AA,N)*CNEG
          !    ENDDO
          !    COLTEL2_WF(CM,Q) = - L_BEAM - L_HOM
          !  ENDDO
          !ENDDO
          DO Q = 1, N_COLUMNWFS
            DO I = 1, NSTREAMS
              CM = C0 + I
              I1 = I + NSTREAMS
              L_BEAM = L_WLOWER(I1,N,Q)
              HELP1  =   T_DELT_EIGEN(1:NSTREAMS,N)   * L_XPOS(I1,1:NSTREAMS,N,Q) &
                     + L_T_DELT_EIGEN(1:NSTREAMS,N,Q) *   XPOS(I1,1:NSTREAMS,N)
              HELP2  = L_XNEG(I1,1:NSTREAMS,N,Q)
              L_HOM  = DOT_PRODUCT ( LCON(1:NSTREAMS,N), HELP1 ) &
                     + DOT_PRODUCT ( MCON(1:NSTREAMS,N), HELP2 )
              COLTEL2_WF(CM,Q) = - L_BEAM - L_HOM
            ENDDO
          ENDDO
        ENDIF

!  Add direct beam solution. This is new code, R. Spurr 02/16/15. 5/24/16
!    --- If necessary, DB term is attenuated upwards from surface to lowest active-layer.

        IF ( DO_INCLUDE_SURFACE.and.DO_INCLUDE_DIRECTBEAM ) THEN
!mick eff 3/22/2017
          !DO I = 1, NSTREAMS
          !  CM = C0 + I
          !  BEAM = DIRECT_BEAM(I,IBEAM)
          !  DO Q = 1, N_COLUMNWFS
          !    L_BEAM = ZERO
          !    DO K = 1, NLAYERS
          !      FAC3 = - BEAM * DELTAU_SLANT(NLAYERS,K,IBEAM)
          !      L_BEAM = L_BEAM + L_DELTAU_VERT(Q,K) * FAC3
          !    ENDDO
          !    COLTEL2_WF(CM,Q) = COLTEL2_WF(CM,Q) + L_BEAM * CUMTRANS(I) + L_CUMTRANS(I,Q) * BEAM
          !  ENDDO
          !ENDDO
          DO Q = 1, N_COLUMNWFS
            DO I = 1, NSTREAMS
              CM = C0 + I
              BEAM = DIRECT_BEAM(I,IBEAM)
              L_BEAM = - BEAM * DOT_PRODUCT(L_DELTAU_VERT(Q,1:NLAYERS),DELTAU_SLANT(NLAYERS,1:NLAYERS,IBEAM))
              COLTEL2_WF(CM,Q) = COLTEL2_WF(CM,Q) + L_BEAM * CUMTRANS(I) + L_CUMTRANS(I,Q) * BEAM
            ENDDO
          ENDDO
        ENDIF

!  Continuation point for the single-active-layer case

      ELSE

!  Zero column vector

        SCOL2_WF = ZERO

!  Top of active layer
!  -------------------

        NS = 1
        N = ACTIVE_LAYERS(NS)

!  layer that is varying,
!    require homogeneous and beam solution linearizations

!mick eff 3/22/2017
        !DO I = 1, NSTREAMS
        !  DO Q = 1, N_COLUMNWFS
        !    L_BEAM = - L_WUPPER(I,N,Q)
        !    L_HOM  = ZERO
        !    DO AA = 1, NSTREAMS
        !      CPOS = L_XPOS(I,AA,N,Q)
        !      CNEG = T_DELT_EIGEN(AA,N)   * L_XNEG(I,AA,N,Q) + &
        !           L_T_DELT_EIGEN(AA,N,Q) *   XNEG(I,AA,N)
        !      L_HOM = L_HOM + LCON(AA,N) * CPOS + MCON(AA,N) * CNEG
        !    ENDDO
        !    SCOL2_WF(I,Q) = L_BEAM - L_HOM
        !  ENDDO
        !ENDDO
        DO Q = 1, N_COLUMNWFS
          DO I = 1, NSTREAMS
            L_BEAM = - L_WUPPER(I,N,Q)
            HELP1  = L_XPOS(I,1:NSTREAMS,N,Q)
            HELP2  =   T_DELT_EIGEN(1:NSTREAMS,N)   * L_XNEG(I,1:NSTREAMS,N,Q) &
                   + L_T_DELT_EIGEN(1:NSTREAMS,N,Q) *   XNEG(I,1:NSTREAMS,N)
            L_HOM  = DOT_PRODUCT ( LCON(1:NSTREAMS,N), HELP1 ) &
                   + DOT_PRODUCT ( MCON(1:NSTREAMS,N), HELP2 )
            SCOL2_WF(I,Q) = L_BEAM - L_HOM
          ENDDO
        ENDDO

!  Bottom of active layer
!  ----------------------

!  Active layer is varying layer

        MBOUNDARY = .true.

!  Get the linearized downward-reflected term. New 5/24/16 Rob Fix

        CALL LC_BVPTEL_SURFACE_SETUP                                           &
          ( DO_INCLUDE_SURFACE, NSTREAMS, NLAYERS, N, MBOUNDARY,               & ! Input
            FOURIER, N_COLUMNWFS, SURFACE_FACTOR, QUAD_STRMWTS, BRDF_F,        & ! Input
            T_DELT_EIGEN, L_T_DELT_EIGEN, T_DELT_DISORDS, L_T_DELT_DISORDS,    & ! Input
            XPOS, XNEG, WLOWER, L_XPOS, L_XNEG, L_WLOWER,                      & ! Input
            R2_L_BEAM, R2_L_HOMP, R2_L_HOMM, CUMTRANS, L_CUMTRANS )              ! Output 

!  Active layer is varying layer

        IF ( DO_INCLUDE_SURFACE ) THEN
!mick eff 3/22/2017
          !DO I = 1, NSTREAMS
          !  I1 = I + NSTREAMS
          !  DO Q = 1, N_COLUMNWFS
          !    L_BEAM = L_WLOWER(I1,N,Q) - R2_L_BEAM(I,Q)
          !    L_HOM  = ZERO
          !    DO AA = 1, NSTREAMS
          !      CPOS = T_DELT_EIGEN(AA,N)   * L_XPOS(I1,AA,N,Q) +  &
          !           L_T_DELT_EIGEN(AA,N,Q) *   XPOS(I1,AA,N)
          !      CPOS =        CPOS       - R2_L_HOMP(I,AA,Q)
          !      CNEG = L_XNEG(I1,AA,N,Q) - R2_L_HOMM(I,AA,Q)
          !      L_HOM = L_HOM + LCON(AA,N)*CPOS + MCON(AA,N)*CNEG
          !    ENDDO
          !    SCOL2_WF(I1,Q) = - L_BEAM - L_HOM
          !  ENDDO
          !ENDDO
          DO Q = 1, N_COLUMNWFS
            DO I = 1, NSTREAMS
              I1 = I + NSTREAMS
              L_BEAM = L_WLOWER(I1,N,Q) - R2_L_BEAM(I,Q)
              HELP1  =   T_DELT_EIGEN(1:NSTREAMS,N)   * L_XPOS(I1,1:NSTREAMS,N,Q) &
                     + L_T_DELT_EIGEN(1:NSTREAMS,N,Q) *   XPOS(I1,1:NSTREAMS,N) &
                     - R2_L_HOMP(I,1:NSTREAMS,Q)
              HELP2  = L_XNEG(I1,1:NSTREAMS,N,Q) - R2_L_HOMM(I,1:NSTREAMS,Q)
              L_HOM  = DOT_PRODUCT ( LCON(1:NSTREAMS,N), HELP1 ) &
                     + DOT_PRODUCT ( MCON(1:NSTREAMS,N), HELP2 )
              SCOL2_WF(I1,Q) = - L_BEAM - L_HOM
            ENDDO
          ENDDO
        ELSE
!mick eff 3/22/2017
          !DO I = 1, NSTREAMS
          !  I1 = I + NSTREAMS
          !  DO Q = 1, N_COLUMNWFS
          !    L_BEAM = L_WLOWER(I1,N,Q)
          !    L_HOM  = ZERO
          !    DO AA = 1, NSTREAMS
          !      CPOS = T_DELT_EIGEN(AA,N)   * L_XPOS(I1,AA,N,Q) +  &
          !           L_T_DELT_EIGEN(AA,N,Q) *   XPOS(I1,AA,N)
          !      CNEG = L_XNEG(I1,AA,N,Q)
          !      L_HOM = L_HOM + LCON(AA,N)*CPOS + MCON(AA,N)*CNEG
          !    ENDDO
          !    SCOL2_WF(I1,Q) = - L_BEAM - L_HOM
          !  ENDDO
          !ENDDO
          DO Q = 1, N_COLUMNWFS
            DO I = 1, NSTREAMS
              I1 = I + NSTREAMS
              L_BEAM = L_WLOWER(I1,N,Q)
              HELP1  =   T_DELT_EIGEN(1:NSTREAMS,N)   * L_XPOS(I1,1:NSTREAMS,N,Q) &
                     + L_T_DELT_EIGEN(1:NSTREAMS,N,Q) *   XPOS(I1,1:NSTREAMS,N)
              HELP2  = L_XNEG(I1,1:NSTREAMS,N,Q)
              L_HOM  = DOT_PRODUCT ( LCON(1:NSTREAMS,N), HELP1 ) &
                     + DOT_PRODUCT ( MCON(1:NSTREAMS,N), HELP2 )
              SCOL2_WF(I1,Q) = - L_BEAM - L_HOM
            ENDDO
          ENDDO
        ENDIF

!  Add direct beam solution. This is new code, R. Spurr 02/16/15. 5/24/16
!    --- If necessary, DB term is attenuated upwards from surface to lowest active-layer.

        IF ( DO_INCLUDE_SURFACE.and.DO_INCLUDE_DIRECTBEAM ) THEN
!mick eff 3/22/2017
          !DO I = 1, NSTREAMS
          !  I1 = I + NSTREAMS
          !  BEAM = DIRECT_BEAM(I,IBEAM)
          !  DO Q = 1, N_COLUMNWFS
          !    L_BEAM = ZERO
          !    DO K = 1, NLAYERS
          !      !FAC3 = - BEAM * DELTAU_SLANT(N,K,IBEAM)
          !      FAC3 = - BEAM * DELTAU_SLANT(NLAYERS,K,IBEAM)
          !      L_BEAM = L_BEAM + L_DELTAU_VERT(Q,K) * FAC3
          !    ENDDO
          !    SCOL2_WF(I1,Q) = SCOL2_WF(I1,Q) + L_BEAM * CUMTRANS(I) + L_CUMTRANS(I,Q) * BEAM
          !  ENDDO
          !ENDDO
          DO Q = 1, N_COLUMNWFS
            DO I = 1, NSTREAMS
              I1 = I + NSTREAMS
              BEAM = DIRECT_BEAM(I,IBEAM)
              L_BEAM = - BEAM * DOT_PRODUCT(L_DELTAU_VERT(Q,1:NLAYERS),DELTAU_SLANT(NLAYERS,1:NLAYERS,IBEAM))
              SCOL2_WF(I1,Q) = SCOL2_WF(I1,Q) + L_BEAM * CUMTRANS(I) + L_CUMTRANS(I,Q) * BEAM
            ENDDO
          ENDDO
        ENDIF

!  End layer clause

      ENDIF

!  finish

      RETURN
END SUBROUTINE LC_BVPTEL_COLUMN_SETUP

!  End

end module lidort_lc_bvproblem_m
