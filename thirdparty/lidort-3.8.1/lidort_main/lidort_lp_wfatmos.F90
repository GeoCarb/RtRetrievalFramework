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

! ##########################################################
! #                                                        #
! # Subroutines in this Module                             #
! #                                                        #
! #   -- Master routines for Profile Jacobians             #
! #                                                        #
! #            UPUSER_PROFILEWF                            #
! #            DNUSER_PROFILEWF                            #
! #            MIFLUX_PROFILEWF                            #
! #                                                        #
! #   -- Master routines for post-processed TOA/BOA fields #
! #                                                        #
! #            GET_LP_TOASOURCE                            #
! #            GET_LP_BOASOURCE                            #
! #                                                        #
! #   -- Master routines for azimuth Convergence           #
! #                                                        #
! #            LIDORT_LP_CONVERGE                          #
! #            LIDORT_LP_CONVERGE_OBSGEO                   #
! #                                                        #
! ##########################################################

!  Notes for Version 3.7 (Taylor series expansions)
!  ------------------------------------------------

!     Rob Fix  05/06/13  - Introduce TAYLOR_LIMIT parameters in module LIDORT_PARS
!     Rob  Fix 05/06/13  - Introduce limiting case scenarios (Taylor series)
!     Mick Fix 09/11/13  - Redefined ZETAs to straight sum and difference
!     Rob  Fix 10/09/13  - Small numbers analysis finalized using Taylor_Order parameter
!     Rob  Fix 01/05/14  - Use N_PPSTREAMS and PPSTREAM_MASK, to cover observation/lattice options

!  Upgrade, Version 3.8.1. June 2019
!    --- Changes to the GET_LP_BOASOURCE to deal with WLeaving

module lidort_lp_wfatmos_m

!  Parameter types

   USE LIDORT_PARS_m, only : fpk, TAYLOR_SMALL

!  Dependencies

   USE lidort_lp_PostProcessing_m

!  Everything public

public 

contains

SUBROUTINE UPUSER_PROFILEWF &
            ( DO_USER_STREAMS, DO_SOLAR_SOURCES, DO_MSMODE_LIDORT,                 & ! Input
              DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY, DO_INCLUDE_MVOUTPUT,    & ! Input
              NSTREAMS, NLAYERS, N_USER_LEVELS, VARIATION_INDEX, K_PARAMETERS,     & ! Input
              N_PPSTREAMS, PPSTREAM_MASK, DO_LAYER_SCATTERING, UTAU_LEVEL_MASK_UP, & ! Input
              PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,        & ! Input
              FOURIER, IBEAM, TAYLOR_ORDER, FLUX_MULTIPLIER, QUAD_STREAMS,         & ! Input
              PARTAU_VERT, DELTAU_VERT, L_DELTAU_VERT,                             & ! Input
              USER_SECANTS, SOLARBEAM_CUTOFF, INITIAL_TRANS, AVERAGE_SECANT,       & ! Input
              T_DELT_MUBAR,   T_UTDN_MUBAR,                      & ! Input
              T_DELT_USERM,   T_UTUP_USERM,                      & ! Input
              T_DELT_EIGEN,   T_UTUP_EIGEN,   T_UTDN_EIGEN,      & ! Input
              T_DELT_DISORDS, T_DISORDS_UTUP,                    & ! Input
              GAMMA_M, GAMMA_P, SIGMA_P, ATERM_SAVE, BTERM_SAVE, & ! Input
              XPOS, LCON_XVEC, MCON_XVEC, LCON, MCON,            & ! Input
              T_WUPPER, BOA_THTONLY_SOURCE,                      & ! Input
              U_XPOS, U_XNEG, U_WPOS, HMULT_1, HMULT_2, EMULT_UP,        & ! Input
              UT_HMULT_UU,  UT_HMULT_UD, UT_EMULT_UP,                    & ! Input
              PMULT_UU, PMULT_UD, UT_PMULT_UU, UT_PMULT_UD,              & ! Input
              UT_GMULT_UP, UT_GMULT_DN, CUMSOURCE_UP,                    & ! Input
              LP_INITIAL_TRANS, LP_AVERAGE_SECANT,                & ! Input
              LP_T_DELT_MUBAR,  LP_T_UTDN_MUBAR,                  & ! Input
              L_T_DELT_USERM,   L_T_UTUP_USERM,                   & ! Input
              L_T_DELT_EIGEN,   L_T_UTUP_EIGEN, L_T_UTDN_EIGEN,   & ! Input
              L_T_DELT_DISORDS, L_T_DISORDS_UTUP,                 & ! Input
              L_KEIGEN, L_ATERM_SAVE, L_BTERM_SAVE,               & ! Input
              L_U_XPOS, L_U_XNEG, LP_U_WPOS, NCON, PCON,                 & ! Input
              L_XPOS, L_XNEG, L_WUPPER, L_WLOWER, NCON_XVEC, PCON_XVEC,  & ! Input
              L_HMULT_1, L_HMULT_2, LP_EMULT_UP,                         & ! Input
              L_UT_HMULT_UU, L_UT_HMULT_UD, LP_UT_EMULT_UP,              & ! Input
              L_LAYER_TSUP_UP, L_LAYER_TSUP_UTUP, L_UT_T_PARTIC, L_T_WUPPER, & ! Input
              L_BOA_THTONLY_SOURCE, LP_BOA_MSSOURCE, LP_BOA_DBSOURCE,        & ! Input
              FLAGS_LP_GMULT, LP_UT_GMULT_UP, LP_UT_GMULT_DN,                & ! Output
              PROFILEWF_F, QUADPROFILEWF )                                     ! Output

!  Upwelling post-processed Profile Jacobian Fourier component

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAXSTREAMS, MAXSTREAMS_2, MAXMOMENTS, MAXLAYERS, MAX_PARTLAYERS, &
                                MAXBEAMS, MAX_USER_LEVELS, MAX_USER_STREAMS, MAX_ATMOSWFS,       &
                                MAX_DIRECTIONS, ZERO, UPIDX

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  Control flags

      LOGICAL  , intent(in)  :: DO_USER_STREAMS
      LOGICAL  , intent(in)  :: DO_MSMODE_LIDORT

!  Solar source term flag

      LOGICAL  , intent(in)  :: DO_SOLAR_SOURCES

!  Thermal emission flag

      LOGICAL  , intent(in)  :: DO_INCLUDE_THERMEMISS

!  Thermal transmittance-only flag

      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY

!  Local control flags

      LOGICAL  , intent(in)  :: DO_INCLUDE_MVOUTPUT

!  Local flags for the solution saving option

      LOGICAL  , intent(in)  :: DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  4/9/19. post-processing control now input

      INTEGER  , intent(in)  :: N_PPSTREAMS, PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)
      
!  Control integers

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: NLAYERS
      INTEGER  , intent(in)  :: N_USER_LEVELS

!  Partial layer bookkeeping

      LOGICAL  , intent(in)  :: PARTLAYERS_OUTFLAG  (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_OUTINDEX (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: UTAU_LEVEL_MASK_UP  (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_LAYERIDX (MAX_PARTLAYERS)

!  Input Fourier number and beam index
!  surface factor (2 for m = 0, 1 otherwise). Not required.

      INTEGER  , intent(in)  :: FOURIER, IBEAM

!  Linearization control

      INTEGER  , intent(in)  :: K_PARAMETERS
      INTEGER  , intent(in)  :: VARIATION_INDEX

!  Flux multiplier = F/4.pi

      REAL(fpk), intent(in)  :: FLUX_MULTIPLIER

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER  , intent(in)  :: TAYLOR_ORDER

!  Regular Inputs
!  --------------

!  Input optical properties after delta-M scaling

      REAL(fpk), intent(in)  :: PARTAU_VERT ( MAX_PARTLAYERS )
      REAL(fpk), intent(in)  :: DELTAU_VERT ( MAXLAYERS )

!  Holding arrays for Multiplier coefficients

      REAL(fpk), intent(in)  :: SIGMA_P ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      REAL(fpk), intent(in)  :: GAMMA_M ( MAXSTREAMS, MAXLAYERS )
      REAL(fpk), intent(in)  :: GAMMA_P ( MAXSTREAMS, MAXLAYERS )

!  User stream cosines

      REAL(fpk), intent(in)  :: USER_SECANTS ( MAX_USER_STREAMS )

!  Green's function particular integral arrays

      REAL(fpk), intent(in)  :: ATERM_SAVE ( MAXSTREAMS, MAXLAYERS )
      REAL(fpk), intent(in)  :: BTERM_SAVE ( MAXSTREAMS, MAXLAYERS )

!  Initial tramsittance factors for solar beams.

      REAL(fpk), intent(in)  :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: SOLARBEAM_CUTOFF(MAXBEAMS)

!  Quadrature streams, only required for Transonly

      REAL(fpk), intent(in)  :: QUAD_STREAMS(MAXSTREAMS)

!  Transmittance factors for +/- eigenvalues
!     User optical depths (UTUP and UTDN)

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: T_UTUP_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: T_UTDN_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)

!  Transmittance factors for average secant stream
!    Computed in the initial setup stage for Fourier m = 0

      REAL(fpk), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS,      MAXBEAMS )
      REAL(fpk), intent(in)  :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

!  Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage for Fourier m = 0

      REAL(fpk), intent(in)  :: T_DELT_USERM ( MAXLAYERS,      MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  Eigenvector solutions

      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Solution constants of integration

      REAL(fpk), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)

!  Solution constants of integration multiplied by homogeneous solutions

      REAL(fpk), intent(in)  :: LCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Eigenvectors defined at user-defined stream angles

      REAL(fpk), intent(in)  :: U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)

!  Particular beam solutions at user-defined stream angles

      REAL(fpk), intent(in)  :: U_WPOS(MAX_USER_STREAMS,MAXLAYERS)

!  Solution multipliers

      REAL(fpk), intent(in)  :: HMULT_1(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: HMULT_2(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: EMULT_UP (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

!  Multipliers

      REAL(fpk), intent(in)  :: UT_EMULT_UP (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS)
      REAL(fpk), intent(in)  :: UT_HMULT_UU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_HMULT_UD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Cumulative source terms

      REAL(fpk), intent(in)  :: CUMSOURCE_UP(MAX_USER_STREAMS,0:MAXLAYERS)

!  Source function integrated Green function multipliers (whole layer)

      REAL(fpk), intent(in)  :: PMULT_UU(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: PMULT_UD(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)

!  Source function integrated Green function multipliers (part layer)

      REAL(fpk), intent(in)  :: UT_PMULT_UU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_PMULT_UD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Green functions multipliers for off-grid optical depths

      REAL(fpk), intent(in)  :: UT_GMULT_UP(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_GMULT_DN(MAXSTREAMS,MAX_PARTLAYERS)

!  Linearized Inputs
!  -----------------

!  Linearized Eigensolutions

      REAL(fpk), intent(in)  :: L_XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Solution constants of integration, and related quantities

      REAL(fpk), intent(in)  :: NCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: PCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: NCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: PCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized General beam solutions at the Lower boundary

      REAL(fpk), intent(in)  :: L_WUPPER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized transmittances, homogeneous solutions

      REAL(fpk), intent(in)  :: L_T_DELT_EIGEN ( MAXSTREAMS, MAXLAYERS,      MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: L_T_UTUP_EIGEN ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: L_T_UTDN_EIGEN ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  Linearized Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage for Fourier m = 0

      REAL(fpk), intent(in)  :: L_T_DELT_USERM(MAXLAYERS,     MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_UTUP_USERM(MAX_PARTLAYERS,MAX_USER_STREAMS,MAX_ATMOSWFS)

!  Linearized input optical properties after delta-M scaling

      REAL(fpk), intent(in)  :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  Linearized (Positive) Eigenvalues

      REAL(fpk), intent(in)  :: L_KEIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearization of pseudo-spherical approximation

      REAL(fpk), intent(in)  :: LP_INITIAL_TRANS  ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS  )
      REAL(fpk), intent(in)  :: LP_AVERAGE_SECANT ( MAXLAYERS, MAXLAYERS,MAXBEAMS, MAX_ATMOSWFS  )

!  Linearizations of Saved quantities for the Green function solution

      REAL(fpk), intent(in)  :: L_ATERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_BTERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized transmittances, solar beam

      REAL(fpk), intent(in)  :: LP_T_DELT_MUBAR ( MAXLAYERS,      MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LP_T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized Eigenvectors defined at user-defined stream angles

      REAL(fpk), intent(in)  :: L_U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Particular beam solution (single scatter), user angles

      REAL(fpk), intent(in)  :: LP_U_WPOS(MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Integrated homogeneous solution multipliers, whole layer

      REAL(fpk), intent(in)  :: L_HMULT_1 (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_HMULT_2 (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized forcing term multipliers, whole layer

      REAL(fpk), intent(in)  :: LP_EMULT_UP (MAX_USER_STREAMS,MAXLAYERS,MAXLAYERS,MAXBEAMS,MAX_ATMOSWFS)

!  Linearized Integrated homogeneous solution multipliers, partial layer

      REAL(fpk), intent(in)  :: L_UT_HMULT_UU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_UT_HMULT_UD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)

!  Linearized forcing term multipliers, partial layer

      REAL(fpk), intent(in)  :: LP_UT_EMULT_UP (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXLAYERS,MAXBEAMS,MAX_ATMOSWFS)

!  Linearized direct thermal solution (upwelling)

      REAL(fpk), intent(in)  :: L_LAYER_TSUP_UP    ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: L_LAYER_TSUP_UTUP  ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  Linearized BOA source terms

      REAL(fpk), intent(in)  :: LP_BOA_MSSOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: LP_BOA_DBSOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)

!  Stuff for thermal transonly computations (Quadratures)
!  ------------------------------------------------------

!  Layer discrete ordinate transmittances and linearizations

      REAL(fpk), intent(in)  :: T_DELT_DISORDS   ( MAXSTREAMS, MAXLAYERS )
      REAL(fpk), intent(in)  :: L_T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS, MAX_ATMOSWFS )

!  Part-Layer discrete ordinate transmittances and linearizations

      REAL(fpk), intent(in)  :: T_DISORDS_UTUP   ( MAXSTREAMS, MAX_PARTLAYERS )
      REAL(fpk), intent(in)  :: L_T_DISORDS_UTUP ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  Thermal solution and linearization

      REAL(fpk), intent(in)  :: T_WUPPER   ( MAXSTREAMS_2, MAXLAYERS )
      REAL(fpk), intent(in)  :: L_T_WUPPER ( MAXSTREAMS_2, MAXLAYERS, MAX_ATMOSWFS )

!  Linearized partial-layer Thermal solution

      REAL(fpk), intent(in)  :: L_UT_T_PARTIC ( MAXSTREAMS_2, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  Special case thermal transmittance - BOA source term + linearization

      REAL(fpk), intent(in)  ::   BOA_THTONLY_SOURCE (MAXSTREAMS)
      REAL(fpk), intent(in)  :: L_BOA_THTONLY_SOURCE (MAXSTREAMS,MAX_ATMOSWFS)

!  Outputs
!  -------

!mick fix 6/29/11 - changed outputs from "out" to "inout"

!  Profile weighting functions at quadrature angles

      REAL(fpk), intent(inout) :: QUADPROFILEWF &
         ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAXSTREAMS, MAXBEAMS, MAX_DIRECTIONS )

!  Profile weighting functions at user angles

      REAL(fpk), intent(inout) :: PROFILEWF_F &
         ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAX_USER_STREAMS, MAXBEAMS, MAX_DIRECTIONS)

!  Linearized Green's function multipliers for off-grid optical depths
!   Will only be calculated as output, if the flag has been set

      LOGICAL  , intent(inout) :: FLAGS_LP_GMULT(MAX_PARTLAYERS)
      REAL(fpk), intent(inout) :: LP_UT_GMULT_UP (MAXSTREAMS,MAX_PARTLAYERS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(inout) :: LP_UT_GMULT_DN (MAXSTREAMS,MAX_PARTLAYERS,MAXLAYERS,MAX_ATMOSWFS)

!  local variables
!  ---------------

      LOGICAL    :: SOURCETERM_FLAG
      INTEGER    :: N, NUT, NSTART, NUT_PREV, NLEVEL
      INTEGER    :: UTA, UM, Q, NC, UT, IB, K, LUM, M

      REAL(fpk)  :: L_CUMUL_SOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(fpk)  :: L_LAYER_SOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(fpk)  :: L_FINAL_SOURCE


!  Indices

      IB = IBEAM
      M  = FOURIER
      K  = VARIATION_INDEX

!  Zero all Fourier components - New rule, better for safety

      DO UTA = 1, N_USER_LEVELS
        DO Q = 1, K_PARAMETERS
          DO LUM = 1, N_PPSTREAMS
            PROFILEWF_F(Q,K,UTA,LUM,IB,UPIDX) = ZERO
          ENDDO
        ENDDO
      ENDDO

!  Initialize post-processing recursion
!  ====================================

!  start the recursion

      IF ( DO_USER_STREAMS ) THEN
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            DO Q = 1, K_PARAMETERS
               L_CUMUL_SOURCE(UM,Q) = LP_BOA_MSSOURCE(UM,Q) + LP_BOA_DBSOURCE(UM,Q)
            ENDDO
        ENDDO
      ENDIF

!  Recursion Loop for linearized Post-processing
!  =============================================

!  initialise cumulative source term loop

      NC  = 0
      NUT = 0
      NSTART = NLAYERS
      NUT_PREV = NSTART + 1

!  loop over all output optical depths
!  -----------------------------------

      DO UTA = N_USER_LEVELS, 1, -1

!  Layer index for given optical depth

        NLEVEL = UTAU_LEVEL_MASK_UP(UTA)

!  Cumulative source terms to layer NUT (user-defined stream angles only)
!    1. Get layer source terms
!    2. Find cumulative source term
!    3. Set multiple scatter source term (MSST) output if flagged

        IF ( DO_USER_STREAMS ) THEN
          NUT = NLEVEL + 1
          DO N = NSTART, NUT, -1
            NC = NLAYERS + 1 - N

            SOURCETERM_FLAG = DO_LAYER_SCATTERING(FOURIER,N)
            CALL LP_WHOLELAYER_STERM_UP &
            ( DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY, & ! input
              DO_MSMODE_LIDORT, SOURCETERM_FLAG, M, TAYLOR_ORDER,            & ! input, FLags + Order
              IB, N, K, K_PARAMETERS, NSTREAMS, N_PPSTREAMS, PPSTREAM_MASK,  & ! input, Numbers
              DELTAU_VERT, L_DELTAU_VERT, USER_SECANTS, T_DELT_USERM,        & ! input, Optical + User_streams
              SOLARBEAM_CUTOFF,  INITIAL_TRANS,    T_DELT_MUBAR,             & ! input, Beam stuff
              LP_AVERAGE_SECANT, LP_INITIAL_TRANS, LP_T_DELT_MUBAR,          & ! input, Beam stuff
              U_XPOS, U_XNEG, U_WPOS, LCON, MCON,                            & ! input, RTE Solutions
              L_KEIGEN, L_U_XPOS, L_U_XNEG, LP_U_WPOS, NCON, PCON,           & ! input, Linearized RTE solutios
              GAMMA_M, GAMMA_P, SIGMA_P, L_LAYER_TSUP_UP,                    & ! input, solutions
              ATERM_SAVE, BTERM_SAVE, L_ATERM_SAVE, L_BTERM_SAVE,            & ! input, Greens Function Saved values
              HMULT_1, HMULT_2, EMULT_UP, PMULT_UU, PMULT_UD,                & ! input, Multipliers
              L_HMULT_1, L_HMULT_2, LP_EMULT_UP,                             & ! input, Linearized HMULT/EMULT
              L_LAYER_SOURCE )                                                 ! output

!  Cumulative sourceterm

            IF ( N.EQ.K ) THEN
               DO LUM = 1, N_PPSTREAMS
                  UM = PPSTREAM_MASK(LUM,IB)
                  DO Q = 1, K_PARAMETERS
                     L_CUMUL_SOURCE(UM,Q) = L_LAYER_SOURCE(UM,Q)         + &
                           T_DELT_USERM(N,UM)   * L_CUMUL_SOURCE(UM,Q)   + &
                         L_T_DELT_USERM(N,UM,Q) *   CUMSOURCE_UP(UM,NC-1)
                  ENDDO
               ENDDO
            ELSE
               DO LUM = 1, N_PPSTREAMS
                  UM = PPSTREAM_MASK(LUM,IB)
                  DO Q = 1, K_PARAMETERS
                     L_CUMUL_SOURCE(UM,Q) = L_LAYER_SOURCE(UM,Q)  + &
                          T_DELT_USERM(N,UM)*L_CUMUL_SOURCE(UM,Q)
                  ENDDO
               ENDDO
            ENDIF

!  End layer recursion

          ENDDO
        ENDIF

!  Offgrid output
!  --------------

        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

          UT = PARTLAYERS_OUTINDEX(UTA)
          N  = PARTLAYERS_LAYERIDX(UT)
          SOURCETERM_FLAG = DO_LAYER_SCATTERING(FOURIER,N)

!  Quadrature output at offgrid optical depths
!  ( Required if mean-value calculations are to be done)
!    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )

          IF ( DO_INCLUDE_MVOUTPUT ) THEN
            FLAGS_LP_GMULT(UT) = .TRUE.
            CALL QUADPROFILEWF_OFFGRID_UP  &
           ( DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,                 & ! Input
             TAYLOR_ORDER, NSTREAMS, NLAYERS, IB, UTA, UT, N, K, K_PARAMETERS,              & ! Input
             QUAD_STREAMS, FLUX_MULTIPLIER, PARTAU_VERT, L_DELTAU_VERT, SOLARBEAM_CUTOFF,   & ! Input
             INITIAL_TRANS, AVERAGE_SECANT,  T_DELT_MUBAR, T_UTDN_MUBAR,                    & ! Input
             T_UTUP_EIGEN, T_UTDN_EIGEN, T_DELT_DISORDS, T_DISORDS_UTUP,                    & ! Input
             GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE, UT_GMULT_UP, UT_GMULT_DN,            & ! Input
             XPOS, LCON_XVEC, MCON_XVEC, LCON, MCON, T_WUPPER, BOA_THTONLY_SOURCE,          & ! Input
             LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR, LP_T_UTDN_MUBAR,         & ! Input
             L_KEIGEN, L_T_UTUP_EIGEN, L_T_UTDN_EIGEN, L_T_DELT_DISORDS, L_T_DISORDS_UTUP,  & ! Input
             L_ATERM_SAVE, L_BTERM_SAVE, L_XPOS, L_XNEG, NCON_XVEC, PCON_XVEC,              & ! Input
             L_UT_T_PARTIC, L_T_WUPPER, L_BOA_THTONLY_SOURCE,                               & ! Input
             QUADPROFILEWF, LP_UT_GMULT_UP, LP_UT_GMULT_DN )                                  ! Output
          ENDIF

!  User-defined stream output, add additional partial layer source term

          IF ( DO_USER_STREAMS ) THEN

            CALL LP_PARTLAYER_STERM_UP  &
           ( DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,       & ! input, Flags
             DO_MSMODE_LIDORT, SOURCETERM_FLAG, M, TAYLOR_ORDER,                  & ! input, FLags + Order
             IB, UT, N, K, K_PARAMETERS, NSTREAMS, N_PPSTREAMS, PPSTREAM_MASK,    & ! input, Numbers
             PARTAU_VERT, DELTAU_VERT, L_DELTAU_VERT,                             & ! input, Optical
             USER_SECANTS, T_UTUP_USERM, L_T_UTUP_USERM,                          & ! input, Optical
             SOLARBEAM_CUTOFF, INITIAL_TRANS, T_DELT_MUBAR, T_UTDN_MUBAR,         & ! input, Beam stuff
             LP_AVERAGE_SECANT, LP_INITIAL_TRANS, LP_T_DELT_MUBAR,                & ! input, Beam stuff
             U_XPOS, U_XNEG, U_WPOS, LCON, MCON,                                  & ! input, RTE Solutions
             L_KEIGEN, L_U_XPOS, L_U_XNEG, LP_U_WPOS, NCON, PCON,                 & ! input, Linearized RTE solutios
             GAMMA_M, GAMMA_P, SIGMA_P, L_LAYER_TSUP_UTUP,                        & ! input, solutions
             ATERM_SAVE, BTERM_SAVE, L_ATERM_SAVE, L_BTERM_SAVE,                  & ! input, Greens Function values
             UT_HMULT_UU, UT_HMULT_UD, UT_EMULT_UP, UT_PMULT_UU, UT_PMULT_UD,     & ! input
             L_UT_HMULT_UU, L_UT_HMULT_UD, LP_UT_EMULT_UP,                        & ! input
             L_LAYER_SOURCE )                                                       ! output

!  Cumulative and final

            IF ( N.EQ.K ) THEN
               DO LUM = 1, N_PPSTREAMS
                  UM = PPSTREAM_MASK(LUM,IB)
                  DO Q = 1, K_PARAMETERS
                    L_FINAL_SOURCE = L_LAYER_SOURCE(UM,Q)             + &
                          T_UTUP_USERM(UT,UM)  * L_CUMUL_SOURCE(UM,Q) + &
                        L_T_UTUP_USERM(UT,UM,Q)*   CUMSOURCE_UP(UM,NC)
                    PROFILEWF_F(Q,K,UTA,LUM,IB,UPIDX) = FLUX_MULTIPLIER * L_FINAL_SOURCE
                  ENDDO
               ENDDO
            ELSE
               DO LUM = 1, N_PPSTREAMS
                  UM = PPSTREAM_MASK(LUM,IB)
                  DO Q = 1, K_PARAMETERS
                    L_FINAL_SOURCE = L_LAYER_SOURCE(UM,Q) + T_UTUP_USERM(UT,UM) * L_CUMUL_SOURCE(UM,Q)
                    PROFILEWF_F(Q,K,UTA,LUM,IB,UPIDX) = FLUX_MULTIPLIER * L_FINAL_SOURCE
                  ENDDO
                ENDDO
            ENDIF

          ENDIF

!  Ongrid output
!  -------------

        ELSE

!  Quadrature output at offgrid optical depths
!  ( Required if mean-value calculations are to be done)
!    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )

          IF ( DO_INCLUDE_MVOUTPUT ) THEN
            CALL QUADPROFILEWF_LEVEL_UP                              &
            ( DO_THERMAL_TRANSONLY, FLUX_MULTIPLIER,                 & ! Input
              NLAYERS, NSTREAMS, IB, UTA, NLEVEL, K, K_PARAMETERS,   & ! Input
              QUAD_STREAMS, L_XPOS, L_XNEG, L_WLOWER, L_WUPPER,      & ! Input
              LCON, LCON_XVEC, NCON_XVEC,   T_DELT_EIGEN,            & ! Input
              MCON, MCON_XVEC, PCON_XVEC, L_T_DELT_EIGEN,            & ! Input
              T_DELT_DISORDS, L_T_DELT_DISORDS, T_WUPPER, L_T_WUPPER,& ! Input
              BOA_THTONLY_SOURCE, L_BOA_THTONLY_SOURCE,              & ! Input
              QUADPROFILEWF )                                          ! Output
          ENDIF

!  User-defined stream output, just set to the cumulative source term

          IF ( DO_USER_STREAMS ) THEN
             DO LUM = 1, N_PPSTREAMS
                UM = PPSTREAM_MASK(LUM,IB)
                DO Q = 1, K_PARAMETERS
                   PROFILEWF_F(Q,K,UTA,LUM,IB,UPIDX) = FLUX_MULTIPLIER * L_CUMUL_SOURCE(UM,Q)
                ENDDO
             ENDDO
          ENDIF

        ENDIF

!  Check for updating the recursion

        IF ( DO_USER_STREAMS ) THEN
          IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
          NUT_PREV = NUT
        ENDIF

!  end loop over optical depth

      ENDDO

!  Finish

      RETURN
END SUBROUTINE UPUSER_PROFILEWF

!

SUBROUTINE DNUSER_PROFILEWF &
            ( DO_USER_STREAMS, DO_SOLAR_SOURCES, DO_MSMODE_LIDORT,                 & ! Input
              DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY, DO_INCLUDE_MVOUTPUT,    & ! Input
              NSTREAMS, NLAYERS, N_USER_LEVELS, VARIATION_INDEX, K_PARAMETERS,     & ! Input
              N_PPSTREAMS, PPSTREAM_MASK, DO_LAYER_SCATTERING, UTAU_LEVEL_MASK_DN, & ! Input
              PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,        & ! Input
              FOURIER, IBEAM, TAYLOR_ORDER, FLUX_MULTIPLIER, QUAD_STREAMS,         & ! Input
              PARTAU_VERT, DELTAU_VERT, L_DELTAU_VERT,                             & ! Input
              USER_SECANTS, SOLARBEAM_CUTOFF, INITIAL_TRANS, AVERAGE_SECANT,       & ! Input              
              T_DELT_MUBAR,   T_UTDN_MUBAR,                      & ! Input
              T_DELT_USERM,   T_UTDN_USERM,                      & ! Input
              T_DELT_EIGEN,   T_UTUP_EIGEN,   T_UTDN_EIGEN,      & ! Input
              T_DELT_DISORDS, T_DISORDS_UTDN,                    & ! Input
              GAMMA_M, GAMMA_P, SIGMA_M, ATERM_SAVE, BTERM_SAVE, & ! Input
              XPOS, LCON_XVEC, MCON_XVEC, LCON, MCON, T_WLOWER,          & ! Input
              U_XPOS, U_XNEG, U_WNEG, HMULT_1, HMULT_2, EMULT_DN,        & ! Input
              UT_HMULT_DU,  UT_HMULT_DD, UT_EMULT_DN,                    & ! Input
              PMULT_DU, PMULT_DD, UT_PMULT_DU, UT_PMULT_DD,              & ! Input
              UT_GMULT_UP, UT_GMULT_DN, CUMSOURCE_DN,                    & ! Input
              LP_INITIAL_TRANS, LP_AVERAGE_SECANT,                & ! Input
              LP_T_DELT_MUBAR,  LP_T_UTDN_MUBAR,                  & ! Input
              L_T_DELT_USERM,   L_T_UTDN_USERM,                   & ! Input
              L_T_DELT_EIGEN,   L_T_UTUP_EIGEN, L_T_UTDN_EIGEN,   & ! Input
              L_T_DELT_DISORDS, L_T_DISORDS_UTDN,                 & ! Input
              L_KEIGEN, L_ATERM_SAVE, L_BTERM_SAVE,               & ! Input
              L_U_XPOS, L_U_XNEG, LP_U_WNEG, NCON, PCON,                 & ! Input
              L_XPOS, L_XNEG, L_WLOWER, NCON_XVEC, PCON_XVEC,            & ! Input
              L_HMULT_1, L_HMULT_2, LP_EMULT_DN,                         & ! Input
              L_UT_HMULT_DU, L_UT_HMULT_DD, LP_UT_EMULT_DN,              & ! Input
              L_LAYER_TSUP_DN, L_LAYER_TSUP_UTDN, L_UT_T_PARTIC,         & ! Input
              L_T_WLOWER, LP_TOA_SOURCE,                                 & ! Input
              FLAGS_LP_GMULT, LP_UT_GMULT_UP, LP_UT_GMULT_DN,            & ! Output
              PROFILEWF_F, QUADPROFILEWF )                                 ! Output

!  Downwelling post-processed Profile Jacobian Fourier component

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAXSTREAMS, MAXSTREAMS_2, MAXMOMENTS, MAXLAYERS, MAX_PARTLAYERS, &
                                MAXBEAMS, MAX_USER_LEVELS, MAX_USER_STREAMS, MAX_ATMOSWFS,       &
                                MAX_DIRECTIONS, ZERO, DNIDX

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  Control flags

      LOGICAL  , intent(in)  :: DO_USER_STREAMS
      LOGICAL  , intent(in)  :: DO_MSMODE_LIDORT

!  Solar source term flag

      LOGICAL  , intent(in)  :: DO_SOLAR_SOURCES

!  Thermal emission flag

      LOGICAL  , intent(in)  :: DO_INCLUDE_THERMEMISS

!  Thermal transmittance-only flag

      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY

!  Local control flags

      LOGICAL  , intent(in)  :: DO_INCLUDE_MVOUTPUT

!  Control integers

      INTEGER  , intent(in)  :: NSTREAMS, NLAYERS
      INTEGER  , intent(in)  :: N_USER_LEVELS

!  Local flags for the solution saving option

      LOGICAL  , intent(in)  :: DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  4/9/19. post-processing control now input

      INTEGER  , intent(in)  :: N_PPSTREAMS, PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)
      
!  Partial layer bookkeeping

      INTEGER  , intent(in)  :: UTAU_LEVEL_MASK_DN  (MAX_USER_LEVELS)
      LOGICAL  , intent(in)  :: PARTLAYERS_OUTFLAG  (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_OUTINDEX (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_LAYERIDX (MAX_PARTLAYERS)

!  Input Fourier number and beam index

      INTEGER  , intent(in)  :: FOURIER, IBEAM

!  Linearization control

      INTEGER  , intent(in)  :: K_PARAMETERS
      INTEGER  , intent(in)  :: VARIATION_INDEX

!  Flux multiplier = F/4.pi

      REAL(fpk), intent(in)  :: FLUX_MULTIPLIER

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER  , intent(in)  :: TAYLOR_ORDER

!  Regular Inputs
!  --------------

!  Input optical properties after delta-M scaling

      REAL(fpk), intent(in)  :: PARTAU_VERT ( MAX_PARTLAYERS )
      REAL(fpk), intent(in)  :: DELTAU_VERT ( MAXLAYERS )

!  Holding arrays for Multiplier coefficients

      REAL(fpk), intent(in)  :: SIGMA_M ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      REAL(fpk), intent(in)  :: GAMMA_M ( MAXSTREAMS, MAXLAYERS )
      REAL(fpk), intent(in)  :: GAMMA_P ( MAXSTREAMS, MAXLAYERS )

!  User stream cosines

      REAL(fpk), intent(in)  :: USER_SECANTS ( MAX_USER_STREAMS )

!  Green's function particular integral arrays

      REAL(fpk), intent(in)  :: ATERM_SAVE ( MAXSTREAMS, MAXLAYERS )
      REAL(fpk), intent(in)  :: BTERM_SAVE ( MAXSTREAMS, MAXLAYERS )

!  initial tramsittance factors for solar beams.

      REAL(fpk), intent(in)  :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: SOLARBEAM_CUTOFF(MAXBEAMS)

!  Quadrature streams, only required for Transonly

      REAL(fpk), intent(in)  :: QUAD_STREAMS(MAXSTREAMS)

!  transmittance factors for +/- eigenvalues
!     User optical depths (UTUP and UTDN)

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: T_UTUP_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: T_UTDN_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)

!  Transmittance factors for average secant stream
!    Computed in the initial setup stage for Fourier m = 0

      REAL(fpk), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS,      MAXBEAMS )
      REAL(fpk), intent(in)  :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

!  Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage for Fourier m = 0

      REAL(fpk), intent(in)  :: T_DELT_USERM ( MAXLAYERS,      MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  Eigenvector solutions

      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Solution constants of integration

      REAL(fpk), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)

!  Solution constants of integration multiplied by homogeneous solutions

      REAL(fpk), intent(in)  :: LCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Eigenvectors defined at user-defined stream angles
!     EP for the positive KEIGEN values, EM for -ve KEIGEN

      REAL(fpk), intent(in)  :: U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)

!  Particular beam solutions at user-defined stream angles

      REAL(fpk), intent(in)  :: U_WNEG(MAX_USER_STREAMS,MAXLAYERS)

!  solution multipliers

      REAL(fpk), intent(in)  :: HMULT_1(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: HMULT_2(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: EMULT_DN (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

!  Multipliers

      REAL(fpk), intent(in)  :: UT_EMULT_DN (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS)
      REAL(fpk), intent(in)  :: UT_HMULT_DU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_HMULT_DD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Cumulative source terms

      REAL(fpk), intent(in)  :: CUMSOURCE_DN(MAX_USER_STREAMS,0:MAXLAYERS)

!  Source function integrated Green function multipliers (whole layer)

      REAL(fpk), intent(in)  :: PMULT_DU(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: PMULT_DD(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)

!   Source function integrated Green function multipliers (part layer)

      REAL(fpk), intent(in)  :: UT_PMULT_DU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_PMULT_DD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Green functions multipliers for off-grid optical depths

      REAL(fpk), intent(in)  :: UT_GMULT_UP(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_GMULT_DN(MAXSTREAMS,MAX_PARTLAYERS)

!  Linearized Inputs
!  -----------------

!  Linearized Eigensolutions

      REAL(fpk), intent(in)  :: L_XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Solution constants of integration, and related quantities

      REAL(fpk), intent(in)  :: NCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: PCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: NCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: PCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized General beam solutions at the Lower boundary

      REAL(fpk), intent(in)  :: L_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized transmittances, homogeneous solutions

      REAL(fpk), intent(in)  :: L_T_DELT_EIGEN ( MAXSTREAMS, MAXLAYERS,      MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: L_T_UTUP_EIGEN ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: L_T_UTDN_EIGEN ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  Linearized Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage for Fourier m = 0

      REAL(fpk), intent(in)  :: L_T_DELT_USERM(MAXLAYERS,     MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_UTDN_USERM(MAX_PARTLAYERS,MAX_USER_STREAMS,MAX_ATMOSWFS)

!  Linearized input optical properties after delta-M scaling

      REAL(fpk), intent(in)  :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  Linearized (Positive) Eigenvalues

      REAL(fpk), intent(in)  :: L_KEIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearizations of Saved quantities for the Green function solution

      REAL(fpk), intent(in)  :: L_ATERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_BTERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearization of pseudo-spherical approximation

      REAL(fpk), intent(in)  :: LP_INITIAL_TRANS ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS  )
      REAL(fpk), intent(in)  :: LP_AVERAGE_SECANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized transmittances, solar beam

      REAL(fpk), intent(in)  :: LP_T_DELT_MUBAR ( MAXLAYERS,      MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LP_T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized Eigenvectors defined at user-defined stream angles

      REAL(fpk), intent(in)  :: L_U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Particular beam solution (single scatter), user angles

      REAL(fpk), intent(in)  :: LP_U_WNEG(MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Integrated homogeneous solution multipliers, whole layer

      REAL(fpk), intent(in)  :: L_HMULT_1 (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_HMULT_2 (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Integrated homogeneous solution multipliers, partial layer

      REAL(fpk), intent(in)  :: L_UT_HMULT_DU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_UT_HMULT_DD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)

!  Linearized forcing term multipliers, whole layer

      REAL(fpk), intent(in)  :: LP_EMULT_DN (MAX_USER_STREAMS,MAXLAYERS,MAXLAYERS,MAXBEAMS,MAX_ATMOSWFS)

!  Linearized forcing term multipliers, partial layer

      REAL(fpk), intent(in)  :: LP_UT_EMULT_DN (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXLAYERS,MAXBEAMS,MAX_ATMOSWFS)

!  Thermal solutions for the Trans-only special case

      REAL(fpk), intent(in)  :: L_LAYER_TSUP_DN   (MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_LAYER_TSUP_UTDN (MAX_USER_STREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)

!  Linearized TOA source terms

      REAL(fpk), intent(in)  :: LP_TOA_SOURCE (MAX_USER_STREAMS,MAX_ATMOSWFS)

!  Stuff for thermal transonly computations (Quadratures)
!  ------------------------------------------------------

!  Layer discrete ordinate transmittances and linearizations

      REAL(fpk), intent(in)  :: T_DELT_DISORDS   (MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: L_T_DELT_DISORDS (MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Part-Layer discrete ordinate transmittances and linearizations

      REAL(fpk), intent(in)  :: T_DISORDS_UTDN   (MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: L_T_DISORDS_UTDN (MAXSTREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)

!  Thermal solutions at the Lower boundary

      REAL(fpk), intent(in)  :: T_WLOWER   (MAXSTREAMS_2,MAXLAYERS)
      REAL(fpk), intent(in)  :: L_T_WLOWER (MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Particular solutions

      REAL(fpk), intent(in)  :: L_UT_T_PARTIC (MAXSTREAMS_2,MAX_PARTLAYERS,MAX_ATMOSWFS)

!  Outputs
!  -------

!mick fix 6/29/11 - changed outputs from "out" to "inout"

!  Profile weighting functions at quadrature angles

      REAL(fpk), intent(inout) :: QUADPROFILEWF &
         ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAXSTREAMS, MAXBEAMS, MAX_DIRECTIONS )

!  Profile weighting functions at user angles

      REAL(fpk), intent(inout) :: PROFILEWF_F &
         ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAX_USER_STREAMS, MAXBEAMS, MAX_DIRECTIONS)

!  Linearized Green functions multipliers for off-grid optical depths

      LOGICAL  , intent(inout) :: FLAGS_LP_GMULT(MAX_PARTLAYERS)
      REAL(fpk), intent(inout) :: LP_UT_GMULT_UP (MAXSTREAMS,MAX_PARTLAYERS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(inout) :: LP_UT_GMULT_DN (MAXSTREAMS,MAX_PARTLAYERS,MAXLAYERS,MAX_ATMOSWFS)

!  local variables
!  ---------------

      LOGICAL    :: LOCAL_LP_GMULT
      LOGICAL    :: SOURCETERM_FLAG
      INTEGER    :: N, NUT, NSTART, NUT_PREV, NLEVEL
      INTEGER    :: UTA, UM, Q, NC, UT, IB, K, LUM, M

      REAL(fpk)  :: L_CUMUL_SOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(fpk)  :: L_LAYER_SOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(fpk)  :: L_FINAL_SOURCE

!  Initialise local indices

      IB = IBEAM
      M  = FOURIER
      K  = VARIATION_INDEX

!  Zero all Fourier component output

      DO UTA = 1, N_USER_LEVELS
        DO Q = 1, K_PARAMETERS
          DO LUM = 1, N_PPSTREAMS
            PROFILEWF_F(Q,K,UTA,LUM,IB,DNIDX) = ZERO
          ENDDO
        ENDDO
      ENDDO

!  Initialize post-processing recursion
!  ====================================

! Get the linearized TOA source terms

      IF ( DO_USER_STREAMS ) THEN
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IBEAM)
            DO Q = 1, K_PARAMETERS
              L_CUMUL_SOURCE(UM,Q) = LP_TOA_SOURCE(UM,Q)
            ENDDO
         ENDDO
      ENDIF

!  Recursion Loop for linearized Post-processing
!  =============================================

!  initialise cumulative source term loop

      NC  = 0
      NUT = 0
      NSTART = 1
      NUT_PREV = NSTART - 1

!  loop over all output optical depths
!  -----------------------------------

      DO UTA = 1, N_USER_LEVELS

!  Layer index for given optical depth

        NLEVEL = UTAU_LEVEL_MASK_DN(UTA)

!  Cumulative source terms to layer NUT (user-defined stream angles only)
!    1. Get layer source terms
!    2. Find cumulative source term
!    3. Set multiple scatter source term output if flagged

        IF ( DO_USER_STREAMS ) THEN
          NUT = NLEVEL
          DO N = NSTART, NUT
            NC = N

!  Layer source term

            SOURCETERM_FLAG = DO_LAYER_SCATTERING(FOURIER,N)
            CALL LP_WHOLELAYER_STERM_DN  &
            ( DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY, & ! input
              DO_MSMODE_LIDORT, SOURCETERM_FLAG, M, TAYLOR_ORDER,            & ! input, FLags + Order
              IB, N, K, K_PARAMETERS, NSTREAMS, N_PPSTREAMS, PPSTREAM_MASK,  & ! input, Numbers
              DELTAU_VERT, L_DELTAU_VERT, USER_SECANTS, T_DELT_USERM,        & ! input, Optical + User_streams
              SOLARBEAM_CUTOFF, AVERAGE_SECANT, INITIAL_TRANS, T_DELT_MUBAR, & ! input, Beam stuff
              LP_AVERAGE_SECANT, LP_INITIAL_TRANS, LP_T_DELT_MUBAR,          & ! input, Beam stuff
              U_XPOS, U_XNEG, U_WNEG, LCON, MCON,                            & ! input, RTE Solutions
              L_KEIGEN, L_U_XPOS, L_U_XNEG, LP_U_WNEG, NCON, PCON,           & ! input, Linearized RTE solutios
              GAMMA_M, GAMMA_P, SIGMA_M, L_LAYER_TSUP_DN,                    & ! input, solutions
              ATERM_SAVE, BTERM_SAVE, L_ATERM_SAVE, L_BTERM_SAVE,            & ! input, Greens Function Saved values
              HMULT_1, HMULT_2, EMULT_DN, PMULT_DU, PMULT_DD,                & ! input, Multipliers
              L_HMULT_1, L_HMULT_2, LP_EMULT_DN,                             & ! input, Linearized HMULT/EMULT
              L_LAYER_SOURCE )                                                 ! output

!  Cumulative source

            IF ( N.EQ.K ) THEN
               DO LUM = 1, N_PPSTREAMS
                  UM = PPSTREAM_MASK(LUM,IBEAM)
                  DO Q = 1, K_PARAMETERS
                    L_CUMUL_SOURCE(UM,Q) = L_LAYER_SOURCE(UM,Q)  + &
                         T_DELT_USERM(N,UM)*L_CUMUL_SOURCE(UM,Q) + &
                        L_T_DELT_USERM(N,UM,Q)*CUMSOURCE_DN(UM,NC-1)
                  ENDDO
                ENDDO
            ELSE
               DO LUM = 1, N_PPSTREAMS
                  UM = PPSTREAM_MASK(LUM,IBEAM)
                  DO Q = 1, K_PARAMETERS
                    L_CUMUL_SOURCE(UM,Q) = L_LAYER_SOURCE(UM,Q) + &
                         T_DELT_USERM(N,UM)*L_CUMUL_SOURCE(UM,Q)
                  ENDDO
                ENDDO
            ENDIF

          ENDDO
        ENDIF

!  Offgrid output
!  --------------

        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

          UT = PARTLAYERS_OUTINDEX(UTA)
          N  = PARTLAYERS_LAYERIDX(UT)
          SOURCETERM_FLAG = DO_LAYER_SCATTERING(FOURIER,N)

!  Quadrature output at offgrid optical depths
!  ( Required if mean-value calculations are to be done)
!    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )

          IF ( DO_INCLUDE_MVOUTPUT ) THEN
            LOCAL_LP_GMULT = FLAGS_LP_GMULT(UT)
            CALL QUADPROFILEWF_OFFGRID_DN                                   &
            ( DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,                 & ! Input
              TAYLOR_ORDER, NSTREAMS, IB, UTA, UT, N, K, K_PARAMETERS,                       & ! Input
              QUAD_STREAMS, FLUX_MULTIPLIER, PARTAU_VERT, L_DELTAU_VERT, SOLARBEAM_CUTOFF,   & ! Input
              LOCAL_LP_GMULT, INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR, T_UTDN_MUBAR,     & ! Input
              T_UTUP_EIGEN, T_UTDN_EIGEN, T_DELT_DISORDS, T_DISORDS_UTDN,                    & ! Input
              GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE, UT_GMULT_UP, UT_GMULT_DN,            & ! Input
              XPOS, LCON_XVEC, MCON_XVEC, LCON, MCON, T_WLOWER,                              & ! Input
              LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR, LP_T_UTDN_MUBAR,         & ! Input
              L_KEIGEN, L_T_UTUP_EIGEN, L_T_UTDN_EIGEN, L_T_DELT_DISORDS, L_T_DISORDS_UTDN,  & ! Input
              L_ATERM_SAVE, L_BTERM_SAVE, L_XPOS, L_XNEG, NCON_XVEC, PCON_XVEC,              & ! Input
              L_UT_T_PARTIC, L_T_WLOWER,                                                     & ! Input
              QUADPROFILEWF, LP_UT_GMULT_UP, LP_UT_GMULT_DN )                                  ! Output
            FLAGS_LP_GMULT(UT) = .FALSE.
          ENDIF

!  User-defined stream output, add additional partial layer source term

          IF ( DO_USER_STREAMS ) THEN
            CALL LP_PARTLAYER_STERM_DN &
           ( DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,    & ! input, Flags
             DO_MSMODE_LIDORT, SOURCETERM_FLAG, M, TAYLOR_ORDER,               & ! input, FLags + Order
             IB, UT, N, K, K_PARAMETERS, NSTREAMS, N_PPSTREAMS, PPSTREAM_MASK, & ! input, Numbers
             PARTAU_VERT, L_DELTAU_VERT, USER_SECANTS, T_UTDN_USERM,           & ! input, Optical
             SOLARBEAM_CUTOFF, INITIAL_TRANS, T_DELT_MUBAR, T_UTDN_MUBAR,      & ! input, Beam stuff
             LP_AVERAGE_SECANT, LP_INITIAL_TRANS, LP_T_DELT_MUBAR,             & ! input, Beam stuff
             U_XPOS, U_XNEG, U_WNEG, LCON, MCON,                               & ! input, RTE Solutions
             L_KEIGEN, L_U_XPOS, L_U_XNEG, LP_U_WNEG, NCON, PCON,              & ! input, Linearized RTE solutios
             GAMMA_M, GAMMA_P, SIGMA_M, L_LAYER_TSUP_UTDN,                     & ! input, solutions
             ATERM_SAVE, BTERM_SAVE, L_ATERM_SAVE, L_BTERM_SAVE,               & ! input, Greens Function Saved values
             UT_HMULT_DU, UT_HMULT_DD, UT_EMULT_DN, UT_PMULT_DU, UT_PMULT_DD,  & ! Input, Multipliers
             L_UT_HMULT_DU, L_UT_HMULT_DD, LP_UT_EMULT_DN,                     & ! Input, Multipliers
             L_LAYER_SOURCE )                                                    ! output

!  Cumulative and final

            IF ( N.EQ.K ) THEN
               DO LUM = 1, N_PPSTREAMS
                  UM = PPSTREAM_MASK(LUM,IBEAM)
                  DO Q = 1, K_PARAMETERS
                    L_FINAL_SOURCE = L_LAYER_SOURCE(UM,Q)            + &
                        T_UTDN_USERM(UT,UM)   * L_CUMUL_SOURCE(UM,Q) + &
                      L_T_UTDN_USERM(UT,UM,Q) *   CUMSOURCE_DN(UM,NC)
                    PROFILEWF_F(Q,K,UTA,LUM,IB,DNIDX) = FLUX_MULTIPLIER * L_FINAL_SOURCE
                  ENDDO
               ENDDO
            ELSE
               DO LUM = 1, N_PPSTREAMS
                  UM = PPSTREAM_MASK(LUM,IBEAM)
                  DO Q = 1, K_PARAMETERS
                    L_FINAL_SOURCE = L_LAYER_SOURCE(UM,Q) + &
                        T_UTDN_USERM(UT,UM) * L_CUMUL_SOURCE(UM,Q)
                    PROFILEWF_F(Q,K,UTA,LUM,IB,DNIDX) = FLUX_MULTIPLIER * L_FINAL_SOURCE
                  ENDDO
               ENDDO
            ENDIF

          ENDIF

!  Ongrid output
!  -------------

        ELSE

!  Quadrature output at offgrid optical depths
!  ( Required if mean-value calculations are to be done)
!    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )

          IF ( DO_INCLUDE_MVOUTPUT ) THEN
            CALL QUADPROFILEWF_LEVEL_DN                  &
            ( DO_THERMAL_TRANSONLY, FLUX_MULTIPLIER,     & ! Input
              NSTREAMS, IB, UTA, NLEVEL, K, K_PARAMETERS,& ! Input
              QUAD_STREAMS, L_XPOS, L_XNEG, L_WLOWER,    & ! Input
              LCON, LCON_XVEC, NCON_XVEC, T_DELT_EIGEN,  & ! Input
              MCON, PCON_XVEC, L_T_DELT_EIGEN,           & ! Input
              T_DELT_DISORDS, L_T_DELT_DISORDS, T_WLOWER,& ! Input
              L_T_WLOWER,                                & ! Input
              QUADPROFILEWF )                              ! Output
          ENDIF

!  User-defined stream output, just set to the cumulative source term

          IF ( DO_USER_STREAMS ) THEN
             DO LUM = 1, N_PPSTREAMS
                UM = PPSTREAM_MASK(LUM,IBEAM)
                DO Q = 1, K_PARAMETERS
                   PROFILEWF_F(Q,K,UTA,LUM,IB,DNIDX) = FLUX_MULTIPLIER * L_CUMUL_SOURCE(UM,Q)
                ENDDO
             ENDDO
          ENDIF

        ENDIF

!  Check for updating the recursion

        IF ( DO_USER_STREAMS ) THEN
          IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
          NUT_PREV = NUT
        ENDIF

!  end loop over optical depth

      ENDDO

!  Finish

      RETURN
END SUBROUTINE DNUSER_PROFILEWF

!

SUBROUTINE MIFLUX_PROFILEWF &
           ( DO_UPWELLING, DO_DNWELLING, DO_INCLUDE_DIRECTBEAM,         & ! Input
             IB, K, K_PARAMETERS, NSTREAMS, N_USER_LEVELS, FLUX_FACTOR, & ! Input
             PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,                   & ! Input
             UTAU_LEVEL_MASK_DN, PARTLAYERS_LAYERIDX,                   & ! Input
             QUAD_WEIGHTS, QUAD_STRMWTS,                                & ! Input
             LP_LEVELS_SOLARTRANS, LP_PARTIALS_SOLARTRANS,              & ! Input ! Added Partials, 1/9/18
             INITIAL_TRANS, SOLARBEAM_CUTOFF, LOCAL_CSZA,               & ! Input
             T_DELT_MUBAR, T_UTDN_MUBAR, QUADPROFILEWF,                 & ! Input
             LP_INITIAL_TRANS, LP_T_DELT_MUBAR, LP_T_UTDN_MUBAR,        & ! Input
             MEANI_DIFFUSE_PROFWF,  FLUX_DIFFUSE_PROFWF,                & ! Output
             DNMEANI_DIRECT_PROFWF, DNFLUX_DIRECT_PROFWF )                ! Output

!  Profile Jacobians for the hemispherically integrated fields

!  Variables renamed, 7/18/17. Distinguish between Diffuse and Direct.
!  1/9/18, Upgrade to fill the LC_PARTIALS_SOLARTRANS hole.

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
                                MAXSTREAMS,   MAXBEAMS,  MAX_PARTLAYERS,  &
                                MAX_DIRECTIONS, ZERO, HALF, PI2, PI4, UPIDX, DNIDX

      IMPLICIT NONE

!  input arguments
!  ---------------

!  Flags

      LOGICAL  , intent(in)  :: DO_UPWELLING
      LOGICAL  , intent(in)  :: DO_DNWELLING
      LOGICAL  , intent(in)  :: DO_INCLUDE_DIRECTBEAM

!  Thread
!      INTEGER  , intent(in)  :: THREAD

!  Index

      INTEGER  , intent(in)  :: IB

!  linearization control

      INTEGER  , intent(in)  :: K, K_PARAMETERS

!  Control integers

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: N_USER_LEVELS

!  Partial layer bookkeeping

      LOGICAL  , intent(in)  :: PARTLAYERS_OUTFLAG  (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_OUTINDEX (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: UTAU_LEVEL_MASK_DN  (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_LAYERIDX (MAX_PARTLAYERS)

!  Flux factor

      REAL(fpk), intent(in)  :: FLUX_FACTOR

!  Quadrature

      REAL(fpk), intent(in)  :: QUAD_WEIGHTS ( MAXSTREAMS )
      REAL(fpk), intent(in)  :: QUAD_STRMWTS ( MAXSTREAMS )

!  Rob fix 7/18/17 Added Arguments

      REAL(fpk), intent(in)  :: LP_LEVELS_SOLARTRANS ( 0:MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Partials added, 1/9/18

      REAL(fpk), intent(in)  :: LP_PARTIALS_SOLARTRANS ( MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  initial tramsittance factors for solar beams.

      REAL(fpk), intent(in)  :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: SOLARBEAM_CUTOFF(MAXBEAMS)

!  local solar zenith angle cosine

      REAL(fpk), intent(in)  :: LOCAL_CSZA ( MAXLAYERS, MAXBEAMS )

!  Transmittance factors for average secant stream

      REAL(fpk), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

!  Linearized initial tramsittance factors for solar beams.

      REAL(fpk), intent(in)  :: LP_INITIAL_TRANS ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS  )

!  Linearized transmittances

      REAL(fpk), intent(in)  :: LP_T_DELT_MUBAR ( MAXLAYERS,      MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LP_T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Quadrature-defined weighting functions

      REAL(fpk), intent(in)  :: QUADPROFILEWF ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
                                                MAXSTREAMS,   MAXBEAMS,  MAX_DIRECTIONS )

!  Output arguments
!  ----------------

!  Mean intensity (actinic flux), Regular Flux, Diffuse values

      REAL(fpk), intent(inout) :: MEANI_DIFFUSE_PROFWF( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS )
      REAL(fpk), intent(inout) :: FLUX_DIFFUSE_PROFWF ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS )

! Direct-beam contributions output separately, added 7/18/17

      REAL(fpk), intent(inout) :: DNMEANI_DIRECT_PROFWF( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAXBEAMS )
      REAL(fpk), intent(inout) :: DNFLUX_DIRECT_PROFWF ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAXBEAMS )

!  local variables
!  ----------------

      INTEGER    :: I, UTA, UT, Q, N
      REAL(fpk)  :: SM, SF, FMU0, HELP
      REAL(fpk)  :: L_TRANS, L_DIRECT_FLUX, L_DIRECT_MEANI
      REAL(fpk)  :: L_TRANS_SCALED, L_DIRECT_FLUX_SCALED, L_DIRECT_MEANI_SCALED

!  mean intensity and flux
!  -----------------------

!  Upwelling

      IF ( DO_UPWELLING ) THEN
        DO UTA = 1, N_USER_LEVELS
         DO Q = 1, K_PARAMETERS
          SM = ZERO
          SF = ZERO
          DO I = 1, NSTREAMS
           SM = SM + QUAD_WEIGHTS(I)*QUADPROFILEWF(Q,K,UTA,I,IB,UPIDX)
           SF = SF + QUAD_STRMWTS(I)*QUADPROFILEWF(Q,K,UTA,I,IB,UPIDX)
          ENDDO
          MEANI_DIFFUSE_PROFWF(Q,K,UTA,IB,UPIDX) = SM * HALF
          FLUX_DIFFUSE_PROFWF (Q,K,UTA,IB,UPIDX) = SF * PI2
         ENDDO
        ENDDO
      ENDIF

!  Downwelling

      IF ( DO_DNWELLING ) THEN

!  Diffuse term contribution

        DO UTA = 1, N_USER_LEVELS
         DO Q = 1, K_PARAMETERS
          SM = ZERO
          SF = ZERO
          DO I = 1, NSTREAMS
           SM = SM + QUAD_WEIGHTS(I)*QUADPROFILEWF(Q,K,UTA,I,IB,DNIDX)
           SF = SF + QUAD_STRMWTS(I)*QUADPROFILEWF(Q,K,UTA,I,IB,DNIDX)
          ENDDO
          MEANI_DIFFUSE_PROFWF(Q,K,UTA,IB,DNIDX) = SM * HALF
          FLUX_DIFFUSE_PROFWF (Q,K,UTA,IB,DNIDX) = SF * PI2
         ENDDO
        ENDDO

!  nothing to do if no reflection of solar source direct beam contribution.
!  NOT RELEVANT HERE.--> Bug 18 july 2017.
!mick fix 3/22/2017 - turned back on
        IF ( .NOT. DO_INCLUDE_DIRECTBEAM ) RETURN

!  For the downward direction, add the direct beam contributions

        DO UTA = 1, N_USER_LEVELS

!  For the offgrid values

          IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
            UT = PARTLAYERS_OUTINDEX(UTA)
            N  = PARTLAYERS_LAYERIDX(UT)

!  Only contributions for layers above the PI cutoff
!    L_INITIAL_TRANS is a logarithmic derivative

            IF ( N .LE. SOLARBEAM_CUTOFF(IB) ) THEN

              FMU0 = LOCAL_CSZA(N,IB) * FLUX_FACTOR
              DO Q = 1, K_PARAMETERS

!  direct transmittances, scaled and unscaled.
!     LTRANS, Placeholder filled, 1/9/18

                HELP = LP_T_UTDN_MUBAR(UT,K,IB,Q) + LP_INITIAL_TRANS(N,K,IB,Q) * T_UTDN_MUBAR(UT,IB)
                L_TRANS_SCALED = HELP * INITIAL_TRANS(N,IB)
                L_TRANS        = LP_PARTIALS_SOLARTRANS(UT,K,IB,Q)  ! Filled PLACEHOLDER, 1/9/18

!  direct calculation with non-scaled linearized transmittances

                L_DIRECT_MEANI = FLUX_FACTOR * L_TRANS / PI4
                L_DIRECT_FLUX  = FMU0  * L_TRANS
                DNMEANI_DIRECT_PROFWF(Q,K,UTA,IB) = L_DIRECT_MEANI
                DNFLUX_DIRECT_PROFWF( Q,K,UTA,IB) = L_DIRECT_FLUX

!  Diffuse calculation

                L_DIRECT_MEANI_SCALED = FLUX_FACTOR * L_TRANS_SCALED / PI4
                L_DIRECT_FLUX_SCALED  = FMU0  * L_TRANS_SCALED
                MEANI_DIFFUSE_PROFWF(Q,K,UTA,IB,DNIDX) = &
                     MEANI_DIFFUSE_PROFWF(Q,K,UTA,IB,DNIDX) + ( L_DIRECT_MEANI_SCALED - L_DIRECT_MEANI )
                FLUX_DIFFUSE_PROFWF (Q,K,UTA,IB,DNIDX) = &
                     FLUX_DIFFUSE_PROFWF (Q,K,UTA,IB,DNIDX) + ( L_DIRECT_FLUX_SCALED  - L_DIRECT_FLUX  )

              ENDDO

            ENDIF

!  For the on-grid balues
!    Direct-beam contributions output separately, 26 May 11
!      Bug Fix, 18 July 2017. for Downwelling. 

          ELSE

            N = UTAU_LEVEL_MASK_DN(UTA)
            IF ( N .LE. SOLARBEAM_CUTOFF(IB) ) THEN

              IF ( N.GT.0 ) THEN
                FMU0 = LOCAL_CSZA(N,IB) * FLUX_FACTOR
                DO Q = 1, K_PARAMETERS

!  direct transmittances, scaled and unscaled.

                  HELP = LP_T_DELT_MUBAR(N,K,IB,Q) + LP_INITIAL_TRANS(N,K,IB,Q) * T_DELT_MUBAR(N,IB)
                  L_TRANS_SCALED = HELP * INITIAL_TRANS(N,IB)
                  L_TRANS        = LP_LEVELS_SOLARTRANS(N,K,IB,Q)

!  direct calculation with non-scaled linearized transmittances

                  L_DIRECT_MEANI = FLUX_FACTOR * L_TRANS / PI4
                  L_DIRECT_FLUX  = FMU0 * L_TRANS
                  DNMEANI_DIRECT_PROFWF(Q,K,UTA,IB) = L_DIRECT_MEANI
                  DNFLUX_DIRECT_PROFWF( Q,K,UTA,IB) = L_DIRECT_FLUX

!  Diffuse calculation

                  L_DIRECT_MEANI_SCALED = FLUX_FACTOR * L_TRANS_SCALED / PI4
                  L_DIRECT_FLUX_SCALED  = FMU0 * L_TRANS_SCALED
                  MEANI_DIFFUSE_PROFWF(Q,K,UTA,IB,DNIDX) = &
                     MEANI_DIFFUSE_PROFWF(Q,K,UTA,IB,DNIDX) + ( L_DIRECT_MEANI_SCALED - L_DIRECT_MEANI )
                  FLUX_DIFFUSE_PROFWF (Q,K,UTA,IB,DNIDX) = &
                     FLUX_DIFFUSE_PROFWF (Q,K,UTA,IB,DNIDX) + ( L_DIRECT_FLUX_SCALED  - L_DIRECT_FLUX  )

                ENDDO
              ENDIF

!  End N>0 and Whole-layer clause

            ENDIF
          ENDIF

!  End loop over optical depth output values

        ENDDO

!  Finish Downwelling

      ENDIF

!  Finish

      RETURN
END SUBROUTINE MIFLUX_PROFILEWF

!

SUBROUTINE GET_LP_TOASOURCE &
     ( N_PPSTREAMS, PPSTREAM_MASK, IBEAM, K_PARAMETERS, LP_TOA_SOURCE )

!  Linearized Top of the atmosphere source term

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAX_USER_STREAMS, MAXBEAMS, MAX_ATMOSWFS, ZERO

      IMPLICIT NONE

!  Subroutine arguments
!  --------------------

!  masking (introduced, 4/9/19)

      INTEGER  , intent(in)  :: N_PPSTREAMS, PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  Beam

      INTEGER  , intent(in)  :: IBEAM

!  Linearization

      INTEGER  , intent(in)  :: K_PARAMETERS

!  Output

      REAL(fpk), intent(out) :: LP_TOA_SOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)

!  local variables
!  ---------------

      INTEGER       :: LUM, UM, Q

!  initialise TOA source function
!  ------------------------------

      DO LUM = 1, N_PPSTREAMS
         UM = PPSTREAM_MASK(LUM,IBEAM)
         DO Q = 1, K_PARAMETERS
            LP_TOA_SOURCE(UM,Q) = ZERO
         ENDDO
      ENDDO

!  Finish

END SUBROUTINE GET_LP_TOASOURCE

!

SUBROUTINE GET_LP_BOASOURCE &
      ( DO_USER_STREAMS, DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS,       & ! Input flags sources
        DO_INCLUDE_MVOUTPUT, DO_THERMAL_TRANSONLY, DO_INCLUDE_SURFACE,  & ! Input flags control
        DO_BRDF_SURFACE, DO_INCLUDE_DIRECTRF, DO_INCLUDE_DIRECTSL,      & ! Input flags inclusion
        NSTREAMS, NLAYERS, IBEAM, FOURIER, K, K_PARAMETERS,             & ! Input
        N_PPSTREAMS, PPSTREAM_MASK, QUAD_STRMWTS, QUAD_WEIGHTS,         & ! input bookkeeping
        DELTAU_SLANT, SURFACE_FACTOR, ALBEDO, BRDF_F, USER_BRDF_F,      & ! Input
        RF_USER_DIRECT_BEAM, SL_USERTERM, LP_TRANS_ATMOS_FINAL,         & ! Input surface radiance
        LCON, MCON, LCON_XVEC, T_DELT_EIGEN, T_DELT_DISORDS, T_WLOWER,  & ! Input RTS solutions
        L_DELTAU_VERT, L_XPOS, L_XNEG, L_T_WLOWER, L_WLOWER,            & ! Input Linearized solutions
        NCON_XVEC, PCON_XVEC, L_T_DELT_EIGEN, L_T_DELT_DISORDS,         & ! Input Linearized solutions
        L_BOA_MSSOURCE, L_BOA_DBSOURCE, L_BOA_THTONLY_SOURCE )            ! output

!  Linearized Bottom of the atmosphere source term

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAXMOMENTS, MAX_USER_STREAMS, MAXSTREAMS, MAXSTREAMS_2, &
                                MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS, ZERO

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  Local control flags

      LOGICAL  , intent(in)  :: DO_USER_STREAMS
      LOGICAL  , intent(in)  :: DO_SOLAR_SOURCES
      LOGICAL  , intent(in)  :: DO_INCLUDE_THERMEMISS
      
      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY
      LOGICAL  , intent(in)  :: DO_INCLUDE_SURFACE
      LOGICAL  , intent(in)  :: DO_INCLUDE_MVOUTPUT

      LOGICAL  , intent(in)  :: DO_BRDF_SURFACE
      LOGICAL  , intent(in)  :: DO_INCLUDE_DIRECTRF
      LOGICAL  , intent(in)  :: DO_INCLUDE_DIRECTSL

!  control integers

      INTEGER  , intent(in)  :: NLAYERS, NSTREAMS

!  Fourier/beam indices

      INTEGER  , intent(in)  :: FOURIER
      INTEGER  , intent(in)  :: IBEAM

!  linearization control

      INTEGER  , intent(in)  :: K, K_PARAMETERS

!  masking (introduced, 4/9/19)

      INTEGER  , intent(in)  :: N_PPSTREAMS, PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  Quadrature

      REAL(fpk), intent(in)  :: QUAD_WEIGHTS ( MAXSTREAMS )
      REAL(fpk), intent(in)  :: QUAD_STRMWTS ( MAXSTREAMS )

!  Slant optical thickness values

      REAL(fpk), intent(in)  :: DELTAU_SLANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS )

!  surface multiplier, albedo

      REAL(fpk), intent(in)  :: SURFACE_FACTOR, ALBEDO

!  Fourier components of BRDF, in the following order (same all threads)
!    ( New code, 23 March 2010 )

!  incident quadrature streams, reflected quadrature streams

      REAL(fpk), intent(in)  :: BRDF_F ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTREAMS )

!  incident quadrature streams, reflected user streams

      REAL(fpk), intent(in)  :: USER_BRDF_F ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTREAMS )

!  Reflected Direct beam and surface-leaving solutions

      REAL(fpk), intent(in)  :: RF_USER_DIRECT_BEAM ( MAX_USER_STREAMS, MAXBEAMS )
      REAL(fpk), intent(in)  :: SL_USERTERM         ( MAX_USER_STREAMS, MAXBEAMS )

!  Linearized transmittance flux for water-leaving

      REAL(fpk), intent(in)  :: LP_TRANS_ATMOS_FINAL (MAXBEAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Solution constants of integration, and related quantities

      REAL(fpk), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: LCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  transmittance factors for +/- eigenvalues

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)

!  Linearized delta taus

      REAL(fpk), intent(in)  :: L_DELTAU_VERT(MAX_ATMOSWFS, MAXLAYERS )

!  Linearized Discrete ordinate transmittances

      REAL(fpk), intent(in)  :: T_DELT_DISORDS   (MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: L_T_DELT_DISORDS (MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Thermal solutions at the Lower boundary

      REAL(fpk), intent(in)  :: T_WLOWER   (MAXSTREAMS_2,MAXLAYERS)
      REAL(fpk), intent(in)  :: L_T_WLOWER (MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Eigensolutions

      REAL(fpk), intent(in)  :: L_XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Solution constants of integration, and related quantities

      REAL(fpk), intent(in)  :: NCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: PCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized General beam solutions at the Lower boundary

      REAL(fpk), intent(in)  :: L_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized transmittance factors for +/- eigenvalues

      REAL(fpk), intent(in)  :: L_T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Subroutine output arguments
!  ---------------------------

      REAL(fpk), intent(out) :: L_BOA_MSSOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: L_BOA_DBSOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: L_BOA_THTONLY_SOURCE(MAXSTREAMS,MAX_ATMOSWFS)

!  local variables
!  ---------------

      LOGICAL    :: DO_QTHTONLY
      INTEGER    :: M, N, J, I, UM, AA, Q, IB, NN, LUM
      REAL(fpk)  :: DOWN (MAXSTREAMS)
      REAL(fpk)  :: L_DOWN (MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk)  :: REFLEC, L_BEAM, FAC, KMULT, SPAR
      REAL(fpk)  :: SHOM, HOM1, HOM2, HOM3, HOM4, HOM5

!  Starting section
!  ----------------

!  Local indices

      M   = FOURIER
      N   = NLAYERS
      IB  = IBEAM

!  Special flag

      DO_QTHTONLY = ( DO_THERMAL_TRANSONLY .and. DO_INCLUDE_MVOUTPUT )

!  initialise linearized BOA source functions

      IF ( DO_USER_STREAMS ) THEN
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IBEAM)
            L_BOA_MSSOURCE(UM,1:K_PARAMETERS) = ZERO
            L_BOA_DBSOURCE(UM,1:K_PARAMETERS) = ZERO
         ENDDO
      ENDIF

!  Thermal tranmsittance only, special term

      IF ( DO_QTHTONLY ) THEN
         DO I = 1, NSTREAMS
            L_BOA_THTONLY_SOURCE(I,1:K_PARAMETERS) = ZERO
         ENDDO
      ENDIF

!  Add reflected direct beam if flagged

      IF ( DO_INCLUDE_DIRECTRF .AND. DO_USER_STREAMS ) THEN
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            FAC = - RF_USER_DIRECT_BEAM(UM,IB)* DELTAU_SLANT(N,K,IB)
            DO Q = 1, K_PARAMETERS
               L_BOA_DBSOURCE(UM,Q) = L_BOA_DBSOURCE(UM,Q) + L_DELTAU_VERT(Q,K) * FAC
            ENDDO
         ENDDO
      ENDIF

!  Add surface-term linearization if flagged
      
      IF ( DO_INCLUDE_DIRECTSL .AND. DO_USER_STREAMS ) THEN
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(UM,IB)
            FAC = SL_USERTERM(UM,IB)
            DO Q = 1, K_PARAMETERS
               L_BOA_DBSOURCE(UM,Q) = L_BOA_DBSOURCE(UM,Q) + LP_TRANS_ATMOS_FINAL(IB,K,Q) * FAC
            ENDDO
         ENDDO
      ENDIF

!  Exit if no surface

      IF ( .NOT. DO_INCLUDE_SURFACE ) RETURN

!  1. Thermal Transmittance only
!     %%%%%%%%%%%%%%%%%%%%%%%%%%

!  Linearization of downwelling quadrature field at surface
!   ---Thermal transmittance solution, build from TOA downwards

      IF ( DO_THERMAL_TRANSONLY ) THEN

!  Initialise

        DO I = 1, NSTREAMS
          DOWN(I) = ZERO
          DO Q = 1, K_PARAMETERS
            L_DOWN(I,Q) = ZERO
          ENDDO
        ENDDO

!  Build

        DO NN = 1, NLAYERS
          IF ( K.EQ.NN ) THEN
            DO I = 1, NSTREAMS
              DO Q = 1, K_PARAMETERS
                L_DOWN(I,Q) = L_DOWN(I,Q) *   T_DELT_DISORDS(I,NN) &
                              + DOWN(I)   * L_T_DELT_DISORDS(I,NN,Q) &
                              + L_T_WLOWER(I,NN,Q)
              ENDDO
            ENDDO
          ELSE
            DO I = 1, NSTREAMS
              DO Q = 1, K_PARAMETERS
                L_DOWN(I,Q) = L_DOWN(I,Q) * T_DELT_DISORDS(I,NN) 
              ENDDO
            ENDDO
          ENDIF
          DO I = 1, NSTREAMS
            DOWN(I) = DOWN(I)*T_DELT_DISORDS(I,NN) + T_WLOWER(I,NN)
          ENDDO
        ENDDO

!  2. Scattering solutions
!     %%%%%%%%%%%%%%%%%%%%

!  Linearization of downwelling quadrature field at surface
!    Set reflectance integrand  a(j).x(j).L_DOWN(-j)

      ELSE

!  Two cases:
!  If  K = N, this is also the layer that is varying --> Extras!
!  If  N > K with variations in layer K above N
!    Add Particular integral linearization

        IF ( K.EQ.N ) THEN
          DO I = 1, NSTREAMS
            DO Q = 1, K_PARAMETERS
              SHOM = ZERO
              DO AA = 1, NSTREAMS
                HOM1 = NCON_XVEC(I,AA,N,Q) * T_DELT_EIGEN(AA,N)
                HOM2 = LCON_XVEC(I,AA,N) * L_T_DELT_EIGEN(AA,N,Q)
                HOM3 = LCON(AA,N)*T_DELT_EIGEN(AA,N)*L_XPOS(I,AA,N,Q)
                HOM4 = PCON_XVEC(I,AA,N,Q)
                HOM5 = MCON(AA,N) * L_XNEG(I,AA,N,Q)
                SHOM = SHOM + HOM1 + HOM2 + HOM3 + HOM4 + HOM5
              ENDDO
              SPAR = ZERO
              IF ( DO_SOLAR_SOURCES .OR. DO_INCLUDE_THERMEMISS ) THEN
                SPAR = L_WLOWER(I,N,Q)
              ENDIF
              L_DOWN(I,Q) = SPAR + SHOM
            ENDDO
          ENDDO
        ELSE IF (K.LT.N.AND.K.NE.0) THEN
          DO I = 1, NSTREAMS
            DO Q = 1, K_PARAMETERS
              SHOM = ZERO
              DO AA = 1, NSTREAMS
                HOM1 = NCON_XVEC(I,AA,N,Q) * T_DELT_EIGEN(AA,N)
                HOM2 = PCON_XVEC(I,AA,N,Q)
                SHOM = SHOM + HOM1 + HOM2
              ENDDO
              SPAR = ZERO
              IF ( DO_SOLAR_SOURCES ) THEN
                SPAR = L_WLOWER(I,N,Q)
              ENDIF
              L_DOWN(I,Q) = SPAR + SHOM
             ENDDO
          ENDDO
        ENDIF

      ENDIF

!    Set reflectance integrand  a(j).x(j).L_DOWN(-j)  Scattering solutions
!    Set reflectance integrand       a(j).L_DOWN(-j)  Thermal tranmsittance

      IF ( DO_THERMAL_TRANSONLY ) THEN
        DO Q = 1, K_PARAMETERS
          DO I = 1, NSTREAMS
            L_DOWN(I,Q) = L_DOWN(I,Q) * QUAD_WEIGHTS(I)
          ENDDO
        ENDDO
      ELSE
        DO Q = 1, K_PARAMETERS
          DO I = 1, NSTREAMS
            L_DOWN(I,Q) = L_DOWN(I,Q) * QUAD_STRMWTS(I)
          ENDDO
        ENDDO
      ENDIF

!  reflected multiple scatter intensity at user defined-angles
!  -----------------------------------------------------------

!  .. integrate reflectance, same for all user-streams in Lambertian case

      IF ( .not. DO_BRDF_SURFACE ) THEN

         KMULT = SURFACE_FACTOR * ALBEDO
         IF ( FOURIER .EQ. 0 ) THEN
           DO Q = 1, K_PARAMETERS
               REFLEC = KMULT * SUM(L_DOWN(1:NSTREAMS,Q))
               IF ( DO_USER_STREAMS ) THEN
                  DO LUM = 1, N_PPSTREAMS
                     UM = PPSTREAM_MASK(LUM,IBEAM)
                     L_BOA_MSSOURCE(UM,Q) = L_BOA_MSSOURCE(UM,Q) + REFLEC
                  ENDDO
               ENDIF
               IF ( DO_QTHTONLY ) THEN
                  DO I = 1, NSTREAMS
                     L_BOA_THTONLY_SOURCE(I,Q) = L_BOA_THTONLY_SOURCE(I,Q) + REFLEC
                  ENDDO
               ENDIF
            ENDDO
         ENDIF

!  .. integrate reflectance, BRDF case

      ELSE IF ( DO_BRDF_SURFACE ) THEN
         DO Q = 1, K_PARAMETERS
            IF ( DO_USER_STREAMS ) THEN
               DO LUM = 1, N_PPSTREAMS
                  UM = PPSTREAM_MASK(LUM,IBEAM)
                  REFLEC = SURFACE_FACTOR * DOT_PRODUCT(L_DOWN(1:NSTREAMS,Q),USER_BRDF_F(M,UM,1:NSTREAMS))
                  L_BOA_MSSOURCE(UM,Q) = L_BOA_MSSOURCE(UM,Q) + REFLEC
               ENDDO
            ENDIF
            IF ( DO_QTHTONLY ) THEN
               DO I = 1, NSTREAMS
                  REFLEC = SURFACE_FACTOR * DOT_PRODUCT(L_DOWN(1:NSTREAMS,Q),BRDF_F(M,I,1:NSTREAMS))
                  L_BOA_THTONLY_SOURCE(I,Q) = L_BOA_THTONLY_SOURCE(I,Q) + REFLEC
               ENDDO
            ENDIF
         ENDDO
      ENDIF

!  Finish

      RETURN
END SUBROUTINE GET_LP_BOASOURCE

!

SUBROUTINE LIDORT_LP_CONVERGE &
      ( DO_FOCORR, DO_FOCORR_ALONE,                      & ! Input flags
        DO_UPWELLING, DO_NO_AZIMUTH,                     & ! Input flags
        N_USER_STREAMS, N_USER_LEVELS, LOCAL_N_USERAZM,  & ! Input numbers
        IBEAM, FOURIER, NLAYERS, LV_FLAG, LV_NUMBER,     & ! Input numbers
        UMOFF, N_DIRECTIONS, WHICH_DIRECTIONS,           & ! Input bookkeeping
        AZMFAC, PROFILEWF_F, PROFILEWF_SS, PROFILEWF_DB, & ! Input Azm, fields
        PROFILEWF )                                        ! Output

!  Just upgrades the weighting function Fourier cosine series
!     Version 3.8, 3/1/17. Logic for FOCORR variables changed
!     Version 3.8, 3/1/17. Argument list rearranged

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_USER_LEVELS, &
                                MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS, MAX_GEOMETRIES,   &
                                MAX_DIRECTIONS, ZERO, UPIDX

      IMPLICIT NONE

!  input variables
!  ---------------

!  FO flags

!      LOGICAL  , intent(in)  :: DO_SSFULL
!      LOGICAL  , intent(in)  :: DO_SSCORR_NADIR
!      LOGICAL  , intent(in)  :: DO_SSCORR_OUTGOING
!      LOGICAL  , intent(in)  :: DO_SS_EXTERNAL

      LOGICAL  , intent(in)  :: DO_FOCORR
      LOGICAL  , intent(in)  :: DO_FOCORR_ALONE

!  Other Flags

      LOGICAL  , intent(in)  :: DO_UPWELLING
      LOGICAL  , intent(in)  :: DO_NO_AZIMUTH

!  Control integers, Local number of azimuths

      INTEGER  , intent(in)  :: NLAYERS
      INTEGER  , intent(in)  :: N_USER_LEVELS
      INTEGER  , intent(in)  :: N_USER_STREAMS
      INTEGER  , intent(in)  :: LOCAL_N_USERAZM

!  Fourier component and beam

      INTEGER  , intent(in)  :: FOURIER, IBEAM

!  Linearization control

      LOGICAL  , intent(in)  :: LV_FLAG   ( MAXLAYERS )
      INTEGER  , intent(in)  :: LV_NUMBER ( MAXLAYERS )

!  Bookkeeping: Offsets for geometry indexing

      INTEGER  , intent(in)  :: UMOFF(MAXBEAMS,MAX_USER_STREAMS)

!  Directional control

      INTEGER  , intent(in)  :: N_DIRECTIONS
      INTEGER  , intent(in)  :: WHICH_DIRECTIONS(2)

!  Azimuth factors

      REAL(fpk), intent(in)  :: AZMFAC (MAX_USER_STREAMS,MAXBEAMS,MAX_USER_RELAZMS)

!  Fourier-component Profile weighting functions at user angles

      REAL(fpk), intent(in)  :: PROFILEWF_F ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
                                              MAX_USER_STREAMS,MAXBEAMS,MAX_DIRECTIONS)

!  Single scatter Profile weighting functions at user angles

      REAL(fpk), intent(in)  :: PROFILEWF_SS ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
                                               MAX_GEOMETRIES, MAX_DIRECTIONS )

!  Direct-bounce Profile weighting functions at user angles

      REAL(fpk), intent(in)  :: PROFILEWF_DB ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
                                               MAX_GEOMETRIES )

!  output
!  ------

!mick fix 6/29/11 - changed output from "out" to "inout"

      REAL(fpk), intent(inout) :: PROFILEWF ( MAX_ATMOSWFS,   MAXLAYERS,  MAX_USER_LEVELS, &
                                              MAX_GEOMETRIES, MAX_DIRECTIONS )

!  local variables

      INTEGER       :: I, IDIR, UT, UA, Q, W, V, N

!  ###################
!  Fourier 0 component
!  ###################

      IF ( FOURIER.EQ.0 ) THEN

!  Copy DIFFUSE Fourier component at all output angles and optical depths
!    If no SSCORR and no DBCORR, then two options apply:
!     (a) Convergence on JACOBIAN = DIFFUSE + SSTRUNCATED + DBTRUNCATED
!              (full radiance, no SS correction, no DB correction)
!     (b) Convergence on JACOBIAN = DIFFUSE alone (MS only mode)
!              (SSTRUNCATED + DBTRUNCATED do not get calculated)

        IF ( .NOT. DO_FOCORR_ALONE ) THEN
         DO N = 1, NLAYERS
          IF ( LV_FLAG(N) ) THEN
           DO Q = 1, LV_NUMBER(N)
            DO IDIR = 1, N_DIRECTIONS
              W = WHICH_DIRECTIONS(IDIR)
              DO UT = 1, N_USER_LEVELS
                DO I = 1, N_USER_STREAMS
                  DO UA = 1, LOCAL_N_USERAZM
                    V = UMOFF(IBEAM,I) + UA
                    PROFILEWF(Q,N,UT,V,W) = PROFILEWF_F(Q,N,UT,I,IBEAM,W)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
           ENDDO
          ENDIF
         ENDDO
        ELSE
         DO N = 1, NLAYERS
          IF ( LV_FLAG(N) ) THEN
           DO Q = 1, LV_NUMBER(N)
            DO IDIR = 1, N_DIRECTIONS
              W = WHICH_DIRECTIONS(IDIR)
              DO UT = 1, N_USER_LEVELS
                DO I = 1, N_USER_STREAMS
                  DO UA = 1, LOCAL_N_USERAZM
                    V = UMOFF(IBEAM,I) + UA
                    PROFILEWF(Q,N,UT,V,W) = ZERO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
           ENDDO
          ENDIF
         ENDDO
        ENDIF

!    Add the single scatter component if flagged
!     If set, now looking at convergence on RADIANCE = DIFFUSE + SSEXACT

!     Version 3.2.   Added outgoing correction flag to this.....
!     Version 3.3    Added Full single scatter flag
!     New 15 March 2012, Introduced DO_SS_EXTERNAL flag
!     Version 3.8    Much simpler condition.

        !IF ( DO_SSFULL.OR.DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING ) THEN
!        IF ( DO_SSFULL .OR. &
!             ((DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING).OR.DO_SS_EXTERNAL) ) THEN

        IF ( DO_FOCORR ) THEN
         DO N = 1, NLAYERS
          IF ( LV_FLAG(N) ) THEN
           DO Q = 1, LV_NUMBER(N)
            DO IDIR = 1, N_DIRECTIONS
              W = WHICH_DIRECTIONS(IDIR)
              DO UT = 1, N_USER_LEVELS
                DO I = 1, N_USER_STREAMS
                  DO UA = 1, LOCAL_N_USERAZM
                    V = UMOFF(IBEAM,I) + UA
                    !CALL TP23E (FOURIER,Q,N,UT,V,W,PROFILEWF,PROFILEWF_SS)
                    PROFILEWF(Q,N,UT,V,W) = &
                      PROFILEWF(Q,N,UT,V,W) + PROFILEWF_SS(Q,N,UT,V,W)
                 ENDDO
                ENDDO
              ENDDO
            ENDDO
           ENDDO
          ENDIF
         ENDDO
        ENDIF

!  Add the Direct bounce to the upwelling
!     New 15 March 2012, Introduced DO_SS_EXTERNAL flag
!     Version 3.8        Changed Logic for SS terms

!        IF ( DO_UPWELLING ) THEN
         !!IF ( DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING ) THEN
         !IF ((DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING).OR.DO_SS_EXTERNAL) THEN

        IF ( DO_FOCORR .AND.DO_UPWELLING ) THEN
          DO N = 1, NLAYERS
            IF ( LV_FLAG(N) ) THEN
              DO Q = 1, LV_NUMBER(N)
                DO UT = 1, N_USER_LEVELS
                  DO I = 1, N_USER_STREAMS
                    DO UA = 1, LOCAL_N_USERAZM
                      V = UMOFF(IBEAM,I) + UA
                      !CALL TP23F (FOURIER,Q,N,UT,V,PROFILEWF,PROFILEWF_DB)
                      PROFILEWF(Q,N,UT,V,UPIDX) = &
                        PROFILEWF(Q,N,UT,V,UPIDX) + PROFILEWF_DB(Q,N,UT,V)
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
          ENDDO
        ENDIF

!  If no_azimuth, then exit

        IF ( DO_NO_AZIMUTH ) RETURN

!  ######################
!  Fourier components > 0
!  ######################

      ELSE

!  Add next Fourier component to output

        DO N = 1, NLAYERS
         IF ( LV_FLAG(N) ) THEN
          DO Q = 1, LV_NUMBER(N)
           DO UA = 1, LOCAL_N_USERAZM
            DO IDIR = 1, N_DIRECTIONS
             W = WHICH_DIRECTIONS(IDIR)
             DO UT = 1, N_USER_LEVELS
              DO I = 1, N_USER_STREAMS
               V = UMOFF(IBEAM,I) + UA
               !CALL TP23G (FOURIER,Q,N,UT,I,IBEAM,V,W,PROFILEWF,PROFILEWF_F)
               PROFILEWF(Q,N,UT,V,W) = PROFILEWF(Q,N,UT,V,W) + &
                     PROFILEWF_F(Q,N,UT,I,IBEAM,W)*AZMFAC(I,IBEAM,UA)
              ENDDO
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDIF
        ENDDO

!  end Fourier clause

      ENDIF

!  Finish

      RETURN
END SUBROUTINE LIDORT_LP_CONVERGE

!

SUBROUTINE LIDORT_LP_CONVERGE_OBSGEO &
      ( DO_FOCORR, DO_FOCORR_ALONE,                       & ! Input flags
        DO_UPWELLING, DO_NO_AZIMUTH,                      & ! Input flags
        N_USER_LEVELS, IBEAM, FOURIER, NLAYERS, LV_FLAG,  & ! Input numbers
        LV_NUMBER, N_DIRECTIONS, WHICH_DIRECTIONS,        & ! Input bookkeeping
        AZMFAC, PROFILEWF_F, PROFILEWF_SS, PROFILEWF_DB,  & ! Input Azm, fields
        PROFILEWF )                                         ! Output

!  Just upgrades the weighting function Fourier cosine series
!     Version 3.8, 3/1/17. Logic for FOCORR variables changed
!     Version 3.8, 3/1/17. Argument list rearranged

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_USER_LEVELS, &
                                MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS, MAX_GEOMETRIES,   &
                                MAX_DIRECTIONS, ZERO, UPIDX

      IMPLICIT NONE

!  input variables
!  ---------------

!  FO flags

!      LOGICAL  , intent(in)  :: DO_SSFULL
!      LOGICAL  , intent(in)  :: DO_SSCORR_NADIR
!      LOGICAL  , intent(in)  :: DO_SSCORR_OUTGOING
!      LOGICAL  , intent(in)  :: DO_SS_EXTERNAL

      LOGICAL  , intent(in)  :: DO_FOCORR
      LOGICAL  , intent(in)  :: DO_FOCORR_ALONE

!  Other Flags

      LOGICAL  , intent(in)  :: DO_UPWELLING
      LOGICAL  , intent(in)  :: DO_NO_AZIMUTH

!  Control integers

      INTEGER  , intent(in)  :: N_USER_LEVELS
      INTEGER  , intent(in)  :: NLAYERS

!  Fourier component and beam

      INTEGER  , intent(in)  :: FOURIER, IBEAM

!  Linearization control

      LOGICAL  , intent(in)  :: LV_FLAG   ( MAXLAYERS )
      INTEGER  , intent(in)  :: LV_NUMBER ( MAXLAYERS )

!  Directional control

      INTEGER  , intent(in)  :: N_DIRECTIONS
      INTEGER  , intent(in)  :: WHICH_DIRECTIONS(2)

!  Azimuth factors

      REAL(fpk), intent(in)  :: AZMFAC (MAX_USER_STREAMS,MAXBEAMS,MAX_USER_RELAZMS)

!  Fourier-component Profile weighting functions at user angles

      REAL(fpk), intent(in)  :: PROFILEWF_F ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
                                              MAX_USER_STREAMS,MAXBEAMS,MAX_DIRECTIONS)

!  Single scatter Profile weighting functions at user angles

      REAL(fpk), intent(in)  :: PROFILEWF_SS ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
                                               MAX_GEOMETRIES, MAX_DIRECTIONS )

!  Direct-bounce Profile weighting functions at user angles

      REAL(fpk), intent(in)  :: PROFILEWF_DB ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
                                               MAX_GEOMETRIES )

!  output
!  ------

!mick fix 6/29/11 - changed output from "out" to "inout"

      REAL(fpk), intent(inout) :: PROFILEWF ( MAX_ATMOSWFS,   MAXLAYERS,  MAX_USER_LEVELS, &
                                              MAX_GEOMETRIES, MAX_DIRECTIONS )

!  local variables

      INTEGER       :: IDIR, UT, Q, W, N, LUM, LUA

!  Local user indices

      LUM = 1
      LUA = 1

!  ###################
!  Fourier 0 component
!  ###################

      IF ( FOURIER.EQ.0 ) THEN

!  Copy DIFFUSE Fourier component at all output angles and optical depths
!    If no SSCORR and no DBCORR, then two options apply:
!     (a) Convergence on JACOBIAN = DIFFUSE + SSTRUNCATED + DBTRUNCATED
!              (full radiance, no SS correction, no DB correction)
!     (b) Convergence on JACOBIAN = DIFFUSE alone (MS only mode)
!              (SSTRUNCATED + DBTRUNCATED do not get calculated)

        IF ( .NOT. DO_FOCORR_ALONE ) THEN
         DO N = 1, NLAYERS
          IF ( LV_FLAG(N) ) THEN
           DO Q = 1, LV_NUMBER(N)
            DO IDIR = 1, N_DIRECTIONS
              W = WHICH_DIRECTIONS(IDIR)
              DO UT = 1, N_USER_LEVELS
                PROFILEWF(Q,N,UT,IBEAM,W) = PROFILEWF_F(Q,N,UT,LUM,IBEAM,W)
              ENDDO
            ENDDO
           ENDDO
          ENDIF
         ENDDO
        ELSE
         DO N = 1, NLAYERS
          IF ( LV_FLAG(N) ) THEN
           DO Q = 1, LV_NUMBER(N)
            DO IDIR = 1, N_DIRECTIONS
              W = WHICH_DIRECTIONS(IDIR)
              DO UT = 1, N_USER_LEVELS
                PROFILEWF(Q,N,UT,IBEAM,W) = ZERO
              ENDDO
            ENDDO
           ENDDO
          ENDIF
         ENDDO
        ENDIF

!    Add the single scatter component if flagged
!     If set, now looking at convergence on RADIANCE = DIFFUSE + SSEXACT

!     Version 3.2.   Added outgoing correction flag to this.....
!     Version 3.3    Added Full single scatter flag
!     New 15 March 2012, Introduced DO_SS_EXTERNAL flag
!     Version 3.8    Much simpler condition.

        !IF ( DO_SSFULL.OR.DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING ) THEN
!        IF ( DO_SSFULL .OR. &
!             ((DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING).OR.DO_SS_EXTERNAL) ) THEN

        IF ( DO_FOCORR ) THEN
         DO N = 1, NLAYERS
          IF ( LV_FLAG(N) ) THEN
           DO Q = 1, LV_NUMBER(N)
            DO IDIR = 1, N_DIRECTIONS
              W = WHICH_DIRECTIONS(IDIR)
              DO UT = 1, N_USER_LEVELS
                !CALL TP23E2 (FOURIER,Q,N,UT,IBEAM,W,PROFILEWF,PROFILEWF_SS)
                PROFILEWF(Q,N,UT,IBEAM,W) = &
                  PROFILEWF(Q,N,UT,IBEAM,W) + PROFILEWF_SS(Q,N,UT,IBEAM,W)
              ENDDO
            ENDDO
           ENDDO
          ENDIF
         ENDDO
        ENDIF

!  Add the Direct bounce to the upwelling
!     New 15 March 2012, Introduced DO_SS_EXTERNAL flag
!     Version 3.8        Changed Logic for SS terms

!        IF ( DO_UPWELLING ) THEN
         !!IF ( DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING ) THEN
         !IF ((DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING).OR.DO_SS_EXTERNAL) THEN

        IF ( DO_FOCORR .AND.DO_UPWELLING ) THEN
          DO N = 1, NLAYERS
            IF ( LV_FLAG(N) ) THEN
              DO Q = 1, LV_NUMBER(N)
                DO UT = 1, N_USER_LEVELS
                  !CALL TP23F2 (FOURIER,Q,N,UT,IBEAM,PROFILEWF,PROFILEWF_DB)
                  PROFILEWF(Q,N,UT,IBEAM,UPIDX) = &
                    PROFILEWF(Q,N,UT,IBEAM,UPIDX) + PROFILEWF_DB(Q,N,UT,IBEAM)
                ENDDO
              ENDDO
            ENDIF
          ENDDO
        ENDIF

!  If no_azimuth, then exit

        IF ( DO_NO_AZIMUTH ) RETURN

!  ######################
!  Fourier components > 0
!  ######################

      ELSE

!  Add next Fourier component to output

        DO N = 1, NLAYERS
         IF ( LV_FLAG(N) ) THEN
          DO Q = 1, LV_NUMBER(N)
           DO IDIR = 1, N_DIRECTIONS
            W = WHICH_DIRECTIONS(IDIR)
            DO UT = 1, N_USER_LEVELS
              !CALL TP23G2 (FOURIER,Q,N,UT,LUM,IBEAM,W,PROFILEWF,PROFILEWF_F)
              PROFILEWF(Q,N,UT,IBEAM,W) = PROFILEWF(Q,N,UT,IBEAM,W) + &
                    PROFILEWF_F(Q,N,UT,LUM,IBEAM,W)*AZMFAC(LUM,IBEAM,LUA)
            ENDDO
           ENDDO
          ENDDO
         ENDIF
        ENDDO

!  end Fourier clause

      ENDIF

!  Finish

      RETURN
END SUBROUTINE LIDORT_LP_CONVERGE_OBSGEO

!  End

end module lidort_lp_wfatmos_m
