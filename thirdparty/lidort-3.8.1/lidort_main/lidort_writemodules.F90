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
! #              LIDORT_WRITE_STD_INPUT                         #
! #              LIDORT_WRITE_SUP_BRDF_INPUT                    #
! #              LIDORT_WRITE_SUP_SS_INPUT                      #
! #              LIDORT_WRITE_SUP_SLEAVE_INPUT                  #
! #                                                             #
! ###############################################################

!  Upgrade, Version 3.8.1, June 2019. Write new STD input variables
!    -- (TOA/BOA illumination, Planetary, WLEAVING

MODULE lidort_writemodules_m

      PRIVATE
      PUBLIC :: LIDORT_WRITE_STD_INPUT, &
                LIDORT_WRITE_SUP_BRDF_INPUT, &
                LIDORT_WRITE_SUP_SS_INPUT, &
                LIDORT_WRITE_SUP_SLEAVE_INPUT

      CONTAINS

      SUBROUTINE LIDORT_WRITE_STD_INPUT ( &
        DO_SOLAR_SOURCES,   DO_THERMAL_EMISSION, DO_THERMAL_TRANSONLY, DO_SURFACE_EMISSION,     & ! Sources
        DO_FULLRAD_MODE,    DO_FOCORR,           DO_FOCORR_EXTERNAL,                            & ! FOCORR
        DO_FOCORR_NADIR,    DO_FOCORR_OUTGOING,  DO_SSCORR_TRUNCATION, DO_SSCORR_USEPHASFUNC,   & ! FOCORR
        DO_RAYLEIGH_ONLY,   DO_ISOTROPIC_ONLY,   DO_NO_AZIMUTH,        DO_ALL_FOURIER,          & ! Performance
        DO_DELTAM_SCALING,  DO_DOUBLE_CONVTEST,  DO_SOLUTION_SAVING,   DO_BVP_TELESCOPING,      & ! Performance
        DO_UPWELLING,       DO_DNWELLING,        DO_PLANE_PARALLEL,    DO_CHAPMAN_FUNCTION,     & ! RT Model
        DO_USER_STREAMS,    DO_OBSERV_GEOMETRY,  DO_ADDITIONAL_MVOUT,  DO_MVOUT_ONLY,           & ! RT model
        DO_REFRAC_GEOMETRY, DO_ALBTRN_MEDIA,     DO_PLANETARY_PROBLEM, DO_TOAFLUX, DO_BOAFLUX,  & ! Media/Planetary/Flux
        DO_BRDF_SURFACE,    DO_SURFACE_LEAVING,  DO_SL_ISOTROPIC,      DO_EXTERNAL_WLEAVE,      & ! Surface
        DO_WATER_LEAVING,   DO_FLUORESCENCE,     DO_WLADJUSTED_OUTPUT, DO_TF_ITERATION,         & ! Surface/Media-props
        TAYLOR_ORDER, TF_MAXITER, NSTREAMS, NLAYERS, NFINELAYERS, N_THERMAL_COEFFS, NMOMENTS_INPUT, & ! Main numbers
        NBEAMS, N_USER_RELAZMS, N_USER_STREAMS, N_USER_OBSGEOMS, N_USER_LEVELS,                 & ! Geometry and level numbers
        BEAM_SZAS, USER_RELAZMS, USER_ANGLES_INPUT, USER_OBSGEOMS, USER_LEVELS,                 & ! Geometry and level control
        LIDORT_ACCURACY, FLUX_FACTOR, EARTH_RADIUS, RFINDEX_PARAMETER, TF_CRITERION,            & ! Flux/Acc/Radius
        HEIGHT_GRID, PRESSURE_GRID, TEMPERATURE_GRID, FINEGRID, GEOMETRY_SPECHEIGHT,            & ! Grids
        DELTAU_VERT_INPUT, OMEGA_TOTAL_INPUT, PHASMOMS_TOTAL_INPUT,                             & ! Optical
        PHASFUNC_INPUT_UP, PHASFUNC_INPUT_DN, LAMBERTIAN_ALBEDO, TOAFLUX, BOAFLUX,              & ! Optical & albedo & Fluxes
        THERMAL_BB_INPUT, SURFACE_BB_INPUT, ATMOS_WAVELENGTH )                                    ! Thermal & spectral

      USE LIDORT_PARS_m, Only : fpk, MAXMOMENTS_INPUT, MAXLAYERS, MAXBEAMS, MAX_USER_RELAZMS,         &
                                MAX_USER_STREAMS, MAX_USER_LEVELS, MAX_USER_OBSGEOMS, MAX_GEOMETRIES, &
                                DWFI, DWFI1, DWFL, DWFR, DWFR1, DWFR2, DWFR1_3

      IMPLICIT NONE
      
!  4/26/19. Record the Media-problem inputs (DO_ALBTRN_MEDIA,DO_PLANETARY_PROBLEM) . Version 3.8.1.
!  4/26/19  Additional Control for SLEAVE (DO_WLADJUSTED_OUTPUT,DO_EXTERNAL_WLEAVE). Version 3.8.1.
!  4/22/19  Additional Control for TO and BOA illumination options.                  Version 3.8.1.

!  -----------------------
!  Standard Inputs - Fixed
!  -----------------------

      LOGICAL, INTENT(IN) :: DO_FULLRAD_MODE

      LOGICAL, INTENT(IN) :: DO_THERMAL_EMISSION
      LOGICAL, INTENT(IN) :: DO_SURFACE_EMISSION
      LOGICAL, INTENT(IN) :: DO_PLANE_PARALLEL

      LOGICAL, INTENT(IN) :: DO_UPWELLING
      LOGICAL, INTENT(IN) :: DO_DNWELLING

      LOGICAL, INTENT(IN) :: DO_BRDF_SURFACE
      LOGICAL, INTENT(IN) :: DO_SURFACE_LEAVING
      LOGICAL, INTENT(IN) :: DO_SL_ISOTROPIC
      LOGICAL, INTENT(IN) :: DO_WATER_LEAVING
      LOGICAL, INTENT(IN) :: DO_FLUORESCENCE
      LOGICAL, INTENT(IN) :: DO_TF_ITERATION

      LOGICAL, INTENT(IN) :: DO_TOAFLUX ! New 4/22/19
      LOGICAL, INTENT(IN) :: DO_BOAFLUX ! New 4/22/19

      INTEGER, INTENT(IN) :: TF_MAXITER
      INTEGER, INTENT(IN) :: TAYLOR_ORDER
      INTEGER, INTENT(IN) :: NSTREAMS
      INTEGER, INTENT(IN) :: NLAYERS
      INTEGER, INTENT(IN) :: NFINELAYERS
      INTEGER, INTENT(IN) :: N_THERMAL_COEFFS
      INTEGER, INTENT(IN) :: N_USER_LEVELS

      REAL(fpk), INTENT(IN) :: LIDORT_ACCURACY

      REAL(fpk), INTENT(IN) :: FLUX_FACTOR

      REAL(fpk), INTENT(IN) :: HEIGHT_GRID ( 0:MAXLAYERS )
      REAL(fpk), INTENT(IN) :: PRESSURE_GRID ( 0:MAXLAYERS )
      REAL(fpk), INTENT(IN) :: TEMPERATURE_GRID ( 0:MAXLAYERS )
      INTEGER  , INTENT(IN) :: FINEGRID ( MAXLAYERS )
      REAL(fpk), INTENT(IN) :: RFINDEX_PARAMETER
      REAL(fpk), INTENT(IN) :: TF_CRITERION

      REAL(fpk), INTENT(IN) :: DELTAU_VERT_INPUT ( MAXLAYERS )
      REAL(fpk), INTENT(IN) :: PHASMOMS_TOTAL_INPUT ( 0:MAXMOMENTS_INPUT, MAXLAYERS )
      REAL(fpk), INTENT(IN) :: PHASFUNC_INPUT_UP ( MAXLAYERS, MAX_GEOMETRIES )
      REAL(fpk), INTENT(IN) :: PHASFUNC_INPUT_DN ( MAXLAYERS, MAX_GEOMETRIES )

      REAL(fpk), INTENT(IN) :: LAMBERTIAN_ALBEDO
      REAL(fpk), INTENT(IN) :: THERMAL_BB_INPUT ( 0:MAXLAYERS )
      REAL(fpk), INTENT(IN) :: SURFACE_BB_INPUT
      REAL(fpk), INTENT(IN) :: ATMOS_WAVELENGTH

!  TOA/BOA Fluxes (new 4/22/19). Version 3.8.1

      REAL(fpk), INTENT (IN) :: TOAFLUX
      REAL(fpk), INTENT (IN) :: BOAFLUX

!  These two added for Version 3.8.1, 4/26/19, 4/28/19

      LOGICAL  , INTENT(IN) :: DO_ALBTRN_MEDIA(2)
      LOGICAL  , INTENT(IN) :: DO_PLANETARY_PROBLEM
      
!  New flags for Version 3.8.1. Introduced 4/22/19.
      
      LOGICAL, INTENT(IN) ::          DO_WLADJUSTED_OUTPUT      
      LOGICAL, INTENT(IN) ::          DO_EXTERNAL_WLEAVE

!  --------------------------
!  Standard Inputs - Modified
!  --------------------------
!mick mod 3/22/2017 - DO_FOCORR_ALONE now defined internally

      LOGICAL, INTENT(IN) :: DO_FOCORR
      LOGICAL, INTENT(IN) :: DO_FOCORR_EXTERNAL
      !LOGICAL, INTENT(IN) :: DO_FOCORR_ALONE
      LOGICAL, INTENT(IN) :: DO_FOCORR_NADIR
      LOGICAL, INTENT(IN) :: DO_FOCORR_OUTGOING

      LOGICAL, INTENT(IN) :: DO_SSCORR_TRUNCATION
      LOGICAL, INTENT(IN) :: DO_SSCORR_USEPHASFUNC

      LOGICAL, INTENT(IN) :: DO_SOLAR_SOURCES
      LOGICAL, INTENT(IN) :: DO_REFRAC_GEOMETRY
      LOGICAL, INTENT(IN) :: DO_CHAPMAN_FUNCTION

      LOGICAL, INTENT(IN) :: DO_RAYLEIGH_ONLY
      LOGICAL, INTENT(IN) :: DO_ISOTROPIC_ONLY
      LOGICAL, INTENT(IN) :: DO_NO_AZIMUTH
      LOGICAL, INTENT(IN) :: DO_ALL_FOURIER

      LOGICAL, INTENT(IN) :: DO_DOUBLE_CONVTEST
      LOGICAL, INTENT(IN) :: DO_DELTAM_SCALING
      LOGICAL, INTENT(IN) :: DO_SOLUTION_SAVING
      LOGICAL, INTENT(IN) :: DO_BVP_TELESCOPING

      LOGICAL, INTENT(IN) :: DO_USER_STREAMS
      LOGICAL, INTENT(IN) :: DO_OBSERV_GEOMETRY

      LOGICAL, INTENT(IN) :: DO_ADDITIONAL_MVOUT
      LOGICAL, INTENT(IN) :: DO_MVOUT_ONLY

      LOGICAL, INTENT(IN) :: DO_THERMAL_TRANSONLY

      INTEGER, INTENT(IN) :: NMOMENTS_INPUT

      INTEGER  , INTENT(IN) :: NBEAMS
      REAL(fpk), INTENT(IN) :: BEAM_SZAS ( MAXBEAMS )

      INTEGER  , INTENT(IN) :: N_USER_RELAZMS
      REAL(fpk), INTENT(IN) :: USER_RELAZMS  ( MAX_USER_RELAZMS )

      INTEGER  , INTENT(IN) :: N_USER_STREAMS
      REAL(fpk), INTENT(IN) :: USER_ANGLES_INPUT ( MAX_USER_STREAMS )

      INTEGER  , INTENT(IN) :: N_USER_OBSGEOMS
      REAL(fpk), INTENT(IN) :: USER_OBSGEOMS (MAX_USER_OBSGEOMS,3)

      REAL(fpk), INTENT(IN) :: USER_LEVELS ( MAX_USER_LEVELS )
      REAL(fpk), INTENT(IN) :: GEOMETRY_SPECHEIGHT

      !REAL(fpk), INTENT(IN) :: CHAPMAN_FACTORS  ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      REAL(fpk), INTENT(IN) :: EARTH_RADIUS

      REAL(fpk), INTENT(IN) :: OMEGA_TOTAL_INPUT ( MAXLAYERS )

!  Local variables

      INTEGER :: OUTUNIT, NGEOMS
      INTEGER :: GEO,SZA,ULEV,LAY,MOM,URA,UVA,UOG

      NGEOMS = NBEAMS
      IF ( .NOT.DO_OBSERV_GEOMETRY ) NGEOMS = NBEAMS * N_USER_STREAMS * N_USER_RELAZMS

!  Open output file

      OUTUNIT = 101
      OPEN (OUTUNIT,file = 'LIDORT_WRITE_STD_INPUT.dbg',status = 'replace')

!  Write all input to file

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,'(A)') '-----------------------'
      WRITE(OUTUNIT,'(A)') 'Standard Inputs - Fixed'
      WRITE(OUTUNIT,'(A)') '-----------------------'

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFL)  'DO_FULLRAD_MODE        = ',DO_FULLRAD_MODE
      WRITE(OUTUNIT,DWFL)  'DO_THERMAL_EMISSION    = ',DO_THERMAL_EMISSION
      WRITE(OUTUNIT,DWFL)  'DO_SURFACE_EMISSION    = ',DO_SURFACE_EMISSION
      WRITE(OUTUNIT,DWFL)  'DO_PLANE_PARALLEL      = ',DO_PLANE_PARALLEL
      WRITE(OUTUNIT,DWFL)  'DO_BRDF_SURFACE        = ',DO_BRDF_SURFACE
      WRITE(OUTUNIT,DWFL)  'DO_UPWELLING           = ',DO_UPWELLING
      WRITE(OUTUNIT,DWFL)  'DO_DNWELLING           = ',DO_DNWELLING
      WRITE(OUTUNIT,DWFL)  'DO_SURFACE_LEAVING     = ',DO_SURFACE_LEAVING
      WRITE(OUTUNIT,DWFL)  'DO_SL_ISOTROPIC        = ',DO_SL_ISOTROPIC
      WRITE(OUTUNIT,DWFL)  'DO_WATER_LEAVING       = ',DO_WATER_LEAVING
      WRITE(OUTUNIT,DWFL)  'DO_FLUORESCENCE        = ',DO_FLUORESCENCE
      WRITE(OUTUNIT,DWFL)  'DO_TF_ITERATION        = ',DO_TF_ITERATION
      WRITE(OUTUNIT,DWFL)  'DO_WLADJUSTED_OUTPUT   = ',DO_WLADJUSTED_OUTPUT   !    Introduced 4/22/19 for Version 3.8.1
      WRITE(OUTUNIT,DWFL)  'DO_TOA_ILLUMINATION    = ',DO_TOAFLUX             !    Introduced 4/22/19 for Version 3.8.1
      WRITE(OUTUNIT,DWFL)  'DO_BOA_ILLUMINATION    = ',DO_BOAFLUX             !    Introduced 4/22/19 for Version 3.8.1

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFI)  'TAYLOR_ORDER     = ',TAYLOR_ORDER
      WRITE(OUTUNIT,DWFI)  'NSTREAMS         = ',NSTREAMS
      WRITE(OUTUNIT,DWFI)  'NLAYERS          = ',NLAYERS
      WRITE(OUTUNIT,DWFI)  'NFINELAYERS      = ',NFINELAYERS
      WRITE(OUTUNIT,DWFI)  'N_THERMAL_COEFFS = ',N_THERMAL_COEFFS

      WRITE(OUTUNIT,DWFR)  'LIDORT_ACCURACY  = ',LIDORT_ACCURACY

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFI)  'TF_MAXITER       = ',TF_MAXITER
      WRITE(OUTUNIT,DWFR)  'TF_CRITERION     = ',TF_CRITERION
      
!  4/26/19. Record the Media-problem and Planetary-problem inputs

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFL)  'DO_ALBTRN_MEDIA(1)   = ',DO_ALBTRN_MEDIA(1)
      WRITE(OUTUNIT,DWFL)  'DO_ALBTRN_MEDIA(2)   = ',DO_ALBTRN_MEDIA(2)
      WRITE(OUTUNIT,DWFL)  'DO_PLANETARY_PROBLEM = ',DO_PLANETARY_PROBLEM

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFR)  'FLUX_FACTOR      = ',FLUX_FACTOR

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFR)  'TOA_ILLUMINATION = ',TOAFLUX          !    Introduced 4/22/19 for Version 3.8.1
      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFR)  'BOA_ILLUMINATION = ',BOAFLUX          !    Introduced 4/22/19 for Version 3.8.1

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFI)  'N_USER_LEVELS    = ',N_USER_LEVELS

      WRITE(OUTUNIT,*)
      DO LAY=0,NLAYERS
        WRITE(OUTUNIT,DWFR1)  'LAY = ',LAY,' HEIGHT_GRID(LAY)      = ',HEIGHT_GRID(LAY)
      END DO

      WRITE(OUTUNIT,*)
      DO LAY=0,NLAYERS
        WRITE(OUTUNIT,DWFR1)  'LAY = ',LAY,' PRESSURE_GRID(LAY)    = ',PRESSURE_GRID(LAY)
      END DO

      WRITE(OUTUNIT,*)
      DO LAY=0,NLAYERS
        WRITE(OUTUNIT,DWFR1)  'LAY = ',LAY,' TEMPERATURE_GRID(LAY) = ',TEMPERATURE_GRID(LAY)
      END DO

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        WRITE(OUTUNIT,DWFI1)  'LAY = ',LAY,' FINEGRID (LAY)        = ',FINEGRID(LAY)
      END DO

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFR)  'RFINDEX_PARAMETER = ',RFINDEX_PARAMETER

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        WRITE(OUTUNIT,DWFR1)  'LAY = ',LAY,' DELTAU_VERT_INPUT(LAY) = ',DELTAU_VERT_INPUT(LAY)
      END DO

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        DO MOM=0,NMOMENTS_INPUT
          WRITE(OUTUNIT,DWFR2)  'LAY = ',LAY,' MOM = ',MOM,&
            ' PHASMOMS_TOTAL_INPUT(MOM,LAY) = ',PHASMOMS_TOTAL_INPUT(MOM,LAY)
        END DO
      END DO

!  New section for 3.8 code. Phase Functions

      WRITE(OUTUNIT,*)
      DO GEO=1,NGEOMS
        DO LAY=1,NLAYERS
          WRITE(OUTUNIT,DWFR2)  'GEO = ',GEO,' LAY = ',LAY,&
                ' PHASFUNC_INPUT_UP(LAY,GEO) = ',PHASFUNC_INPUT_UP(LAY,GEO)
        END DO
      END DO

      WRITE(OUTUNIT,*)
      DO GEO=1,NGEOMS
        DO LAY=1,NLAYERS
          WRITE(OUTUNIT,DWFR2)  'GEO = ',GEO,' LAY = ',LAY,&
                ' PHASFUNC_INPUT_DN(LAY,GEO) = ',PHASFUNC_INPUT_DN(LAY,GEO)
        END DO
      END DO

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFR)   'LAMBERTIAN_ALBEDO = ',LAMBERTIAN_ALBEDO

      WRITE(OUTUNIT,*)
      DO LAY=0,NLAYERS
        WRITE(OUTUNIT,DWFR1)  'LAY = ',LAY,&
          ' THERMAL_BB_INPUT(LAY) = ',THERMAL_BB_INPUT(LAY)
      END DO

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFR)   'SURFACE_BB_INPUT  = ',SURFACE_BB_INPUT

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFR)   'ATMOS_WAVELENGTH  = ',ATMOS_WAVELENGTH

!      WRITE(OUTUNIT,*)
!      WRITE(OUTUNIT,DWFL)  'DO_DEBUG_WRITE          = ',DO_DEBUG_WRITE
!      WRITE(OUTUNIT,DWFL)  'DO_WRITE_INPUT          = ',DO_WRITE_INPUT
!      WRITE(OUTUNIT,DWFL)  'DO_WRITE_SCENARIO       = ',DO_WRITE_SCENARIO
!      WRITE(OUTUNIT,DWFL)  'DO_WRITE_FOURIER        = ',DO_WRITE_FOURIER
!      WRITE(OUTUNIT,DWFL)  'DO_WRITE_RESULTS        = ',DO_WRITE_RESULTS
!      WRITE(OUTUNIT,DWFC)  'INPUT_WRITE_FILENAME    = ',INPUT_WRITE_FILENAME
!      WRITE(OUTUNIT,DWFC)  'SCENARIO_WRITE_FILENAME = ',SCENARIO_WRITE_FILENAME
!      WRITE(OUTUNIT,DWFC)  'FOURIER_WRITE_FILENAME  = ',FOURIER_WRITE_FILENAME
!      WRITE(OUTUNIT,DWFC)  'RESULTS_WRITE_FILENAME  = ',RESULTS_WRITE_FILENAME

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,'(A)') '--------------------------'
      WRITE(OUTUNIT,'(A)') 'Standard Inputs - Modified'
      WRITE(OUTUNIT,'(A)') '--------------------------'

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFL)  'DO_FOCORR               = ',DO_FOCORR
      WRITE(OUTUNIT,DWFL)  'DO_FOCORR_EXTERNAL      = ',DO_FOCORR_EXTERNAL
      !WRITE(OUTUNIT,DWFL)  'DO_FOCORR_ALONE         = ',DO_FOCORR_ALONE
      WRITE(OUTUNIT,DWFL)  'DO_FOCORR_NADIR         = ',DO_FOCORR_NADIR
      WRITE(OUTUNIT,DWFL)  'DO_FOCORR_OUTGOING      = ',DO_FOCORR_OUTGOING
      WRITE(OUTUNIT,DWFL)  'DO_SSCORR_TRUNCATION    = ',DO_SSCORR_TRUNCATION
      WRITE(OUTUNIT,DWFL)  'DO_SSCORR_USEPHASFUNC   = ',DO_SSCORR_USEPHASFUNC
      WRITE(OUTUNIT,DWFL)  'DO_DOUBLE_CONVTEST      = ',DO_DOUBLE_CONVTEST
      WRITE(OUTUNIT,DWFL)  'DO_EXTERNAL_WLEAVE      = ',DO_EXTERNAL_WLEAVE   !    Introduced 4/22/19 for Version 3.8.1
      WRITE(OUTUNIT,DWFL)  'DO_SOLAR_SOURCES        = ',DO_SOLAR_SOURCES
      WRITE(OUTUNIT,DWFL)  'DO_REFRACTIVE_GEOMETRY  = ',DO_REFRAC_GEOMETRY
      WRITE(OUTUNIT,DWFL)  'DO_CHAPMAN_FUNCTION     = ',DO_CHAPMAN_FUNCTION
      WRITE(OUTUNIT,DWFL)  'DO_RAYLEIGH_ONLY        = ',DO_RAYLEIGH_ONLY
      WRITE(OUTUNIT,DWFL)  'DO_ISOTROPIC_ONLY       = ',DO_ISOTROPIC_ONLY
      WRITE(OUTUNIT,DWFL)  'DO_NO_AZIMUTH           = ',DO_NO_AZIMUTH
      WRITE(OUTUNIT,DWFL)  'DO_ALL_FOURIER          = ',DO_ALL_FOURIER
      WRITE(OUTUNIT,DWFL)  'DO_DELTAM_SCALING       = ',DO_DELTAM_SCALING
      WRITE(OUTUNIT,DWFL)  'DO_SOLUTION_SAVING      = ',DO_SOLUTION_SAVING
      WRITE(OUTUNIT,DWFL)  'DO_BVP_TELESCOPING      = ',DO_BVP_TELESCOPING
      WRITE(OUTUNIT,DWFL)  'DO_USER_STREAMS         = ',DO_USER_STREAMS
      WRITE(OUTUNIT,DWFL)  'DO_ADDITIONAL_MVOUT     = ',DO_ADDITIONAL_MVOUT
      WRITE(OUTUNIT,DWFL)  'DO_MVOUT_ONLY           = ',DO_MVOUT_ONLY
      WRITE(OUTUNIT,DWFL)  'DO_THERMAL_TRANSONLY    = ',DO_THERMAL_TRANSONLY
      WRITE(OUTUNIT,DWFL)  'DO_OBSERVATION_GEOMETRY = ',DO_OBSERV_GEOMETRY

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFI)  'NMOMENTS_INPUT = ',NMOMENTS_INPUT

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFI)  'NBEAMS         = ',NBEAMS
      DO SZA=1,NBEAMS
        WRITE(OUTUNIT,DWFR1)  'SZA = ',SZA,' BEAM_SZAS(SZA) = ',BEAM_SZAS(SZA)
      END DO

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFI)  'N_USER_RELAZMS = ',N_USER_RELAZMS
      DO URA=1,N_USER_RELAZMS
        WRITE(OUTUNIT,DWFR1)  'URA = ',URA,' USER_RELAZMS(URA) = ',USER_RELAZMS(URA)
      END DO

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFI)  'N_USER_STREAMS = ',N_USER_STREAMS
      DO UVA=1,N_USER_STREAMS
        WRITE(OUTUNIT,DWFR1)  'UVA = ',UVA,' USER_ANGLES_INPUT(UVA) = ',USER_ANGLES_INPUT(UVA)
      END DO

      WRITE(OUTUNIT,*)
      DO ULEV=1,N_USER_LEVELS
        WRITE(OUTUNIT,DWFR1)  'ULEV = ',ULEV,' USER_LEVELS(ULEV)  = ',USER_LEVELS(ULEV)
      END DO

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFR)  'GEOMETRY_SPECHEIGHT = ',GEOMETRY_SPECHEIGHT
      WRITE(OUTUNIT,DWFI)  'N_USER_OBSGEOMS     = ',N_USER_OBSGEOMS

      IF (N_USER_OBSGEOMS > 0) WRITE(OUTUNIT,*)
      DO UOG=1,N_USER_OBSGEOMS
        WRITE(OUTUNIT,DWFR1_3)  'UOG = ',UOG,&
          ' USER_OBSGEOMS(UOG,1:3) = ',USER_OBSGEOMS(UOG,1),&
                                   ',',USER_OBSGEOMS(UOG,2),&
                                   ',',USER_OBSGEOMS(UOG,3)
      END DO

! Note: Not used as external input yet
!      WRITE(OUTUNIT,*)
!      DO SZA=1,NBEAMS
!        DO LAYJ=1,NLAYERS
!          DO LAYI=1,NLAYERS
!            WRITE(OUTUNIT,DWFR3)  'SZA = ',SZA,' LAYJ = ',LAYJ,' LAYI = ',LAYI,&
!              ' CHAPMAN_FACTORS(LAYI,LAYJ,SZA) = ',&
!                CHAPMAN_FACTORS(LAYI,LAYJ,SZA)
!          END DO
!        END DO
!      END DO

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFR)  'EARTH_RADIUS        = ',EARTH_RADIUS

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        WRITE(OUTUNIT,DWFR1) 'LAY = ',LAY,' OMEGA_TOTAL_INPUT(LAY) = ',OMEGA_TOTAL_INPUT(LAY)
      END DO

!  Close output file

      CLOSE (OUTUNIT)

      END SUBROUTINE LIDORT_WRITE_STD_INPUT

!

      SUBROUTINE LIDORT_WRITE_SUP_BRDF_INPUT ( &
        DO_USER_STREAMS, DO_SURFACE_EMISSION, &
        NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, NSTREAMS, NMOMS, &
        EXACTDB_BRDFUNC, BRDF_F_0, BRDF_F, USER_BRDF_F_0, &
        USER_BRDF_F, EMISSIVITY, USER_EMISSIVITY )

      USE LIDORT_PARS_m, Only : fpk, MAXSTREAMS, MAX_USER_STREAMS, MAX_USER_RELAZMS, &
                                MAXBEAMS, MAXMOMENTS, DWFR1, DWFR2, DWFR3, DWFR4

      IMPLICIT NONE

      LOGICAL, INTENT(IN) :: DO_USER_STREAMS
      LOGICAL, INTENT(IN) :: DO_SURFACE_EMISSION

      INTEGER, INTENT(IN) :: NBEAMS
      INTEGER, INTENT(IN) :: N_USER_STREAMS
      INTEGER, INTENT(IN) :: N_USER_RELAZMS
      INTEGER, INTENT(IN) :: NSTREAMS
      INTEGER, INTENT(IN) :: NMOMS

      REAL(fpk), INTENT(IN) :: EXACTDB_BRDFUNC ( MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

      REAL(fpk), INTENT(IN) :: BRDF_F_0        ( 0:MAXMOMENTS, MAXSTREAMS, MAXBEAMS )
      REAL(fpk), INTENT(IN) :: BRDF_F          ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTREAMS )
      REAL(fpk), INTENT(IN) :: USER_BRDF_F_0   ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXBEAMS )
      REAL(fpk), INTENT(IN) :: USER_BRDF_F     ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTREAMS )

      REAL(fpk), INTENT(IN) :: EMISSIVITY      ( MAXSTREAMS )
      REAL(fpk), INTENT(IN) :: USER_EMISSIVITY ( MAX_USER_STREAMS )

!  Local variables

      INTEGER :: OUTUNIT
      INTEGER :: MOM,STRM,USTRM,IB,URA,STRMI,STRMJ

!  Open output file

      OUTUNIT = 102
      OPEN (OUTUNIT,file = 'LIDORT_WRITE_SUP_BRDF_INPUT.dbg',status = 'replace')

!  Write all BRDF input to file

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,'(A)') '----------------------'
      WRITE(OUTUNIT,'(A)') 'BRDF Supplement Inputs'
      WRITE(OUTUNIT,'(A)') '----------------------'

      WRITE(OUTUNIT,*)
      DO IB=1,NBEAMS
        DO URA=1,N_USER_RELAZMS
          DO USTRM=1,N_USER_STREAMS
              WRITE(OUTUNIT,DWFR3) &
                'IB = ',IB,' URA = ',URA,' USTRM = ',USTRM,&
                ' EXACTDB_BRDFUNC(USTRM,URA,IB) = ',&
                  EXACTDB_BRDFUNC(USTRM,URA,IB)
          END DO
        END DO
      END DO

      WRITE(OUTUNIT,*)
      DO IB=1,NBEAMS
        DO STRM=1,NSTREAMS
          DO MOM=0,NMOMS
            WRITE(OUTUNIT,DWFR3) &
              'IB = ',IB,' STRM = ',STRM,' MOM = ',MOM,&
              ' BRDF_F_0(MOM,STRM,IB) = ',BRDF_F_0(MOM,STRM,IB)
          END DO
        END DO
      END DO

      WRITE(OUTUNIT,*)
      DO STRMJ=1,NSTREAMS
        DO STRMI=1,NSTREAMS
          DO MOM=0,NMOMS
            WRITE(OUTUNIT,DWFR3) &
              'STRMJ = ',STRMJ,' STRMI = ',STRMI,' MOM = ',MOM,&
              ' BRDF_F(MOM,STRMI,STRMJ) = ',BRDF_F(MOM,STRMI,STRMJ)
          END DO
        END DO
      END DO

      IF ( DO_USER_STREAMS ) THEN
        WRITE(OUTUNIT,*)
        DO IB=1,NBEAMS
          DO USTRM=1,N_USER_STREAMS
            DO MOM=0,NMOMS
              WRITE(OUTUNIT,DWFR3) &
                'IB = ',IB,' USTRM = ',USTRM,' MOM = ',MOM,&
                ' USER_BRDF_F_0(MOM,USTRM,IB) = ',&
                  USER_BRDF_F_0(MOM,USTRM,IB)
            END DO
          END DO
        END DO

        WRITE(OUTUNIT,*)
        DO STRM=1,NSTREAMS
          DO USTRM=1,N_USER_STREAMS
            DO MOM=0,NMOMS
              WRITE(OUTUNIT,DWFR3) &
                'STRM = ',STRM,' USTRM = ',USTRM,' MOM = ',MOM,&
                ' USER_BRDF_F(MOM,USTRM,STRM) = ',&
                  USER_BRDF_F(MOM,USTRM,STRM)
            END DO
          END DO
        END DO
      END IF

      IF ( DO_SURFACE_EMISSION ) THEN
        WRITE(OUTUNIT,*)
        DO STRM=1,NSTREAMS
          WRITE(OUTUNIT,DWFR1) ' STRM = ',STRM,&
            ' EMISSIVITY (STRM) = ',EMISSIVITY (STRM)
        END DO

        IF ( DO_USER_STREAMS ) THEN
          WRITE(OUTUNIT,*)
          DO USTRM=1,N_USER_STREAMS
            WRITE(OUTUNIT,DWFR1)  ' USTRM = ',USTRM,&
              ' USER_EMISSIVITY (USTRM) = ',USER_EMISSIVITY(USTRM)
          END DO
        END IF
      END IF

!  Close output file

      CLOSE (OUTUNIT)

      END SUBROUTINE LIDORT_WRITE_SUP_BRDF_INPUT

!

      SUBROUTINE LIDORT_WRITE_SUP_SS_INPUT ( &
        N_USER_LEVELS, INTENSITY_SS, INTENSITY_DB )

      USE LIDORT_PARS_m, Only : fpk, MAX_USER_LEVELS, MAX_GEOMETRIES, MAX_DIRECTIONS, DWFR2, DWFR3, DWFR4

      IMPLICIT NONE

      INTEGER, INTENT(IN) ::  N_USER_LEVELS

      REAL(fpk), INTENT(IN) ::  INTENSITY_SS ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAX_DIRECTIONS )
      REAL(fpk), INTENT(IN) ::  INTENSITY_DB ( MAX_USER_LEVELS, MAX_GEOMETRIES )

!  Local variables

      INTEGER :: OUTUNIT
      INTEGER :: DIR,GEO,ULEV

!  Open output file

      OUTUNIT = 103
      OPEN (OUTUNIT,file = 'LIDORT_WRITE_SUP_SS_INPUT.dbg',&
            status = 'replace')

!  Write all single-scatter input to file

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,'(A)') '--------------------'
      WRITE(OUTUNIT,'(A)') 'SS Supplement Inputs'
      WRITE(OUTUNIT,'(A)') '--------------------'

      WRITE(OUTUNIT,*)
      DO DIR=1,MAX_DIRECTIONS
        DO GEO=1,MAX_GEOMETRIES
          DO ULEV=1,N_USER_LEVELS
          WRITE(OUTUNIT,DWFR3) &
            'DIR = ',DIR,' GEO = ',GEO,' ULEV = ',ULEV,&
            ' INTENSITY_SS(ULEV,GEO,DIR) = ',INTENSITY_SS(ULEV,GEO,DIR)
          END DO
        END DO
      END DO

      WRITE(OUTUNIT,*)
      DO GEO=1,MAX_GEOMETRIES
        DO ULEV=1,N_USER_LEVELS
        WRITE(OUTUNIT,DWFR2) &
          'GEO = ',GEO,' ULEV = ',ULEV,&
          ' INTENSITY_DB(ULEV,GEO) = ',INTENSITY_DB(ULEV,GEO)
        END DO
      END DO

!  Close output file

      CLOSE (OUTUNIT)

      END SUBROUTINE LIDORT_WRITE_SUP_SS_INPUT

!

      SUBROUTINE LIDORT_WRITE_SUP_SLEAVE_INPUT ( &
        DO_USER_STREAMS, &
        NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, NSTREAMS, NMOMS, &
        SLTERM_ISOTROPIC, SLTERM_USERANGLES, SLTERM_F_0, USER_SLTERM_F_0 )

      USE LIDORT_PARS_m, Only : fpk, MAXMOMENTS, MAXSTREAMS, MAX_USER_STREAMS, &
                                MAX_USER_RELAZMS, MAXBEAMS, DWFR2, DWFR4, DWFR1, DWFR3

      IMPLICIT NONE

      LOGICAL, INTENT(IN) :: DO_USER_STREAMS

      INTEGER, INTENT(IN) :: NBEAMS
      INTEGER, INTENT(IN) :: N_USER_STREAMS
      INTEGER, INTENT(IN) :: N_USER_RELAZMS
      INTEGER, INTENT(IN) :: NSTREAMS
      INTEGER, INTENT(IN) :: NMOMS

      REAL(fpk), INTENT(IN) :: SLTERM_ISOTROPIC  ( MAXBEAMS )
      REAL(fpk), INTENT(IN) :: SLTERM_USERANGLES ( MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS  )

      REAL(fpk), INTENT(IN) :: SLTERM_F_0      ( 0:MAXMOMENTS, MAXSTREAMS,       MAXBEAMS )
      REAL(fpk), INTENT(IN) :: USER_SLTERM_F_0 ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXBEAMS )

!  Local variables

      INTEGER :: OUTUNIT
      INTEGER :: MOM,STRM,USTRM,IB,URA

!  Open output file

      OUTUNIT = 104
      OPEN (OUTUNIT,file = 'LIDORT_WRITE_SUP_SLEAVE_INPUT.dbg',status = 'replace')

!  Write all surface-leaving input to file

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,'(A)') '------------------------'
      WRITE(OUTUNIT,'(A)') 'SLEAVE Supplement Inputs'
      WRITE(OUTUNIT,'(A)') '------------------------'

      WRITE(OUTUNIT,*)
      DO IB=1,NBEAMS
        WRITE(OUTUNIT,DWFR1) &
          'IB = ',IB,&
          ' SLTERM_ISOTROPIC(IB) = ',SLTERM_ISOTROPIC(IB)
      END DO

      IF ( DO_USER_STREAMS ) THEN
        WRITE(OUTUNIT,*)
        DO IB=1,NBEAMS
          DO URA=1,N_USER_RELAZMS
            DO USTRM=1,N_USER_STREAMS
              WRITE(OUTUNIT,DWFR3) &
                'IB = ',IB,' URA = ',URA,' USTRM = ',USTRM,&
                ' SLTERM_USERANGLES(USTRM,URA,IB) = ',&
                  SLTERM_USERANGLES(USTRM,URA,IB)
            END DO
          END DO
        END DO
      END IF

      WRITE(OUTUNIT,*)
      DO IB=1,NBEAMS
        DO STRM=1,NSTREAMS
          DO MOM=0,NMOMS
            WRITE(OUTUNIT,DWFR3) &
              'IB = ',IB,' STRM = ',STRM,' MOM = ',MOM,&
              ' SLTERM_F_0(MOM,STRM,IB) = ',SLTERM_F_0(MOM,STRM,IB)
          END DO
        END DO
      END DO

      IF ( DO_USER_STREAMS ) THEN
        WRITE(OUTUNIT,*)
        DO IB=1,NBEAMS
          DO USTRM=1,N_USER_STREAMS
            DO MOM=0,NMOMS
              WRITE(OUTUNIT,DWFR3) &
                'IB = ',IB,' USTRM = ',USTRM,' MOM = ',MOM,&
                ' USER_SLTERM_F_0(MOM,USTRM,IB) = ',&
                  USER_SLTERM_F_0(MOM,USTRM,IB)
            END DO
          END DO
        END DO
      END IF

!  Close output file

      CLOSE (OUTUNIT)

      END SUBROUTINE LIDORT_WRITE_SUP_SLEAVE_INPUT

      END MODULE lidort_writemodules_m

