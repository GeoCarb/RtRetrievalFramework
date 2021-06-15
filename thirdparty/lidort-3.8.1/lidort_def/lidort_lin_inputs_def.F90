
      module LIDORT_Lin_Inputs_def_m

!  Version 3.7, Internal threading removed  , 02 May   2014
!  Version 3.8, Inclusion of phase functions, 03 March 2017
!  Version 3.8.1, Upgrade, June 2019

!  This Module contains the following LIDORT Input Structures, with Intents :

!          LIDORT_Fixed_LinControl    nested in LIDORT_Fixed_LinInputs
!          LIDORT_Fixed_LinOptical    nested in LIDORT_Fixed_LinInputs
!           LIDORT_Fixed_LinInputs    Intent(In)

      use LIDORT_PARS_m, only : fpk, MAXMOMENTS_INPUT, MAXLAYERS, MAX_GEOMETRIES, MAX_ATMOSWFS

      implicit none

! #####################################################################
! #####################################################################

      TYPE LIDORT_Fixed_LinControl

!  Control for atmospheric linearizations, layer by layer

      LOGICAL, dimension(MAXLAYERS)  :: TS_LAYER_VARY_FLAG
      INTEGER, dimension(MAXLAYERS)  :: TS_LAYER_VARY_NUMBER

!  Total number of column Jacobians

      INTEGER  :: TS_N_TOTALCOLUMN_WFS

!  Total number of Surface Jacobians

      INTEGER  :: TS_N_SURFACE_WFS

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@ Addition of SLEAVE WF stuff, R. Spurr, 22 August 2012 @@@@@@@@@
!  Total number of Sleave Jacobians
      INTEGER  :: TS_N_SLEAVE_WFS
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Jacobian names

      CHARACTER (LEN=31), dimension(MAX_ATMOSWFS) :: TS_COLUMNWF_NAMES
      CHARACTER (LEN=31), dimension(MAX_ATMOSWFS) :: TS_PROFILEWF_NAMES

      END TYPE LIDORT_Fixed_LinControl

! #####################################################################
! #####################################################################

      TYPE LIDORT_Fixed_LinOptical

!  Optical property linearizations
!  Layer linearization (bulk property variation) input
!  Layer linearization (phase function variation) input

      REAL(fpk), dimension(MAX_ATMOSWFS,MAXLAYERS) :: TS_L_DELTAU_VERT_INPUT
      REAL(fpk), dimension(MAX_ATMOSWFS,MAXLAYERS) :: TS_L_OMEGA_TOTAL_INPUT
      REAL(fpk), dimension(MAX_ATMOSWFS,0:MAXMOMENTS_INPUT,MAXLAYERS) :: TS_L_PHASMOMS_TOTAL_INPUT

!  Linearized Phase function inputs for FO calculations. Version 3.8. 3/3/17

      REAL(fpk) :: TS_L_PHASFUNC_INPUT_UP ( MAX_ATMOSWFS, MAXLAYERS, MAX_GEOMETRIES )
      REAL(fpk) :: TS_L_PHASFUNC_INPUT_DN ( MAX_ATMOSWFS, MAXLAYERS, MAX_GEOMETRIES )

      END TYPE LIDORT_Fixed_LinOptical

! #####################################################################
! #####################################################################

      TYPE LIDORT_Fixed_LinInputs

      TYPE(LIDORT_Fixed_LinControl)    :: Cont
      TYPE(LIDORT_Fixed_LinOptical)    :: Optical

      END TYPE LIDORT_Fixed_LinInputs

! #####################################################################
! #####################################################################

      TYPE LIDORT_Modified_LinControl

!  Structure renamed in line with VLIDORT 2.8.
!  Control linearization. (3/3/17) Revised Version 3.8, in line with VLIDORT 2.8
!    Some of thiese were moved from the Fixed Type

      LOGICAL  :: TS_DO_COLUMN_LINEARIZATION
      LOGICAL  :: TS_DO_PROFILE_LINEARIZATION
      LOGICAL  :: TS_DO_ATMOS_LINEARIZATION

      LOGICAL  :: TS_DO_SURFACE_LINEARIZATION
      LOGICAL  :: TS_DO_LINEARIZATION

      LOGICAL  :: TS_DO_SIMULATION_ONLY

!  BlackBody Jacobian Flags, Introduced March 18th 2014

      LOGICAL  :: TS_DO_ATMOS_LBBF
      LOGICAL  :: TS_DO_SURFACE_LBBF

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@ Addition of SLEAVE WF stuff, R. Spurr, 22 August 2012 @@@@@@@@@
!  Total number of Sleave Jacobians
      LOGICAL  :: TS_DO_SLEAVE_WFS
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Old 3.7 code was just nothing
!      INTEGER :: Dummy

      END TYPE LIDORT_Modified_LinControl

! #####################################################################
! #####################################################################

      TYPE LIDORT_Modified_LinInputs


      TYPE(LIDORT_Modified_LinControl)    :: MCont


      END TYPE LIDORT_Modified_LinInputs

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

      PRIVATE
      PUBLIC :: LIDORT_Fixed_LinControl,    &
                LIDORT_Fixed_LinOptical,    &
                LIDORT_Fixed_LinInputs,     &
                LIDORT_Modified_LinControl, &
                LIDORT_Modified_LinInputs

      end module LIDORT_Lin_Inputs_def_m
