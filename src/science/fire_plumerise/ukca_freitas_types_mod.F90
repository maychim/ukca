! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Dimensions and state type for the Freitas et al. (2007, 2010) 1D
!   plume-rise model.
!
!   All variables of the plume model for one plume calculation are held
!   in a variable of type freitas_column_type. A new local variable of this
!   type is created for every calculation (all components are default-
!   initialised to zero), so no state is kept between calculations and
!   the model has no module-level variables.
!
! Part of the UKCA model, a community model supported by
! The Met Office and NCAS, with components provided initially
! by The University of Cambridge, University of Leeds and
! The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code Description:
!   Language:  FORTRAN 2003
!   This code is written to UMDP3 programming standards.
!
! ---------------------------------------------------------------------------
MODULE ukca_freitas_types_mod

IMPLICIT NONE

PRIVATE

! Number of levels of the plume-model grid
INTEGER, PARAMETER, PUBLIC :: nkp_prm = 200

! Maximum burn duration (minutes); dimension of the heating time series
INTEGER, PARAMETER, PUBLIC :: ntime_prm = 200

! Vertical grid spacing of the plume-model grid (m)
REAL, PARAMETER, PUBLIC :: dz_prm = 100.0

! ---------------------------------------------------------------------------
! State of one plume calculation.
! Units follow the original model: mixing ratios in g/g, plume density in
! g m-3, pressure in kPa, heights in m.
! ---------------------------------------------------------------------------
TYPE, PUBLIC :: freitas_column_type

  ! -- Plume-model grid --
  REAL :: zt(nkp_prm)  = 0.0     ! height of thermodynamic levels (m)
  REAL :: zm(nkp_prm)  = 0.0     ! height of dynamical (w) levels (m)
  REAL :: dzt(nkp_prm) = 0.0     ! 1 / thickness around thermo. levels
  REAL :: dzm(nkp_prm) = 0.0     ! 1 / thickness around dyn. levels

  ! -- Plume (blob) properties --
  REAL :: w(nkp_prm)      = 0.0  ! vertical velocity (m s-1)
  REAL :: t(nkp_prm)      = 0.0  ! temperature (K)
  REAL :: qv(nkp_prm)     = 0.0  ! water vapour mixing ratio (g/g)
  REAL :: qc(nkp_prm)     = 0.0  ! cloud water (g/g)
  REAL :: qh(nkp_prm)     = 0.0  ! rain / hail (g/g)
  REAL :: qi(nkp_prm)     = 0.0  ! ice (g/g)
  REAL :: sc(nkp_prm)     = 0.0  ! passive tracer
  REAL :: vth(nkp_prm)    = 0.0  ! rain fall speed (m s-1)
  REAL :: vti(nkp_prm)    = 0.0  ! ice fall speed (m s-1)
  REAL :: rho(nkp_prm)    = 0.0  ! air density (g m-3)
  REAL :: txs(nkp_prm)    = 0.0  ! temperature excess over environment (K)
  REAL :: est(nkp_prm)    = 0.0  ! saturation vapour pressure (kPa)
  REAL :: qsat(nkp_prm)   = 0.0  ! saturation mixing ratio (g/g)
  REAL :: vel_p(nkp_prm)  = 0.0  ! horizontal plume velocity (m s-1)
  REAL :: rad_p(nkp_prm)  = 0.0  ! plume radius, prognostic (m)
  REAL :: radius(nkp_prm) = 0.0  ! plume (entrainment) radius (m)
  REAL :: visc(nkp_prm)   = 0.0  ! eddy viscosity (m2 s-1)
  REAL :: cvh(nkp_prm)    = 0.0  ! rain ventilation coefficient
  REAL :: cvi(nkp_prm)    = 0.0  ! ice ventilation coefficient
  REAL :: rbuoy(nkp_prm)  = 0.0  ! buoyancy acceleration (m s-2)

  ! -- Tendencies and work arrays --
  REAL :: wc(nkp_prm)    = 0.0   ! w at the current time level
  REAL :: wt(nkp_prm)    = 0.0   ! tendency of w
  REAL :: tt(nkp_prm)    = 0.0   ! tendency of t
  REAL :: qvt(nkp_prm)   = 0.0   ! tendency of qv
  REAL :: qct(nkp_prm)   = 0.0   ! tendency of qc
  REAL :: qht(nkp_prm)   = 0.0   ! tendency of qh
  REAL :: qit(nkp_prm)   = 0.0   ! tendency of qi
  REAL :: sct(nkp_prm)   = 0.0   ! tendency of sc
  REAL :: vel_t(nkp_prm) = 0.0   ! tendency of vel_p
  REAL :: rad_t(nkp_prm) = 0.0   ! tendency of rad_p

  ! -- Environment on the plume-model grid --
  REAL :: te(nkp_prm)    = 0.0   ! temperature (K)
  REAL :: pe(nkp_prm)    = 0.0   ! pressure (kPa)
  REAL :: qvenv(nkp_prm) = 0.0   ! water vapour mixing ratio (g/g)
  REAL :: the(nkp_prm)   = 0.0   ! potential temperature (K)
  REAL :: thve(nkp_prm)  = 0.0   ! virtual potential temperature (K)
  REAL :: dne(nkp_prm)   = 0.0   ! dry air density (kg m-3)
  REAL :: upe(nkp_prm)   = 0.0   ! zonal wind (m s-1)
  REAL :: vpe(nkp_prm)   = 0.0   ! meridional wind (m s-1)
  REAL :: vel_e(nkp_prm) = 0.0   ! wind speed (m s-1)

  ! -- Fire --
  REAL    :: heating(ntime_prm) = 0.0  ! convective heat flux (W m-2)
  REAL    :: area      = 0.0     ! burnt area (m2)
  REAL    :: rsurf     = 0.0     ! plume radius at the surface (m)
  REAL    :: alpha     = 0.0     ! entrainment constant
  REAL    :: fmoist    = 0.0     ! fuel moisture fraction
  REAL    :: tdur      = 0.0     ! burn duration (s)
  INTEGER :: mdur      = 0       ! burn duration (min)
  INTEGER :: maxtime   = 0       ! end of the model run (s)

  ! -- Time stepping and diagnostics --
  REAL    :: dt        = 0.0     ! time step (s)
  REAL    :: time      = 0.0     ! elapsed time (s)
  REAL    :: viscosity = 0.0     ! viscosity constant
  REAL    :: wbar      = 0.0     ! layer-mean w used by the microphysics
  REAL    :: dqsdz     = 0.0     ! d(qsat)/dz used by the microphysics
  REAL    :: ztop      = 0.0     ! current plume top (m)
  REAL    :: ztop_(ntime_prm) = 0.0  ! plume top for each minute (m)
  INTEGER :: mintime   = 0       ! elapsed time (min)
  INTEGER :: nm1       = 0       ! top level of the integration
  INTEGER :: n         = 0       ! top level of the environment profile

END TYPE freitas_column_type

END MODULE ukca_freitas_types_mod
