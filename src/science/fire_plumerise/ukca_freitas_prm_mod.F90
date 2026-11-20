! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Freitas et al. 1D plume-rise model for vegetation fires
!   (CPTEC/INPE 2005-2006, 2009), adapted for UKCA.
!
!   freitas_column_plumerise runs the model in one host-model column, once
!   for each of a set of fire heat fluxes, and returns the plume top for
!   each. Following Freitas et al. (2007), runs with the lower and upper
!   heat flux of a fire give the bottom and top of its injection layer.
!
!   Error codes returned in ierr (0 = OK):
!     1  bad input column (too few levels, or burnt area <= 0)
!     2  host column does not cover the plume-model grid (interpolation)
!     4  burn duration longer than ntime_prm minutes
!    10  numerically unstable (dt < 0.25 s); plume top returned as 0
!
! References:
!   Freitas, S. R., et al.: Including the sub-grid scale plume rise of
!     vegetation fires in low resolution atmospheric transport models,
!     Atmos. Chem. Phys., 7, 3385-3398, 2007.
!   Freitas, S. R., et al.: Technical Note: Sensitivity of 1-D smoke plume
!     rise models to the inclusion of environmental wind drag,
!     Atmos. Chem. Phys., 10, 585-594, 2010.
!   Latham, D.: PLUMP: a one-dimensional plume predictor and cloud model
!     for fire and smoke managers, USDA Forest Service, INT-GTR-314, 1994.
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
MODULE ukca_freitas_prm_mod

USE ukca_freitas_types_mod,    ONLY: freitas_column_type, nkp_prm,           &
                                     ntime_prm, dz_prm
USE ukca_freitas_micro_mod,    ONLY: freitas_fallpart, freitas_waterbal
USE ukca_freitas_dyn_mod,      ONLY: freitas_tend0, freitas_vel_advect,      &
                                     freitas_scl_advect, freitas_scl_misc,   &
                                     freitas_scl_dyn_entrain,                &
                                     freitas_damp_grav_wave,                 &
                                     freitas_buoyancy, freitas_entrainment,  &
                                     freitas_visc_w, freitas_update,         &
                                     freitas_hadvance
USE ukca_um_legacy_mod,        ONLY: g, rgas => r, cp
USE ukca_config_constants_mod, ONLY: repsilon, c_virtual
USE ukca_constants,            ONLY: pi, zerodegc

USE parkind1,                  ONLY: jpim, jprb      ! DrHook
USE yomhook,                   ONLY: lhook, dr_hook  ! DrHook

IMPLICIT NONE

PRIVATE

PUBLIC :: freitas_column_plumerise

CHARACTER(LEN=*), PARAMETER :: ModuleName = 'UKCA_FREITAS_PRM_MOD'

! Error codes (see header)
INTEGER, PARAMETER :: ierr_bad_input   = 1
INTEGER, PARAMETER :: ierr_interp      = 2
INTEGER, PARAMETER :: ierr_burn_length = 4
INTEGER, PARAMETER :: ierr_unstable    = 10

! Fire parameters (Freitas et al., 2007)
REAL,    PARAMETER :: conv_frac = 0.55     ! convective fraction of heat flux
INTEGER, PARAMETER :: mdur_fire = 53       ! burn duration (min)
REAL,    PARAMETER :: heat_comb = 19.3e6   ! heat of combustion (J kg-1),
                                           ! tropical forest (Alta Floresta)
INTEGER, PARAMETER :: moist_fire = 10      ! fuel moisture (%)
REAL,    PARAMETER :: alpha_fire = 0.05 * 0.2  ! entrainment constant,
                                               ! x 0.2 as in Ma et al.

! Model set-up
REAL,    PARAMETER :: viscosity_prm = 500.0  ! viscosity constant
REAL,    PARAMETER :: tstpf  = 2.0       ! timestep factor
INTEGER, PARAMETER :: kkmax0 = 10        ! initial top level of integration
INTEGER, PARAMETER :: deltak = 20        ! levels above the plume integrated
REAL,    PARAMETER :: wtop   = 1.0       ! w threshold defining the top (m/s)

CONTAINS

! ---------------------------------------------------------------------------
SUBROUTINE freitas_column_plumerise(nk, z_t, theta_in, qv_in, p_in,     &
                                    t_in, u_in, v_in, burnt_area,            &
                                    nflux, heat_flux, ztop, ierr,            &
                                    w_prof, txs_prof)
! ---------------------------------------------------------------------------
! Description:
!   Runs the plume-rise model in one column, once for each heat flux.
!
!   Input profiles are on host-model levels 1..nk (bottom to top):
!     z_t      height of the theta levels above the surface (m)
!     theta_in potential temperature (K)
!     qv_in    water vapour mixing ratio (kg/kg)
!     p_in     pressure (Pa)
!     t_in     temperature (K)
!     u_in     zonal wind (m s-1)
!     v_in     meridional wind (m s-1)
!   burnt_area   fire area (m2)
!   heat_flux    total fire heat flux (W m-2) for each run
!
!   Output:
!     ztop       plume top (m above surface) for each run; 0 if no plume
!     ierr       0 if OK, otherwise the largest error code (see header)
!     w_prof     (optional) final w on the plume-model grid (m s-1)
!     txs_prof   (optional) final plume temperature excess (K)
! ---------------------------------------------------------------------------
IMPLICIT NONE

INTEGER, INTENT(IN)  :: nk
REAL,    INTENT(IN)  :: z_t(nk)
REAL,    INTENT(IN)  :: theta_in(nk)
REAL,    INTENT(IN)  :: qv_in(nk)
REAL,    INTENT(IN)  :: p_in(nk)
REAL,    INTENT(IN)  :: t_in(nk)
REAL,    INTENT(IN)  :: u_in(nk)
REAL,    INTENT(IN)  :: v_in(nk)
REAL,    INTENT(IN)  :: burnt_area
INTEGER, INTENT(IN)  :: nflux
REAL,    INTENT(IN)  :: heat_flux(nflux)
REAL,    INTENT(OUT) :: ztop(nflux)
INTEGER, INTENT(OUT) :: ierr
REAL, OPTIONAL, INTENT(OUT) :: w_prof(nkp_prm, nflux)
REAL, OPTIONAL, INTENT(OUT) :: txs_prof(nkp_prm, nflux)

INTEGER :: iflux, ierr_run

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'FREITAS_COLUMN_PLUMERISE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

ztop(:) = 0.0
ierr    = 0
IF (PRESENT(w_prof))   w_prof(:,:)   = 0.0
IF (PRESENT(txs_prof)) txs_prof(:,:) = 0.0

IF (nk < 2 .OR. burnt_area <= 0.0) THEN
  ierr = ierr_bad_input
ELSE
  DO iflux = 1, nflux
    IF (heat_flux(iflux) <= 0.0) CYCLE
    IF (PRESENT(w_prof) .AND. PRESENT(txs_prof)) THEN
      CALL run_plume(nk, z_t, theta_in, qv_in, p_in, t_in, u_in, v_in, &
                     burnt_area, heat_flux(iflux), ztop(iflux), ierr_run,    &
                     w_prof(:,iflux), txs_prof(:,iflux))
    ELSE
      CALL run_plume(nk, z_t, theta_in, qv_in, p_in, t_in, u_in, v_in, &
                     burnt_area, heat_flux(iflux), ztop(iflux), ierr_run)
    END IF
    ierr = MAX(ierr, ierr_run)
    ! Problems with the environment affect every run: stop here
    IF (ierr_run /= 0 .AND. ierr_run /= ierr_unstable) EXIT
  END DO
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

END SUBROUTINE freitas_column_plumerise

! ---------------------------------------------------------------------------
SUBROUTINE run_plume(nk, z_t, theta_in, qv_in, p_in, t_in, u_in, v_in, &
                     burnt_area, heat_flux, ztop, ierr, w_prof, txs_prof)
! ---------------------------------------------------------------------------
! Description:
!   One plume calculation. All model state is held in the local variable
!   prm, so every calculation starts from a clean state.
! ---------------------------------------------------------------------------
IMPLICIT NONE

INTEGER, INTENT(IN)  :: nk
REAL,    INTENT(IN)  :: z_t(nk), theta_in(nk), qv_in(nk)
REAL,    INTENT(IN)  :: p_in(nk), t_in(nk), u_in(nk), v_in(nk)
REAL,    INTENT(IN)  :: burnt_area, heat_flux
REAL,    INTENT(OUT) :: ztop
INTEGER, INTENT(OUT) :: ierr
REAL, OPTIONAL, INTENT(OUT) :: w_prof(nkp_prm), txs_prof(nkp_prm)

! Model state. Allocated on the heap (it is too large for the stack) and
! default-initialised to zero on allocation, so every calculation starts
! from a clean state and calls from different threads are independent.
TYPE(freitas_column_type), ALLOCATABLE :: prm
INTEGER :: kmt

ALLOCATE(prm)
ztop = 0.0
ierr = 0

CALL set_grid(prm)

CALL get_env_condition(prm, nk, z_t, theta_in, qv_in, p_in, t_in,      &
                       u_in, v_in, kmt, ierr)

IF (ierr == 0) CALL fire_properties(prm, heat_flux, burnt_area, ierr)

IF (ierr == 0) CALL makeplume(prm, kmt, ztop, ierr)

IF (PRESENT(w_prof))   w_prof(:)   = prm%w(:)
IF (PRESENT(txs_prof)) txs_prof(:) = prm%t(:) - prm%te(:)

DEALLOCATE(prm)

END SUBROUTINE run_plume

! ---------------------------------------------------------------------------
SUBROUTINE set_grid(prm)
! ---------------------------------------------------------------------------
! Description:
!   Plume-model grid: constant spacing dz_prm from the surface.
!   zt = thermodynamic and water levels, zm = dynamical levels.
! ---------------------------------------------------------------------------
IMPLICIT NONE

TYPE(freitas_column_type), INTENT(IN OUT) :: prm

INTEGER :: k

prm%zt(1) = 0.0
prm%zm(1) = 0.0
prm%zt(2) = prm%zt(1) + 0.5 * dz_prm
prm%zm(2) = prm%zm(1) + dz_prm

DO k = 3, nkp_prm
  prm%zt(k) = prm%zt(k-1) + dz_prm
  prm%zm(k) = prm%zm(k-1) + dz_prm
END DO

DO k = 1, nkp_prm - 1
  prm%dzm(k) = 1.0 / (prm%zt(k+1) - prm%zt(k))
END DO
prm%dzm(nkp_prm) = prm%dzm(nkp_prm-1)

DO k = 2, nkp_prm
  prm%dzt(k) = 1.0 / (prm%zm(k) - prm%zm(k-1))
END DO
prm%dzt(1) = prm%dzt(2) * prm%dzt(2) / prm%dzt(3)

END SUBROUTINE set_grid

! ---------------------------------------------------------------------------
SUBROUTINE get_env_condition(prm, nk, z_t, theta_in, qv_in, p_in,       &
                             t_in, u_in, v_in, kmt, ierr)
! ---------------------------------------------------------------------------
! Description:
!   Interpolates the host-model column onto the plume-model grid and
!   derives the other environmental fields. ln(p) is interpolated linearly
!   in height. kmt is the top plume-model level inside the host column
!   (at most nkp_prm - 1, so that level kmt + 1 exists).
! ---------------------------------------------------------------------------
IMPLICIT NONE

TYPE(freitas_column_type), INTENT(IN OUT) :: prm
INTEGER, INTENT(IN)  :: nk
REAL,    INTENT(IN)  :: z_t(nk), theta_in(nk), qv_in(nk)
REAL,    INTENT(IN)  :: p_in(nk), t_in(nk), u_in(nk), v_in(nk)
INTEGER, INTENT(OUT) :: kmt
INTEGER, INTENT(OUT) :: ierr

INTEGER :: k
REAL    :: lnp(nk)

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'GET_ENV_CONDITION'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

ierr = 0

! Top plume-model level inside the host column
kmt = 0
DO k = nkp_prm - 1, 1, -1
  IF (prm%zt(k) < z_t(nk)) THEN
    kmt = k
    EXIT
  END IF
END DO

IF (kmt < 3) THEN
  ierr = ierr_bad_input                    ! host column too shallow
ELSE
  lnp(:) = LOG(p_in(:))

  CALL htint(nk, u_in,     z_t, kmt, prm%upe,   prm%zt, ierr)
  CALL htint(nk, v_in,     z_t, kmt, prm%vpe,   prm%zt, ierr)
  CALL htint(nk, theta_in, z_t, kmt, prm%the,   prm%zt, ierr)
  CALL htint(nk, qv_in,    z_t, kmt, prm%qvenv, prm%zt, ierr)
  CALL htint(nk, t_in,     z_t, kmt, prm%te,    prm%zt, ierr)
  CALL htint(nk, lnp,      z_t, kmt, prm%pe,    prm%zt, ierr)
END IF

IF (ierr == 0) THEN
  DO k = 1, kmt
    prm%pe(k)    = EXP(prm%pe(k))                          ! Pa
    prm%qvenv(k) = MAX(prm%qvenv(k), 1.0e-8)
    prm%thve(k)  = prm%the(k) * (1.0 + c_virtual * prm%qvenv(k))
    prm%dne(k)   = prm%pe(k) /                                               &
                   (rgas * prm%te(k) * (1.0 + c_virtual * prm%qvenv(k)))
    prm%vel_e(k) = SQRT(prm%upe(k)**2 + prm%vpe(k)**2)
    prm%pe(k)    = prm%pe(k) * 1.0e-3                      ! Pa to kPa
  END DO
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

END SUBROUTINE get_env_condition

! ---------------------------------------------------------------------------
SUBROUTINE htint(nzz1, vctra, eleva, nzz2, vctrb, elevb, ierr)
! ---------------------------------------------------------------------------
! Description:
!   Linear interpolation of vctra (at heights eleva(1:nzz1)) to the heights
!   elevb(1:nzz2), with linear extrapolation outside eleva.
!   Sets ierr = ierr_interp if a target height cannot be bracketed.
! ---------------------------------------------------------------------------
IMPLICIT NONE

INTEGER, INTENT(IN)     :: nzz1, nzz2
REAL,    INTENT(IN)     :: vctra(nzz1), eleva(nzz1)
REAL,    INTENT(IN)     :: elevb(nzz2)
REAL,    INTENT(OUT)    :: vctrb(nzz2)
INTEGER, INTENT(IN OUT) :: ierr

INTEGER :: l, k
REAL    :: wt

l = 1

DO k = 1, nzz2
  DO
    IF ((elevb(k) < eleva(1)) .OR.                                           &
        ((elevb(k) >= eleva(l)) .AND. (elevb(k) <= eleva(l+1)))) THEN
      wt       = (elevb(k) - eleva(l)) / (eleva(l+1) - eleva(l))
      vctrb(k) = vctra(l) + (vctra(l+1) - vctra(l)) * wt
      EXIT
    ELSE IF (elevb(k) > eleva(nzz1)) THEN
      wt       = (elevb(k) - eleva(nzz1)) / (eleva(nzz1-1) - eleva(nzz1))
      vctrb(k) = vctra(nzz1) + (vctra(nzz1-1) - vctra(nzz1)) * wt
      EXIT
    END IF

    l = l + 1
    IF (l == nzz1) THEN
      ierr = ierr_interp
      vctrb(k:nzz2) = vctra(nzz1)
      RETURN
    END IF
  END DO
END DO

END SUBROUTINE htint

! ---------------------------------------------------------------------------
SUBROUTINE fire_properties(prm, heat_flux, burnt_area, ierr)
! ---------------------------------------------------------------------------
! Description:
!   Fire size and surface heating time series: constant convective heat
!   flux (conv_frac * heat_flux) for mdur_fire minutes, ramped up over the
!   first minutes for stability.
! ---------------------------------------------------------------------------
IMPLICIT NONE

TYPE(freitas_column_type), INTENT(IN OUT) :: prm
REAL,    INTENT(IN)  :: heat_flux      ! total fire heat flux (W m-2)
REAL,    INTENT(IN)  :: burnt_area     ! m2
INTEGER, INTENT(OUT) :: ierr

INTEGER :: i
REAL    :: hinc

ierr = 0

prm%area    = burnt_area
prm%mdur    = mdur_fire
prm%maxtime = prm%mdur - 1
prm%alpha   = alpha_fire

IF (prm%mdur > ntime_prm) THEN
  ierr = ierr_burn_length
  RETURN
END IF

IF (MOD(prm%maxtime, 2) /= 0) prm%maxtime = prm%maxtime + 1  ! make it even
prm%maxtime = prm%maxtime * 60                             ! in seconds
prm%rsurf   = SQRT(prm%area / pi)     ! entrainment radius at surface (m)
prm%fmoist  = moist_fire / 100.0      ! fuel moisture fraction
prm%tdur    = prm%mdur * 60.0         ! duration of the burn (s)

prm%heating(:) = 0.0001               ! avoid possible divide by 0
DO i = 1, prm%mdur
  prm%heating(i) = heat_flux * conv_frac     ! convective part (W m-2)
END DO

hinc = prm%heating(1) / 4.0           ! ramp up over the first minutes
prm%heating(1) = 0.1
prm%heating(2) = hinc
prm%heating(3) = 2.0 * hinc
prm%heating(4) = 3.0 * hinc

END SUBROUTINE fire_properties

! ---------------------------------------------------------------------------
SUBROUTINE burn(prm, eflux, water)
! ---------------------------------------------------------------------------
! Description:
!   Energy flux (W m-2) and fire water (g m-2) at the lower boundary for
!   the current time step. The emission factor for water is 0.5, so the
!   water produced is fuel mass * (0.5 + fuel moisture fraction).
! ---------------------------------------------------------------------------
IMPLICIT NONE

TYPE(freitas_column_type), INTENT(IN) :: prm
REAL, INTENT(OUT) :: eflux, water

IF (prm%time > prm%tdur) THEN          ! burn is over
  eflux = 0.000001                     ! prevent a potential divide by zero
  water = 0.0
ELSE
  eflux = prm%heating(prm%mintime)                                ! W m-2
  water = eflux * (prm%dt / heat_comb) * (0.5 + prm%fmoist)                  &
          / conv_frac                                             ! kg m-2
  water = water * 1000.0                                          ! g m-2
END IF

END SUBROUTINE burn

! ---------------------------------------------------------------------------
SUBROUTINE initial(prm, kmt)
! ---------------------------------------------------------------------------
! Description:
!   Initial conditions: plume equal to the environment, Turner-style
!   entrainment radius, viscosity profile, and the lower boundary.
! ---------------------------------------------------------------------------
IMPLICIT NONE

TYPE(freitas_column_type), INTENT(IN OUT) :: prm
INTEGER,                   INTENT(IN)     :: kmt

INTEGER :: k
REAL    :: es

prm%n = kmt

DO k = 1, prm%n
  prm%txs(k)   = 0.0
  prm%w(k)     = 0.0
  prm%t(k)     = prm%te(k)              ! plume set to environment
  prm%wc(k)    = 0.0
  prm%wt(k)    = 0.0
  prm%qv(k)    = prm%qvenv(k)
  prm%vth(k)   = 0.0                    ! initial rain velocity
  prm%vti(k)   = 0.0                    ! initial ice velocity
  prm%qh(k)    = 0.0                    ! no rain
  prm%qi(k)    = 0.0                    ! no ice
  prm%qc(k)    = 0.0                    ! no cloud drops
  es           = 0.1 * esat(prm%t(k))   ! saturation vapour pressure (kPa)
  prm%est(k)   = es
  prm%qsat(k)  = (repsilon * es) / (prm%pe(k) - es)
  prm%rho(k)   = 1.0e6 / rgas * prm%pe(k) / prm%t(k)   ! g m-3
  prm%vel_p(k) = 0.0
  prm%rad_p(k) = 0.0
END DO

! Entrainment radius, Turner-style plume
prm%radius(1) = prm%rsurf
DO k = 2, prm%n
  prm%radius(k) = prm%radius(k-1)                                            &
                + (6.0 / 5.0) * prm%alpha * (prm%zt(k) - prm%zt(k-1))
  prm%rad_p(k)  = prm%radius(k)
END DO
prm%rad_p(1) = prm%rsurf

! Viscosity
prm%visc(1) = prm%viscosity
DO k = 2, prm%n
  prm%visc(k) = MAX(1.0e-3,                                                  &
                    prm%visc(k-1) - 1.0 * prm%viscosity / REAL(nkp_prm))
END DO

CALL lbound_mtt(prm)

END SUBROUTINE initial

! ---------------------------------------------------------------------------
SUBROUTINE lbound_mtt(prm)
! ---------------------------------------------------------------------------
! Description:
!   Boundary conditions at the surface for plume and cloud, from a virtual
!   buoyant point source below ground (Turner, 1973, p.172; Briggs, 1969,
!   p.28). Gives w(1) and t(1), which vary with time, and the fire water,
!   which first saturates level 1 and then goes into qc(1).
! ---------------------------------------------------------------------------
IMPLICIT NONE

TYPE(freitas_column_type), INTENT(IN OUT) :: prm

REAL, PARAMETER :: e1 = 1.0 / 3.0
REAL, PARAMETER :: e2 = 5.0 / 3.0
REAL :: es, eflux, water, pres, c1, c2, f, zv, denscor, xwater

prm%qh(1) = prm%qh(2)                  ! soak up hydrometeors
prm%qi(1) = prm%qi(2)
prm%qc(1) = 0.0                        ! no cloud here

CALL burn(prm, eflux, water)

pres = prm%pe(1) * 1000.0              ! Pa
c1   = 5.0 / (6.0 * prm%alpha)
c2   = 0.9 * prm%alpha

! Buoyancy flux
f = eflux / (pres * cp * pi)
f = g * rgas * f * prm%area

zv = c1 * prm%rsurf                    ! virtual source height

prm%w(1) = c1 * ((c2 * f)**e1) / zv**e1
denscor  = c1 * f / g / (c2 * f)**e1 / zv**e2
prm%t(1) = prm%te(1) / (1.0 - denscor)

prm%wc(1)    = prm%w(1)
prm%vel_p(1) = 0.0
prm%rad_p(1) = prm%rsurf
prm%sc(1)    = 1000.0
prm%vth(1)   = -4.0
prm%vti(1)   = -3.0
prm%txs(1)   = prm%t(1) - prm%te(1)
prm%visc(1)  = prm%viscosity
prm%rho(1)   = 1.0e6 / rgas * prm%pe(1) / prm%t(1)     ! g m-3

xwater    = water / (prm%w(1) * prm%dt * prm%rho(1))   ! fire water (g/g)
prm%qv(1) = xwater + prm%qvenv(1)

es          = 0.1 * esat(prm%t(1))     ! kPa
prm%est(1)  = es
prm%qsat(1) = (repsilon * es) / (prm%pe(1) - es)

IF (prm%qv(1) > prm%qsat(1)) THEN      ! excess goes into cloud drops
  prm%qc(1) = prm%qv(1) - prm%qsat(1) + prm%qc(1)
  prm%qv(1) = prm%qsat(1)
END IF

CALL freitas_waterbal(prm, 1)

END SUBROUTINE lbound_mtt

! ---------------------------------------------------------------------------
SUBROUTINE makeplume(prm, kmt, ztopmax, ierr)
! ---------------------------------------------------------------------------
! Description:
!   Time integration of the plume model. The plume top is the highest level
!   where w > wtop; ztopmax is its maximum over the run (m above surface).
!
!   Equation sources: Kessler, Met. Monograph 32 (K); Weinstein, JAS 27,
!   246-255 (W); Ogura and Takahashi, MWR 99, 895-911 (OT); Pielke,
!   Mesoscale Meteorological Modeling, 1984. Originally developed by
!   D. Latham (USFS).
! ---------------------------------------------------------------------------
IMPLICIT NONE

TYPE(freitas_column_type), INTENT(IN OUT) :: prm
INTEGER, INTENT(IN)  :: kmt
REAL,    INTENT(OUT) :: ztopmax
INTEGER, INTENT(OUT) :: ierr

INTEGER, PARAMETER :: n_sub_step = 3   ! microphysics sub-steps
INTEGER :: k, kk, kkmax, i_micro, l
REAL    :: wmax, rmaxtime, es, dt_save, dvel

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'MAKEPLUME'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

ierr = 0

! Problem set-up and initial conditions
prm%mintime   = 1
ztopmax       = 0.0
prm%ztop      = 0.0
prm%time      = 0.0
prm%dt        = 1.0
wmax          = 1.0
kkmax         = kkmax0
prm%viscosity = viscosity_prm

CALL initial(prm, kmt)

rmaxtime = REAL(prm%maxtime)

DO WHILE (prm%time <= rmaxtime)

  ! Top of the integration
  prm%nm1 = MIN(kmt, kkmax + deltak)

  ! Time step (CFL, and PNNL recommendation for the horizontal motion)
  prm%dt = MIN(5.0, (prm%zm(2) - prm%zm(1)) / (tstpf * wmax))
  dvel   = MAX(ABS(prm%vel_p(2) - prm%vel_e(2)), 1.0e-6)
  prm%dt = MIN(prm%dt, 0.8 * pi * prm%radius(2) / (2.0 * dvel))

  IF (prm%dt < 0.25 .OR. MAXVAL(prm%t(:)) < 0.0) THEN
    ! Numerically unstable: no plume top
    prm%ztop = 0.0
    prm%ztop_(prm%mintime) = prm%ztop
    ztopmax = 0.0
    ierr    = ierr_unstable
    EXIT
  END IF

  prm%time    = prm%time + prm%dt                 ! elapsed time (s)
  prm%mintime = 1 + INT(prm%time) / 60            ! elapsed time (min)
  wmax = 1.0                                      ! no zeroes allowed

  CALL freitas_tend0(prm)

  ! Surface boundary conditions (k = 1)
  CALL lbound_mtt(prm)

  ! Advection of w and scalars
  CALL freitas_vel_advect(prm)
  CALL freitas_scl_advect(prm)

  ! Scalars: dry adiabat and lateral entrainment, dynamic entrainment
  CALL freitas_scl_misc(prm)
  CALL freitas_scl_dyn_entrain(prm)

  ! Sponge layer for t
  CALL freitas_damp_grav_wave(prm, 1, deltak)

  ! Microphysics, with sub-steps
  dt_save = prm%dt
  prm%dt  = prm%dt / REAL(n_sub_step)

  DO i_micro = 1, n_sub_step
    CALL freitas_fallpart(prm)
    DO l = 2, prm%nm1 - 1
      prm%wbar    = 0.5 * (prm%w(l) + prm%w(l-1))
      es          = 0.1 * esat(prm%t(l))          ! kPa
      prm%qsat(l) = (repsilon * es) / (prm%pe(l) - es)
      prm%est(l)  = es
      prm%rho(l)  = 1.0e6 / rgas * prm%pe(l) / prm%t(l)   ! g m-3
      prm%dqsdz   = (prm%qsat(l+1) - prm%qsat(l-1))                          &
                  / (prm%zt(l+1) - prm%zt(l-1))
      CALL freitas_waterbal(prm, l)
    END DO
  END DO
  prm%dt = dt_save

  ! Viscosity for stability
  CALL freitas_visc_w(prm, kmt)

  ! Update scalars
  CALL freitas_update(prm, 'S')
  CALL freitas_hadvance(prm, 1)

  ! Buoyancy and entrainment
  CALL freitas_buoyancy(prm)
  CALL freitas_entrainment(prm)

  ! Update w
  CALL freitas_update(prm, 'W')
  CALL freitas_hadvance(prm, 2)

  ! Diagnostics
  DO k = 2, prm%nm1
    es          = 0.1 * esat(prm%t(k))            ! kPa
    prm%qsat(k) = (repsilon * es) / (prm%pe(k) - es)
    prm%est(k)  = es
    prm%txs(k)  = prm%t(k) - prm%te(k)
    prm%rho(k)  = 1.0e6 / rgas * prm%pe(k) / prm%t(k)     ! g m-3
    IF (ABS(prm%wc(k)) > wmax) wmax = ABS(prm%wc(k))      ! largest w
  END DO

  ! Sponge layer for w
  CALL freitas_damp_grav_wave(prm, 2, deltak)

  ! Update radius
  DO k = 2, prm%nm1
    prm%radius(k) = prm%rad_p(k)
  END DO

  ! Plume top (above the surface)
  kk = 1
  DO WHILE (prm%w(kk) > wtop .AND. kk < nkp_prm)
    kk = kk + 1
    prm%ztop = prm%zm(kk)
  END DO

  prm%ztop_(prm%mintime) = prm%ztop
  ztopmax = MAX(prm%ztop, ztopmax)
  kkmax   = MAX(kk, kkmax)

  ! Exit if the solution has become stationary
  IF (prm%mintime > 40) THEN
    IF (ABS(prm%ztop_(prm%mintime) - prm%ztop_(prm%mintime-10)) < dz_prm)    &
      EXIT
  END IF

END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

END SUBROUTINE makeplume

! ---------------------------------------------------------------------------
REAL FUNCTION esat(t)
! ---------------------------------------------------------------------------
! Description:
!   Saturation vapour pressure over water (hPa) at temperature t (K).
! ---------------------------------------------------------------------------
IMPLICIT NONE

REAL, INTENT(IN) :: t

REAL :: tc

tc   = t - zerodegc
esat = 6.1078 * EXP((17.2693882 * tc) / (tc + 237.3))

END FUNCTION esat

END MODULE ukca_freitas_prm_mod
