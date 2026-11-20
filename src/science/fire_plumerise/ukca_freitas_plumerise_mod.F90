! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Couples the Freitas et al. (2007, 2010) 1D plume-rise model to the
!   'stratos_inj' biomass-burning emissions in UKCA.
!
!   ukca_freitas_plumerise : runs the plume model (ukca_freitas_prm_mod)
!                            in every column with
!                            stratos_inj emissions and returns the bottom
!                            and top of the injection layer
!   freitas_vert_profile   : converts these heights into a vertical
!                            emission profile (fraction per model level)
!
!   Intended for case studies of a single large fire: the fire size and
!   heat-flux bounds are module parameters, applied to every fire column.
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
MODULE ukca_freitas_plumerise_mod

USE ukca_freitas_prm_mod,    ONLY: freitas_column_plumerise
USE ereport_mod,             ONLY: ereport
USE errormessagelength_mod,  ONLY: errormessagelength
USE parkind1,                ONLY: jpim, jprb      ! DrHook
USE yomhook,                 ONLY: lhook, dr_hook  ! DrHook

IMPLICIT NONE

PRIVATE
PUBLIC :: ukca_freitas_plumerise, freitas_vert_profile

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName =                           &
                                        'UKCA_FREITAS_PLUMERISE_MOD'

! ---------------------------------------------------------------------------
! Case-study fire settings. Edit these (or move them to ukca_config) for
! your event. The plume top is very sensitive to all three.
! ---------------------------------------------------------------------------
! Fire area per fire column (m2)
REAL, PARAMETER :: freitas_burnt_area   = 600.0e4
! Total fire heat flux (W m-2) for the bottom and the top of the injection
! layer (Freitas et al., 2007). Forest fires: 30e3 - 80e3.
! From fire radiative power (W): heat_flux = beta * FRP / area / 0.55,
! beta = 0.88 (Paugam et al., 2015) or 5.0 (Wooster et al., 2005).
REAL, PARAMETER :: freitas_heat_flux_lo = 30.0e3
REAL, PARAMETER :: freitas_heat_flux_hi = 80.0e3
! Columns with a stratos_inj flux above this (kg m-2 s-1) are fire columns
REAL, PARAMETER :: freitas_bb_threshold = 0.0
! .TRUE.: stop the run if the plume model fails in a fire column.
! .FALSE.: print a warning and return ztop = 0 for that column, so that
!          the caller can fall back to the static stratos_inj profile.
LOGICAL, PARAMETER :: l_freitas_stop_on_failure = .TRUE.

CONTAINS

! ---------------------------------------------------------------------------
SUBROUTINE ukca_freitas_plumerise(row_length, rows, model_levels, bb_sfc,      &
                                  r_theta_levels, r_rho_levels, theta, q,      &
                                  p_theta_levels, t_theta_levels,              &
                                  u_rho_levels, v_rho_levels,                  &
                                  ztop_lo, ztop_hi)
! ---------------------------------------------------------------------------
! Description:
!   Runs the plume model in each column where bb_sfc > freitas_bb_threshold,
!   once with freitas_heat_flux_lo and once with freitas_heat_flux_hi.
!   Returns the two plume tops (m above the surface); 0 in other columns.
! ---------------------------------------------------------------------------
IMPLICIT NONE

INTEGER, INTENT(IN) :: row_length, rows, model_levels

! Total stratos_inj biomass-burning flux in each column (kg m-2 s-1)
REAL, INTENT(IN) :: bb_sfc(row_length, rows)
! Distance from the centre of the Earth (m); level 0 is the surface
REAL, INTENT(IN) :: r_theta_levels(row_length, rows, 0:model_levels)
REAL, INTENT(IN) :: r_rho_levels  (row_length, rows, model_levels)
REAL, INTENT(IN) :: theta         (row_length, rows, model_levels)  ! K
REAL, INTENT(IN) :: q             (row_length, rows, model_levels)  ! kg/kg
REAL, INTENT(IN) :: p_theta_levels(row_length, rows, model_levels)  ! Pa
REAL, INTENT(IN) :: t_theta_levels(row_length, rows, model_levels)  ! K
! Winds on rho levels (environment fields; assumed-shape, as allocated there)
REAL, INTENT(IN) :: u_rho_levels(:,:,:)                             ! m/s
REAL, INTENT(IN) :: v_rho_levels(:,:,:)                             ! m/s

REAL, INTENT(OUT) :: ztop_lo(row_length, rows)  ! plume top, low flux (m)
REAL, INTENT(OUT) :: ztop_hi(row_length, rows)  ! plume top, high flux (m)

! Local variables
INTEGER :: i, j, k, ierr
! Column profiles passed to the plume model
REAL :: z_t(model_levels)      ! height of theta levels above surface (m)
REAL :: z_r(model_levels)      ! height of rho levels above surface (m)
REAL :: qv(model_levels)       ! water vapour mixing ratio (kg/kg)
REAL :: th_col(model_levels)   ! potential temperature (K)
REAL :: p_col(model_levels)    ! pressure (Pa)
REAL :: t_col(model_levels)    ! temperature (K)
REAL :: u_col(model_levels)    ! zonal wind (m s-1)
REAL :: v_col(model_levels)    ! meridional wind (m s-1)
REAL :: ztop(2)                ! plume tops for the two heat fluxes (m)

INTEGER                            :: errcode
CHARACTER (LEN=errormessagelength) :: cmessage

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'UKCA_FREITAS_PLUMERISE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

ztop_lo(:,:) = 0.0
ztop_hi(:,:) = 0.0

IF (SIZE(u_rho_levels, 3) < model_levels .OR.                                 &
    SIZE(v_rho_levels, 3) < model_levels) THEN
  errcode  = 1
  cmessage = 'u_rho_levels/v_rho_levels not available on all model levels'
  CALL ereport(ModuleName//':'//RoutineName, errcode, cmessage)
END IF

DO j = 1, rows
  DO i = 1, row_length

    IF (bb_sfc(i,j) <= freitas_bb_threshold) CYCLE

    ! Column profiles for the plume model
    DO k = 1, model_levels
      z_t(k)    = r_theta_levels(i,j,k) - r_theta_levels(i,j,0)
      z_r(k)    = r_rho_levels(i,j,k)   - r_theta_levels(i,j,0)
      qv(k)     = q(i,j,k) / (1.0 - q(i,j,k))   ! specific hum. to mix. ratio
      th_col(k) = theta(i,j,k)
      p_col(k)  = p_theta_levels(i,j,k)
      t_col(k)  = t_theta_levels(i,j,k)
      u_col(k)  = u_rho_levels(i,j,k)
      v_col(k)  = v_rho_levels(i,j,k)
    END DO

    CALL freitas_column_plumerise(model_levels, z_t, z_r, th_col, qv,        &
           p_col, t_col, u_col, v_col, freitas_burnt_area, 2,                &
           [freitas_heat_flux_lo, freitas_heat_flux_hi], ztop, ierr)

    ! A failed run returns ztop = 0 for that heat flux
    IF (ierr /= 0 .OR. ztop(2) <= 0.0) THEN
      WRITE(cmessage, '(A,I0,A,2I6)') 'Freitas plume model failed (code ',    &
            ierr, ') in fire column i,j =', i, j
      IF (l_freitas_stop_on_failure) THEN
        errcode = MAX(ierr, 1)          ! positive: fatal
      ELSE
        errcode = -MAX(ierr, 1)         ! negative: warning only
        ztop(:) = 0.0                   ! caller falls back to static profile
      END IF
      CALL ereport(ModuleName//':'//RoutineName, errcode, cmessage)
    END IF

    ztop_lo(i,j) = ztop(1)
    ztop_hi(i,j) = ztop(2)

  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

END SUBROUTINE ukca_freitas_plumerise

! ---------------------------------------------------------------------------
SUBROUTINE freitas_vert_profile(row_length, rows, model_levels,                &
                                r_theta_levels, ztop_lo, ztop_hi, vf_prm)
! ---------------------------------------------------------------------------
! Description:
!   Fraction of the column emission to put in each model level: uniform in
!   height between ztop_lo and ztop_hi (weighted by the overlap of each
!   level with that layer, as for 'stratos_inj'), summing to 1. If the
!   layer is thinner than a model level, everything goes in the level that
!   contains ztop_hi. Columns with ztop_hi <= 0 get vf_prm = 0.
! ---------------------------------------------------------------------------
IMPLICIT NONE

INTEGER, INTENT(IN) :: row_length, rows, model_levels
REAL, INTENT(IN)  :: r_theta_levels(row_length, rows, 0:model_levels)
REAL, INTENT(IN)  :: ztop_lo(row_length, rows), ztop_hi(row_length, rows)
REAL, INTENT(OUT) :: vf_prm(row_length, rows, model_levels)

INTEGER :: i, j, k
REAL    :: zbot, ztop, z0, z1, total

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'FREITAS_VERT_PROFILE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

vf_prm(:,:,:) = 0.0

DO j = 1, rows
  DO i = 1, row_length

    IF (ztop_hi(i,j) <= 0.0) CYCLE              ! no plume: leave zeros

    zbot = MIN(ztop_lo(i,j), ztop_hi(i,j))
    ztop = MAX(ztop_lo(i,j), ztop_hi(i,j))

    ! UM mass layer k spans theta levels k-1 .. k (heights above surface).
    ! Weight = thickness of the overlap of layer k with [zbot, ztop].
    DO k = 1, model_levels
      z0 = r_theta_levels(i,j,k-1) - r_theta_levels(i,j,0)
      z1 = r_theta_levels(i,j,k)   - r_theta_levels(i,j,0)
      vf_prm(i,j,k) = MAX(0.0, MIN(z1, ztop) - MAX(z0, zbot))
    END DO

    total = SUM(vf_prm(i,j,:))
    IF (total > 0.0) THEN
      vf_prm(i,j,:) = vf_prm(i,j,:) / total
    ELSE
      ! Zero-thickness layer: put everything in the level containing ztop
      DO k = 1, model_levels
        IF (r_theta_levels(i,j,k) - r_theta_levels(i,j,0) >= ztop) THEN
          vf_prm(i,j,k) = 1.0
          EXIT
        END IF
      END DO
    END IF

  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

END SUBROUTINE freitas_vert_profile

END MODULE ukca_freitas_plumerise_mod
