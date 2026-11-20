! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Dynamics and advection of the Freitas et al. (2007, 2010) 1D
!   plume-rise model: advection of w and scalars, lateral and dynamic
!   entrainment, buoyancy, eddy viscosity, sponge layer and the
!   leapfrog/Asselin time stepping (after RAMS).
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
MODULE ukca_freitas_dyn_mod

USE ukca_freitas_types_mod,    ONLY: freitas_column_type, nkp_prm
USE ukca_um_legacy_mod,        ONLY: g, cp
USE ukca_config_constants_mod, ONLY: repsilon
USE ukca_constants,            ONLY: pi

IMPLICIT NONE

PRIVATE

PUBLIC :: freitas_tend0, freitas_vel_advect, freitas_scl_advect,             &
          freitas_scl_misc, freitas_scl_dyn_entrain, freitas_damp_grav_wave,  &
          freitas_buoyancy, freitas_entrainment, freitas_visc_w,              &
          freitas_update, freitas_hadvance

! Virtual-mass coefficient: compensates for the missing acceleration from
! non-hydrostatic pressure perturbations
REAL, PARAMETER :: gama_vir = 0.5

CONTAINS

! ---------------------------------------------------------------------------
SUBROUTINE freitas_tend0(prm)
! ---------------------------------------------------------------------------
! Description:
!   Sets all tendencies to zero.
! ---------------------------------------------------------------------------
IMPLICIT NONE

TYPE(freitas_column_type), INTENT(IN OUT) :: prm

prm%wt   (1:prm%nm1) = 0.0
prm%tt   (1:prm%nm1) = 0.0
prm%qvt  (1:prm%nm1) = 0.0
prm%qct  (1:prm%nm1) = 0.0
prm%qht  (1:prm%nm1) = 0.0
prm%qit  (1:prm%nm1) = 0.0
prm%vel_t(1:prm%nm1) = 0.0
prm%rad_t(1:prm%nm1) = 0.0
prm%sct  (1:prm%nm1) = 0.0

END SUBROUTINE freitas_tend0

! ---------------------------------------------------------------------------
SUBROUTINE freitas_vel_advect(prm)
! ---------------------------------------------------------------------------
! Description:
!   Advection tendency of w.
! ---------------------------------------------------------------------------
IMPLICIT NONE

TYPE(freitas_column_type), INTENT(IN OUT) :: prm

INTEGER :: k, m1
REAL    :: flxw(nkp_prm), dn0(nkp_prm)
REAL, PARAMETER :: c1z = 0.5

m1 = prm%nm1

dn0(1:m1) = prm%rho(1:m1) * 1.0e-3         ! g m-3 to kg m-3

flxw(1) = prm%wc(1) * dn0(1)
DO k = 2, m1 - 1
  flxw(k) = prm%wc(k) * 0.5 * (dn0(k) + dn0(k+1))
END DO

DO k = 2, m1 - 2
  prm%wt(k) = prm%wt(k) + c1z * prm%dzm(k) / (dn0(k) + dn0(k+1)) *           &
              ((flxw(k) + flxw(k-1)) * (prm%wc(k) + prm%wc(k-1))              &
             - (flxw(k) + flxw(k+1)) * (prm%wc(k) + prm%wc(k+1))              &
             + (flxw(k+1) - flxw(k-1)) * 2.0 * prm%wc(k))
END DO

END SUBROUTINE freitas_vel_advect

! ---------------------------------------------------------------------------
SUBROUTINE freitas_scl_advect(prm)
! ---------------------------------------------------------------------------
! Description:
!   Advection tendencies of the scalars (t, qv, qc, qi, qh, vel_p, rad_p,
!   sc), using a positive-definite flux scheme (after RAMS).
! ---------------------------------------------------------------------------
IMPLICIT NONE

TYPE(freitas_column_type), INTENT(IN OUT) :: prm

INTEGER :: k, mzp
REAL    :: dtlto2
REAL    :: vt3dc(nkp_prm), vt3df(nkp_prm), vt3dg(nkp_prm), vt3dk(nkp_prm)
REAL    :: vctr1(nkp_prm), vctr2(nkp_prm)

mzp = prm%nm1

dtlto2   = 0.5 * prm%dt
vt3dc(1) = (prm%w(1) + prm%wc(1)) * dtlto2 * prm%rho(1) * 1.0e-3
vt3df(1) = 0.5 * (prm%w(1) + prm%wc(1)) * dtlto2 * prm%dzm(1)

DO k = 2, mzp
  vt3dc(k) = (prm%w(k) + prm%wc(k)) * dtlto2 * 0.5                           &
             * (prm%rho(k) + prm%rho(k+1)) * 1.0e-3
  vt3df(k) = (prm%w(k) + prm%wc(k)) * dtlto2 * 0.5 * prm%dzm(k)
END DO

DO k = 1, mzp
  vctr1(k) = (prm%zt(k+1) - prm%zm(k)) * prm%dzm(k)
  vctr2(k) = (prm%zm(k)   - prm%zt(k)) * prm%dzm(k)
  vt3dk(k) = prm%dzt(k) / (prm%rho(k) * 1.0e-3)
END DO

CALL advect_one(mzp, prm%t,     prm%tt)
CALL advect_one(mzp, prm%qv,    prm%qvt)
CALL advect_one(mzp, prm%qc,    prm%qct)
CALL advect_one(mzp, prm%qi,    prm%qit)
CALL advect_one(mzp, prm%qh,    prm%qht)
CALL advect_one(mzp, prm%vel_p, prm%vel_t)
CALL advect_one(mzp, prm%rad_p, prm%rad_t)
CALL advect_one(mzp, prm%sc,    prm%sct)

CONTAINS

  ! Advection tendency of one scalar
  SUBROUTINE advect_one(m1, scp, sct)
  IMPLICIT NONE
  INTEGER, INTENT(IN)     :: m1
  REAL,    INTENT(IN)     :: scp(nkp_prm)  ! scalar
  REAL,    INTENT(IN OUT) :: sct(nkp_prm)  ! its tendency
  REAL :: scr1(nkp_prm), dtli
  INTEGER :: kk
  REAL, PARAMETER :: dfact = 0.5

  scr1(:) = scp(:)

  ! Scalar flux
  DO kk = 1, m1 - 1
    vt3dg(kk) = vt3dc(kk) * (vctr1(kk) * scr1(kk) + vctr2(kk) * scr1(kk+1)   &
              + vt3df(kk) * (scr1(kk) - scr1(kk+1)))
  END DO

  ! Keep scalars positive-definite: if a flux would remove more than half
  ! the quantity in a timestep, reduce it to first order.
  DO kk = 1, m1 - 1
    IF (vt3dc(kk) > 0.0) THEN
      IF (vt3dg(kk) * vt3dk(kk) > dfact * scr1(kk)) THEN
        vt3dg(kk) = vt3dc(kk) * scr1(kk)
      END IF
    ELSE IF (vt3dc(kk) < 0.0) THEN
      IF (-vt3dg(kk) * vt3dk(kk+1) > dfact * scr1(kk+1)) THEN
        vt3dg(kk) = vt3dc(kk) * scr1(kk+1)
      END IF
    END IF
  END DO

  ! Flux divergence
  DO kk = 2, m1 - 1
    scr1(kk) = scr1(kk) + vt3dk(kk) * (vt3dg(kk-1) - vt3dg(kk)               &
             + scp(kk) * (vt3dc(kk) - vt3dc(kk-1)))
  END DO

  ! Tendency
  dtli = 1.0 / prm%dt
  DO kk = 2, m1 - 1
    sct(kk) = sct(kk) + (scr1(kk) - scp(kk)) * dtli
  END DO

  END SUBROUTINE advect_one

END SUBROUTINE freitas_scl_advect

! ---------------------------------------------------------------------------
SUBROUTINE freitas_scl_misc(prm)
! ---------------------------------------------------------------------------
! Description:
!   Scalar tendencies from the dry adiabat and from lateral (turbulent)
!   entrainment. Also sets prm%wbar (used later by the microphysics at the
!   lower boundary, as in the original model).
! ---------------------------------------------------------------------------
IMPLICIT NONE

TYPE(freitas_column_type), INTENT(IN OUT) :: prm

INTEGER :: k
REAL    :: adiabat, dmdtm

DO k = 2, prm%nm1 - 1
  prm%wbar = 0.5 * (prm%w(k) + prm%w(k-1))

  adiabat = -prm%wbar * g / cp                                   ! dry adiabat
  dmdtm   = 2.0 * prm%alpha * ABS(prm%wbar) / prm%radius(k)      ! (1/M)dM/dt

  prm%tt(k)    = prm%tt(k) + adiabat - dmdtm * (prm%t(k) - prm%te(k))
  prm%qvt(k)   = prm%qvt(k) - dmdtm * (prm%qv(k) - prm%qvenv(k))
  prm%qct(k)   = prm%qct(k) - dmdtm * prm%qc(k)
  prm%qht(k)   = prm%qht(k) - dmdtm * prm%qh(k)
  prm%qit(k)   = prm%qit(k) - dmdtm * prm%qi(k)
  prm%vel_t(k) = prm%vel_t(k) - dmdtm * (prm%vel_p(k) - prm%vel_e(k))
  prm%rad_t(k) = prm%rad_t(k) + 0.5 * dmdtm * (6.0 / 5.0) * prm%radius(k)
END DO

END SUBROUTINE freitas_scl_misc

! ---------------------------------------------------------------------------
SUBROUTINE freitas_scl_dyn_entrain(prm)
! ---------------------------------------------------------------------------
! Description:
!   Scalar tendencies from dynamic entrainment, driven by the horizontal
!   velocity difference between plume and environment.
! ---------------------------------------------------------------------------
IMPLICIT NONE

TYPE(freitas_column_type), INTENT(IN OUT) :: prm

INTEGER :: k
REAL    :: dmdtm

DO k = 2, prm%nm1 - 1
  prm%rad_t(k) = prm%rad_t(k) + ABS(prm%vel_e(k) - prm%vel_p(k)) / pi

  dmdtm = (2.0 / pi) * ABS(prm%vel_e(k) - prm%vel_p(k)) / prm%radius(k)

  prm%vel_t(k) = prm%vel_t(k) - dmdtm * (prm%vel_p(k) - prm%vel_e(k))
  prm%tt(k)    = prm%tt(k)    - dmdtm * (prm%t(k) - prm%te(k))
  prm%qvt(k)   = prm%qvt(k)   - dmdtm * (prm%qv(k) - prm%qvenv(k))
  prm%qct(k)   = prm%qct(k)   - dmdtm * prm%qc(k)
  prm%qht(k)   = prm%qht(k)   - dmdtm * prm%qh(k)
  prm%qit(k)   = prm%qit(k)   - dmdtm * prm%qi(k)
END DO

END SUBROUTINE freitas_scl_dyn_entrain

! ---------------------------------------------------------------------------
SUBROUTINE freitas_damp_grav_wave(prm, ifrom, deltak)
! ---------------------------------------------------------------------------
! Description:
!   Rayleigh-friction (sponge) layer at the top of the integration domain,
!   so that gravity waves are absorbed instead of reflected.
!     ifrom = 1: relax t towards te (through the tendency tt)
!     ifrom = 2: relax w towards 0 (direct update)
! ---------------------------------------------------------------------------
IMPLICIT NONE

TYPE(freitas_column_type), INTENT(IN OUT) :: prm
INTEGER,                   INTENT(IN)     :: ifrom
INTEGER,                   INTENT(IN)     :: deltak

INTEGER :: k, kf
REAL    :: zmkf, ztop, distim, c1, c2

kf = prm%nm1 - deltak
IF (kf < 1) RETURN

zmkf = prm%zm(kf)                          ! height of sponge base
ztop = prm%zm(prm%nm1)                     ! height of model top

distim = MIN(3.0 * prm%dt, 60.0)           ! relaxation time scale at the top
c1 = 1.0 / (distim * (ztop - zmkf))        ! tendency coefficient
c2 = prm%dt * c1                           ! direct-update coefficient

SELECT CASE (ifrom)
CASE (1)
  DO k = prm%nm1, 2, -1
    IF (prm%zt(k) <= zmkf) EXIT
    prm%tt(k) = prm%tt(k) + c1 * (prm%zt(k) - zmkf) * (prm%te(k) - prm%t(k))
  END DO
CASE (2)
  DO k = prm%nm1, 2, -1
    IF (prm%zt(k) <= zmkf) EXIT
    prm%w(k) = prm%w(k) + c2 * (prm%zt(k) - zmkf) * (0.0 - prm%w(k))
  END DO
END SELECT

END SUBROUTINE freitas_damp_grav_wave

! ---------------------------------------------------------------------------
SUBROUTINE freitas_buoyancy(prm)
! ---------------------------------------------------------------------------
! Description:
!   Buoyancy acceleration of the plume relative to the environment, added
!   to the w tendency. Hydrometeor loading acts as drag.
! ---------------------------------------------------------------------------
IMPLICIT NONE

TYPE(freitas_column_type), INTENT(IN OUT) :: prm

INTEGER :: k, m1
REAL    :: tv, tve, qwtotl, umgamai
REAL    :: scr1(nkp_prm)

m1 = prm%nm1
umgamai = 1.0 / (1.0 + gama_vir)

DO k = 2, m1 - 1
  tv  = prm%t(k)  * (1.0 + (prm%qv(k)    / repsilon)) / (1.0 + prm%qv(k))
  tve = prm%te(k) * (1.0 + (prm%qvenv(k) / repsilon)) / (1.0 + prm%qvenv(k))
  qwtotl  = prm%qh(k) + prm%qi(k) + prm%qc(k)          ! qwtotl*g is drag
  scr1(k) = g * umgamai * ((tv - tve) / tve - qwtotl)
END DO

DO k = 2, m1 - 2
  prm%rbuoy(k) = 0.5 * (scr1(k) + scr1(k+1))
  prm%wt(k)    = prm%wt(k) + 0.5 * (scr1(k) + scr1(k+1))
END DO
prm%rbuoy(1) = prm%rbuoy(2)

END SUBROUTINE freitas_buoyancy

! ---------------------------------------------------------------------------
SUBROUTINE freitas_entrainment(prm)
! ---------------------------------------------------------------------------
! Description:
!   Deceleration of w by entrainment of environmental air:
!   1) lateral (turbulent) entrainment,
!        (1/M) dM/dt = 1/(1+gama_vir) * 2 * alpha * |w| / R
!   2) dynamic entrainment, driven by the horizontal velocity difference
!      between plume and environment,
!        (2/pi) * 0.5 * |(vel_p - vel_e)(k) + (vel_p - vel_e)(k-1)| / R
!   R is the plume radius at the w level.
! ---------------------------------------------------------------------------
IMPLICIT NONE

TYPE(freitas_column_type), INTENT(IN OUT) :: prm

INTEGER :: k
REAL    :: dmdtm, wbar, radius_bar, umgamai, dyn_entr

umgamai = 1.0 / (1.0 + gama_vir)

DO k = 2, prm%nm1 - 1
  wbar       = prm%w(k)
  radius_bar = 0.5 * (prm%radius(k) + prm%radius(k-1))

  ! Lateral entrainment
  dmdtm = umgamai * 2.0 * prm%alpha * ABS(wbar) / radius_bar
  prm%wt(k) = prm%wt(k) - dmdtm * wbar

  ! Dynamic entrainment
  dyn_entr = (2.0 / pi) * 0.5 * ABS(prm%vel_p(k) - prm%vel_e(k)              &
           + prm%vel_p(k-1) - prm%vel_e(k-1)) / radius_bar
  prm%wt(k) = prm%wt(k) - dyn_entr * wbar
END DO

END SUBROUTINE freitas_entrainment

! ---------------------------------------------------------------------------
SUBROUTINE freitas_visc_w(prm, kmt)
! ---------------------------------------------------------------------------
! Description:
!   Vertical diffusion (eddy viscosity) tendencies of w, t, water species,
!   horizontal plume velocity and plume radius, for numerical stability.
! ---------------------------------------------------------------------------
IMPLICIT NONE

TYPE(freitas_column_type), INTENT(IN OUT) :: prm
INTEGER,                   INTENT(IN)     :: kmt

INTEGER :: k, kf
REAL    :: dz1t, dz1m, dz2t, dz2m

kf = MIN(prm%nm1, kmt - 1)

DO k = 2, kf
  dz1t = 0.5 * (prm%zt(k+1) - prm%zt(k-1))
  dz2t = prm%visc(k) / (dz1t * dz1t)
  dz1m = 0.5 * (prm%zm(k+1) - prm%zm(k-1))
  dz2m = prm%visc(k) / (dz1m * dz1m)

  prm%wt(k)    = prm%wt(k)    + d2dz(prm%w,     dz2m)
  prm%tt(k)    = prm%tt(k)    + d2dz(prm%t,     dz2t)
  prm%qvt(k)   = prm%qvt(k)   + d2dz(prm%qv,    dz2t)
  prm%qct(k)   = prm%qct(k)   + d2dz(prm%qc,    dz2t)
  prm%qht(k)   = prm%qht(k)   + d2dz(prm%qh,    dz2t)
  prm%qit(k)   = prm%qit(k)   + d2dz(prm%qi,    dz2t)
  prm%vel_t(k) = prm%vel_t(k) + d2dz(prm%vel_p, dz2t)
  prm%rad_t(k) = prm%rad_t(k) + d2dz(prm%rad_p, dz2t)
END DO

CONTAINS

  ! Second derivative of x at level k, times the factor fac
  REAL FUNCTION d2dz(x, fac)
  IMPLICIT NONE
  REAL, INTENT(IN) :: x(nkp_prm), fac
  d2dz = (x(k+1) - 2 * x(k) + x(k-1)) * fac
  END FUNCTION d2dz

END SUBROUTINE freitas_visc_w

! ---------------------------------------------------------------------------
SUBROUTINE freitas_update(prm, varn)
! ---------------------------------------------------------------------------
! Description:
!   Applies the tendencies: varn = 'W' for w, anything else for the scalars.
! ---------------------------------------------------------------------------
IMPLICIT NONE

TYPE(freitas_column_type), INTENT(IN OUT) :: prm
CHARACTER(LEN=*),          INTENT(IN)     :: varn

INTEGER :: k

IF (varn == 'W') THEN
  DO k = 2, prm%nm1 - 1
    prm%w(k) = prm%w(k) + prm%wt(k) * prm%dt
  END DO
ELSE
  DO k = 2, prm%nm1 - 1
    prm%t(k)  = prm%t(k)  + prm%tt(k)  * prm%dt
    prm%qv(k) = prm%qv(k) + prm%qvt(k) * prm%dt
    prm%qc(k) = prm%qc(k) + prm%qct(k) * prm%dt   ! cloud drops travel with air
    prm%qh(k) = prm%qh(k) + prm%qht(k) * prm%dt
    prm%qi(k) = prm%qi(k) + prm%qit(k) * prm%dt

    prm%qv(k) = MAX(0.0, prm%qv(k))
    prm%qc(k) = MAX(0.0, prm%qc(k))
    prm%qh(k) = MAX(0.0, prm%qh(k))
    prm%qi(k) = MAX(0.0, prm%qi(k))

    prm%vel_p(k) = prm%vel_p(k) + prm%vel_t(k) * prm%dt
    prm%rad_p(k) = prm%rad_p(k) + prm%rad_t(k) * prm%dt
    prm%sc(k)    = prm%sc(k)    + prm%sct(k)   * prm%dt
  END DO
END IF

END SUBROUTINE freitas_update

! ---------------------------------------------------------------------------
SUBROUTINE freitas_hadvance(prm, iac)
! ---------------------------------------------------------------------------
! Description:
!   Asselin filter for w (after RAMS), applied in two stages:
!     iac = 1: {wc} = wc + eps * (w - 2 * wc)
!     iac = 2: completes the filter and advances wc to the new time level.
! ---------------------------------------------------------------------------
IMPLICIT NONE

TYPE(freitas_column_type), INTENT(IN OUT) :: prm
INTEGER,                   INTENT(IN)     :: iac

INTEGER :: m
REAL    :: eps2, af

eps2 = 0.2
IF (prm%mintime == 1) eps2 = 0.5

IF (iac == 1) THEN
  DO m = 1, prm%nm1
    prm%wc(m) = prm%wc(m) + eps2 * (prm%w(m) - 2.0 * prm%wc(m))
  END DO
ELSE IF (iac == 2) THEN
  DO m = 1, prm%nm1
    af        = prm%w(m)
    prm%w(m)  = prm%wc(m) + eps2 * af
    prm%wc(m) = af
  END DO
END IF

END SUBROUTINE freitas_hadvance

END MODULE ukca_freitas_dyn_mod
