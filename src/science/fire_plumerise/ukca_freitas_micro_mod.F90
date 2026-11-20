! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Microphysics of the Freitas et al. (2007, 2010) 1D plume-rise model:
!   sedimentation of rain and ice, and the water balance (evaporation,
!   condensation, sublimation, glaciation, melting, autoconversion and
!   accretion) after Kessler (1969), Ogura and Takahashi (1971) and
!   Berry (1968).
!
!   Units follow the original model: mixing ratios in g/g, density in
!   g m-3, latent heats in J g-1.
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
MODULE ukca_freitas_micro_mod

USE ukca_freitas_types_mod,    ONLY: freitas_column_type
USE ukca_um_legacy_mod,        ONLY: cp
USE ukca_config_constants_mod, ONLY: lc

IMPLICIT NONE

PRIVATE

PUBLIC :: freitas_fallpart, freitas_waterbal

! Plume-model specific constants
REAL, PARAMETER :: tfreeze  = 269.3    ! ice formation temperature (K)
REAL, PARAMETER :: tmelt    = 273.3    ! melting temperature (K)
REAL, PARAMETER :: heatsubl = 2834.0   ! latent heat of sublimation (J g-1)
REAL, PARAMETER :: heatfus  = 334.0    ! latent heat of fusion (J g-1)
REAL, PARAMETER :: herc     = 5.44e-4  ! rain evaporation constant
REAL, PARAMETER :: f0       = 0.75     ! ice velocity factor
REAL, PARAMETER :: qmin     = 1.0e-10  ! threshold to avoid underflow (g/g)

CONTAINS

! ---------------------------------------------------------------------------
SUBROUTINE freitas_fallpart(prm)
! ---------------------------------------------------------------------------
! Description:
!   Sedimentation (fallout) of rain/hail (qh) and ice (qi), and the
!   ventilation coefficients used for evaporation and sublimation.
! ---------------------------------------------------------------------------
IMPLICIT NONE

TYPE(freitas_column_type), INTENT(IN OUT) :: prm

INTEGER :: k
REAL    :: vhrel, virel, dfhz, dfiz, dz1

DO k = 2, prm%nm1 - 1

  ! Hydrometeor fall speeds (Kessler, table 4): small variation with q
  prm%vth(k) = -4.0
  vhrel = prm%w(k) + prm%vth(k)            ! relative to surrounding cloud

  ! Rain ventilation coefficient for evaporation
  prm%cvh(k) = 1.6 + 0.57e-3 * (ABS(vhrel))**1.5

  prm%vti(k) = -3.0
  virel = prm%w(k) + prm%vti(k)            ! relative to surrounding cloud

  ! Ice ventilation coefficient for sublimation
  prm%cvi(k) = 1.6 + 0.57e-3 * (ABS(virel))**1.5 / f0

  IF (vhrel >= 0.0) THEN
    dfhz = prm%qh(k) * (prm%rho(k)   * prm%vth(k)                            &
                      - prm%rho(k-1) * prm%vth(k-1)) / prm%rho(k-1)
  ELSE
    dfhz = prm%qh(k) * (prm%rho(k+1) * prm%vth(k+1)                          &
                      - prm%rho(k)   * prm%vth(k))   / prm%rho(k)
  END IF

  IF (virel >= 0.0) THEN
    dfiz = prm%qi(k) * (prm%rho(k)   * prm%vti(k)                            &
                      - prm%rho(k-1) * prm%vti(k-1)) / prm%rho(k-1)
  ELSE
    dfiz = prm%qi(k) * (prm%rho(k+1) * prm%vti(k+1)                          &
                      - prm%rho(k)   * prm%vti(k))   / prm%rho(k)
  END IF

  dz1 = prm%zm(k) - prm%zm(k-1)
  prm%qht(k) = prm%qht(k) - dfhz / dz1
  prm%qit(k) = prm%qit(k) - dfiz / dz1

END DO

END SUBROUTINE freitas_fallpart

! ---------------------------------------------------------------------------
SUBROUTINE freitas_waterbal(prm, l)
! ---------------------------------------------------------------------------
! Description:
!   Water balance at level l. Uses prm%wbar and prm%dqsdz, which are set
!   by the caller.
! ---------------------------------------------------------------------------
IMPLICIT NONE

TYPE(freitas_column_type), INTENT(IN OUT) :: prm
INTEGER,                   INTENT(IN)     :: l

IF (prm%qc(l) <= qmin) prm%qc(l) = 0.0      ! defeat underflow problem
IF (prm%qh(l) <= qmin) prm%qh(l) = 0.0
IF (prm%qi(l) <= qmin) prm%qi(l) = 0.0

CALL evaporate(prm, l)     ! vapour to cloud, cloud to vapour
CALL sublimate(prm, l)     ! vapour to ice
CALL glaciate(prm, l)      ! rain to ice
CALL melt(prm, l)          ! ice to rain
CALL convert(prm, l)       ! (auto)conversion and accretion

END SUBROUTINE freitas_waterbal

! ---------------------------------------------------------------------------
SUBROUTINE evaporate(prm, l)
! ---------------------------------------------------------------------------
! Description:
!   Condenses vapour, or evaporates cloud, rain and ice to saturation.
! ---------------------------------------------------------------------------
IMPLICIT NONE

TYPE(freitas_column_type), INTENT(IN OUT) :: prm
INTEGER,                   INTENT(IN)     :: l

REAL :: frc, src, heatcond
REAL :: evhdt, evidt, evrate, evap, sd, quant, dividend, divisor, devidt

heatcond = lc * 1.0e-3                     ! J g-1
frc = heatcond / (cp / 1000.0)
src = heatsubl / (cp / 1000.0)

sd = prm%qsat(l) - prm%qv(l)               ! vapour deficit
IF (sd == 0.0) RETURN

evhdt = 0.0
evidt = 0.0

evrate = ABS(prm%wbar * prm%dqsdz)         ! evaporation rate (Kessler 8.32)
evap   = evrate * prm%dt                   ! what we can get in dt

IF (sd <= 0.0) THEN
  ! Condense. sd is negative
  IF (evap >= ABS(sd)) THEN                ! we get it all
    prm%qc(l) = prm%qc(l) - sd
    prm%qv(l) = prm%qsat(l)                ! set the vapour to saturation
    prm%t(l)  = prm%t(l) - sd * frc        ! heat gained through condensation
  ELSE
    prm%qc(l) = prm%qc(l) + evap           ! get what we can in dt
    prm%qv(l) = prm%qv(l) - evap           ! remove it from the vapour
    prm%t(l)  = prm%t(l) + evap * frc      ! get some heat
  END IF
  RETURN
END IF

! Not saturated: saturate if possible, using cloud, then rain, then ice.
! sd is positive.

! -- Cloud --
IF (evap <= prm%qc(l)) THEN                ! enough cloud to last dt
  IF (sd <= evap) THEN                     ! enough time to saturate
    prm%qc(l) = prm%qc(l) - sd             ! remove cloud
    prm%qv(l) = prm%qsat(l)                ! saturate
    prm%t(l)  = prm%t(l) - sd * frc        ! cool the parcel
    RETURN
  ELSE                                     ! not enough time
    sd = sd - evap                         ! use what there is
    prm%qv(l) = prm%qv(l) + evap           ! add vapour
    prm%t(l)  = prm%t(l) - evap * frc      ! lose heat
    prm%qc(l) = prm%qc(l) - evap           ! lose cloud
  END IF
ELSE                                       ! not enough cloud to last dt
  IF (sd <= prm%qc(l)) THEN                ! but there is enough to saturate
    prm%qv(l) = prm%qsat(l)
    prm%qc(l) = prm%qc(l) - sd
    prm%t(l)  = prm%t(l) - sd * frc
    RETURN
  ELSE                                     ! not enough to saturate
    sd = sd - prm%qc(l)
    prm%qv(l) = prm%qv(l) + prm%qc(l)
    prm%t(l)  = prm%t(l) - prm%qc(l) * frc
    prm%qc(l) = 0.0                        ! all gone
  END IF
END IF

! -- Rain -- (W12, OT25, K table 4). Only dt to evaporate; rain first,
! then ice. Saturation may not be possible in one dt.
IF (prm%qh(l) > qmin) THEN
  quant = (prm%qsat(l) - prm%qc(l) - prm%qv(l)) * prm%rho(l)   ! g m-3
  evhdt = (prm%dt * herc * quant * (prm%qh(l) * prm%rho(l))**0.65)            &
          / prm%rho(l)                     ! rain evaporation in dt

  IF (evhdt <= prm%qh(l)) THEN             ! enough rain to last dt
    IF (sd <= evhdt) THEN                  ! enough time to saturate
      prm%qh(l) = prm%qh(l) - sd
      prm%qv(l) = prm%qsat(l)
      prm%t(l)  = prm%t(l) - sd * frc
      RETURN
    ELSE                                   ! not enough time
      sd = sd - evhdt
      prm%qv(l) = prm%qv(l) + evhdt
      prm%t(l)  = prm%t(l) - evhdt * frc
      prm%qh(l) = prm%qh(l) - evhdt
    END IF
  ELSE                                     ! not enough rain to last dt
    IF (sd <= prm%qh(l)) THEN              ! but there is enough to saturate
      prm%qv(l) = prm%qsat(l)
      prm%qh(l) = prm%qh(l) - sd
      prm%t(l)  = prm%t(l) - sd * frc
      RETURN
    ELSE                                   ! not enough to saturate
      sd = sd - prm%qh(l)
      prm%qv(l) = prm%qv(l) + prm%qh(l)
      prm%t(l)  = prm%t(l) - prm%qh(l) * frc
      prm%qh(l) = 0.0                      ! all gone
    END IF
  END IF
END IF

! -- Ice -- (OT, with unit correction factors)
IF (prm%qi(l) <= qmin) RETURN              ! no ice there

dividend = ((1.0e6 / prm%rho(l))**0.475) * (sd / prm%qsat(l) - 1.0)          &
           * (prm%qi(l)**0.525) * 1.13
divisor  = 7.0e5 + 4.1e6 / (10.0 * prm%est(l))
devidt   = -prm%cvi(l) * dividend / divisor   ! rate of change
evidt    = devidt * prm%dt                    ! what we could get

IF (evidt <= prm%qi(l)) THEN               ! enough ice to last dt
  IF (sd <= evidt) THEN                    ! enough time to saturate
    prm%qi(l) = prm%qi(l) - sd
    prm%qv(l) = prm%qsat(l)
    prm%t(l)  = prm%t(l) - sd * src
  ELSE                                     ! not enough time
    sd = sd - evidt
    prm%qv(l) = prm%qv(l) + evidt
    prm%t(l)  = prm%t(l) - evidt * src
    prm%qi(l) = prm%qi(l) - evidt
  END IF
ELSE                                       ! not enough ice to last dt
  IF (sd <= prm%qi(l)) THEN                ! but there is enough to saturate
    prm%qv(l) = prm%qsat(l)
    prm%qi(l) = prm%qi(l) - sd
    prm%t(l)  = prm%t(l) - sd * src
  ELSE                                     ! not enough to saturate
    sd = sd - prm%qi(l)
    prm%qv(l) = prm%qv(l) + prm%qi(l)
    prm%t(l)  = prm%t(l) - prm%qi(l) * src
    prm%qi(l) = 0.0                        ! all gone
  END IF
END IF

END SUBROUTINE evaporate

! ---------------------------------------------------------------------------
SUBROUTINE sublimate(prm, l)
! ---------------------------------------------------------------------------
! Description:
!   Vapour to ice (Ogura and Takahashi, eq. 22).
! ---------------------------------------------------------------------------
IMPLICIT NONE

TYPE(freitas_column_type), INTENT(IN OUT) :: prm
INTEGER,                   INTENT(IN)     :: l

REAL :: src, dtsubh, dividend, divisor, subl

src = heatsubl / (cp / 1000.0)

! Selection criteria for sublimation
IF (prm%t(l)  >  tfreeze)      RETURN
IF (prm%qv(l) <= prm%qsat(l))  RETURN

dividend = ((1.0e6 / prm%rho(l))**0.475) * (prm%qv(l) / prm%qsat(l) - 1.0)   &
           * (prm%qi(l)**0.525) * 1.13
divisor  = 7.0e5 + 4.1e6 / (10.0 * prm%est(l))

dtsubh = ABS(dividend / divisor)           ! sublimation rate
subl   = dtsubh * prm%dt                   ! and amount possible

IF (subl < prm%qv(l)) THEN
  prm%qv(l) = prm%qv(l) - subl             ! lose vapour
  prm%qi(l) = prm%qi(l) + subl             ! gain ice
  prm%t(l)  = prm%t(l) + subl * src        ! energy change, warms air
ELSE
  ! Use what there is. Note: as in the original model, qi is set to qv
  ! (not incremented) here.
  prm%qi(l) = prm%qv(l)
  prm%t(l)  = prm%t(l) + prm%qv(l) * src   ! warm the air
  prm%qv(l) = 0.0
END IF

END SUBROUTINE sublimate

! ---------------------------------------------------------------------------
SUBROUTINE glaciate(prm, l)
! ---------------------------------------------------------------------------
! Description:
!   Conversion of rain to ice (Ogura and Takahashi, eq. 16).
! ---------------------------------------------------------------------------
IMPLICIT NONE

TYPE(freitas_column_type), INTENT(IN OUT) :: prm
INTEGER,                   INTENT(IN)     :: l

REAL, PARAMETER :: glconst = 0.025         ! glaciation time constant (s-1)
REAL :: frc, dfrzh

frc = heatfus / (cp / 1000.0)

! Selection rules for glaciation
IF (prm%qh(l) <= 0.0)          RETURN
IF (prm%qv(l) <  prm%qsat(l))  RETURN
IF (prm%t(l)  >  tfreeze)      RETURN

dfrzh = prm%dt * glconst * prm%qh(l)       ! rate of mass gain in ice

IF (dfrzh < prm%qh(l)) THEN
  prm%qi(l) = prm%qi(l) + dfrzh
  prm%qh(l) = prm%qh(l) - dfrzh
  prm%t(l)  = prm%t(l) + frc * dfrzh       ! warms air
ELSE
  prm%qi(l) = prm%qi(l) + prm%qh(l)
  prm%t(l)  = prm%t(l) + frc * prm%qh(l)
  prm%qh(l) = 0.0
END IF

END SUBROUTINE glaciate

! ---------------------------------------------------------------------------
SUBROUTINE melt(prm, l)
! ---------------------------------------------------------------------------
! Description:
!   Melting of ice to rain (after Mason, 1956).
! ---------------------------------------------------------------------------
IMPLICIT NONE

TYPE(freitas_column_type), INTENT(IN OUT) :: prm
INTEGER,                   INTENT(IN)     :: l

REAL, PARAMETER :: frc = 332.27            ! heatfus / cp (K per g/g)
REAL :: dtmelt

! Selection rules
IF (prm%qi(l) <= 0.0)   RETURN
IF (prm%t(l)  <  tmelt) RETURN

dtmelt = prm%dt * (2.27 / prm%rho(l)) * prm%cvi(l) * (prm%t(l) - tmelt)       &
         * ((prm%rho(l) * prm%qi(l) * 1.0e-6)**0.525) * (f0**(-0.42))

IF (dtmelt < prm%qi(l)) THEN
  prm%qh(l) = prm%qh(l) + dtmelt
  prm%qi(l) = prm%qi(l) - dtmelt
  prm%t(l)  = prm%t(l) - frc * dtmelt      ! cools air
ELSE
  prm%qh(l) = prm%qh(l) + prm%qi(l)        ! get all there is to get
  prm%t(l)  = prm%t(l) - frc * prm%qi(l)
  prm%qi(l) = 0.0
END IF

END SUBROUTINE melt

! ---------------------------------------------------------------------------
SUBROUTINE convert(prm, l)
! ---------------------------------------------------------------------------
! Description:
!   Autoconversion (Berry) and accretion (Kessler) of cloud water to rain.
! ---------------------------------------------------------------------------
IMPLICIT NONE

TYPE(freitas_column_type), INTENT(IN OUT) :: prm
INTEGER,                   INTENT(IN)     :: l

REAL, PARAMETER :: ak2    = 0.0052    ! collection (accretion) rate
REAL, PARAMETER :: anbase = 100000.0  ! Berry number at cloud base
                                      ! (continental, # cm-3 in this formula)
REAL, PARAMETER :: bdisp  = 0.146     ! Berry size dispersion (continental)

REAL :: accrete, con, q, h, total

! Selection rules
IF (prm%t(l)  <= tfreeze) RETURN           ! process not allowed above ice
IF (prm%qc(l) == 0.0)     RETURN

accrete = 0.0
q = prm%rho(l) * prm%qc(l)
h = prm%rho(l) * prm%qh(l)

IF (prm%qh(l) > 0.0) accrete = ak2 * q * (h**0.875)   ! accretion, Kessler

! Autoconversion, Berry
con = q * q * q * bdisp / (60.0 * (5.0 * q * bdisp + 0.0366 * anbase))

total = (con + accrete) * prm%dt / prm%rho(l)

IF (total < prm%qc(l)) THEN
  prm%qc(l) = prm%qc(l) - total
  prm%qh(l) = prm%qh(l) + total            ! no phase change involved
ELSE
  prm%qh(l) = prm%qh(l) + prm%qc(l)        ! uses all there is
  prm%qc(l) = 0.0
END IF

END SUBROUTINE convert

END MODULE ukca_freitas_micro_mod
