!WRF:MODEL_LAYER:CONSTANTS
!

module module_model_constants

   !  2. Following are constants for use in defining real number bounds.

   !  A really small number.

   real, parameter :: epsilon = 1.e-15

   !  4. Following is information related to the physical constants.

   !  These are the physical constants used within the model.

! JM NOTE -- can we name this grav instead?
   real, parameter :: g = 9.81 ! acceleration due to gravity (m {s}^-2)

   real, parameter :: r_d = 287. ! gas constant of dry air (J deg^-1 kg^-1)
   real, parameter :: cp = 7.*r_d/2. !

   real, parameter :: r_v = 461.6 ! gas constant for water vapor (J deg^-1 kg^-1)
   real, parameter :: cv = cp - r_d ! Specific heat of air at contant volume (J deg^-1 kg^-1)
   real, parameter :: cpv = 4.*r_v
   real, parameter :: cvv = cpv - r_v !
   real, parameter :: cvpm = -cv/cp
   real, parameter :: cliq = 4190. ! specific heat of liquid water at 0^oC
   real, parameter :: cice = 2106. ! specific heat of ice at 0^oC
   real, parameter :: psat = 610.78
   real, parameter :: rcv = r_d/cv !
   real, parameter :: rcp = r_d/cp
   real, parameter :: rovg = r_d/g
   real, parameter :: c2 = cp*rcv
   real, parameter :: mwdry = 28.966 ! molecular weight of dry air (g/mole)

   real, parameter :: p1000mb = 100000. ! pressure at 1000 hPa (pa)
   real, parameter :: t0 = 300. ! base state tempertaure (K)
   real, parameter :: p0 = p1000mb ! base state surface pressure (pa)
   real, parameter :: cpovcv = cp/(cp - r_d)
   real, parameter :: cvovcp = 1./cpovcv
   real, parameter :: rvovrd = r_v/r_d

   real, parameter :: reradius = 1./6370.0e03 ! reciprocal of earth radius (m^-1)

   real, parameter :: asselin = .025
!   REAL    , PARAMETER :: asselin      = .0
   real, parameter :: cb = 25.

   real, parameter :: XLV0 = 3.15e6 !  constant defined for calculation of latent heating
   real, parameter :: XLV1 = 2370. !  constant defined for calculation of latent heating
   real, parameter :: XLS0 = 2.905e6 !  constant defined for calculation of latent heating
   real, parameter :: XLS1 = 259.532 !  constant defined for calculation of latent heating

   real, parameter :: XLS = 2.85e6 ! latent heat of sublimation of water at 0^oC (J kg^-1)
   real, parameter :: XLV = 2.5e6 ! latent heat of vaporization of water at 0^oC (J kg^-1)
   real, parameter :: XLF = 3.50e5 ! latent heat of fusion of water at 0^oC (J kg^-1)

   real, parameter :: rhowater = 1000. ! density of liquid water at 0^oC (kg m^-3)
   real, parameter :: rhosnow = 100. ! density of snow (kg m^-3)
   real, parameter :: rhoair0 = 1.28 ! density of dry air at 0^oC and 1000mb pressure (kg m^-3)

   real, parameter :: RE_QC_BG = 2.49e-6 ! effective radius of cloud for background (m)
   real, parameter :: RE_QI_BG = 4.99e-6 ! effective radius of ice for background (m)
   real, parameter :: RE_QS_BG = 9.99e-6 ! effective radius of snow for background (m)
!
! Now namelist-specified parameter: ccn_conc - RAS
!   REAL    , PARAMETER :: n_ccn0       = 1.0E8
!
   real, parameter :: piconst = 3.1415926535897932384626433 ! constant of PI
   real, parameter :: DEGRAD = piconst/180. ! radians per degree
   real, parameter :: DPD = 360./365.

   real, parameter ::  SVP1 = 0.6112 ! constant for saturation vapor pressure calculation (dimensionless)
   real, parameter ::  SVP2 = 17.67 ! constant for saturation vapor pressure calculation (dimensionless)
   real, parameter ::  SVP3 = 29.65 ! constant for saturation vapor pressure calculation (K)
   real, parameter ::  SVPT0 = 273.15 ! constant for saturation vapor pressure calculation (K)
   real, parameter ::  EP_1 = R_v/R_d - 1. !  constant for virtual temperature (r_v/r_d - 1) (dimensionless)
   real, parameter ::  EP_2 = R_d/R_v ! constant for specific humidity calculation (dimensionless)
   real, parameter ::  KARMAN = 0.4 ! von Karman constant
   real, parameter ::  EOMEG = 7.2921e-5 ! angular velocity of rotation (rad^-1)
   real, parameter ::  STBOLT = 5.67051e-8 ! Stefan-Boltzmann constant (W m^-2 deg^-4)

   real, parameter ::  prandtl = 1./3.0 ! prandtl's mixing length (m)
   ! constants for w-damping option
   real, parameter ::  w_alpha = 0.3 ! strength m/s/s
   real, parameter ::  w_beta = 1.0 ! activation cfl number

   real, parameter ::  pq0 = 379.90516 !
   real, parameter ::  epsq2 = 0.2 ! initial TKE for camuw PBL scheme (m2 s^-2)
   real, parameter ::  a2 = 17.2693882
   real, parameter ::  a3 = 273.16
   real, parameter ::  a4 = 35.86
   real, parameter ::  epsq = 1.e-12 ! threshold specified for SPECIFIC HUMIDITY calculation in BMJ cumulus scheme (kg kg^-1)
   real, parameter ::  p608 = rvovrd - 1.
   real, parameter ::  climit = 1.e-20
   real, parameter ::  cm1 = 2937.4
   real, parameter ::  cm2 = 4.9283
   real, parameter ::  cm3 = 23.5518
!       REAL , PARAMETER ::  defc=8.0
!       REAL , PARAMETER ::  defm=32.0
   real, parameter ::  defc = 0.0
   real, parameter ::  defm = 99999.0
   real, parameter ::  epsfc = 1./1.05
   real, parameter ::  epswet = 0.0
   real, parameter ::  fcdif = 1./3.
   real, parameter ::  fcm = 0.00003
   real, parameter ::  gma = -r_d*(1.-rcp)*0.5
   real, parameter ::  p400 = 40000.0
   real, parameter ::  phitp = 15000.0
   real, parameter ::  pi2 = 2.*3.1415926, pi1 = 3.1415926
   real, parameter ::  plbtm = 105000.0
   real, parameter ::  plomd = 64200.0
   real, parameter ::  pmdhi = 35000.0
   real, parameter ::  q2ini = 0.50
   real, parameter ::  rfcp = 0.25/cp
   real, parameter ::  rhcrit_land = 0.75
   real, parameter ::  rhcrit_sea = 0.80
   real, parameter ::  rlag = 14.8125
   real, parameter ::  rlx = 0.90
   real, parameter ::  scq2 = 50.0
   real, parameter ::  slopht = 0.001
   real, parameter ::  tlc = 2.*0.703972477
   real, parameter ::  wa = 0.15
   real, parameter ::  wght = 0.35
   real, parameter ::  wpc = 0.075
   real, parameter ::  z0land = 0.10 ! surface roughness length over land (m)
   real, parameter ::  z0max = 0.008 !  maximum roughness length (m)
   real, parameter ::  z0sea = 0.001 ! roughness length over ocean (m)

   !  Earth

   !  The value for P2SI *must* be set to 1.0 for Earth
   !  Although, now we may not need this declaration here (see above)
   !REAL    , PARAMETER :: P2SI         = 1.0

   !  Orbital constants:

   integer, parameter :: PLANET_YEAR = 365 ! number of days in a calendar year
   real, parameter :: OBLIQUITY = 23.5 ! solar obliquity (degree)
   real, parameter :: ECCENTRICITY = 0.014 ! Orbital eccentricity
   real, parameter :: SEMIMAJORAXIS = 1.0 ! Ratio of semi-major axis of planet / semi-major axis of earth
   real, parameter :: zero_date = 0.0 ! Time of perihelion passage
   real, parameter :: EQUINOX_FRACTION = 0.0 ! Fraction into the year (from perhelion) of the occurrence of the Northern Spring Equinox

! 2012103
#if (EM_CORE == 1)
! for calls to set_tiles
   integer, parameter :: ZONE_SOLVE_EM = 1
   integer, parameter :: ZONE_SFS = 2
#endif

contains
   subroutine init_module_model_constants
   end subroutine init_module_model_constants
end module module_model_constants
