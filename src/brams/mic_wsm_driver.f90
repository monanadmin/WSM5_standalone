!==========================================================================================
!- This code is a 3, 5, 6 and 7-class species microphysics scheme of the 
!  Single-Moment MicroPhyiscs (WSMMP).
!- Adapted to BRAMS 6.0+ by Saulo Freitas Apr/2022
!- version from WRF 4.3.3
!==========================================================================================

SUBROUTINE micro_wsm( )

  use mem_basic, only : &
       basic_g            ! INTENT(INOUT)

  use mem_micro, only:  &
       micro_g            ! INTENT(INOUT)

  use mem_grid, only:   &
       ngrids,          & ! INTENT(IN)
       ngrid,           & ! INTENT(IN)
       zm,              & ! INTENT(IN)
       dzt,             & ! INTENT(IN)
       dtlt,            & ! INTENT(IN)
       jdim,            & ! INTENT(IN)
       maxnzp,          & ! INTENT(IN)
       time,            & ! INTENT(IN)
       zt,              & ! INTENT(IN)
       itime1,          & ! INTENT(IN)
       if_adap,         & ! INTENT(IN)
       grid_g,          & ! INTENT(IN)
       nnzp,npatch,imonth1,dtlongn,timmax! INTENT(IN)
!
  use node_mod, only :  &
       mzp,             & ! INTENT(IN)
       mxp,             & ! INTENT(IN)
       myp,             & ! INTENT(IN)
       ja,              & ! INTENT(IN)
       jz,              & ! INTENT(IN)
       ia,              & ! INTENT(IN)
       iz,              & ! INTENT(IN)
       mynum,i0,j0        ! INTENT(IN)


  use io_params, only : frqanl !INTENT(IN)
  
  use micphys,   only:  &
       mcphys_type        ! INTENT(IN)

  use mem_radiate, ONLY: ilwrtyp, iswrtyp ! INTENT(IN)

  use mem_leaf, only: leaf_g
   IMPLICIT NONE

   INTEGER,PARAMETER :: &
         IDS=1, IDE=2, JDS=1, JDE=2, KDS=1, &
         IMS=1, IME=2, JMS=1, JME=2, KMS=1, &
         ITS=1, ITE=1, JTS=1, JTE=1, KTS=1  !- rams 1st level is below surface => kts=2

   INTEGER :: KDE, &
              KME, &
              KTE

   INTEGER :: i,j

   LOGICAL :: diagflag=.false. 
   REAL :: ocean_fraction

  

   !- converting WRF setting to BRAMS
   !ids=1   ;ide=mxp ;jds=1   ;jde=myp ;kds=1; kde=mzp                 
   !ims=1   ;ime=mxp ;jms=1   ;jme=myp ;kms=1; kme=mzp                           
   !its=ia  ;ite=iz  ;jts=ja  ;jte=jz  ;kts=2; kte=mzp-1  
   !- converting WRF setting to BRAMS
   kde=mzp             
   kme=mzp                       
   kte=mzp-1  
   !
   !- flag for diagnostic time
   diagflag = .false.
   if(mod(time,frqanl)<dtlongn(1).or.time>=timmax - 0.01*dtlongn(1) ) then 
      diagflag = .true. 
   endif

   do j = ja,jz
      do i = ia,iz
       ocean_fraction = leaf_g(ngrid)%patch_area(i,j,1)
              
       call brams_to_mic_wsm(ia,ja,iz,jz,mzp,mxp,myp &
            ,mcphys_type &
            ,ilwrtyp     &
            ,iswrtyp     &
            ,j           &
            ,i           &
            ,IDS, IDE, JDS, JDE, KDS, KDE   &
            ,IMS, IME, JMS, JME, KMS, KME   &
            ,ITS, ITE, JTS, JTE, KTS, KTE   &
            ,ngrid    &
            ,mynum    &
            ,if_adap  &
            !
            ,diagflag &
            !
            ,dtlt     &
            ,time     &
            ,zm       &
            ,dzt      &            
            ,zt       &
            ,basic_g(ngrid) &
            ,grid_g (ngrid) &
            ,micro_g(ngrid) &
	         ,ocean_fraction &
            )
      enddo
  enddo
  !- for consistency with surface and radiation schemes, the total
  !- precip will be also stored in the pcpg array
  micro_g(ngrid)%pcpg(:,:)=micro_g(ngrid)%pcprr(:,:)
   
 
END SUBROUTINE micro_wsm
!=======================================================================================
!
  SUBROUTINE brams_to_mic_wsm(ia,ja,iz,jz,m1,m2,m3 &
            ,mcphys_type &
            ,ilwrtyp     &
            ,iswrtyp     &
            ,j &
            ,i &
            ,IDS, IDE, JDS, JDE, KDS, KDE   &
            ,IMS, IME, JMS, JME, KMS, KME   &
            ,ITS, ITE, JTS, JTE, KTS, KTE   &
            ,ngrid &
            ,mynum &
            ,if_adap  &
            !
            ,diagflag &
            !
            ,dtlt  &
            ,time  &
            ,zm    &
            ,dzt   &
            ,zt    &
            !
            ,basic &
            !
            ,grd &
            !
            ,mic &
            !
            ,ocean_fraction &
            )
         
   USE module_mp_wsm3, only: wsm3init, wsm3
   USE module_mp_wsm5, only: wsm5init, wsm5
   USE module_mp_wsm6, only: wsm6init, wsm6
   USE module_mp_wsm7, only: wsm7init, wsm7
   USE rconstants, only: p00,cpor,alvl,alvi,cpi,cpi4,cp253i
   use mem_basic, only : &
       basic_vars            ! INTENT(INOUT)
   use mem_micro, only:  &
       micro_vars            ! INTENT(INOUT)

   use mem_grid, only:   &
       grid_vars         &   ! INTENT(IN)
      ,grid_g             ! INTENT(IN)           
    
   
   IMPLICIT NONE

   type(basic_vars) ::basic
   type(grid_vars ) ::grd
   type(micro_vars) ::mic

   INTEGER, INTENT(IN) ::  &          
             mcphys_type   &
            ,ilwrtyp       &
            ,iswrtyp       &
            ,j &
            ,i &
            ,IDS, IDE, JDS, JDE, KDS, KDE   &
            ,IMS, IME, JMS, JME, KMS, KME   &
            ,ITS, ITE, JTS, JTE, KTS, KTE   &
            ,m1, m2, m3                     &
            ,ngrid &
            ,mynum &
            ,if_adap ,ia,ja,iz,jz

   REAL, INTENT(IN) ::   &
             dtlt  &
            ,time  &
	         ,ocean_fraction 

   REAL,  INTENT(IN)   ,DIMENSION(m1) :: &
             zm    &
            ,dzt   &
            ,zt    
              
   LOGICAL, INTENT(IN) ::  diagflag

!- in the context of BRAMS, the variables below are "local"
   REAL, DIMENSION( ims:ime , kms:kme , jms:jme )  ::   &
                  th       &
                 ,dz8w     &
                 ,pi_phy   &
                 ,p        &
                 ,air_dens 

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ) ::                &
                  w                                               &
                 ,qv_curr,qc_curr,qr_curr,qi_curr,qs_curr,qg_curr &
                 ,qh_curr &
                 ,re_cloud, re_ice, re_snow, orho
                 !,qnc_curr, qnr_curr,  qni_curr          &
                 !,qnwfa_curr,qnifa_curr                  &
                     
   REAL, DIMENSION( ims:ime , jms:jme )  ::  &
                   RAINNC     &
                  ,RAINNCV    &
                  ,SNOWNC     &
                  ,SNOWNCV    &
                  ,GRAUPELNC  &
                  ,GRAUPELNCV &
                  ,HAIL       & 
                  ,HAILNCV    &
                  ,SR
                                                                   
!----------------------------------------------------------------------
! qv              water vapor    mixing ratio (kg/kg)
! qc              cloud water    mixing ratio (kg/kg)
! qr              rain water     mixing ratio (kg/kg)
! qi              cloud ice      mixing ratio (kg/kg)
! qs              snow            mixing ratio (kg/kg)
! qg              graupel            mixing ratio (kg/kg)
!
! qnc             cloud water number concentration (#/kg)
! qni             cloud ice   number concentration (#/kg)
! qnr             rain        number concentration (#/kg)
! qnwfa      water friendly aerosol number concentration (#/kg) - CCN
! qnifa      ice   friendly aerosol number concentration (#/kg) - IFN
!
!-- th            potential temperature    (K)
!-- w             vertical velocity (cartesian) (m/s)
!-- rho           density of air           (kg/m^3)
!-- pi_phy        exner function           (dimensionless)
!-- p             pressure                 (Pa)
!-- RAINNC        grid scale precipitation (mm)
!-- RAINNCV       one time step grid scale precipitation (mm/step)
!-- SNOWNC        grid scale snow and ice (mm)
!-- SNOWNCV       one time step grid scale snow and ice (mm/step)
!-- GRAUPELNC     grid scale graupel (mm)
!-- GRAUPELNCV    one time step grid scale graupel (mm/step)
!-- HAILNC        grid scale hail (mm)
!-- HAILNCV       one time step grid scale hail (mm/step)
!-- SR            one time step mass ratio of snow to total precip
!-- z             Height above sea level   (m)
!-- dt            Time step              (s)
!-- G             acceleration due to gravity  (m/s^2)
!-- CP            heat capacity at constant pressure for dry air (J/kg/K)
!-- R_d           gas constant for dry air (J/kg/K)
!-- R_v           gas constant for water vapor (J/kg/K)
!-- XLS           latent heat of sublimation   (J/kg)
!-- XLV           latent heat of vaporization  (J/kg)
!-- XLF           latent heat of melting       (J/kg)
!-- id            grid id number
!-- ids           start index for i in domain
!-- ide           end index for i in domain
!-- jds           start index for j in domain
!-- jde           end index for j in domain
!-- kds           start index for k in domain
!-- kde           end index for k in domain
!-- ims           start index for i in memory
!-- ime           end index for i in memory
!-- jms           start index for j in memory
!-- jme           end index for j in memory
!-- kms           start index for k in memory
!-- kme           end index for k in memory
!-- i_start       start indices for i in tile
!-- i_end         end indices for i in tile
!-- j_start       start indices for j in tile
!-- j_end         end indices for j in tile
!-- its           start index for i in tile
!-- ite           end index for i in tile
!-- jts           start index for j in tile
!-- jte           end index for j in tile
!-- kts           start index for k in tile
!-- kte           end index for k in tile
!-- diagflag      Logical to tell us when to produce diagnostics for history or restart
!-- rainprod      total tendency of conversion of cloud water/ice and graupel to rain (kg kg-1 s-1)

     real, dimension(ims:ime,kms:kme,jms:jme) :: hgt,refl_10cm 
     real, dimension(ims:ime, kms:kme, jms:jme):: rainprod,evapprod
     
     real :: tempk,rliq,rice,til,qhydm,tairstr
    ! integer :: itimestep = 1    ! not used in mp_wsm6
     integer :: do_radar_ref = 0 ! flag to compute radar reflectivity  
     integer :: has_reqc, has_reqi, has_reqs ! flags to calculate effec radius
                                             ! for radiation (1=on, 0=off
     real :: dx, dy              ! grid spacing (m)
     real :: dt                  ! model timestep (s)
     logical :: start_of_simulation =.true. 
     integer, save ::it=0
     integer  :: ke_diag              ! check this latter
     logical  :: wetscav_on = .false. ! check this latter

     REAL,  DIMENSION(m1) :: &
      thp    &
     ,theta  &
     ,pp     &
     ,rtp    &
     ,rv     &
     ,wp     &
     ,dn0    &
     ,pi0

     REAL :: rtgt

     REAL,DIMENSION(m1) :: &
      rcp     &
     ,rrp     &
     ,rpp     &
     ,rsp     &
     ,rgp     &
     ,rhp     &
     ,rei     &
     ,rel     &
     ,crp     &
     ,cpp     &
     ,ccp     &
     ,cccnp   &
     ,cifnp   


     REAL ::  &
      accpr   &! kg/m2 - rain+ice+snow+graupel
     ,pcprr   &! kg/m2 - rain+ice+snow+graupel
     ,accps   &! kg/m2 - ice+snow
     ,pcprs   &! kg/m2 - ice+snow
     ,accpg   &! kg/m2 - graupel
     ,pcprg   & ! kg/m2 - graupel
     ,accph   &
     ,pcprh 

     REAL, PARAMETER :: nt_c_ocean=100.E6 &
                       ,nt_c_land =200.E6 

     REAL                :: nt_c_var
     REAL    , PARAMETER :: epsilon         = 1.E-15
     REAL    , PARAMETER :: g = 9.81  ! acceleration due to gravity (m {s}^-2)
  
     REAL    , PARAMETER :: r_d          = 287.       ! gas constant of dry air (J deg^-1 kg^-1)
     REAL    , PARAMETER :: cp           = 7.*r_d/2.  ! 
     REAL    , PARAMETER :: rhowater     = 1000.      ! density of liquid water at 0^oC (kg m^-3)
     REAL    , PARAMETER :: rhosnow      = 100.       ! density of snow (kg m^-3)
     REAL    , PARAMETER :: rhoair0      = 1.28       ! density of dry air at 0^oC and 1000mb pressure (kg m^-3)
     REAL    , PARAMETER :: r_v          = 461.6      ! gas constant for water vapor (J deg^-1 kg^-1)
     REAL    , PARAMETER :: cv           = cp-r_d     ! Specific heat of air at contant volume (J deg^-1 kg^-1)
     REAL    , PARAMETER :: cpv          = 4.*r_v
     REAL    , PARAMETER :: cvv          = cpv-r_v    ! 
     REAL    , PARAMETER :: cvpm         = -cv/cp
     REAL    , PARAMETER :: cliq         = 4190.      ! specific heat of liquid water at 0^oC
     REAL    , PARAMETER :: cice         = 2106.      ! specific heat of ice at 0^oC
     REAL    , PARAMETER :: psat         = 610.78
     REAL    , PARAMETER :: XLV0         = 3.15E6     !  constant defined for calculation of latent heating 
     REAL    , PARAMETER :: XLV1         = 2370.      !  constant defined for calculation of latent heating 
     REAL    , PARAMETER :: XLS0         = 2.905E6    !  constant defined for calculation of latent heating
     REAL    , PARAMETER :: XLS1         = 259.532    !  constant defined for calculation of latent heating
     REAL    , PARAMETER :: SVP1         = 0.6112     !  constant for saturation vapor pressure calculation (dimensionless)
     REAL    , PARAMETER :: SVP2         = 17.67      ! constant for saturation vapor pressure calculation (dimensionless)
     REAL    , PARAMETER :: SVP3         = 29.65      ! constant for saturation vapor pressure calculation (K)
     REAL    , PARAMETER :: SVPT0        = 273.15     ! constant for saturation vapor pressure calculation (K)

     REAL    , PARAMETER :: XLS          = 2.85E6     ! latent heat of sublimation of water at 0^oC (J kg^-1) 
     REAL    , PARAMETER :: XLV          = 2.5E6      ! latent heat of vaporization of water at 0^oC (J kg^-1)
     REAL    , PARAMETER :: XLF          = 3.50E5     ! latent heat of fusion of water at 0^oC (J kg^-1)
     REAL    , PARAMETER :: EP_1         = R_v/R_d-1. !  constant for virtual temperature (r_v/r_d - 1) (dimensionless)
     REAL    , PARAMETER :: EP_2         = R_d/R_v    ! constant for specific humidity calculation (dimensionless)

     INTEGER :: k,kr,ii,jj
     
     integer, parameter :: output = 1
     integer :: l_unit = 4
     character(len=5) :: ctime,cmynum
     character(len=2) :: cmcphys_type
     integer :: rec_size, irec, nz,nlev,nz1,nz2,int_byte_size
     real, dimension(m2,m3) :: lons,lats, glat, glon
     real   :: real_byte_size
        !===============================================
        if(OUTPUT==1 .and. i==ia .and. j==ja .and. mod(time,3600.)< dtlt) then 
            do jj=1,m3
             do ii=1,m2
                lons (ii,jj) = grid_g(ngrid)%glon(ii,jj)*3.14159/180. !- convert to rad
                lats (ii,jj) = grid_g(ngrid)%glat(ii,jj)*3.14159/180. !- convert to rad
            enddo; enddo
            cmynum = '00000'
            cmcphys_type = '00'
            write(ctime ,fmt="(I5.5)") int(time)
            write(cmynum,fmt="(I5.5)") mynum
            write(cmcphys_type,fmt="(I2.2)") mcphys_type
            
            write(45,*) "lons: ", lons
            write(45,*) "lats: ", lats

            open(newunit = l_unit,file = "WSM"//trim(cmcphys_type)//"_dataIn-"//ctime//"-"//trim(cmynum)//".bin",ACCESS = "stream", action="write", status="replace")
            write(l_unit) m1,m2,m3, mynum
            write(l_unit) ia,iz,ja,jz
            write(l_unit) ids,ide, jds,jde, kds,kde 
            write(l_unit) ims,ime, jms,jme, kms,kme 
            write(l_unit) its,ite, jts,jte, kts,kte 
            write(l_unit) mcphys_type , ilwrtyp, iswrtyp
            write(l_unit) time, dtlt
            write(l_unit)  g      
            write(l_unit) cp      
            write(l_unit) cpv     
            write(l_unit) r_d     
            write(l_unit) r_v     
            write(l_unit) svpt0   
            write(l_unit) ep_1    
            write(l_unit) ep_2    
            write(l_unit) epsilon 
            write(l_unit) xls     
            write(l_unit) xlv     
            write(l_unit) xlf     
            write(l_unit) rhoair0 
            write(l_unit) rhowater
            write(l_unit) rhosnow   
            write(l_unit) cliq    
            write(l_unit) cice    
            write(l_unit) psat
            write(l_unit) dzt     
            write(l_unit) zt
            write(l_unit) grid_g(ngrid)%glon
            write(l_unit) grid_g(ngrid)%glat
            !
            write(l_unit) basic%thp  
            write(l_unit) basic%theta
            write(l_unit) basic%pp   
            write(l_unit) basic%rtp  
            write(l_unit) basic%rv   
            write(l_unit) basic%wp   
            write(l_unit) basic%dn0  
            write(l_unit) basic%pi0  
            write(l_unit) mic%rcp    
            write(l_unit) mic%rrp    
             
            write(l_unit) grd%rtgt 
            write(l_unit) mic%accpr
            write(l_unit) mic%pcprr
                  
            if(mcphys_type == 7 )  then 
               write(l_unit)  mic%rhp  
               write(l_unit)  mic%accph
               write(l_unit)  mic%pcprh
            endif
            
            if(mcphys_type == 6 .or. mcphys_type == 7) then
               write(l_unit)  mic%rgp  
               write(l_unit)  mic%accpg
               write(l_unit)  mic%pcprg
            endif

            if(mcphys_type == 5 .or. mcphys_type == 6 .or. mcphys_type == 7 ) then
               write(l_unit)  mic%rpp  
               write(l_unit)  mic%rsp  
               write(l_unit)  mic%accps
               write(l_unit)  mic%pcprs
            endif

            close(l_unit)
        endif
        !===============================================

        nt_c_var =    ocean_fraction *nt_c_ocean + &
	               (1.-ocean_fraction)*nt_c_land

        !- column quantities
	     thp  (1:m1)= basic%thp  (1:m1,i,j)
        theta(1:m1)= basic%theta(1:m1,i,j)
        pp   (1:m1)= basic%pp   (1:m1,i,j)
        rtp  (1:m1)= basic%rtp  (1:m1,i,j)
        rv   (1:m1)= basic%rv   (1:m1,i,j)
        wp   (1:m1)= basic%wp   (1:m1,i,j)
        dn0  (1:m1)= basic%dn0  (1:m1,i,j)
        pi0  (1:m1)= basic%pi0  (1:m1,i,j)
        !--- mass mixing ratio
        rcp  (1:m1)= mic%rcp    (1:m1,i,j)
        rrp  (1:m1)= mic%rrp    (1:m1,i,j)
        
	     !- surface quantities
	     rtgt = grd%rtgt (i,j)
	     accpr= mic%accpr(i,j)
        pcprr= mic%pcprr(i,j)

        if(mcphys_type == 7 )  then 
          rhp  (1:m1)= mic%rhp  (1:m1,i,j)
          accph      = mic%accph(i,j)
          pcprh      = mic%pcprh(i,j)
        else
          rhp  (1:m1)= 0.
          accph      = 0.
          pcprh      = 0.
        endif

        if(mcphys_type == 6 .or. mcphys_type == 7) then
          rgp  (1:m1)= mic%rgp  (1:m1,i,j)
          accpg      = mic%accpg(i,j)
          pcprg      = mic%pcprg(i,j)
        else
          rgp  (1:m1)= 0.
          accpg      = 0.
          pcprg      = 0.
        endif

        if(mcphys_type == 5 .or. mcphys_type == 6 .or. mcphys_type == 7 ) then
          rpp  (1:m1)= mic%rpp  (1:m1,i,j)
          rsp  (1:m1)= mic%rsp  (1:m1,i,j)
          accps      = mic%accps(i,j)
          pcprs      = mic%pcprs(i,j)      
        else
          rpp  (1:m1)= 0.
          rsp  (1:m1)= 0.
          accps      = 0.
          pcprs      = 0.  
        endif
!
!
!- for coupling with brams
!       !- converting WRF setting to BRAMS
!       ids=1   ;ide=mxp ;jds=1   ;jde=myp ;kds=1; kde=mzp              
!       ims=1   ;ime=mxp ;jms=1   ;jme=myp ;kms=1; kme=mzp                        
!       its=ia  ;ite=iz  ;jts=ja  ;jte=jz  ;kts=1; kte=mzp-1  

        ! flags to calculate effec radius
        IF( (ilwrtyp==6 .or. iswrtyp==6)) then           
           has_reqc= 1 ; has_reqi= 1 ; has_reqs= 1 
        ELSE
           has_reqc= 0 ; has_reqi= 0 ; has_reqs= 0 
        ENDIF
        
        dt= dtlt        ! time step            (s)
        rainprod  =0.0  ! for scaveging   aerosols/gases
        evapprod  =0.0  ! for evaporation aerosols/gases
        SR        =0.0  ! fraction of snow of the total water
                        ! ( for land surface models)
        refl_10cm =0.0  ! 
        ke_diag   = kte
       
        !- surface precipitation (total accumulated)
        RAINNC    (1,1)=  accpr !- rain+ice+snow+graupel+hail
        SNOWNC    (1,1)=  accps !- ice+snow
        GRAUPELNC (1,1)=  accpg !- graupel
        HAIL      (1,1)=  accph !- hail
        
        DO k=1,kme-1
          kr = k + 1
          qv_curr (1,k,1)= max(1.e-12,rtp(kr) - &    ! QV
                               (rcp(kr)+rrp(kr)+rpp(kr)+rsp(kr)+rgp(kr)+rhp(kr)))                    
          qc_curr (1,k,1)= max(0.0,rcp(kr))          ! QC     
          qr_curr (1,k,1)= max(0.0,rrp(kr))          ! QR   
          qi_curr (1,k,1)= max(0.0,rpp(kr))          ! QI   
          qs_curr (1,k,1)= max(0.0,rsp(kr))          ! QS   
          qg_curr (1,k,1)= max(0.0,rgp(kr))          ! QG

          qh_curr (1,k,1)= max(0.0,rhp(kr))          ! QH   
         
          pi_phy  (1,k,1)= (pp(kr)+pi0(kr))*cpi ! Exner function/cp (dimensionless)
          
          P       (1,k,1)= ( (pp(kr)+pi0(kr))*cpi )** cpor * p00        ! pressure(Pa)
          W       (1,k,1)= wp(kr)       ! vertical velocity (m/s) ! must be at center or face? ASK
          
          dz8w    (1,k,1)= rtgt/dzt(kr) ! layer thickness (m) 

        ENDDO

        !- get potential temperature (theta) from theta_il (thp) and condensates
        DO k=1,kme -1 
           kr = k + 1
           tempK    = theta(kr)* (pp(kr)+pi0(kr))*cpi 
              til   = thp  (kr)* (pp(kr)+pi0(kr))*cpi 
         
           rliq     =  qc_curr(1,k,1) + qr_curr(1,k,1)                
           rice     =  qi_curr(1,k,1) + qs_curr(1,k,1) + qg_curr(1,k,1) + qh_curr(1,k,1)
           qhydm    =  alvl * rliq + alvi * rice
           
           if (tempK .gt. 253.) then
              tairstr = 0.5 * (til + sqrt(til * (til + cpi4 * qhydm)))
           else
              tairstr = til * (1. + qhydm * cp253i)
           endif
           !- updated potential temperature TH in Kelvin (adv+dif+rad+conv+)
           TH (1,k,1) = tairstr / pi_phy(1,k,1)

           !- air density
           air_dens(1,k,1) = P(1,k,1)/(287.04*tempK*(1.+0.608*qv_curr(1,k,1)))
          !air_dens(1,k,1)= dn0(kr) 

        ENDDO
        
        !
        IF(start_of_simulation) THEN !.or.restart.)   
           IF(mcphys_type == 30 ) &   
              CALL wsm3init(rhoair0,rhowater,rhosnow,cliq,cpv,   .false. )
           IF(mcphys_type ==  5 ) &   
              CALL wsm5init(rhoair0,rhowater,rhosnow,cliq,cpv,   .false. )
           IF(mcphys_type ==  6 ) &   
              CALL wsm6init(rhoair0,rhowater,rhosnow,cliq,cpv, 0,.false. )
           IF(mcphys_type ==  7 ) &   
              CALL wsm7init(rhoair0,rhowater,rhosnow,cliq,cpv,   .false. )

           start_of_simulation =.false.
        ENDIF

        IF(mcphys_type == 30)                   &   
            
             CALL wsm3(                         &
                     TH,                        &! potential temperature    (K)
                     qv_curr,                   &! QV=qv_curr,     
                     qc_curr,                   &! QC=qc_curr,     
                     qr_curr,                   &! QR=qr_curr,     
                     !
                     w,                         &! check if used
                     air_dens,                  &          
                     pi_phy,                    &! exner function (dimensionless)
                     P,                         &! pressure(Pa)
                     dz8w,                      &! deltaz
                     !
                     dt,      &                  ! time step              (s)
                      g,      &
                     cp,      &
                     cpv,     &
                     r_d,     &
                     r_v,     &
                     svpt0,   &
                     ep_1,    &
                     ep_2,    &
                     epsilon, &
                     xls,     &
                     xlv,     &
                     xlf,     &
                     rhoair0, &
                     rhowater,&  
                     cliq,    &
                     cice,    &
                     psat,    & 
                     !
                     RAINNC,                    &
                     RAINNCV,                   &
                     SNOWNC,                    &
                     SNOWNCV,                   &
                     SR,                        &
                     !
                     ! 
                     has_reqc,                  & 
                     has_reqi,                  &  
                     has_reqs,                  & 
                     !
                     re_cloud,                  & 
                     re_ice,                    &
                     re_snow,                   &
                     IDS,IDE, JDS,JDE, KDS,KDE, &
                     IMS,IME, JMS,JME, KMS,KME, &
                     ITS,ITE, JTS,JTE, KTS,KTE  &
                     )
                   

        
        IF(mcphys_type == 5)                    &   
               
             CALL wsm5(                         &
                     TH,                        &! potential temperature    (K)
                     qv_curr,                   &! QV=qv_curr,     
                     qc_curr,                   &! QC=qc_curr,     
                     qr_curr,                   &! QR=qr_curr,     
                     qi_curr,                   &! QI=qi_curr,     
                     qs_curr,                   &! QS=qs_curr,     
                     !
                     air_dens,                  &          
                     pi_phy,                    &! exner function (dimensionless)
                     P,                         &! pressure(Pa)
                     dz8w,                      &! deltaz
                     !
                     dt,      &                  ! time step              (s)
                     g,       &
                     cp,      &
                     cpv,     &
                     r_d,     &
                     r_v,     &

                     svpt0,   &
                     ep_1,    &
                     ep_2,    &
                     epsilon, &
                     xls,     &
                     xlv,     &
                     xlf,     &
                     rhoair0, &
                     rhowater,&  
                     cliq,    &
                     cice,    &
                     psat,    & 
                    
                     RAINNC,                    &
                     RAINNCV,                   &
                     SNOWNC,                    &
                     SNOWNCV,                   &
                     SR,                        &
                     !
                     refl_10cm,                 &
                     diagflag,                  &
                     do_radar_ref,              &
                     ! 
                     has_reqc,                  & 
                     has_reqi,                  &  
                     has_reqs,                  & 
                     !
                     re_cloud,                  & 
                     re_ice,                    &
                     re_snow,                   &
                     IDS,IDE, JDS,JDE, KDS,KDE, &
                     IMS,IME, JMS,JME, KMS,KME, &
                     ITS,ITE, JTS,JTE, KTS,KTE  &
                     )
        
        IF(mcphys_type == 6)                    &   
             
             CALL wsm6(                         &
                     TH,                        &! potential temperature    (K)
                     qv_curr,                   &! QV=qv_curr,     
                     qc_curr,                   &! QC=qc_curr,     
                     qr_curr,                   &! QR=qr_curr,     
                     qi_curr,                   &! QI=qi_curr,     
                     qs_curr,                   &! QS=qs_curr,     
                     qg_curr,                   &! QG=qg_curr,     
                     
                     air_dens,                  &          
                     pi_phy,                    &! exner function (dimensionless)
                     P,                         &! pressure(Pa)
                     dz8w,                      &! deltaz

                     dt,      &                  ! time step              (s)
                     g,       &
                     cp,      &
                     cpv,     &
                     r_d,     &
                     r_v,     &
                     svpt0,   &
                     ep_1,    &
                     ep_2,    &
                     epsilon, &
                     xls,     &
                     xlv,     &
                     xlf,     &
                     rhoair0, &
                     rhowater,&  
                     cliq,    &
                     cice,    &
                     psat,    & 
                    
                     RAINNC,                    &
                     RAINNCV,                   &
                     SNOWNC,                    &
                     SNOWNCV,                   &
                     SR,                        &
                     !
                     refl_10cm,                 &
                     diagflag,                  &
                     do_radar_ref,              &
                     GRAUPELNC,                 & 
                     GRAUPELNCV,                &
                     !
                     !ke_diag,                  &
                     ! 
                     has_reqc,                  & 
                     has_reqi,                  &  
                     has_reqs,                  & 
                     !
                     re_cloud,                  & 
                     re_ice,                    &
                     re_snow,                   &
                     IDS,IDE, JDS,JDE, KDS,KDE, &
                     IMS,IME, JMS,JME, KMS,KME, &
                     ITS,ITE, JTS,JTE, KTS,KTE, &
                     !
                     wetscav_on,                &
                     evapprod,                  &
                     rainprod                   &
                     )

           IF(mcphys_type == 7)                 &   
             
             CALL wsm7(                         &
                     TH,                        &! potential temperature    (K)
                     qv_curr,                   &! QV=qv_curr,     
                     qc_curr,                   &! QC=qc_curr,     
                     qr_curr,                   &! QR=qr_curr,     
                     qi_curr,                   &! QI=qi_curr,     
                     qs_curr,                   &! QS=qs_curr,     
                     qg_curr,                   &! QG=qg_curr,     
                     qh_curr,                   &! Qh=qg_curr,     
                     
                     air_dens,                  &          
                     pi_phy,                    &! exner function (dimensionless)
                     P,                         &! pressure(Pa)
                     dz8w,                      &! deltaz

                     dt,      &                  ! time step              (s)
                     g,       &
                     cp,      &
                     cpv,     &
                     r_d,     &
                     r_v,     &
                     svpt0,   &
                     ep_1,    &
                     ep_2,    &
                     epsilon, &
                     xls,     &
                     xlv,     &
                     xlf,     &
                     rhoair0, &
                     rhowater,&  
                     cliq,    &
                     cice,    &
                     psat,    & 
                    
                     RAINNC,                    &
                     RAINNCV,                   &
                     SNOWNC,                    &
                     SNOWNCV,                   &
                     SR,                        &
                     !
                     refl_10cm,                 &
                     diagflag,                  &
                     do_radar_ref,              &
                     GRAUPELNC,                 & 
                     GRAUPELNCV,                &
                     HAIL,                      & 
                     HAILNCV,                   &
                     !
                     !ke_diag,                  &
                     ! 
                     has_reqc,                  & 
                     has_reqi,                  &  
                     has_reqs,                  & 
                     !
                     re_cloud,                  & 
                     re_ice,                    &
                     re_snow,                   &
                     IDS,IDE, JDS,JDE, KDS,KDE, &
                     IMS,IME, JMS,JME, KMS,KME, &
                     ITS,ITE, JTS,JTE, KTS,KTE  &
                     )


                                                                    
        !- updated variables after microphysics processes (from WSM to BRAMS)
        DO k=1,kme-1
         kr=k+1
         rtp(kr)= qv_curr(1,k,1) + &
                  qc_curr(1,k,1) + &     
                  qr_curr(1,k,1) + &    
                  qi_curr(1,k,1) + &    
                  qs_curr(1,k,1) + &    
                  qg_curr(1,k,1) + &    
                  qh_curr(1,k,1)  

         rcp(kr)= qc_curr(1,k,1)
         rrp(kr)= qr_curr(1,k,1)
         rpp(kr)= qi_curr(1,k,1)
         rsp(kr)= qs_curr(1,k,1)
         rgp(kr)= qg_curr(1,k,1)
         rhp(kr)= qh_curr(1,k,1)  
         
         rv (kr)= max(1.0e-12, rtp(kr) -(rcp(kr)+rrp(kr)+rpp(kr)+rsp(kr)+rgp(kr)+rhp(kr)))
      
         theta(kr) =  TH(1,k,1)
         tempK     =  TH(1,k,1)*pi_phy(1,k,1)
         
         rliq     =  qc_curr(1,k,1) + qr_curr(1,k,1)                 
         rice     =  qi_curr(1,k,1) + qs_curr(1,k,1) + qg_curr(1,k,1) + qh_curr(1,k,1)
         
         !- update liq-ice potential temperature THP in Kelvin including microphysics processes
         thp(kr)  =   TH(1,k,1)*(1. + alvl * rliq/(cp * max(tempK,253.))  &
                                    + alvi * rice/(cp * max(tempK,253.)) ) **(-1.0)      
        ENDDO
        !- definition for k=1
          rtp(1)  = rtp(2)  
          rcp(1)  = rcp(2)  
          rrp(1)  = rrp(2)  
          rpp(1)  = rpp(2)  
          rsp(1)  = rsp(2)  
          rgp(1)  = rgp(2)  
          rhp(1)  = rhp(2)
          rv (1)  = rv (2)  
          theta(1)= theta(2)
          thp(1)  = thp(2)  
        !
        IF( (ilwrtyp==6 .or. iswrtyp==6)) then           
          DO k=1,kme-1
            kr=k+1
            rel (kr) = re_cloud (1,k,1) * 1.e+6 ! RRTM requires in micrometer
            rei (kr) = re_ice   (1,k,1) * 1.e+6 ! RRTM requires in micrometer
          ENDDO
          rel (1) =rel (2) !;  rel (kme) =rel (kme-1) 
          rei (1) =rei (2) !;  rei (kme) =rei (kme-1)
        ENDIF        

        !- surface precipitation (units are kg/m^2 = mm)
        !- RAINNC and RAINNCV constains all precipitation hidrometeors (rain, graupel, snow, ...)
        accpr = RAINNC    (1,1) ! a = accum
        pcprr = RAINNCV   (1,1) ! p = for each dt  (or per time step)
        accps = SNOWNC    (1,1) 
        pcprs = SNOWNCV   (1,1) 
        accpg = GRAUPELNC (1,1) 
        pcprg = GRAUPELNCV(1,1) 
        accph = HAIL      (1,1) 
        pcprh = HAILNCV   (1,1) 

        !- column quantities
        basic%thp  (1:m1,i,j) =thp  (1:m1) 
        basic%theta(1:m1,i,j) =theta(1:m1)
        basic%rtp  (1:m1,i,j) =rtp  (1:m1)
        basic%rv   (1:m1,i,j) =rv   (1:m1)  
              
        mic%rcp    (1:m1,i,j) =rcp  (1:m1)   
        mic%rrp    (1:m1,i,j) =rrp  (1:m1)

        if(mcphys_type == 5 .or. mcphys_type == 6 .or. mcphys_type == 7 ) then   
          mic%rpp    (1:m1,i,j) =rpp  (1:m1)   
          mic%rsp    (1:m1,i,j) =rsp  (1:m1)   
        endif
        
        if(mcphys_type == 6 .or. mcphys_type == 7) mic%rgp(1:m1,i,j) =rgp  (1:m1)  
        if(mcphys_type == 7                      ) mic%rhp(1:m1,i,j) =rhp  (1:m1)
        
        if( ilwrtyp==6 .or. iswrtyp==6 ) then
          mic%rei   (1:m1,i,j) =rei  (1:m1)
          mic%rel   (1:m1,i,j) =rel  (1:m1)
        endif 

	!- surface quantities
        mic%accpr(i,j) = accpr !constains all precipitation hidrometeors (rain, graupel, snow, ...)
        mic%pcprr(i,j) = pcprr !constains all precipitation hidrometeors (rain, graupel, snow, ...)
        if(mcphys_type == 5 .or. mcphys_type == 6 .or. mcphys_type == 7 ) then
          mic%accps(i,j) = accps
          mic%pcprs(i,j) = pcprs
        endif
        if(mcphys_type == 6 .or. mcphys_type == 7) then
          mic%accpg(i,j) = accpg
          mic%pcprg(i,j) = pcprg
        endif
        if(mcphys_type == 7) then
          mic%accph(i,j) = accph
          mic%pcprh(i,j) = pcprh
        endif

        !===============================================
        if(OUTPUT==1 .and. i==iz .and. j==jz .and. mod(time,3600.)< dtlt) then 
if(maxval(3600*mic%pcprr ) > 10. .or. maxval(basic%wp) > 5.) then
print*,"================================"
print*,"mynum",time,mynum,maxval(3600*mic%pcprr ),maxval(basic%wp) ; call flush(6)
print*,"================================"; call flush(6)
         !glon = grid_g(ngrid)%glon
         !glat = grid_g(ngrid)%glat
         inquire (iolength=int_byte_size) real_byte_size  ! inquire by output list
         rec_size = m2*m3*int_byte_size
         nz1 = 2
         nz2 = m1-1
         open(newunit = l_unit, file="WSM"//trim(cmcphys_type)//"_dataOut-ref-"//ctime//"-"//trim(cmynum)//".gra", &
                  form='unformatted', access='direct', status='replace', recl=rec_size)
             irec=1
             do nz=nz1,nz2;  write(l_unit,rec=irec) basic%wp         (nz,1:m2,1:m3); irec=irec+1 ;enddo
             do nz=nz1,nz2;  write(l_unit,rec=irec) basic%thp        (nz,1:m2,1:m3); irec=irec+1 ;enddo
             do nz=nz1,nz2;  write(l_unit,rec=irec) basic%theta      (nz,1:m2,1:m3); irec=irec+1 ;enddo
             do nz=nz1,nz2;  write(l_unit,rec=irec) 1000.*basic%rtp  (nz,1:m2,1:m3); irec=irec+1 ;enddo
             do nz=nz1,nz2;  write(l_unit,rec=irec) 1000.*basic%rv   (nz,1:m2,1:m3); irec=irec+1 ;enddo
             do nz=nz1,nz2;  write(l_unit,rec=irec) 1000.*mic%rrp    (nz,1:m2,1:m3); irec=irec+1 ;enddo
             do nz=nz1,nz2;  write(l_unit,rec=irec) 1000.*mic%rcp    (nz,1:m2,1:m3); irec=irec+1 ;enddo

           !-- 2d 
             write(l_unit,rec=irec) 3600*mic%accpr(1:m2,1:m3); irec=irec+1 
             write(l_unit,rec=irec) 3600*mic%pcprr(1:m2,1:m3); irec=irec+1 

         close(l_unit)
         !== number of levels 2: m1 -1 
         nlev = nz2-nz1+1

         open(newunit = l_unit, file="WSM"//trim(cmcphys_type)//"_dataOut-ref-"//ctime//"-"//trim(cmynum)//".ctl", &
                action='write', status='replace')

            write(l_unit,*) 'dset ^'//"WSM"//trim(cmcphys_type)//"_dataOut-ref-"//ctime//"-"//trim(cmynum)//".gra"
            !writing others infos to ctl
            write(l_unit,*) 'undef -0.9990000E+34'
            write(l_unit,*) 'title WSM '
            write(l_unit,*) 'xdef ',m2,' linear ',grid_g(ngrid)%glon(1,1),grid_g(ngrid)%glon(2,1)-grid_g(ngrid)%glon(1,1)
            write(l_unit,*) 'ydef ',m3,' linear ',grid_g(ngrid)%glat(1,1),grid_g(ngrid)%glat(1,2)-grid_g(ngrid)%glat(1,1)
            write(l_unit,2005) nlev, zt(nz1:nz2)
            write(l_unit,*) 'tdef 1 linear 00:00Z01JAN200 1mo'
            write(l_unit,*) 'vars ',9
            write(l_unit,*) 'wp'    ,nlev,'99 ','K'
            write(l_unit,*) 'thp'   ,nlev,'99 ','K'
            write(l_unit,*) 'theta' ,nlev,'99 ','K'
            write(l_unit,*) 'rtp'   ,nlev,'99 ','g/kg'
            write(l_unit,*) 'rv'    ,nlev,'99 ','g/kg'
            write(l_unit,*) 'rain'  ,nlev,'99 ','g/kg'
            write(l_unit,*) 'cloud' ,nlev,'99 ','g/kg'
            write(l_unit,*) 'apr',' 0 ','99 ','K'
            write(l_unit,*) 'ppr',' 0 ','99 ','K'
            write(l_unit,*) 'endvars'
         close(l_unit)
         2005 format('zdef ',i4,' levels ',100f8.1)
        endif
endif
!
END  SUBROUTINE brams_to_mic_wsm
   