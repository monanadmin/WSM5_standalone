program wsm_test

   use module_mp_wsm3, only: wsm3init, wsm3
   use module_mp_wsm5, only: wsm5init, wsm5
   use module_mp_wsm6, only: wsm6init, wsm6
   use module_mp_wsm7, only: wsm7init, wsm7
   use rconstants, only: p00,cpor,alvl,alvi,cpi,cpi4,cp253i

   implicit none     
     
   integer :: m1,m2,m3, mynum,ia,iz,ja,jz
   integer :: ids,ide, jds,jde, kds,kde 
   integer :: ims,ime, jms,jme, kms,kme 
   integer :: its,ite, jts,jte, kts,kte 
   integer :: mcphys_type , ilwrtyp, iswrtyp
   real ::  time, dtlt
   real ::   g      
   real ::  cp      
   real ::  cpv     
   real ::  r_d     
   real ::  r_v     
   real ::  svpt0   
   real ::  ep_1    
   real ::  ep_2    
   real ::  epsilon 
   real ::  xls     
   real ::  xlv     
   real ::  xlf     
   real ::  rhoair0 
   real ::  rhowater, rhosnow  
   real ::  cliq    
   real ::  cice    
   real ::  psat     
   
   integer ,parameter :: nzpmax=100
   integer :: l_unit = 4, i,j,k,kr
   character(len=5) :: ctime,cmynum
   character(len=2) :: cmcphys_type
   character(len=180) :: filein, fileout
   integer :: rec_size, irec, nz,nlev,nz1,nz2,int_byte_size 
   logical :: first_read = .true.

   real :: tempk,rliq,rice,til,qhydm,tairstr,real_byte_size
   integer :: do_radar_ref = 0 
   integer :: has_reqc, has_reqi, has_reqs 
   real :: dx, dy              ! grid spacing (m)
   real :: dt                  ! model timestep (s)
   integer  :: ke_diag            
   logical  :: wetscav_on = .false.     
   logical :: start_of_simulation =.true. 
   logical :: diagflag=.false. 
   character(len=180) :: prefix
   character(len=180) :: var_datain
   character(len=180) :: var_dataout
!- this for the namelist wsm.inp
   namelist /run/ prefix, var_datain, var_dataout, time, mcphys_type, mynum

   real   , allocatable, dimension(:)  :: &
              dzt, zt 

   real   , allocatable, dimension(:,:)  :: &
              rtgt   &
             ,accpr  &
             ,pcprr  &    
             ,accph  &
             ,pcprh  &
             ,accpg  &
             ,pcprg  &
             ,accps  &
             ,pcprs  &
             ,glon   &
             ,glat

   real   , allocatable, dimension(:,:,:)  :: &
              thp    &
             ,theta  &
             ,pp     & 
             ,rtp    &
             ,rv     &
             ,wp     &
             ,dn0    &
             ,pi0    &
             ,rcp    &
             ,rrp    &
             ,rhp    & 
             ,rgp    &
             ,rpp    &
             ,rsp    &
             ,rei    &
             ,rel
!--- local vars
   real                :: &
              rtgt_1d     &
             ,accpr_1d    &
             ,pcprr_1d    &    
             ,accph_1d    &
             ,pcprh_1d    &
             ,accpg_1d    &
             ,pcprg_1d    &
             ,accps_1d    &
             ,pcprs_1d  

   real, dimension(nzpmax)  :: &
              theta_1d&
             ,pp_1d   & 
             ,wp_1d   &
             ,dn0_1d  &
             ,pi0_1d  &
             ,thp_1d  &
             ,rtp_1d  &
             ,rv_1d   &
             ,rcp_1d  &
             ,rrp_1d  &
             ,rhp_1d  & 
             ,rgp_1d  &
             ,rpp_1d  &
             ,rsp_1d  &
             ,rel_1d  &
             ,rei_1d             


   real, allocatable, dimension(:,:,:)  ::   &
              th       &
             ,dz8w     &
             ,pi_phy   &
             ,p        &
             ,air_dens &
             ,w 		  &
             ,qv_curr,qc_curr,qr_curr,qi_curr,qs_curr,qg_curr &
             ,qh_curr,re_cloud, re_ice, re_snow, orho	      &
             ,hgt, refl_10cm , rainprod,evapprod
 
   real, allocatable, dimension(:,:)  ::   &
              rainnc	    &
             ,rainncv	 &
             ,snownc	    &
             ,snowncv	 &
             ,graupelnc  &
             ,graupelncv &
             ,hail	    & 
             ,hailncv	 &
             ,sr

!- 
   integer :: nvar_out
   
!- read namelist  
   open(15,file='wsm.inp',status='old',form='formatted')    
    read(15,nml=run)
   close(15)

   !time=18000
   !mynum=365
   !mcphys_type=30
   !------------------------ reading input data
   write(ctime ,fmt="(I5.5)") int(time)
   write(cmynum,fmt="(I5.5)") mynum
   write(cmcphys_type,fmt="(I2.2)") mcphys_type

   filein= "WSM"//trim(cmcphys_type)//"_dataIn-"//ctime//"-"//trim(cmynum)//".bin"
   
   print*,"opening file for reading: ",trim(var_datain)//trim(filein)

   open(newunit = l_unit,file =trim(var_datain)//trim(filein),ACCESS = "stream",action="read", status="old")
        
      read(l_unit) m1,m2,m3, mynum
      print*,'dim1=', m1,m2,m3, mynum
      read(l_unit) ia,iz,ja,jz
      
      read(l_unit) ids,ide, jds,jde, kds,kde 
      read(l_unit) ims,ime, jms,jme, kms,kme 
      read(l_unit) its,ite, jts,jte, kts,kte 

      if(first_read) then
         first_read = .false.
         allocate(dzt   (m1))       ;dzt    = 0.0
         allocate(zt    (m1))       ;zt     = 0.0
         allocate(rtgt  (m2,m3))    ;rtgt   = 0.0 
         allocate(accpr (m2,m3))    ;accpr  = 0.0
         allocate(pcprr (m2,m3))    ;pcprr  = 0.0    
         allocate(accph (m2,m3))    ;accph  = 0.0 
         allocate(pcprh (m2,m3))    ;pcprh  = 0.0 
         allocate(accpg (m2,m3))    ;accpg  = 0.0 
         allocate(pcprg (m2,m3))    ;pcprg  = 0.0 
         allocate(accps (m2,m3))    ;accps  = 0.0 
         allocate(pcprs (m2,m3))    ;pcprs  = 0.0
         allocate(glon  (m2,m3))    ;glon   = 0.0
         allocate(glat  (m2,m3))    ;glat   = 0.0
          
                      
         allocate(thp  (m1,m2,m3))   ;thp   = 0.0 
         allocate(theta(m1,m2,m3))   ;theta = 0.0
         allocate(pp   (m1,m2,m3))   ;pp    = 0.0 
         allocate(rtp  (m1,m2,m3))   ;rtp   = 0.0 
         allocate(rv   (m1,m2,m3))   ;rv    = 0.0 
         allocate(wp   (m1,m2,m3))   ;wp    = 0.0 
         allocate(dn0  (m1,m2,m3))   ;dn0   = 0.0 
         allocate(pi0  (m1,m2,m3))   ;pi0   = 0.0 
         allocate(rcp  (m1,m2,m3))   ;rcp   = 0.0 
         allocate(rrp  (m1,m2,m3))   ;rrp   = 0.0 
         allocate(rhp  (m1,m2,m3))   ;rhp   = 0.0 
         allocate(rgp  (m1,m2,m3))   ;rgp   = 0.0 
         allocate(rpp  (m1,m2,m3))   ;rpp   = 0.0 
         allocate(rsp  (m1,m2,m3))   ;rsp   = 0.0 
         allocate(rei  (m1,m2,m3))   ;rei   = 0.0 
         allocate(rel  (m1,m2,m3))   ;rel   = 0.0 

         !---------  Local vars
         allocate(th        ( ims:ime, kms:kme, jms:jme )) ;th       = 0.0 
         allocate(dz8w      ( ims:ime, kms:kme, jms:jme )) ;dz8w     = 0.0 
         allocate(pi_phy    ( ims:ime, kms:kme, jms:jme )) ;pi_phy   = 0.0 
         allocate(p         ( ims:ime, kms:kme, jms:jme )) ;p        = 0.0 
         allocate(air_dens  ( ims:ime, kms:kme, jms:jme )) ;air_dens = 0.0 
         allocate(w         ( ims:ime, kms:kme, jms:jme )) ;w        = 0.0 
         allocate(qv_curr   ( ims:ime, kms:kme, jms:jme )) ;qv_curr  = 0.0 
         allocate(qc_curr   ( ims:ime, kms:kme, jms:jme )) ;qc_curr  = 0.0 
         allocate(qr_curr   ( ims:ime, kms:kme, jms:jme )) ;qr_curr  = 0.0
         allocate(qi_curr   ( ims:ime, kms:kme, jms:jme )) ;qi_curr  = 0.0
         allocate(qs_curr   ( ims:ime, kms:kme, jms:jme )) ;qs_curr  = 0.0
         allocate(qg_curr   ( ims:ime, kms:kme, jms:jme )) ;qg_curr  = 0.0
         allocate(qh_curr   ( ims:ime, kms:kme, jms:jme )) ;qh_curr  = 0.0
         allocate(re_cloud  ( ims:ime, kms:kme, jms:jme )) ;re_cloud = 0.0
         allocate(re_ice    ( ims:ime, kms:kme, jms:jme )) ;re_ice   = 0.0
         allocate(re_snow   ( ims:ime, kms:kme, jms:jme )) ;re_snow  = 0.0
         allocate(orho      ( ims:ime, kms:kme, jms:jme )) ;orho     = 0.0   
         allocate(hgt       ( ims:ime, kms:kme, jms:jme )) ;hgt      = 0.0   
         allocate(refl_10cm ( ims:ime, kms:kme, jms:jme )) ;refl_10cm= 0.0   
         allocate(rainprod  ( ims:ime, kms:kme, jms:jme )) ;rainprod = 0.0   
         allocate(evapprod  ( ims:ime, kms:kme, jms:jme )) ;evapprod = 0.0   
      
         allocate(rainnc    ( ims:ime,jms:jme )) ;rainnc       = 0.0
         allocate(rainncv   ( ims:ime,jms:jme )) ;rainncv      = 0.0
         allocate(snownc    ( ims:ime,jms:jme )) ;snownc       = 0.0
         allocate(snowncv   ( ims:ime,jms:jme )) ;snowncv      = 0.0
         allocate(graupelnc ( ims:ime,jms:jme )) ;graupelnc    = 0.0
         allocate(graupelncv( ims:ime,jms:jme )) ;graupelncv   = 0.0
         allocate(hail      ( ims:ime,jms:jme )) ;hail         = 0.0
         allocate(hailncv   ( ims:ime,jms:jme )) ;hailncv      = 0.0
         allocate(sr        ( ims:ime,jms:jme )) ;sr           = 0.0

      endif

      read(l_unit) mcphys_type , ilwrtyp, iswrtyp
      print*,'mcyphys=', mcphys_type
      read(l_unit) time, dtlt
      read(l_unit)  g      
      read(l_unit) cp      
      read(l_unit) cpv     
      read(l_unit) r_d     
      read(l_unit) r_v     
      read(l_unit) svpt0   
      read(l_unit) ep_1    
      read(l_unit) ep_2    
      read(l_unit) epsilon 
      read(l_unit) xls     
      read(l_unit) xlv     
      read(l_unit) xlf     
      read(l_unit) rhoair0 
      read(l_unit) rhowater  
      read(l_unit) rhosnow      
      read(l_unit) cliq    
      read(l_unit) cice    
      read(l_unit) psat
      read(l_unit) dzt
      read(l_unit) zt
      read(l_unit) glon
      read(l_unit) glat

      !-- 1st section all microphysics
      read(l_unit) thp
      read(l_unit) theta
      read(l_unit) pp   
      read(l_unit) rtp  
      read(l_unit) rv   
      read(l_unit) wp   
      read(l_unit) dn0  
      read(l_unit) pi0  
      read(l_unit) rcp    
      read(l_unit) rrp    

 
      read(l_unit) rtgt 
      read(l_unit) accpr
      read(l_unit) pcprr

      if(mcphys_type == 7 )  then 
         read(l_unit)  rhp  
         read(l_unit)  accph
         read(l_unit)  pcprh
      endif
      
      if(mcphys_type == 6 .or. mcphys_type == 7) then
         read(l_unit)  rgp  
         read(l_unit)  accpg
         read(l_unit)  pcprg
      endif

      if(mcphys_type == 5 .or. mcphys_type == 6 .or. mcphys_type == 7 ) then
         read(l_unit)  rpp  
         read(l_unit)  rsp  
         read(l_unit)  accps
         read(l_unit)  pcprs
      endif
   close(l_unit)
   
   
   !------------------------------ process cloud microphysics -----------------------
   do j = ja,jz
    do i = ia,iz
 
    !- column quantities
      thp_1d   (1:m1)= thp  (1:m1,i,j)
      theta_1d (1:m1)= theta(1:m1,i,j)
      pp_1d    (1:m1)= pp   (1:m1,i,j)
      rtp_1d   (1:m1)= rtp  (1:m1,i,j)
      rv_1d    (1:m1)= rv   (1:m1,i,j)
      wp_1d    (1:m1)= wp   (1:m1,i,j)
      dn0_1d   (1:m1)= dn0  (1:m1,i,j)
      pi0_1d   (1:m1)= pi0  (1:m1,i,j)
           
      !--- mass mixing ratio
      rcp_1d   (1:m1)= rcp    (1:m1,i,j)
      rrp_1d   (1:m1)= rrp    (1:m1,i,j)
           
      !- surface quantities
      rtgt_1d  = rtgt (i,j)
      accpr_1d = accpr(i,j)
      pcprr_1d = pcprr(i,j)
      
      if(mcphys_type == 7 )  then 
        rhp_1d  (1:m1)= rhp  (1:m1,i,j)
        accph_1d      = accph(i,j)
        pcprh_1d      = pcprh(i,j)
      else
        rhp_1d  (1:m1)= 0.
        accph_1d      = 0.
        pcprh_1d      = 0.
      endif

      if(mcphys_type == 6 .or. mcphys_type == 7) then
        rgp_1d  (1:m1)= rgp  (1:m1,i,j)
        accpg_1d      = accpg(i,j)
        pcprg_1d      = pcprg(i,j)
      else
        rgp_1d  (1:m1)= 0.
        accpg_1d      = 0.
        pcprg_1d      = 0.
      endif

      if(mcphys_type == 5 .or. mcphys_type == 6 .or. mcphys_type == 7 ) then
        rpp_1d  (1:m1)= rpp  (1:m1,i,j)
        rsp_1d  (1:m1)= rsp  (1:m1,i,j)
        accps_1d      = accps(i,j)
        pcprs_1d      = pcprs(i,j)      
      else
        rpp_1d  (1:m1)= 0.
        rsp_1d  (1:m1)= 0.
        accps_1d      = 0.
        pcprs_1d      = 0.  
      endif

      !--- coupling with the WSM MPs
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
      RAINNC    (1,1)=  accpr_1d !- rain+ice+snow+graupel+hail
      SNOWNC    (1,1)=  accps_1d !- ice+snow
      GRAUPELNC (1,1)=  accpg_1d !- graupel
      HAIL      (1,1)=  accph_1d !- hail

      do k = 1,kme-1
        
        kr = k + 1
        qv_curr (1,k,1)= max(1.e-12,rtp_1d(kr) - &    ! QV
                 (rcp_1d(kr)+rrp_1d(kr)+rpp_1d(kr)+rsp_1d(kr)+rgp_1d(kr)+rhp_1d(kr)))   
        qc_curr (1,k,1)= max(0.0,rcp_1d(kr))       ! QC     
        qr_curr (1,k,1)= max(0.0,rrp_1d(kr))       ! QR   
        qi_curr (1,k,1)= max(0.0,rpp_1d(kr))       ! QI   
        qs_curr (1,k,1)= max(0.0,rsp_1d(kr))       ! QS   
        qg_curr (1,k,1)= max(0.0,rgp_1d(kr))       ! QG

        qh_curr (1,k,1)= max(0.0,rhp_1d(kr))       ! QH   

        pi_phy  (1,k,1)= (pp_1d(kr)+pi0_1d(kr))*cpi ! Exner function/cp (dimensionless)

        P   (1,k,1)= ( (pp_1d(kr)+pi0_1d(kr))*cpi )** cpor * p00      ! pressure(Pa)
        W   (1,k,1)= wp_1d(kr)    ! vertical velocity (m/s) ! must be at center or face? ASK

        dz8w   (1,k,1)= rtgt_1d/dzt(kr) ! layer thickness (m) 
        !print*,'dz8',k,dz8w   (1,k,1)
      enddo

      !- get potential temperature (theta) from theta_il (thp) and condensates
      DO k=1,kme -1 
         kr = k + 1
         tempK    = theta_1d(kr)* (pp_1d(kr)+pi0_1d(kr))*cpi 
            til   = thp_1d  (kr)* (pp_1d(kr)+pi0_1d(kr))*cpi 
       
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
           

                 
      !- updated variables after microphysics processes (from WSM to BRAMS)
         DO k=1,kme-1
          kr=k+1
          rtp_1d(kr)= qv_curr(1,k,1) + &
                      qc_curr(1,k,1) + &     
                      qr_curr(1,k,1) + &    
                      qi_curr(1,k,1) + &    
                      qs_curr(1,k,1) + &    
                      qg_curr(1,k,1) + &    
                      qh_curr(1,k,1)  

          rcp_1d(kr)= qc_curr(1,k,1)
          rrp_1d(kr)= qr_curr(1,k,1)
          rpp_1d(kr)= qi_curr(1,k,1)
          rsp_1d(kr)= qs_curr(1,k,1)
          rgp_1d(kr)= qg_curr(1,k,1)
          rhp_1d(kr)= qh_curr(1,k,1)  
          
          rv_1d(kr)= max(1.0e-12, rtp_1d(kr) -(rcp_1d(kr)+rrp_1d(kr)+rpp_1d(kr)+rsp_1d(kr)+rgp_1d(kr)+rhp_1d(kr)))
       
          theta_1d(kr) =  TH(1,k,1)
          tempK        =  TH(1,k,1)*pi_phy(1,k,1)
          
          rliq     =  qc_curr(1,k,1) + qr_curr(1,k,1)       
          rice     =  qi_curr(1,k,1) + qs_curr(1,k,1) + qg_curr(1,k,1) + qh_curr(1,k,1)
          
          !- update liq-ice potential temperature THP in Kelvin including microphysics processes
          thp_1d(kr)  =   TH(1,k,1)*(1. + alvl * rliq/(cp * max(tempK,253.))  &
                                        + alvi * rice/(cp * max(tempK,253.)) ) **(-1.0)      
         ENDDO
         !- definition for k=1
           rtp_1d(1)  = rtp_1d(2)  
           rcp_1d(1)  = rcp_1d(2)  
           rrp_1d(1)  = rrp_1d(2)  
           rpp_1d(1)  = rpp_1d(2)  
           rsp_1d(1)  = rsp_1d(2)  
           rgp_1d(1)  = rgp_1d(2)  
           rhp_1d(1)  = rhp_1d(2)
           rv_1d (1)  = rv_1d (2)  
           thp_1d(1)  = thp_1d(2)  
           theta_1d(1)= theta_1d(2)
        
        IF( (ilwrtyp==6 .or. iswrtyp==6)) then           
          DO k=1,kme-1
            kr=k+1
            rel_1d (kr) = re_cloud (1,k,1) * 1.e+6 ! RRTM requires in micrometer
            rei_1d (kr) = re_ice   (1,k,1) * 1.e+6 ! RRTM requires in micrometer
          ENDDO
          rel_1d (1) =rel_1d (2) !;  rel (kme) =rel (kme-1) 
          rei_1d (1) =rei_1d (2) !;  rei (kme) =rei (kme-1)
        ENDIF        
        !- surface precipitation (units are kg/m^2 = mm)
        !- RAINNC and RAINNCV constains all precipitation hidrometeors (rain, graupel, snow, ...)
        accpr_1d = RAINNC    (1,1) ! a = accum
        pcprr_1d = RAINNCV   (1,1) ! p = for each dt  (or per time step)
        accps_1d = SNOWNC    (1,1) 
        pcprs_1d = SNOWNCV   (1,1) 
        accpg_1d = GRAUPELNC (1,1) 
        pcprg_1d = GRAUPELNCV(1,1) 
        accph_1d = HAIL      (1,1) 
        pcprh_1d = HAILNCV   (1,1) 

        !- column quantities
         thp  (1:m1,i,j) =thp_1d  (1:m1) 
         theta(1:m1,i,j) =theta_1d(1:m1)
         rtp  (1:m1,i,j) =rtp_1d  (1:m1)
         rv   (1:m1,i,j) =rv_1d   (1:m1)  
        
         rcp     (1:m1,i,j) =rcp_1d  (1:m1)   
         rrp     (1:m1,i,j) =rrp_1d  (1:m1)

         if(mcphys_type == 5 .or. mcphys_type == 6 .or. mcphys_type == 7 ) then   
           rpp    (1:m1,i,j) =rpp_1d (1:m1)   
           rsp    (1:m1,i,j) =rsp_1d (1:m1)   
         endif
         
         if(mcphys_type == 6 .or. mcphys_type == 7) rgp(1:m1,i,j) =rgp_1d  (1:m1)  
         if(mcphys_type == 7           )            rhp(1:m1,i,j) =rhp_1d  (1:m1)
         
         if( ilwrtyp==6 .or. iswrtyp==6 ) then
           rei   (1:m1,i,j) =rei_1d  (1:m1)
           rel   (1:m1,i,j) =rel_1d  (1:m1)
         endif 

         !- surface quantities
         accpr(i,j) = accpr_1d !constains all precipitation hidrometeors (rain, graupel, snow, ...)
         pcprr(i,j) = pcprr_1d !constains all precipitation hidrometeors (rain, graupel, snow, ...)
         if(mcphys_type == 5 .or. mcphys_type == 6 .or. mcphys_type == 7 ) then
           accps(i,j) = accps_1d
           pcprs(i,j) = pcprs_1d
         endif
         if(mcphys_type == 6 .or. mcphys_type == 7) then
           accpg(i,j) = accpg_1d
           pcprg(i,j) = pcprg_1d
         endif
         if(mcphys_type == 7) then
           accph(i,j) = accph_1d
           pcprh(i,j) = pcprh_1d
         endif
   enddo; enddo ! loop i,j
   
   !=======================================================================
    inquire (iolength=int_byte_size) real_byte_size  ! inquire by output list
    rec_size = m2*m3*int_byte_size    

    nz1 = 2
    nz2 = m1-1
    
    fileout=trim(prefix)//"_WSM"//trim(cmcphys_type)//"_dataOut-"//ctime//"-"//trim(cmynum)
    print*," ===> writing output in file: ",trim(var_dataout)//trim(fileout)


    open(newunit = l_unit, file=trim(var_dataout)//trim(fileout)//".gra", &
             form='unformatted', access='direct', status='replace', recl=rec_size)
        irec=1
        nlev= m1-1
        nvar_out = 0
        !-- 3d updated fields
        do nz=nz1,nz2;  write(l_unit,rec=irec) wp         (nz,1:m2,1:m3); irec=irec+1 ;enddo; nvar_out = nvar_out + 1
        do nz=nz1,nz2;  write(l_unit,rec=irec) thp        (nz,1:m2,1:m3); irec=irec+1 ;enddo; nvar_out = nvar_out + 1
        do nz=nz1,nz2;  write(l_unit,rec=irec) theta      (nz,1:m2,1:m3); irec=irec+1 ;enddo; nvar_out = nvar_out + 1
        do nz=nz1,nz2;  write(l_unit,rec=irec) 1000.*rtp  (nz,1:m2,1:m3); irec=irec+1 ;enddo; nvar_out = nvar_out + 1
        do nz=nz1,nz2;  write(l_unit,rec=irec) 1000.*rv   (nz,1:m2,1:m3); irec=irec+1 ;enddo; nvar_out = nvar_out + 1
        do nz=nz1,nz2;  write(l_unit,rec=irec) 1000.*rrp  (nz,1:m2,1:m3); irec=irec+1 ;enddo; nvar_out = nvar_out + 1
        do nz=nz1,nz2;  write(l_unit,rec=irec) 1000.*rcp  (nz,1:m2,1:m3); irec=irec+1 ;enddo; nvar_out = nvar_out + 1

        !-- 2d 
        write(l_unit,rec=irec) 3600*accpr(1:m2,1:m3); irec=irec+1 ; nvar_out = nvar_out + 1
        write(l_unit,rec=irec) 3600*pcprr(1:m2,1:m3); irec=irec+1 ; nvar_out = nvar_out + 1

    close(l_unit)  
    !== number of levels 2: m1 -1 
    nlev = nz2-nz1+1
   
    open(newunit = l_unit, file=trim(var_dataout)//trim(fileout)//".ctl", action='write', status='replace')

       fileout=trim(fileout)//".gra"
       write(l_unit,*) 'dset ^'//trim(fileout) 
       !writing others infos to ctl
       write(l_unit,*) 'undef -0.9990000E+34'
       write(l_unit,*) 'title WSM '
       write(l_unit,*) 'xdef ',m2,' linear ',glon(1,1),glon(2,1)-glon(1,1)
       write(l_unit,*) 'ydef ',m3,' linear ',glat(1,1),glat(1,2)-glat(1,1)
       write(l_unit,2005) nlev, zt(nz1:nz2)
       write(l_unit,*) 'tdef 1 linear 00:00Z01JAN200 1mo'
       write(l_unit,*) 'vars ' , nvar_out
       write(l_unit,*) 'wp'    ,nlev,'99 ','K'
       write(l_unit,*) 'thp'   ,nlev,'99 ','K'
       write(l_unit,*) 'theta' ,nlev,'99 ','K'
       write(l_unit,*) 'rtp'   ,nlev,'99 ','g/kg'
       write(l_unit,*) 'rv'    ,nlev,'99 ','g/kg'
       write(l_unit,*) 'rrp'   ,nlev,'99 ','g/kg'
       write(l_unit,*) 'rcp'   ,nlev,'99 ','g/kg'
       write(l_unit,*) 'apr',' 0 ','99 ','K'
       write(l_unit,*) 'ppr',' 0 ','99 ','K'
 
       write(l_unit,*) 'endvars'

    close(l_unit)

2005 format('zdef ',i4,' levels ',100f8.1)
2001 format('dset ',A100)

end program wsm_test
