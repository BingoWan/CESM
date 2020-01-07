!!! start of uwshcu

  subroutine compute_uwshcu( mix      , mkx       , iend         , ncnst    , dt        , &
                             ps0_in   , zs0_in    , p0_in        , z0_in    , dp0_in    , &
                             u0_in    , v0_in     , qv0_in       , ql0_in   , qi0_in    , &
                             t0_in    , s0_in     , tr0_in       ,                        &
                             tke_in   , cldfrct_in, concldfrct_in,  pblh_in , cush_inout, &
                             umf_out  , slflx_out , qtflx_out    ,                        &
                             flxprc1_out  , flxsnow1_out  , 			 	  &
                             qvten_out, qlten_out , qiten_out    ,                        &
                             sten_out , uten_out  , vten_out     , trten_out,             &
                             qrten_out, qsten_out , precip_out   , snow_out , evapc_out , &
                             cufrc_out, qcu_out   , qlu_out      , qiu_out  ,             &
                             cbmf_out , qc_out    , rliq_out     ,                        &
                             cnt_out  , cnb_out   , lchnk        , dpdry0_in )

    ! ------------------------------------------------------------ !
    !                                                              !
    !  University of Washington Shallow Convection Scheme          !
    !                                                              !
    !  Described in Park and Bretherton. 2008. J. Climate :        !
    !                                                              !
    ! 'The University of Washington shallow convection and         !
    !  moist turbulent schemes and their impact on climate         !
    !  simulations with the Community Atmosphere Model'            !
    !                                                              !
    !  Coded by Sungsu Park. Oct.2005.                             !
    !                        May.2008.                             !
    !  For questions, send an email to sungsup@ucar.edu or         !
    !                                  sungsu@atmos.washington.edu !
    !                                                              !
    ! ------------------------------------------------------------ !

    use cam_history,     only : outfld, addfld, phys_decomp
    use constituents,    only : qmin, cnst_get_type_byind, cnst_get_ind
    use wv_saturation,   only : findsp_vc

    implicit none

    ! ---------------------- !
    ! Input-Output Variables !
    ! ---------------------- !

!!! start of initialization

    integer , intent(in)    :: lchnk
    ! add for extend qsat
    real(r8) :: t_tmp     ! intermediate temperature for es look-up
    real(r8) :: weight ! Weight for interpolation
!!! end of initialization

    !--------------------------------------------------------------!
    !                                                              !
    ! Start the column i loop where i is a horizontal column index !
    !                                                              !
    !--------------------------------------------------------------!

    ! Compute wet-bulb temperature and specific humidity
    ! for treating evaporation of precipitation.

    ! "True" means ice will be taken into account
    do k = 1, mkx
       call findsp_vc(qv0_in(:iend,k), t0_in(:iend,k), p0_in(:iend,k), .true., &
            tw0_in(:iend,k), qw0_in(:iend,k))
    end do
  do i = 1, iend

      id_exit = .false.
      id_exit_ex = .false.

      ! -------------------------------------------- !
      ! Define 1D input variables at each grid point !
      ! -------------------------------------------- !

      ps0(0:mkx)       = ps0_in(i,0:mkx)
      zs0(0:mkx)       = zs0_in(i,0:mkx)
      p0(:mkx)         = p0_in(i,:mkx)
      z0(:mkx)         = z0_in(i,:mkx)
      dp0(:mkx)        = dp0_in(i,:mkx)
      dpdry0(:mkx)     = dpdry0_in(i,:mkx)
      u0(:mkx)         = u0_in(i,:mkx)
      v0(:mkx)         = v0_in(i,:mkx)
      qv0(:mkx)        = qv0_in(i,:mkx)
      ql0(:mkx)        = ql0_in(i,:mkx)
      qi0(:mkx)        = qi0_in(i,:mkx)
      t0(:mkx)         = t0_in(i,:mkx)
      s0(:mkx)         = s0_in(i,:mkx)
      tke(0:mkx)       = tke_in(i,0:mkx)
      cldfrct(:mkx)    = cldfrct_in(i,:mkx)
      concldfrct(:mkx) = concldfrct_in(i,:mkx)
      pblh             = pblh_in(i)
      cush             = cush_inout(i)
      do m = 1, ncnst
         tr0(:mkx,m)   = tr0_in(i,:mkx,m)
      enddo

!!! start count
      ! --------------------------------------------------------- !
      ! Compute other basic thermodynamic variables directly from !
      ! the input variables at each grid point                    !
      ! --------------------------------------------------------- !

      !----- 1. Compute internal environmental variables

      exn0(:mkx)   = (p0(:mkx)/p00)**rovcp
      exns0(0:mkx) = (ps0(0:mkx)/p00)**rovcp
      qt0(:mkx)    = (qv0(:mkx) + ql0(:mkx) + qi0(:mkx))
      thl0(:mkx)   = (t0(:mkx) - xlv*ql0(:mkx)/cp - xls*qi0(:mkx)/cp)/exn0(:mkx)
      thvl0(:mkx)  = (1._r8 + zvir*qt0(:mkx))*thl0(:mkx)

      !----- 2. Compute slopes of environmental variables in each layer
      !         Dimension of ssthl0(:mkx) is implicit.

      ssthl0       = slope(mkx,thl0,p0)
      ssqt0        = slope(mkx,qt0 ,p0)
      ssu0         = slope(mkx,u0  ,p0)
      ssv0         = slope(mkx,v0  ,p0)
      do m = 1, ncnst
         sstr0(:mkx,m) = slope(mkx,tr0(:mkx,m),p0)
      enddo

      !----- 3. Compute "thv0" and "thvl0" at the top/bottom interfaces in each layer
      !         There are computed from the reconstructed thl, qt at the top/bottom.

      do k = 1, mkx

         thl0bot = thl0(k) + ssthl0(k)*(ps0(k-1) - p0(k))
         qt0bot  = qt0(k)  + ssqt0(k) *(ps0(k-1) - p0(k))
         call conden(ps0(k-1),thl0bot,qt0bot,thj,qvj,qlj,qij,qse,id_check)
         if( id_check .eq. 1 ) then
             exit_conden(i) = 1._r8
             id_exit = .true.
             go to 333
         end if
         thv0bot(k)  = thj*(1._r8 + zvir*qvj - qlj - qij)
         thvl0bot(k) = thl0bot*(1._r8 + zvir*qt0bot)

         thl0top = thl0(k) + ssthl0(k)*(ps0(k) - p0(k))
         qt0top  =  qt0(k) + ssqt0(k) *(ps0(k) - p0(k))
         call conden(ps0(k),thl0top,qt0top,thj,qvj,qlj,qij,qse,id_check)
         if( id_check .eq. 1 ) then
             exit_conden(i) = 1._r8
             id_exit = .true.
             go to 333
         end if
         thv0top(k)  = thj*(1._r8 + zvir*qvj - qlj - qij)
         thvl0top(k) = thl0top*(1._r8 + zvir*qt0top)

      end do

      ! ------------------------------------------------------------ !
      ! Save input and related environmental thermodynamic variables !
      ! for use at "iter_cin=2" when "del_CIN >= 0"                  !
      ! ------------------------------------------------------------ !

      qv0_o(:mkx)          = qv0(:mkx)
      ql0_o(:mkx)          = ql0(:mkx)
      qi0_o(:mkx)          = qi0(:mkx)
      t0_o(:mkx)           = t0(:mkx)
      s0_o(:mkx)           = s0(:mkx)
      u0_o(:mkx)           = u0(:mkx)
      v0_o(:mkx)           = v0(:mkx)
      qt0_o(:mkx)          = qt0(:mkx)
      thl0_o(:mkx)         = thl0(:mkx)
      thvl0_o(:mkx)        = thvl0(:mkx)
      ssthl0_o(:mkx)       = ssthl0(:mkx)
      ssqt0_o(:mkx)        = ssqt0(:mkx)
      thv0bot_o(:mkx)      = thv0bot(:mkx)
      thv0top_o(:mkx)      = thv0top(:mkx)
      thvl0bot_o(:mkx)     = thvl0bot(:mkx)
      thvl0top_o(:mkx)     = thvl0top(:mkx)
      ssu0_o(:mkx)         = ssu0(:mkx)
      ssv0_o(:mkx)         = ssv0(:mkx)
      do m = 1, ncnst
         tr0_o(:mkx,m)     = tr0(:mkx,m)
         sstr0_o(:mkx,m)   = sstr0(:mkx,m)
      enddo

      ! ---------------------------------------------- !
      ! Initialize output variables at each grid point !
      ! ---------------------------------------------- !

      umf(0:mkx)          = 0.0_r8
      emf(0:mkx)          = 0.0_r8
      slflx(0:mkx)        = 0.0_r8
      qtflx(0:mkx)        = 0.0_r8
      uflx(0:mkx)         = 0.0_r8
      vflx(0:mkx)         = 0.0_r8
      qvten(:mkx)         = 0.0_r8
      qlten(:mkx)         = 0.0_r8
      qiten(:mkx)         = 0.0_r8
      sten(:mkx)          = 0.0_r8
      uten(:mkx)          = 0.0_r8
      vten(:mkx)          = 0.0_r8
      qrten(:mkx)         = 0.0_r8
      qsten(:mkx)         = 0.0_r8
      dwten(:mkx)         = 0.0_r8
      diten(:mkx)         = 0.0_r8
      precip              = 0.0_r8
      snow                = 0.0_r8
      evapc(:mkx)         = 0.0_r8
      cufrc(:mkx)         = 0.0_r8
      qcu(:mkx)           = 0.0_r8
      qlu(:mkx)           = 0.0_r8
      qiu(:mkx)           = 0.0_r8
      fer(:mkx)           = 0.0_r8
      fdr(:mkx)           = 0.0_r8
      cin                 = 0.0_r8
      cbmf                = 0.0_r8
      qc(:mkx)            = 0.0_r8
      qc_l(:mkx)          = 0.0_r8
      qc_i(:mkx)          = 0.0_r8
      rliq                = 0.0_r8
      cnt                 = real(mkx, r8)
      cnb                 = 0.0_r8
      qtten(:mkx)         = 0.0_r8
      slten(:mkx)         = 0.0_r8
      ufrc(0:mkx)         = 0.0_r8

      thlu(0:mkx)         = 0.0_r8
      qtu(0:mkx)          = 0.0_r8
      uu(0:mkx)           = 0.0_r8
      vu(0:mkx)           = 0.0_r8
      wu(0:mkx)           = 0.0_r8
      thvu(0:mkx)         = 0.0_r8
      thlu_emf(0:mkx)     = 0.0_r8
      qtu_emf(0:mkx)      = 0.0_r8
      uu_emf(0:mkx)       = 0.0_r8
      vu_emf(0:mkx)       = 0.0_r8

      ufrcinvbase         = 0.0_r8
      ufrclcl             = 0.0_r8
      winvbase            = 0.0_r8
      wlcl                = 0.0_r8
      emfkbup             = 0.0_r8
      cbmflimit           = 0.0_r8
      excessu_arr(:mkx)   = 0.0_r8
      excess0_arr(:mkx)   = 0.0_r8
      xc_arr(:mkx)        = 0.0_r8
      aquad_arr(:mkx)     = 0.0_r8
      bquad_arr(:mkx)     = 0.0_r8
      cquad_arr(:mkx)     = 0.0_r8
      bogbot_arr(:mkx)    = 0.0_r8
      bogtop_arr(:mkx)    = 0.0_r8

      uemf(0:mkx)         = 0.0_r8
      comsub(:mkx)        = 0.0_r8
      qlten_sink(:mkx)    = 0.0_r8
      qiten_sink(:mkx)    = 0.0_r8
      nlten_sink(:mkx)    = 0.0_r8
      niten_sink(:mkx)    = 0.0_r8

      do m = 1, ncnst
         trflx(0:mkx,m)   = 0.0_r8
         trten(:mkx,m)    = 0.0_r8
         tru(0:mkx,m)     = 0.0_r8
         tru_emf(0:mkx,m) = 0.0_r8
      enddo

    !-----------------------------------------------!
    ! Below 'iter' loop is for implicit CIN closure !
    !-----------------------------------------------!

    ! ----------------------------------------------------------------------------- !
    ! It is important to note that this iterative cin loop is located at the outest !
    ! shell of the code. Thus, source air properties can also be changed during the !
    ! iterative cin calculation, because cumulus convection induces non-zero fluxes !
    ! even at interfaces below PBL top height through 'fluxbelowinv' subroutine.    !
    ! ----------------------------------------------------------------------------- !
    do iter = 1, iter_cin

       ! ---------------------------------------------------------------------- !
       ! Cumulus scale height                                                   !
       ! In contrast to the premitive code, cumulus scale height is iteratively !
       ! calculated at each time step, and at each iterative cin step.          !
       ! It is not clear whether I should locate below two lines within or  out !
       ! of the iterative cin loop.                                             !
       ! ---------------------------------------------------------------------- !

       tscaleh = cush
       cush    = -1._r8

       ! ----------------------------------------------------------------------- !
       ! Find PBL top height interface index, 'kinv-1' where 'kinv' is the layer !
       ! index with PBLH in it. When PBLH is exactly at interface, 'kinv' is the !
       ! layer index having PBLH as a lower interface.                           !
       ! In the previous code, I set the lower limit of 'kinv' by 2  in order to !
       ! be consistent with the other parts of the code. However in the modified !
       ! code, I allowed 'kinv' to be 1 & if 'kinv = 1', I just exit the program !
       ! without performing cumulus convection. This new approach seems to be    !
       ! more reasonable: if PBL height is within 'kinv=1' layer, surface is STL !
       ! interface (bflxs <= 0) and interface just above the surface should be   !
       ! either non-turbulent (Ri>0.19) or stably turbulent (0<=Ri<0.19 but this !
       ! interface is identified as a base external interface of upperlying CL.  !
       ! Thus, when 'kinv=1', PBL scheme guarantees 'bflxs <= 0'.  For this case !
       ! it is reasonable to assume that cumulus convection does not happen.     !
       ! When these is SBCL, PBL height from the PBL scheme is likely to be very !
       ! close at 'kinv-1' interface, but not exactly, since 'zi' information is !
       ! changed between two model time steps. In order to ensure correct identi !
       ! fication of 'kinv' for general case including SBCL, I imposed an offset !
       ! of 5 [m] in the below 'kinv' finding block.                             !
       ! ----------------------------------------------------------------------- !

       do k = mkx - 1, 1, -1
          if( (pblh + 5._r8 - zs0(k))*(pblh + 5._r8 - zs0(k+1)) .lt. 0._r8 ) then
               kinv = k + 1
               go to 15
          endif
       end do
       kinv = 1
15     continue

       if( kinv .le. 1 ) then
           exit_kinv1(i) = 1._r8
           id_exit = .true.
           go to 333
       endif
       ! From here, it must be 'kinv >= 2'.

       ! -------------------------------------------------------------------------- !
       ! Find PBL averaged tke ('tkeavg') and minimum 'thvl' ('thvlmin') in the PBL !
       ! In the current code, 'tkeavg' is obtained by averaging all interfacial TKE !
       ! within the PBL. However, in order to be conceptually consistent with   PBL !
       ! scheme, 'tkeavg' should be calculated by considering surface buoyancy flux.!
       ! If surface buoyancy flux is positive ( bflxs >0 ), surface interfacial TKE !
       ! should be included in calculating 'tkeavg', while if bflxs <= 0,   surface !
       ! interfacial TKE should not be included in calculating 'tkeavg'.   I should !
       ! modify the code when 'bflxs' is available as an input of cumulus scheme.   !
       ! 'thvlmin' is a minimum 'thvl' within PBL obtained by comparing top &  base !
       ! interface values of 'thvl' in each layers within the PBL.                  !
       ! -------------------------------------------------------------------------- !

       dpsum    = 0._r8
       tkeavg   = 0._r8
       thvlmin  = 1000._r8
       do k = 0, kinv - 1   ! Here, 'k' is an interfacial layer index.
          if( k .eq. 0 ) then
              dpi = ps0(0) - p0(1)
          elseif( k .eq. (kinv-1) ) then
              dpi = p0(kinv-1) - ps0(kinv-1)
          else
              dpi = p0(k) - p0(k+1)
          endif
          dpsum  = dpsum  + dpi
          tkeavg = tkeavg + dpi*tke(k)
          if( k .ne. 0 ) thvlmin = min(thvlmin,min(thvl0bot(k),thvl0top(k)))
       end do
       tkeavg  = tkeavg/dpsum

       ! ------------------------------------------------------------------ !
       ! Find characteristics of cumulus source air: qtsrc,thlsrc,usrc,vsrc !
       ! Note that 'thlsrc' was con-cocked using 'thvlsrc' and 'qtsrc'.     !
       ! 'qtsrc' is defined as the lowest layer mid-point value;   'thlsrc' !
       ! is from 'qtsrc' and 'thvlmin=thvlsrc'; 'usrc' & 'vsrc' are defined !
       ! as the values just below the PBL top interface.                    !
       ! ------------------------------------------------------------------ !

       qtsrc   = qt0(1)
       thvlsrc = thvlmin
       thlsrc  = thvlsrc / ( 1._r8 + zvir * qtsrc )
       usrc    = u0(kinv-1) + ssu0(kinv-1) * ( ps0(kinv-1) - p0(kinv-1) )
       vsrc    = v0(kinv-1) + ssv0(kinv-1) * ( ps0(kinv-1) - p0(kinv-1) )
       do m = 1, ncnst
          trsrc(m) = tr0(1,m)
       enddo

       ! ------------------------------------------------------------------ !
       ! Find LCL of the source air and a layer index containing LCL (klcl) !
       ! When the LCL is exactly at the interface, 'klcl' is a layer index  !
       ! having 'plcl' as the lower interface similar to the 'kinv' case.   !
       ! In the previous code, I assumed that if LCL is located within the  !
       ! lowest model layer ( 1 ) or the top model layer ( mkx ), then  no  !
       ! convective adjustment is performed and just exited.   However, in  !
       ! the revised code, I relaxed the first constraint and  even though  !
       ! LCL is at the lowest model layer, I allowed cumulus convection to  !
       ! be initiated. For this case, cumulus convection should be started  !
       ! from the PBL top height, as shown in the following code.           !
       ! When source air is already saturated even at the surface, klcl is  !
       ! set to 1.                                                          !
       ! ------------------------------------------------------------------ !

       plcl = qsinvert(qtsrc,thlsrc,ps0(0))
       do k = 0, mkx
          if( ps0(k) .lt. plcl ) then
              klcl = k
              go to 25
          endif
       end do
       klcl = mkx
25     continue
       klcl = max(1,klcl)

       if( plcl .lt. 30000._r8 ) then
     ! if( klcl .eq. mkx ) then
           exit_klclmkx(i) = 1._r8
           id_exit = .true.
           go to 333
       endif

       ! ------------------------------------------------------------- !
       ! Calculate environmental virtual potential temperature at LCL, !
       !'thv0lcl' which is solely used in the 'cin' calculation. Note  !
       ! that 'thv0lcl' is calculated first by calculating  'thl0lcl'  !
       ! and 'qt0lcl' at the LCL, and performing 'conden' afterward,   !
       ! in fully consistent with the other parts of the code.         !
       ! ------------------------------------------------------------- !

       thl0lcl = thl0(klcl) + ssthl0(klcl) * ( plcl - p0(klcl) )
       qt0lcl  = qt0(klcl)  + ssqt0(klcl)  * ( plcl - p0(klcl) )
       call conden(plcl,thl0lcl,qt0lcl,thj,qvj,qlj,qij,qse,id_check)
       if( id_check .eq. 1 ) then
           exit_conden(i) = 1._r8
           id_exit = .true.
           go to 333
       end if
       thv0lcl = thj * ( 1._r8 + zvir * qvj - qlj - qij )

       ! ------------------------------------------------------------------------ !
       ! Compute Convective Inhibition, 'cin' & 'cinlcl' [J/kg]=[m2/s2] TKE unit. !
       !                                                                          !
       ! 'cin' (cinlcl) is computed from the PBL top interface to LFC (LCL) using !
       ! piecewisely reconstructed environmental profiles, assuming environmental !
       ! buoyancy profile within each layer ( or from LCL to upper interface in   !
       ! each layer ) is simply a linear profile. For the purpose of cin (cinlcl) !
       ! calculation, we simply assume that lateral entrainment does not occur in !
       ! updrafting cumulus plume, i.e., cumulus source air property is conserved.!
       ! Below explains some rules used in the calculations of cin (cinlcl).   In !
       ! general, both 'cin' and 'cinlcl' are calculated from a PBL top interface !
       ! to LCL and LFC, respectively :                                           !
       ! 1. If LCL is lower than the PBL height, cinlcl = 0 and cin is calculated !
       !    from PBL height to LFC.                                               !
       ! 2. If LCL is higher than PBL height,   'cinlcl' is calculated by summing !
       !    both positive and negative cloud buoyancy up to LCL using 'single_cin'!
       !    From the LCL to LFC, however, only negative cloud buoyancy is counted !
       !    to calculate final 'cin' upto LFC.                                    !
       ! 3. If either 'cin' or 'cinlcl' is negative, they are set to be zero.     !
       ! In the below code, 'klfc' is the layer index containing 'LFC' similar to !
       ! 'kinv' and 'klcl'.                                                       !
       ! ------------------------------------------------------------------------ !

        cin    = 0._r8
        cinlcl = 0._r8
        plfc   = 0._r8
        klfc   = mkx

        ! ------------------------------------------------------------------------- !
        ! Case 1. LCL height is higher than PBL interface ( 'pLCL <= ps0(kinv-1)' ) !
        ! ------------------------------------------------------------------------- !

        if( klcl .ge. kinv ) then

            do k = kinv, mkx - 1
               if( k .lt. klcl ) then
                   thvubot = thvlsrc
                   thvutop = thvlsrc
                   cin     = cin + single_cin(ps0(k-1),thv0bot(k),ps0(k),thv0top(k),thvubot,thvutop)
               elseif( k .eq. klcl ) then
                   !----- Bottom to LCL
                   thvubot = thvlsrc
                   thvutop = thvlsrc
                   cin     = cin + single_cin(ps0(k-1),thv0bot(k),plcl,thv0lcl,thvubot,thvutop)
                   if( cin .lt. 0._r8 ) limit_cinlcl(i) = 1._r8
                   cinlcl  = max(cin,0._r8)
                   cin     = cinlcl
                   !----- LCL to Top
                   thvubot = thvlsrc
                   call conden(ps0(k),thlsrc,qtsrc,thj,qvj,qlj,qij,qse,id_check)
                   if( id_check .eq. 1 ) then
                       exit_conden(i) = 1._r8
                       id_exit = .true.
                       go to 333
                   end if
                   thvutop = thj * ( 1._r8 + zvir*qvj - qlj - qij )
                   call getbuoy(plcl,thv0lcl,ps0(k),thv0top(k),thvubot,thvutop,plfc,cin)
                   if( plfc .gt. 0._r8 ) then
                       klfc = k
                       go to 35
                   end if
               else
                   thvubot = thvutop
                   call conden(ps0(k),thlsrc,qtsrc,thj,qvj,qlj,qij,qse,id_check)
                   if( id_check .eq. 1 ) then
                       exit_conden(i) = 1._r8
                       id_exit = .true.
                       go to 333
                   end if
                   thvutop = thj * ( 1._r8 + zvir*qvj - qlj - qij )
                   call getbuoy(ps0(k-1),thv0bot(k),ps0(k),thv0top(k),thvubot,thvutop,plfc,cin)
                   if( plfc .gt. 0._r8 ) then
                       klfc = k
                       go to 35
                   end if
               endif
            end do

       ! ----------------------------------------------------------------------- !
       ! Case 2. LCL height is lower than PBL interface ( 'pLCL > ps0(kinv-1)' ) !
       ! ----------------------------------------------------------------------- !

       else
          cinlcl = 0._r8
          do k = kinv, mkx - 1
             call conden(ps0(k-1),thlsrc,qtsrc,thj,qvj,qlj,qij,qse,id_check)
             if( id_check .eq. 1 ) then
                 exit_conden(i) = 1._r8
                 id_exit = .true.
                 go to 333
             end if
             thvubot = thj * ( 1._r8 + zvir*qvj - qlj - qij )
             call conden(ps0(k),thlsrc,qtsrc,thj,qvj,qlj,qij,qse,id_check)
             if( id_check .eq. 1 ) then
                 exit_conden(i) = 1._r8
                 id_exit = .true.
                 go to 333
             end if
             thvutop = thj * ( 1._r8 + zvir*qvj - qlj - qij )
             call getbuoy(ps0(k-1),thv0bot(k),ps0(k),thv0top(k),thvubot,thvutop,plfc,cin)
             if( plfc .gt. 0._r8 ) then
                 klfc = k
                 go to 35
             end if
          end do
       endif  ! End of CIN case selection

 35    continue
       if( cin .lt. 0._r8 ) limit_cin(i) = 1._r8
       cin = max(0._r8,cin)
       if( klfc .ge. mkx ) then
           klfc = mkx
         ! write(iulog,*) 'klfc >= mkx'
           exit_klfcmkx(i) = 1._r8
           id_exit = .true.
           go to 333
       endif

       ! ---------------------------------------------------------------------- !
       ! In order to calculate implicit 'cin' (or 'cinlcl'), save the initially !
       ! calculated 'cin' and 'cinlcl', and other related variables. These will !
       ! be restored after calculating implicit CIN.                            !
       ! ---------------------------------------------------------------------- !

       if( iter .eq. 1 ) then
           cin_i       = cin
           cinlcl_i    = cinlcl
           ke          = rbuoy / ( rkfre * tkeavg + epsvarw )
           kinv_o      = kinv
           klcl_o      = klcl
           klfc_o      = klfc
           plcl_o      = plcl
           plfc_o      = plfc
           tkeavg_o    = tkeavg
           thvlmin_o   = thvlmin
           qtsrc_o     = qtsrc
           thvlsrc_o   = thvlsrc
           thlsrc_o    = thlsrc
           usrc_o      = usrc
           vsrc_o      = vsrc
           thv0lcl_o   = thv0lcl
           do m = 1, ncnst
              trsrc_o(m) = trsrc(m)
           enddo
       endif
!!!end count
    return

  end subroutine compute_uwshcu
!!! end of uwshcu
