!!! start of uwshcu
!!!start count
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
!!!end count
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
     ! Modification : If I impose w = max(0.1_r8, w) up to the top interface of
     !                klfc, I should only use cinlfc.  That is, if I want to
     !                use cinlcl, I should not impose w = max(0.1_r8, w).
     !                Using cinlcl is equivalent to treating only 'saturated'
     !                moist convection. Note that in this sense, I should keep
     !                the functionality of both cinlfc and cinlcl.
     !                However, the treatment of penetrative entrainment level becomes
     !                ambiguous if I choose 'cinlcl'. Thus, the best option is to use
     !                'cinlfc'.

       ! -------------------------------------------------------------------------- !
       ! Calculate implicit 'cin' by averaging initial and final cins.    Note that !
       ! implicit CIN is adopted only when cumulus convection stabilized the system,!
       ! i.e., only when 'del_CIN >0'. If 'del_CIN<=0', just use explicit CIN. Note !
       ! also that since 'cinlcl' is set to zero whenever LCL is below the PBL top, !
       ! (see above CIN calculation part), the use of 'implicit CIN=cinlcl'  is not !
       ! good. Thus, when using implicit CIN, always try to only use 'implicit CIN= !
       ! cin', not 'implicit CIN=cinlcl'. However, both 'CIN=cin' and 'CIN=cinlcl'  !
       ! are good when using explicit CIN.                                          !
       ! -------------------------------------------------------------------------- !
       if( iter .ne. 1 ) then

           cin_f = cin
           cinlcl_f = cinlcl
           if( use_CINcin ) then
               del_CIN = cin_f - cin_i
           else
               del_CIN = cinlcl_f - cinlcl_i
           endif

           if( del_CIN .gt. 0._r8 ) then

               ! -------------------------------------------------------------- !
               ! Calculate implicit 'cin' and 'cinlcl'. Note that when we chose !
               ! to use 'implicit CIN = cin', choose 'cinlcl = cinlcl_i' below: !
               ! because iterative CIN only aims to obtain implicit CIN,  once  !
               ! we obtained 'implicit CIN=cin', it is good to use the original !
               ! profiles information for all the other variables after that.   !
               ! Note 'cinlcl' will be explicitly used in calculating  'wlcl' & !
               ! 'ufrclcl' after calculating 'winv' & 'ufrcinv'  at the PBL top !
               ! interface later, after calculating 'cbmf'.                     !
               ! -------------------------------------------------------------- !

               alpha = compute_alpha( del_CIN, ke )
               cin   = cin_i + alpha * del_CIN
               if( use_CINcin ) then
                   cinlcl = cinlcl_i
               else
                   cinlcl = cinlcl_i + alpha * del_cinlcl
               endif

               ! ----------------------------------------------------------------- !
               ! Restore the original values from the previous 'iter_cin' step (1) !
               ! to compute correct tendencies for (n+1) time step by implicit CIN !
               ! ----------------------------------------------------------------- !

               kinv      = kinv_o
               klcl      = klcl_o
               klfc      = klfc_o
               plcl      = plcl_o
               plfc      = plfc_o
               tkeavg    = tkeavg_o
               thvlmin   = thvlmin_o
               qtsrc     = qtsrc_o
               thvlsrc   = thvlsrc_o
               thlsrc    = thlsrc_o
               usrc      = usrc_o
               vsrc      = vsrc_o
               thv0lcl   = thv0lcl_o
               do m = 1, ncnst
                  trsrc(m) = trsrc_o(m)
               enddo

               qv0(:mkx)            = qv0_o(:mkx)
               ql0(:mkx)            = ql0_o(:mkx)
               qi0(:mkx)            = qi0_o(:mkx)
               t0(:mkx)             = t0_o(:mkx)
               s0(:mkx)             = s0_o(:mkx)
               u0(:mkx)             = u0_o(:mkx)
               v0(:mkx)             = v0_o(:mkx)
               qt0(:mkx)            = qt0_o(:mkx)
               thl0(:mkx)           = thl0_o(:mkx)
               thvl0(:mkx)          = thvl0_o(:mkx)
               ssthl0(:mkx)         = ssthl0_o(:mkx)
               ssqt0(:mkx)          = ssqt0_o(:mkx)
               thv0bot(:mkx)        = thv0bot_o(:mkx)
               thv0top(:mkx)        = thv0top_o(:mkx)
               thvl0bot(:mkx)       = thvl0bot_o(:mkx)
               thvl0top(:mkx)       = thvl0top_o(:mkx)
               ssu0(:mkx)           = ssu0_o(:mkx)
               ssv0(:mkx)           = ssv0_o(:mkx)
               do m = 1, ncnst
                  tr0(:mkx,m)   = tr0_o(:mkx,m)
                  sstr0(:mkx,m) = sstr0_o(:mkx,m)
               enddo

               ! ------------------------------------------------------ !
               ! Initialize all fluxes, tendencies, and other variables !
               ! in association with cumulus convection.                !
               ! ------------------------------------------------------ !

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
               qc(:mkx)            = 0.0_r8
               qc_l(:mkx)          = 0.0_r8
               qc_i(:mkx)          = 0.0_r8
               rliq                = 0.0_r8
               cbmf                = 0.0_r8
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

               do m = 1, ncnst
                  trflx(0:mkx,m)   = 0.0_r8
                  trten(:mkx,m)    = 0.0_r8
                  tru(0:mkx,m)     = 0.0_r8
                  tru_emf(0:mkx,m) = 0.0_r8
               enddo

               ! -------------------------------------------------- !
               ! Below are diagnostic output variables for detailed !
               ! analysis of cumulus scheme.                        !
               ! -------------------------------------------------- !

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

          else ! When 'del_CIN < 0', use explicit CIN instead of implicit CIN.

               ! ----------------------------------------------------------- !
               ! Identifier showing whether explicit or implicit CIN is used !
               ! ----------------------------------------------------------- !

               ind_delcin(i) = 1._r8

               ! --------------------------------------------------------- !
               ! Restore original output values of "iter_cin = 1" and exit !
               ! --------------------------------------------------------- !

               umf_out(i,0:mkx)         = umf_s(0:mkx)
               qvten_out(i,:mkx)        = qvten_s(:mkx)
               qlten_out(i,:mkx)        = qlten_s(:mkx)
               qiten_out(i,:mkx)        = qiten_s(:mkx)
               sten_out(i,:mkx)         = sten_s(:mkx)
               uten_out(i,:mkx)         = uten_s(:mkx)
               vten_out(i,:mkx)         = vten_s(:mkx)
               qrten_out(i,:mkx)        = qrten_s(:mkx)
               qsten_out(i,:mkx)        = qsten_s(:mkx)
               precip_out(i)            = precip_s
               snow_out(i)              = snow_s
               evapc_out(i,:mkx)        = evapc_s(:mkx)
               cush_inout(i)            = cush_s
               cufrc_out(i,:mkx)        = cufrc_s(:mkx)
               slflx_out(i,0:mkx)       = slflx_s(0:mkx)
               qtflx_out(i,0:mkx)       = qtflx_s(0:mkx)
               qcu_out(i,:mkx)          = qcu_s(:mkx)
               qlu_out(i,:mkx)          = qlu_s(:mkx)
               qiu_out(i,:mkx)          = qiu_s(:mkx)
               cbmf_out(i)              = cbmf_s
               qc_out(i,:mkx)           = qc_s(:mkx)
               rliq_out(i)              = rliq_s
               cnt_out(i)               = cnt_s
               cnb_out(i)               = cnb_s
               do m = 1, ncnst
                  trten_out(i,:mkx,m)   = trten_s(:mkx,m)
               enddo

               ! ------------------------------------------------------------------------------ !
               ! Below are diagnostic output variables for detailed analysis of cumulus scheme. !
               ! The order of vertical index is reversed for this internal diagnostic output.   !
               ! ------------------------------------------------------------------------------ !

               fer_out(i,mkx:1:-1)      = fer_s(:mkx)
               fdr_out(i,mkx:1:-1)      = fdr_s(:mkx)
               cinh_out(i)              = cin_s
               cinlclh_out(i)           = cinlcl_s
               qtten_out(i,mkx:1:-1)    = qtten_s(:mkx)
               slten_out(i,mkx:1:-1)    = slten_s(:mkx)
               ufrc_out(i,mkx:0:-1)     = ufrc_s(0:mkx)
               uflx_out(i,mkx:0:-1)     = uflx_s(0:mkx)
               vflx_out(i,mkx:0:-1)     = vflx_s(0:mkx)

               ufrcinvbase_out(i)       = ufrcinvbase_s
               ufrclcl_out(i)           = ufrclcl_s
               winvbase_out(i)          = winvbase_s
               wlcl_out(i)              = wlcl_s
               plcl_out(i)              = plcl_s
               pinv_out(i)              = pinv_s
               plfc_out(i)              = plfc_s
               pbup_out(i)              = pbup_s
               ppen_out(i)              = ppen_s
               qtsrc_out(i)             = qtsrc_s
               thlsrc_out(i)            = thlsrc_s
               thvlsrc_out(i)           = thvlsrc_s
               emfkbup_out(i)           = emfkbup_s
               cbmflimit_out(i)         = cbmflimit_s
               tkeavg_out(i)            = tkeavg_s
               zinv_out(i)              = zinv_s
               rcwp_out(i)              = rcwp_s
               rlwp_out(i)              = rlwp_s
               riwp_out(i)              = riwp_s

               wu_out(i,mkx:0:-1)       = wu_s(0:mkx)
               qtu_out(i,mkx:0:-1)      = qtu_s(0:mkx)
               thlu_out(i,mkx:0:-1)     = thlu_s(0:mkx)
               thvu_out(i,mkx:0:-1)     = thvu_s(0:mkx)
               uu_out(i,mkx:0:-1)       = uu_s(0:mkx)
               vu_out(i,mkx:0:-1)       = vu_s(0:mkx)
               qtu_emf_out(i,mkx:0:-1)  = qtu_emf_s(0:mkx)
               thlu_emf_out(i,mkx:0:-1) = thlu_emf_s(0:mkx)
               uu_emf_out(i,mkx:0:-1)   = uu_emf_s(0:mkx)
               vu_emf_out(i,mkx:0:-1)   = vu_emf_s(0:mkx)
               uemf_out(i,mkx:0:-1)     = uemf_s(0:mkx)

               dwten_out(i,mkx:1:-1)    = dwten_s(:mkx)
               diten_out(i,mkx:1:-1)    = diten_s(:mkx)
               flxrain_out(i,mkx:0:-1)  = flxrain_s(0:mkx)
               flxsnow_out(i,mkx:0:-1)  = flxsnow_s(0:mkx)
               ntraprd_out(i,mkx:1:-1)  = ntraprd_s(:mkx)
               ntsnprd_out(i,mkx:1:-1)  = ntsnprd_s(:mkx)

               excessu_arr_out(i,mkx:1:-1)  = excessu_arr_s(:mkx)
               excess0_arr_out(i,mkx:1:-1)  = excess0_arr_s(:mkx)
               xc_arr_out(i,mkx:1:-1)       = xc_arr_s(:mkx)
               aquad_arr_out(i,mkx:1:-1)    = aquad_arr_s(:mkx)
               bquad_arr_out(i,mkx:1:-1)    = bquad_arr_s(:mkx)
               cquad_arr_out(i,mkx:1:-1)    = cquad_arr_s(:mkx)
               bogbot_arr_out(i,mkx:1:-1)   = bogbot_arr_s(:mkx)
               bogtop_arr_out(i,mkx:1:-1)   = bogtop_arr_s(:mkx)

               do m = 1, ncnst
                  trflx_out(i,mkx:0:-1,m)   = trflx_s(0:mkx,m)
                  tru_out(i,mkx:0:-1,m)     = tru_s(0:mkx,m)
                  tru_emf_out(i,mkx:0:-1,m) = tru_emf_s(0:mkx,m)
               enddo

               id_exit = .false.
               go to 333

          endif

       endif

       ! ------------------------------------------------------------------ !
       ! Define a release level, 'prel' and release layer, 'krel'.          !
       ! 'prel' is the lowest level from which buoyancy sorting occurs, and !
       ! 'krel' is the layer index containing 'prel' in it, similar to  the !
       ! previous definitions of 'kinv', 'klcl', and 'klfc'.    In order to !
       ! ensure that only PBL scheme works within the PBL,  if LCL is below !
       ! PBL top height, then 'krel = kinv', while if LCL is above  PBL top !
       ! height, then 'krel = klcl'.   Note however that regardless of  the !
       ! definition of 'krel', cumulus convection induces fluxes within PBL !
       ! through 'fluxbelowinv'.  We can make cumulus convection start from !
       ! any level, even within the PBL by appropriately defining 'krel'  & !
       ! 'prel' here. Then it must be accompanied by appropriate definition !
       ! of source air properties, CIN, and re-setting of 'fluxbelowinv', & !
       ! many other stuffs.                                                 !
       ! Note that even when 'prel' is located above the PBL top height, we !
       ! still have cumulus convection between PBL top height and 'prel':   !
       ! we simply assume that no lateral mixing occurs in this range.      !
       ! ------------------------------------------------------------------ !

       if( klcl .lt. kinv ) then
           krel    = kinv
           prel    = ps0(krel-1)
           thv0rel = thv0bot(krel)
       else
           krel    = klcl
           prel    = plcl
           thv0rel = thv0lcl
       endif

       ! --------------------------------------------------------------------------- !
       ! Calculate cumulus base mass flux ('cbmf'), fractional area ('ufrcinv'), and !
       ! and mean vertical velocity (winv) of cumulus updraft at PBL top interface.  !
       ! Also, calculate updraft fractional area (ufrclcl) and vertical velocity  at !
       ! the LCL (wlcl). When LCL is below PBLH, cinlcl = 0 and 'ufrclcl = ufrcinv', !
       ! and 'wlcl = winv.                                                           !
       ! Only updrafts strong enough to overcome CIN can rise over PBL top interface.!
       ! Thus,  in order to calculate cumulus mass flux at PBL top interface, 'cbmf',!
       ! we need to know 'CIN' ( the strength of potential energy barrier ) and      !
       ! 'sigmaw' ( a standard deviation of updraft vertical velocity at the PBL top !
       ! interface, a measure of turbulentce strength in the PBL ).   Naturally, the !
       ! ratio of these two variables, 'mu' - normalized CIN by TKE- is key variable !
       ! controlling 'cbmf'.  If 'mu' becomes large, only small fraction of updrafts !
       ! with very strong TKE can rise over the PBL - both 'cbmf' and 'ufrc' becomes !
       ! small, but 'winv' becomes large ( this can be easily understood by PDF of w !
       ! at PBL top ).  If 'mu' becomes small, lots of updraft can rise over the PBL !
       ! top - both 'cbmf' and 'ufrc' becomes large, but 'winv' becomes small. Thus, !
       ! all of the key variables associated with cumulus convection  at the PBL top !
       ! - 'cbmf', 'ufrc', 'winv' where 'cbmf = rho*ufrc*winv' - are a unique functi !
       ! ons of 'mu', normalized CIN. Although these are uniquely determined by 'mu',!
       ! we usually impose two comstraints on 'cbmf' and 'ufrc': (1) because we will !
       ! simply assume that subsidence warming and drying of 'kinv-1' layer in assoc !
       ! iation with 'cbmf' at PBL top interface is confined only in 'kinv-1' layer, !
       ! cbmf must not be larger than the mass within the 'kinv-1' layer. Otherwise, !
       ! instability will occur due to the breaking of stability con. If we consider !
       ! semi-Lagrangian vertical advection scheme and explicitly consider the exten !
       ! t of vertical movement of each layer in association with cumulus mass flux, !
       ! we don't need to impose this constraint. However,  using a  semi-Lagrangian !
       ! scheme is a future research subject. Note that this constraint should be ap !
       ! plied for all interfaces above PBL top as well as PBL top interface.   As a !
       ! result, this 'cbmf' constraint impose a 'lower' limit on mu - 'mumin0'. (2) !
       ! in order for mass flux parameterization - rho*(w'a')= M*(a_c-a_e) - to   be !
       ! valid, cumulus updraft fractional area should be much smaller than 1.    In !
       ! current code, we impose 'rmaxfrac = 0.1 ~ 0.2'   through the whole vertical !
       ! layers where cumulus convection occurs. At the PBL top interface,  the same !
       ! constraint is made by imposing another lower 'lower' limit on mu, 'mumin1'. !
       ! After that, also limit 'ufrclcl' to be smaller than 'rmaxfrac' by 'mumin2'. !
       ! --------------------------------------------------------------------------- !

       ! --------------------------------------------------------------------------- !
       ! Calculate normalized CIN, 'mu' satisfying all the three constraints imposed !
       ! on 'cbmf'('mumin0'), 'ufrc' at the PBL top - 'ufrcinv' - ( by 'mumin1' from !
       ! a parameter sentence), and 'ufrc' at the LCL - 'ufrclcl' ( by 'mumin2').    !
       ! Note that 'cbmf' does not change between PBL top and LCL  because we assume !
       ! that buoyancy sorting does not occur when cumulus updraft is unsaturated.   !
       ! --------------------------------------------------------------------------- !

       if( use_CINcin ) then
           wcrit = sqrt( 2._r8 * cin * rbuoy )
       else
           wcrit = sqrt( 2._r8 * cinlcl * rbuoy )
       endif
       sigmaw = sqrt( rkfre * tkeavg + epsvarw )
       mu = wcrit/sigmaw/1.4142_r8
       if( mu .ge. 3._r8 ) then
         ! write(iulog,*) 'mu >= 3'
           id_exit = .true.
           go to 333
       endif
       rho0inv = ps0(kinv-1)/(r*thv0top(kinv-1)*exns0(kinv-1))
       cbmf = (rho0inv*sigmaw/2.5066_r8)*exp(-mu**2)
       ! 1. 'cbmf' constraint
       cbmflimit = 0.9_r8*dp0(kinv-1)/g/dt
       mumin0 = 0._r8
       if( cbmf .gt. cbmflimit ) mumin0 = sqrt(-log(2.5066_r8*cbmflimit/rho0inv/sigmaw))
       ! 2. 'ufrcinv' constraint
       mu = max(max(mu,mumin0),mumin1)
       ! 3. 'ufrclcl' constraint
       mulcl = sqrt(2._r8*cinlcl*rbuoy)/1.4142_r8/sigmaw
       mulclstar = sqrt(max(0._r8,2._r8*(exp(-mu**2)/2.5066_r8)**2*(1._r8/erfc(mu)**2-0.25_r8/rmaxfrac**2)))
       if( mulcl .gt. 1.e-8_r8 .and. mulcl .gt. mulclstar ) then
           mumin2 = compute_mumin2(mulcl,rmaxfrac,mu)
           if( mu .gt. mumin2 ) then
               write(iulog,*) 'Critical error in mu calculation in UW_ShCu'
               call endrun
           endif
           mu = max(mu,mumin2)
           if( mu .eq. mumin2 ) limit_ufrc(i) = 1._r8
       endif
       if( mu .eq. mumin0 ) limit_cbmf(i) = 1._r8
       if( mu .eq. mumin1 ) limit_ufrc(i) = 1._r8

       ! ------------------------------------------------------------------- !
       ! Calculate final ['cbmf','ufrcinv','winv'] at the PBL top interface. !
       ! Note that final 'cbmf' here is obtained in such that 'ufrcinv' and  !
       ! 'ufrclcl' are smaller than ufrcmax with no instability.             !
       ! ------------------------------------------------------------------- !

       cbmf = (rho0inv*sigmaw/2.5066_r8)*exp(-mu**2)
       winv = sigmaw*(2._r8/2.5066_r8)*exp(-mu**2)/erfc(mu)
       ufrcinv = cbmf/winv/rho0inv

       ! ------------------------------------------------------------------- !
       ! Calculate ['ufrclcl','wlcl'] at the LCL. When LCL is below PBL top, !
       ! it automatically becomes 'ufrclcl = ufrcinv' & 'wlcl = winv', since !
       ! it was already set to 'cinlcl=0' if LCL is below PBL top interface. !
       ! Note 'cbmf' at the PBL top is the same as 'cbmf' at the LCL.  Note  !
       ! also that final 'cbmf' here is obtained in such that 'ufrcinv' and  !
       ! 'ufrclcl' are smaller than ufrcmax and there is no instability.     !
       ! By construction, it must be 'wlcl > 0' but for assurance, I checked !
       ! this again in the below block. If 'ufrclcl < 0.1%', just exit.      !
       ! ------------------------------------------------------------------- !

       wtw = winv * winv - 2._r8 * cinlcl * rbuoy
       if( wtw .le. 0._r8 ) then
         ! write(iulog,*) 'wlcl < 0 at the LCL'
           exit_wtw(i) = 1._r8
           id_exit = .true.
           go to 333
       endif
       wlcl = sqrt(wtw)
       ufrclcl = cbmf/wlcl/rho0inv
       wrel = wlcl
       if( ufrclcl .le. 0.0001_r8 ) then
         ! write(iulog,*) 'ufrclcl <= 0.0001'
           exit_ufrc(i) = 1._r8
           id_exit = .true.
           go to 333
       endif
       ufrc(krel-1) = ufrclcl

       ! ----------------------------------------------------------------------- !
       ! Below is just diagnostic output for detailed analysis of cumulus scheme !
       ! ----------------------------------------------------------------------- !

       ufrcinvbase        = ufrcinv
       winvbase           = winv
       umf(kinv-1:krel-1) = cbmf
       wu(kinv-1:krel-1)  = winv

       ! -------------------------------------------------------------------------- !
       ! Define updraft properties at the level where buoyancy sorting starts to be !
       ! happening, i.e., by definition, at 'prel' level within the release layer.  !
       ! Because no lateral entrainment occurs upto 'prel', conservative scalars of !
       ! cumulus updraft at release level is same as those of source air.  However, !
       ! horizontal momentums of source air are modified by horizontal PGF forcings !
       ! from PBL top interface to 'prel'.  For this case, we should add additional !
       ! horizontal momentum from PBL top interface to 'prel' as will be done below !
       ! to 'usrc' and 'vsrc'. Note that below cumulus updraft properties - umf, wu,!
       ! thlu, qtu, thvu, uu, vu - are defined all interfaces not at the layer mid- !
       ! point. From the index notation of cumulus scheme, wu(k) is the cumulus up- !
       ! draft vertical velocity at the top interface of k layer.                   !
       ! Diabatic horizontal momentum forcing should be treated as a kind of 'body' !
       ! forcing without actual mass exchange between convective updraft and        !
       ! environment, but still taking horizontal momentum from the environment to  !
       ! the convective updrafts. Thus, diabatic convective momentum transport      !
       ! vertically redistributes environmental horizontal momentum.                !
       ! -------------------------------------------------------------------------- !

       emf(krel-1)  = 0._r8
       umf(krel-1)  = cbmf
       wu(krel-1)   = wrel
       thlu(krel-1) = thlsrc
       qtu(krel-1)  = qtsrc
       call conden(prel,thlsrc,qtsrc,thj,qvj,qlj,qij,qse,id_check)
       if( id_check .eq. 1 ) then
           exit_conden(i) = 1._r8
           id_exit = .true.
           go to 333
       endif
       thvu(krel-1) = thj * ( 1._r8 + zvir*qvj - qlj - qij )

       uplus = 0._r8
       vplus = 0._r8
       if( krel .eq. kinv ) then
           uplus = PGFc * ssu0(kinv) * ( prel - ps0(kinv-1) )
           vplus = PGFc * ssv0(kinv) * ( prel - ps0(kinv-1) )
       else
           do k = kinv, max(krel-1,kinv)
              uplus = uplus + PGFc * ssu0(k) * ( ps0(k) - ps0(k-1) )
              vplus = vplus + PGFc * ssv0(k) * ( ps0(k) - ps0(k-1) )
           end do
           uplus = uplus + PGFc * ssu0(krel) * ( prel - ps0(krel-1) )
           vplus = vplus + PGFc * ssv0(krel) * ( prel - ps0(krel-1) )
       end if
       uu(krel-1) = usrc + uplus
       vu(krel-1) = vsrc + vplus

       do m = 1, ncnst
          tru(krel-1,m)  = trsrc(m)
       enddo

       ! -------------------------------------------------------------------------- !
       ! Define environmental properties at the level where buoyancy sorting occurs !
       ! ('pe', normally, layer midpoint except in the 'krel' layer). In the 'krel' !
       ! layer where buoyancy sorting starts to occur, however, 'pe' is defined     !
       ! differently because LCL is regarded as lower interface for mixing purpose. !
       ! -------------------------------------------------------------------------- !

       pe      = 0.5_r8 * ( prel + ps0(krel) )
       dpe     = prel - ps0(krel)
       exne    = exnf(pe)
       thvebot = thv0rel
       thle    = thl0(krel) + ssthl0(krel) * ( pe - p0(krel) )
       qte     = qt0(krel)  + ssqt0(krel)  * ( pe - p0(krel) )
       ue      = u0(krel)   + ssu0(krel)   * ( pe - p0(krel) )
       ve      = v0(krel)   + ssv0(krel)   * ( pe - p0(krel) )
       do m = 1, ncnst
          tre(m) = tr0(krel,m)  + sstr0(krel,m) * ( pe - p0(krel) )
       enddo
!!! end of uwshcu
