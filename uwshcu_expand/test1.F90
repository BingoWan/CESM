!!! start of uwshcu
!!!start count


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
