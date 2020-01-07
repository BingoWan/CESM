  module uwshcu

  use cam_history,    only: outfld, addfld, phys_decomp
  use shr_spfn_mod,   only: erfc => shr_spfn_erfc
  use cam_logfile,    only: iulog
  use ppgrid,         only: pcols, pver, pverp
  use cam_abortutils, only: endrun
  use spmd_utils,     only: masterproc
  use wv_saturation,  only: qsat, qsat_o, qsat_wv, qsat_wv_gam
  use perf_mod, only               : t_startf, t_stopf  ! _EXTERNAL


  implicit none
  private
  save

  public &
     uwshcu_readnl,      &
     init_uwshcu,        &
     compute_uwshcu,     &
     compute_uwshcu_inv

  integer , parameter :: r8 = selected_real_kind(12)    !  8 byte real
  real(r8), parameter :: unset_r8 = huge(1.0_r8)
  real(r8)            :: xlv                            !  Latent heat of vaporization
  real(r8)            :: xlf                            !  Latent heat of fusion
  real(r8)            :: xls                            !  Latent heat of sublimation = xlv + xlf
  real(r8)            :: cp                             !  Specific heat of dry air
  real(r8)            :: zvir                           !  rh2o/rair - 1
  real(r8)            :: r                              !  Gas constant for dry air
  real(r8)            :: g                              !  Gravitational constant
  real(r8)            :: ep2                            !  mol wgt water vapor / mol wgt dry air
  real(r8)            :: p00                            !  Reference pressure for exner function
  real(r8)            :: rovcp                          !  R/cp

  ! Tuning parameters set via namelist
  real(r8) :: rpen          !  For penetrative entrainment efficiency

!===============================================================================
contains
!===============================================================================

  real(r8) function exnf(pressure)
           real(r8), intent(in)              :: pressure
           exnf = (pressure/p00)**rovcp
           return
  end function exnf

!===============================================================================

subroutine uwshcu_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'uwshcu_readnl'

   ! Namelist variables
   real(r8) :: uwshcu_rpen =  unset_r8    !  For penetrative entrainment efficiency

   namelist /uwshcu_nl/ uwshcu_rpen
   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'uwshcu_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, uwshcu_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(uwshcu_rpen,            1, mpir8,  0, mpicom)
#endif

   rpen=uwshcu_rpen


end subroutine uwshcu_readnl

!===============================================================================

  subroutine init_uwshcu( kind, xlv_in, cp_in, xlf_in, zvir_in, r_in, g_in, ep2_in )

    !------------------------------------------------------------- !
    ! Purpose:                                                     !
    ! Initialize key constants for the shallow convection package. !
    !------------------------------------------------------------- !

    use cam_history,   only: outfld, addfld, phys_decomp
    use ppgrid,        only: pcols, pver, pverp
    implicit none
    integer , intent(in) :: kind       !  kind of reals being passed in
    real(r8), intent(in) :: xlv_in     !  Latent heat of vaporization
    real(r8), intent(in) :: xlf_in     !  Latent heat of fusion
    real(r8), intent(in) :: cp_in      !  Specific heat of dry air
    real(r8), intent(in) :: zvir_in    !  rh2o/rair - 1
    real(r8), intent(in) :: r_in       !  Gas constant for dry air
    real(r8), intent(in) :: g_in       !  Gravitational constant
    real(r8), intent(in) :: ep2_in     !  mol wgt water vapor / mol wgt dry air

    character(len=*), parameter :: subname = 'init_uwshcu'

    ! ------------------------- !
    ! Internal Output Variables !
    ! ------------------------- !

    call addfld( 'qtflx_Cu'       , 'kg/m2/s' , pverp , 'A' , 'Convective qt flux'                                  , phys_decomp )
    call addfld( 'slflx_Cu'       , 'J/m2/s'  , pverp , 'A' , 'Convective sl flux'                                  , phys_decomp )
    call addfld( 'uflx_Cu'        , 'kg/m/s2' , pverp , 'A' , 'Convective  u flux'                                  , phys_decomp )
    call addfld( 'vflx_Cu'        , 'kg/m/s2' , pverp , 'A' , 'Convective  v flux'                                  , phys_decomp )

    call addfld( 'qtten_Cu'       , 'kg/kg/s' , pver  , 'A' , 'qt tendency by convection'                           , phys_decomp )
    call addfld( 'slten_Cu'       , 'J/kg/s'  , pver  , 'A' , 'sl tendency by convection'                           , phys_decomp )
    call addfld( 'uten_Cu'        , 'm/s2'    , pver  , 'A' , ' u tendency by convection'                           , phys_decomp )
    call addfld( 'vten_Cu'        , 'm/s2'    , pver  , 'A' , ' v tendency by convection'                           , phys_decomp )
    call addfld( 'qvten_Cu'       , 'kg/kg/s' , pver  , 'A' , 'qv tendency by convection'                           , phys_decomp )
    call addfld( 'qlten_Cu'       , 'kg/kg/s' , pver  , 'A' , 'ql tendency by convection'                           , phys_decomp )
    call addfld( 'qiten_Cu'       , 'kg/kg/s' , pver  , 'A' , 'qi tendency by convection'                           , phys_decomp )

    call addfld( 'cbmf_Cu'        , 'kg/m2/s' , 1     , 'A' , 'Cumulus base mass flux'                              , phys_decomp )
    call addfld( 'ufrcinvbase_Cu' , 'fraction', 1     , 'A' , 'Cumulus fraction at PBL top'                         , phys_decomp )
    call addfld( 'ufrclcl_Cu'     , 'fraction', 1     , 'A' , 'Cumulus fraction at LCL'                             , phys_decomp )
    call addfld( 'winvbase_Cu'    , 'm/s'     , 1     , 'A' , 'Cumulus vertical velocity at PBL top'                , phys_decomp )
    call addfld( 'wlcl_Cu'        , 'm/s'     , 1     , 'A' , 'Cumulus vertical velocity at LCL'                    , phys_decomp )
    call addfld( 'plcl_Cu'        , 'Pa'      , 1     , 'A' , 'LCL of source air'                                   , phys_decomp )
    call addfld( 'pinv_Cu'        , 'Pa'      , 1     , 'A' , 'PBL top pressure'                                    , phys_decomp )
    call addfld( 'plfc_Cu'        , 'Pa'      , 1     , 'A' , 'LFC of source air'                                   , phys_decomp )
    call addfld( 'pbup_Cu'        , 'Pa'      , 1     , 'A' , 'Highest interface level of positive cumulus buoyancy', phys_decomp )
    call addfld( 'ppen_Cu'        , 'Pa'      , 1     , 'A' , 'Highest level where cumulus w is 0'                  , phys_decomp )
    call addfld( 'qtsrc_Cu'       , 'kg/kg'   , 1     , 'A' , 'Cumulus source air qt'                               , phys_decomp )
    call addfld( 'thlsrc_Cu'      , 'K'       , 1     , 'A' , 'Cumulus source air thl'                              , phys_decomp )
    call addfld( 'thvlsrc_Cu'     , 'K'       , 1     , 'A' , 'Cumulus source air thvl'                             , phys_decomp )
    call addfld( 'emfkbup_Cu'     , 'kg/m2/s' , 1     , 'A' , 'Penetrative mass flux at kbup'                       , phys_decomp )
    call addfld( 'cin_Cu'         , 'J/kg'    , 1     , 'A' , 'CIN upto LFC'                                        , phys_decomp )
    call addfld( 'cinlcl_Cu'      , 'J/kg'    , 1     , 'A' , 'CIN upto LCL'                                        , phys_decomp )
    call addfld( 'cbmflimit_Cu'   , 'kg/m2/s' , 1     , 'A' , 'cbmflimiter'                                         , phys_decomp )
    call addfld( 'tkeavg_Cu'      , 'm2/s2'   , 1     , 'A' , 'Average tke within PBL for convection scheme'        , phys_decomp )
    call addfld( 'zinv_Cu'        , 'm'       , 1     , 'A' , 'PBL top height'                                      , phys_decomp )
    call addfld( 'rcwp_Cu'        , 'kg/m2'   , 1     , 'A' , 'Cumulus LWP+IWP'                                     , phys_decomp )
    call addfld( 'rlwp_Cu'        , 'kg/m2'   , 1     , 'A' , 'Cumulus LWP'                                         , phys_decomp )
    call addfld( 'riwp_Cu'        , 'kg/m2'   , 1     , 'A' , 'Cumulus IWP'                                         , phys_decomp )
    call addfld( 'tophgt_Cu'      , 'm'       , 1     , 'A' , 'Cumulus top height'                                  , phys_decomp )

    call addfld( 'wu_Cu'          , 'm/s'     , pverp , 'A' , 'Convective updraft vertical velocity'                , phys_decomp )
    call addfld( 'ufrc_Cu'        , 'fraction', pverp , 'A' , 'Convective updraft fractional area'                  , phys_decomp )
    call addfld( 'qtu_Cu'         , 'kg/kg'   , pverp , 'A' , 'Cumulus updraft qt'                                  , phys_decomp )
    call addfld( 'thlu_Cu'        , 'K'       , pverp , 'A' , 'Cumulus updraft thl'                                 , phys_decomp )
    call addfld( 'thvu_Cu'        , 'K'       , pverp , 'A' , 'Cumulus updraft thv'                                 , phys_decomp )
    call addfld( 'uu_Cu'          , 'm/s'     , pverp , 'A' , 'Cumulus updraft uwnd'                                , phys_decomp )
    call addfld( 'vu_Cu'          , 'm/s'     , pverp , 'A' , 'Cumulus updraft vwnd'                                , phys_decomp )
    call addfld( 'qtu_emf_Cu'     , 'kg/kg'   , pverp , 'A' , 'qt of penatratively entrained air'                   , phys_decomp )
    call addfld( 'thlu_emf_Cu'    , 'K'       , pverp , 'A' , 'thl of penatratively entrained air'                  , phys_decomp )
    call addfld( 'uu_emf_Cu'      , 'm/s'     , pverp , 'A' , 'uwnd of penatratively entrained air'                 , phys_decomp )
    call addfld( 'vu_emf_Cu'      , 'm/s'     , pverp , 'A' , 'vwnd of penatratively entrained air'                 , phys_decomp )
    call addfld( 'umf_Cu'         , 'kg/m2/s' , pverp , 'A' , 'Cumulus updraft mass flux'                           , phys_decomp )
    call addfld( 'uemf_Cu'        , 'kg/m2/s' , pverp , 'A' , 'Cumulus net ( updraft + entrainment ) mass flux'     , phys_decomp )
    call addfld( 'qcu_Cu'         , 'kg/kg'   , pver  , 'A' , 'Cumulus updraft LWC+IWC'                             , phys_decomp )
    call addfld( 'qlu_Cu'         , 'kg/kg'   , pver  , 'A' , 'Cumulus updraft LWC'                                 , phys_decomp )
    call addfld( 'qiu_Cu'         , 'kg/kg'   , pver  , 'A' , 'Cumulus updraft IWC'                                 , phys_decomp )
    call addfld( 'cufrc_Cu'       , 'fraction', pver  , 'A' , 'Cumulus cloud fraction'                              , phys_decomp )
    call addfld( 'fer_Cu'         , '1/m'     , pver  , 'A' , 'Cumulus lateral fractional entrainment rate'         , phys_decomp )
    call addfld( 'fdr_Cu'         , '1/m'     , pver  , 'A' , 'Cumulus lateral fractional detrainment Rate'         , phys_decomp )

    call addfld( 'dwten_Cu'       , 'kg/kg/s' , pver  , 'A' , 'Expellsion rate of cumulus cloud water to env.'      , phys_decomp )
    call addfld( 'diten_Cu'       , 'kg/kg/s' , pver  , 'A' , 'Expellsion rate of cumulus ice water to env.'        , phys_decomp )
    call addfld( 'qrten_Cu'       , 'kg/kg/s' , pver  , 'A' , 'Production rate of rain by cumulus'                  , phys_decomp )
    call addfld( 'qsten_Cu'       , 'kg/kg/s' , pver  , 'A' , 'Production rate of snow by cumulus'                  , phys_decomp )
    call addfld( 'flxrain_Cu'     , 'kg/m2/s' , pverp , 'A' , 'Rain flux induced by Cumulus'                        , phys_decomp )
    call addfld( 'flxsnow_Cu'     , 'kg/m2/s' , pverp , 'A' , 'Snow flux induced by Cumulus'                        , phys_decomp )
    call addfld( 'ntraprd_Cu'     , 'kg/kg/s' , pver  , 'A' , 'Net production rate of rain by Cumulus'              , phys_decomp )
    call addfld( 'ntsnprd_Cu'     , 'kg/kg/s' , pver  , 'A' , 'Net production rate of snow by Cumulus'              , phys_decomp )

    call addfld( 'excessu_Cu'     , 'no'      , pver  , 'A' , 'Updraft saturation excess'                           , phys_decomp )
    call addfld( 'excess0_Cu'     , 'no'      , pver  , 'A' , 'Environmental saturation excess'                     , phys_decomp )
    call addfld( 'xc_Cu'          , 'no'      , pver  , 'A' , 'Critical mixing ratio'                               , phys_decomp )
    call addfld( 'aquad_Cu'       , 'no'      , pver  , 'A' , 'aquad'                                               , phys_decomp )
    call addfld( 'bquad_Cu'       , 'no'      , pver  , 'A' , 'bquad'                                               , phys_decomp )
    call addfld( 'cquad_Cu'       , 'no'      , pver  , 'A' , 'cquad'                                               , phys_decomp )
    call addfld( 'bogbot_Cu'      , 'no'      , pver  , 'A' , 'Cloud buoyancy at the bottom interface'              , phys_decomp )
    call addfld( 'bogtop_Cu'      , 'no'      , pver  , 'A' , 'Cloud buoyancy at the top interface'                 , phys_decomp )

    call addfld('exit_UWCu_Cu'    , 'no'      , 1     , 'A' , 'exit_UWCu'                                           , phys_decomp )
    call addfld('exit_conden_Cu'  , 'no'      , 1     , 'A' , 'exit_conden'                                         , phys_decomp )
    call addfld('exit_klclmkx_Cu' , 'no'      , 1     , 'A' , 'exit_klclmkx'                                        , phys_decomp )
    call addfld('exit_klfcmkx_Cu' , 'no'      , 1     , 'A' , 'exit_klfcmkx'                                        , phys_decomp )
    call addfld('exit_ufrc_Cu'    , 'no'      , 1     , 'A' , 'exit_ufrc'                                           , phys_decomp )
    call addfld('exit_wtw_Cu'     , 'no'      , 1     , 'A' , 'exit_wtw'                                            , phys_decomp )
    call addfld('exit_drycore_Cu' , 'no'      , 1     , 'A' , 'exit_drycore'                                        , phys_decomp )
    call addfld('exit_wu_Cu'      , 'no'      , 1     , 'A' , 'exit_wu'                                             , phys_decomp )
    call addfld('exit_cufilter_Cu', 'no'      , 1     , 'A' , 'exit_cufilter'                                       , phys_decomp )
    call addfld('exit_kinv1_Cu'   , 'no'      , 1     , 'A' , 'exit_kinv1'                                          , phys_decomp )
    call addfld('exit_rei_Cu'     , 'no'      , 1     , 'A' , 'exit_rei'                                            , phys_decomp )

    call addfld('limit_shcu_Cu'   , 'no'      , 1     , 'A' , 'limit_shcu'                                          , phys_decomp )
    call addfld('limit_negcon_Cu' , 'no'      , 1     , 'A' , 'limit_negcon'                                        , phys_decomp )
    call addfld('limit_ufrc_Cu'   , 'no'      , 1     , 'A' , 'limit_ufrc'                                          , phys_decomp )
    call addfld('limit_ppen_Cu'   , 'no'      , 1     , 'A' , 'limit_ppen'                                          , phys_decomp )
    call addfld('limit_emf_Cu'    , 'no'      , 1     , 'A' , 'limit_emf'                                           , phys_decomp )
    call addfld('limit_cinlcl_Cu' , 'no'      , 1     , 'A' , 'limit_cinlcl'                                        , phys_decomp )
    call addfld('limit_cin_Cu'    , 'no'      , 1     , 'A' , 'limit_cin'                                           , phys_decomp )
    call addfld('limit_cbmf_Cu'   , 'no'      , 1     , 'A' , 'limit_cbmf'                                          , phys_decomp )
    call addfld('limit_rei_Cu'    , 'no'      , 1     , 'A' , 'limit_rei'                                           , phys_decomp )
    call addfld('ind_delcin_Cu'   , 'no'      , 1     , 'A' , 'ind_delcin'                                          , phys_decomp )

    if( kind .ne. r8 ) then
        write(iulog,*) subname//': ERROR -- real KIND does not match internal specification.'
        call endrun(subname//': ERROR -- real KIND does not match internal specification.')
    endif

    xlv   = xlv_in
    xlf   = xlf_in
    xls   = xlv + xlf
    cp    = cp_in
    zvir  = zvir_in
    r     = r_in
    g     = g_in
    ep2   = ep2_in
    p00   = 1.e5_r8
    rovcp = r/cp

    if (rpen == unset_r8) then
       call endrun(subname//': uwshcu_rpen must be set in the namelist')
    end if

    if ( masterproc ) then
       write(iulog,*) subname//': tuning parameters: rpen=',rpen
    endif

  end subroutine init_uwshcu

  subroutine compute_uwshcu_inv( mix      , mkx        , iend          , ncnst     , dt       ,  &
                                 ps0_inv  , zs0_inv    , p0_inv        , z0_inv    , dp0_inv  ,  &
                                 u0_inv   , v0_inv     , qv0_inv       , ql0_inv   , qi0_inv  ,  &
                                 t0_inv   , s0_inv     , tr0_inv       ,                         &
                                 tke_inv  , cldfrct_inv, concldfrct_inv, pblh      , cush     ,  &
                                 umf_inv  , slflx_inv  , qtflx_inv     ,                         &
                                 flxprc1_inv, flxsnow1_inv,     				 &
                                 qvten_inv, qlten_inv  , qiten_inv     ,                         &
                                 sten_inv , uten_inv   , vten_inv      , trten_inv ,             &
                                 qrten_inv, qsten_inv  , precip        , snow      , evapc_inv,  &
                                 cufrc_inv, qcu_inv    , qlu_inv       , qiu_inv   ,             &
                                 cbmf     , qc_inv     , rliq          ,                         &
                                 cnt_inv  , cnb_inv    , lchnk         , dpdry0_inv )

    implicit none
    integer , intent(in)    :: lchnk
    integer , intent(in)    :: mix
    integer , intent(in)    :: mkx
    integer , intent(in)    :: iend
    integer , intent(in)    :: ncnst
    real(r8), intent(in)    :: dt                       !  Time step : 2*delta_t [ s ]
    real(r8), intent(in)    :: ps0_inv(mix,mkx+1)       !  Environmental pressure at the interfaces [ Pa ]
    real(r8), intent(in)    :: zs0_inv(mix,mkx+1)       !  Environmental height at the interfaces   [ m ]
    real(r8), intent(in)    :: p0_inv(mix,mkx)          !  Environmental pressure at the layer mid-point [ Pa ]
    real(r8), intent(in)    :: z0_inv(mix,mkx)          !  Environmental height at the layer mid-point [ m ]
    real(r8), intent(in)    :: dp0_inv(mix,mkx)         !  Environmental layer pressure thickness [ Pa ] > 0.
    real(r8), intent(in)    :: dpdry0_inv(mix,mkx)      !  Environmental dry layer pressure thickness [ Pa ]
    real(r8), intent(in)    :: u0_inv(mix,mkx)          !  Environmental zonal wind [ m/s ]
    real(r8), intent(in)    :: v0_inv(mix,mkx)          !  Environmental meridional wind [ m/s ]
    real(r8), intent(in)    :: qv0_inv(mix,mkx)         !  Environmental water vapor specific humidity [ kg/kg ]
    real(r8), intent(in)    :: ql0_inv(mix,mkx)         !  Environmental liquid water specific humidity [ kg/kg ]
    real(r8), intent(in)    :: qi0_inv(mix,mkx)         !  Environmental ice specific humidity [ kg/kg ]
    real(r8), intent(in)    :: t0_inv(mix,mkx)          !  Environmental temperature [ K ]
    real(r8), intent(in)    :: s0_inv(mix,mkx)          !  Environmental dry static energy [ J/kg ]
    real(r8), intent(in)    :: tr0_inv(mix,mkx,ncnst)   !  Environmental tracers [ #, kg/kg ]
    real(r8), intent(in)    :: tke_inv(mix,mkx+1)       !  Turbulent kinetic energy at the interfaces [ m2/s2 ]
    real(r8), intent(in)    :: cldfrct_inv(mix,mkx)     !  Total cloud fraction at the previous time step [ fraction ]
    real(r8), intent(in)    :: concldfrct_inv(mix,mkx)  !  Total convective ( shallow + deep ) cloud fraction
                                                        !  at the previous time step [ fraction ]
    real(r8), intent(in)    :: pblh(mix)                !  Height of PBL [ m ]
    real(r8), intent(inout) :: cush(mix)                !  Convective scale height [ m ]
    real(r8), intent(out)   :: umf_inv(mix,mkx+1)       !  Updraft mass flux at the interfaces [ kg/m2/s ]
    real(r8), intent(out)   :: qvten_inv(mix,mkx)       !  Tendency of water vapor specific humidity [ kg/kg/s ]
    real(r8), intent(out)   :: qlten_inv(mix,mkx)       !  Tendency of liquid water specific humidity [ kg/kg/s ]
    real(r8), intent(out)   :: qiten_inv(mix,mkx)       !  Tendency of ice specific humidity [ kg/kg/s ]
    real(r8), intent(out)   :: sten_inv(mix,mkx)        !  Tendency of dry static energy [ J/kg/s ]
    real(r8), intent(out)   :: uten_inv(mix,mkx)        !  Tendency of zonal wind [ m/s2 ]
    real(r8), intent(out)   :: vten_inv(mix,mkx)        !  Tendency of meridional wind [ m/s2 ]
    real(r8), intent(out)   :: trten_inv(mix,mkx,ncnst) !  Tendency of tracers [ #/s, kg/kg/s ]
    real(r8), intent(out)   :: qrten_inv(mix,mkx)       !  Tendency of rain water specific humidity [ kg/kg/s ]
    real(r8), intent(out)   :: qsten_inv(mix,mkx)       !  Tendency of snow specific humidity [ kg/kg/s ]
    real(r8), intent(out)   :: precip(mix)              !  Precipitation ( rain + snow ) flux at the surface [ m/s ]
    real(r8), intent(out)   :: snow(mix)                !  Snow flux at the surface [ m/s ]
    real(r8), intent(out)   :: evapc_inv(mix,mkx)       !  Evaporation of precipitation [ kg/kg/s ]
    real(r8), intent(out)   :: rliq(mix)                !  Vertical integral of tendency of detrained cloud condensate qc [ m/s ]
    real(r8), intent(out)   :: slflx_inv(mix,mkx+1)     !  Updraft liquid static energy flux [ J/kg * kg/m2/s ]
    real(r8), intent(out)   :: qtflx_inv(mix,mkx+1)     !  Updraft total water flux [ kg/kg * kg/m2/s ]
    real(r8), intent(out)   :: flxprc1_inv(mix,mkx+1)   ! uw grid-box mean rain+snow flux (kg m^-2 s^-1)
                                                        ! for physics buffer calls in convect_shallow.F90
    real(r8), intent(out)   :: flxsnow1_inv(mix,mkx+1)  ! uw grid-box mean snow flux (kg m^-2 s^-1)
                                                        ! for physics buffer calls in convect_shallow.F90

    real(r8), intent(out)   :: cufrc_inv(mix,mkx)       !  Shallow cumulus cloud fraction at the layer mid-point [ fraction ]
    real(r8), intent(out)   :: qcu_inv(mix,mkx)         !  Liquid+ice specific humidity within cumulus updraft [ kg/kg ]
    real(r8), intent(out)   :: qlu_inv(mix,mkx)         !  Liquid water specific humidity within cumulus updraft [ kg/kg ]
    real(r8), intent(out)   :: qiu_inv(mix,mkx)         !  Ice specific humidity within cumulus updraft [ kg/kg ]
    real(r8), intent(out)   :: qc_inv(mix,mkx)          !  Tendency of cumulus condensate detrained into the environment [ kg/kg/s ]
    real(r8), intent(out)   :: cbmf(mix)                !  Cumulus base mass flux [ kg/m2/s ]
    real(r8), intent(out)   :: cnt_inv(mix)             !  Cumulus top  interface index, cnt = kpen [ no ]
    real(r8), intent(out)   :: cnb_inv(mix)             !  Cumulus base interface index, cnb = krel - 1 [ no ]

    real(r8)                :: ps0(mix,0:mkx)           !  Environmental pressure at the interfaces [ Pa ]
    real(r8)                :: zs0(mix,0:mkx)           !  Environmental height at the interfaces   [ m ]
    real(r8)                :: p0(mix,mkx)              !  Environmental pressure at the layer mid-point [ Pa ]
    real(r8)                :: z0(mix,mkx)              !  Environmental height at the layer mid-point [ m ]
    real(r8)                :: dp0(mix,mkx)             !  Environmental layer pressure thickness [ Pa ] > 0.
    real(r8)                :: dpdry0(mix,mkx)          !  Environmental dry layer pressure thickness [ Pa ]
    real(r8)                :: u0(mix,mkx)              !  Environmental zonal wind [ m/s ]
    real(r8)                :: v0(mix,mkx)              !  Environmental meridional wind [ m/s ]
    real(r8)                :: tke(mix,0:mkx)           !  Turbulent kinetic energy at the interfaces [ m2/s2 ]
    real(r8)                :: cldfrct(mix,mkx)         !  Total cloud fraction at the previous time step [ fraction ]
    real(r8)                :: concldfrct(mix,mkx)      !  Total convective ( shallow + deep ) cloud fraction
                                                        ! at the previous time step [ fraction ]
    real(r8)                :: qv0(mix,mkx)             !  Environmental water vapor specific humidity [ kg/kg ]
    real(r8)                :: ql0(mix,mkx)             !  Environmental liquid water specific humidity [ kg/kg ]
    real(r8)                :: qi0(mix,mkx)             !  Environmental ice specific humidity [ kg/kg ]
    real(r8)                :: t0(mix,mkx)              !  Environmental temperature [ K ]
    real(r8)                :: s0(mix,mkx)              !  Environmental dry static energy [ J/kg ]
    real(r8)                :: tr0(mix,mkx,ncnst)       !  Environmental tracers [ #, kg/kg ]
    real(r8)                :: umf(mix,0:mkx)           !  Updraft mass flux at the interfaces [ kg/m2/s ]
    real(r8)                :: qvten(mix,mkx)           !  Tendency of water vapor specific humidity [ kg/kg/s ]
    real(r8)                :: qlten(mix,mkx)           !  Tendency of liquid water specific humidity [ kg/kg/s ]
    real(r8)                :: qiten(mix,mkx)           !  tendency of ice specific humidity [ kg/kg/s ]
    real(r8)                :: sten(mix,mkx)            !  Tendency of static energy [ J/kg/s ]
    real(r8)                :: uten(mix,mkx)            !  Tendency of zonal wind [ m/s2 ]
    real(r8)                :: vten(mix,mkx)            !  Tendency of meridional wind [ m/s2 ]
    real(r8)                :: trten(mix,mkx,ncnst)     !  Tendency of tracers [ #/s, kg/kg/s ]
    real(r8)                :: qrten(mix,mkx)           !  Tendency of rain water specific humidity [ kg/kg/s ]
    real(r8)                :: qsten(mix,mkx)           !  Tendency of snow speficif humidity [ kg/kg/s ]
    real(r8)                :: evapc(mix,mkx)           !  Tendency of evaporation of precipitation [ kg/kg/s ]
    real(r8)                :: slflx(mix,0:mkx)         !  Updraft liquid static energy flux [ J/kg * kg/m2/s ]
    real(r8)                :: qtflx(mix,0:mkx)         !  Updraft total water flux [ kg/kg * kg/m2/s ]
    real(r8)                :: flxprc1(mix,0:mkx)       ! uw grid-box mean rain+snow flux (kg m^-2 s^-1)
                                                        ! for physics buffer calls in convect_shallow.F90
    real(r8)                :: flxsnow1(mix,0:mkx)      ! uw grid-box mean snow flux (kg m^-2 s^-1)
                                                        ! for physics buffer calls in convect_shallow.F90
    real(r8)                :: cufrc(mix,mkx)           !  Shallow cumulus cloud fraction at the layer mid-point [ fraction ]
    real(r8)                :: qcu(mix,mkx)             !  Condensate water specific humidity within cumulus updraft
                                                        ! at the layer mid-point [ kg/kg ]
    real(r8)                :: qlu(mix,mkx)             !  Liquid water specific humidity within cumulus updraft
                                                        ! at the layer mid-point [ kg/kg ]
    real(r8)                :: qiu(mix,mkx)             !  Ice specific humidity within cumulus updraft
                                                        ! at the layer mid-point [ kg/kg ]
    real(r8)                :: qc(mix,mkx)              !  Tendency of cumulus condensate detrained into the environment [ kg/kg/s ]
    real(r8)                :: cnt(mix)                 !  Cumulus top  interface index, cnt = kpen [ no ]
    real(r8)                :: cnb(mix)                 !  Cumulus base interface index, cnb = krel - 1 [ no ]
    integer                 :: k                        !  Vertical index for local fields [ no ]
    integer                 :: k_inv                    !  Vertical index for incoming fields [ no ]
    integer                 :: m                        !  Tracer index [ no ]

    do k = 1, mkx
       k_inv               = mkx + 1 - k
       p0(:iend,k)         = p0_inv(:iend,k_inv)
       u0(:iend,k)         = u0_inv(:iend,k_inv)
       v0(:iend,k)         = v0_inv(:iend,k_inv)
       z0(:iend,k)         = z0_inv(:iend,k_inv)
       dp0(:iend,k)        = dp0_inv(:iend,k_inv)
       dpdry0(:iend,k)     = dpdry0_inv(:iend,k_inv)
       qv0(:iend,k)        = qv0_inv(:iend,k_inv)
       ql0(:iend,k)        = ql0_inv(:iend,k_inv)
       qi0(:iend,k)        = qi0_inv(:iend,k_inv)
       t0(:iend,k)         = t0_inv(:iend,k_inv)
       s0(:iend,k)         = s0_inv(:iend,k_inv)
       cldfrct(:iend,k)    = cldfrct_inv(:iend,k_inv)
       concldfrct(:iend,k) = concldfrct_inv(:iend,k_inv)
       do m = 1, ncnst
          tr0(:iend,k,m)   = tr0_inv(:iend,k_inv,m)
       enddo
    enddo

    do k = 0, mkx
       k_inv               = mkx + 1 - k
       ps0(:iend,k)        = ps0_inv(:iend,k_inv)
       zs0(:iend,k)        = zs0_inv(:iend,k_inv)
       tke(:iend,k)        = tke_inv(:iend,k_inv)
    end do

    call compute_uwshcu( mix  , mkx    , iend      , ncnst , dt   , &
                         ps0  , zs0    , p0        , z0    , dp0  , &
                         u0   , v0     , qv0       , ql0   , qi0  , &
                         t0   , s0     , tr0       ,                &
                         tke  , cldfrct, concldfrct, pblh  , cush , &
                         umf  , slflx  , qtflx     ,                &
                         flxprc1  , flxsnow1  ,		            &
                         qvten, qlten  , qiten     ,                &
                         sten , uten   , vten      , trten ,        &
                         qrten, qsten  , precip    , snow  , evapc, &
                         cufrc, qcu    , qlu       , qiu   ,        &
                         cbmf , qc     , rliq      ,                &
                         cnt  , cnb    , lchnk     , dpdry0 )

    ! Reverse cloud top/base interface indices

       cnt_inv(:iend) = mkx + 1 - cnt(:iend)
       cnb_inv(:iend) = mkx + 1 - cnb(:iend)

    do k = 0, mkx
       k_inv                  = mkx + 1 - k
       umf_inv(:iend,k_inv)   = umf(:iend,k)
       slflx_inv(:iend,k_inv) = slflx(:iend,k)
       qtflx_inv(:iend,k_inv) = qtflx(:iend,k)
       flxprc1_inv(:iend,k_inv) = flxprc1(:iend,k)     ! reversed for output to cam
       flxsnow1_inv(:iend,k_inv) = flxsnow1(:iend,k)   ! ""
    end do

    do k = 1, mkx
       k_inv                         = mkx + 1 - k
       qvten_inv(:iend,k_inv)        = qvten(:iend,k)
       qlten_inv(:iend,k_inv)        = qlten(:iend,k)
       qiten_inv(:iend,k_inv)        = qiten(:iend,k)
       sten_inv(:iend,k_inv)         = sten(:iend,k)
       uten_inv(:iend,k_inv)         = uten(:iend,k)
       vten_inv(:iend,k_inv)         = vten(:iend,k)
       qrten_inv(:iend,k_inv)        = qrten(:iend,k)
       qsten_inv(:iend,k_inv)        = qsten(:iend,k)
       evapc_inv(:iend,k_inv)        = evapc(:iend,k)
       cufrc_inv(:iend,k_inv)        = cufrc(:iend,k)
       qcu_inv(:iend,k_inv)          = qcu(:iend,k)
       qlu_inv(:iend,k_inv)          = qlu(:iend,k)
       qiu_inv(:iend,k_inv)          = qiu(:iend,k)
       qc_inv(:iend,k_inv)           = qc(:iend,k)
       do m = 1, ncnst
          trten_inv(:iend,k_inv,m)   = trten(:iend,k,m)
       enddo

    enddo

  end subroutine compute_uwshcu_inv

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
    integer , intent(in)    :: mix
    integer , intent(in)    :: mkx
    integer , intent(in)    :: iend
    integer , intent(in)    :: ncnst
    real(r8), intent(in)    :: dt                             !  Time step : 2*delta_t [ s ]
    real(r8), intent(in)    :: ps0_in(mix,0:mkx)              !  Environmental pressure at the interfaces [ Pa ]
    real(r8), intent(in)    :: zs0_in(mix,0:mkx)              !  Environmental height at the interfaces [ m ]
    real(r8), intent(in)    :: p0_in(mix,mkx)                 !  Environmental pressure at the layer mid-point [ Pa ]
    real(r8), intent(in)    :: z0_in(mix,mkx)                 !  Environmental height at the layer mid-point [ m ]
    real(r8), intent(in)    :: dp0_in(mix,mkx)                !  Environmental layer pressure thickness [ Pa ] > 0.
    real(r8), intent(in)    :: dpdry0_in(mix,mkx)             !  Environmental dry layer pressure thickness [ Pa ]
    real(r8), intent(in)    :: u0_in(mix,mkx)                 !  Environmental zonal wind [ m/s ]
    real(r8), intent(in)    :: v0_in(mix,mkx)                 !  Environmental meridional wind [ m/s ]
    real(r8), intent(in)    :: qv0_in(mix,mkx)                !  Environmental water vapor specific humidity [ kg/kg ]
    real(r8), intent(in)    :: ql0_in(mix,mkx)                !  Environmental liquid water specific humidity [ kg/kg ]
    real(r8), intent(in)    :: qi0_in(mix,mkx)                !  Environmental ice specific humidity [ kg/kg ]
    real(r8), intent(in)    :: t0_in(mix,mkx)                 !  Environmental temperature [ K ]
    real(r8), intent(in)    :: s0_in(mix,mkx)                 !  Environmental dry static energy [ J/kg ]
    real(r8), intent(in)    :: tr0_in(mix,mkx,ncnst)          !  Environmental tracers [ #, kg/kg ]
    real(r8), intent(in)    :: tke_in(mix,0:mkx)              !  Turbulent kinetic energy at the interfaces [ m2/s2 ]
    real(r8), intent(in)    :: cldfrct_in(mix,mkx)            !  Total cloud fraction at the previous time step [ fraction ]
    real(r8), intent(in)    :: concldfrct_in(mix,mkx)         !  Total convective cloud fraction
                                                              ! at the previous time step [ fraction ]
    real(r8), intent(in)    :: pblh_in(mix)                   !  Height of PBL [ m ]
    real(r8), intent(inout) :: cush_inout(mix)                !  Convective scale height [ m ]

    real(r8)                   tw0_in(mix,mkx)                !  Wet bulb temperature [ K ]
    real(r8)                   qw0_in(mix,mkx)                !  Wet-bulb specific humidity [ kg/kg ]

    real(r8), intent(out)   :: umf_out(mix,0:mkx)             !  Updraft mass flux at the interfaces [ kg/m2/s ]
    real(r8), intent(out)   :: qvten_out(mix,mkx)             !  Tendency of water vapor specific humidity [ kg/kg/s ]
    real(r8), intent(out)   :: qlten_out(mix,mkx)             !  Tendency of liquid water specific humidity [ kg/kg/s ]
    real(r8), intent(out)   :: qiten_out(mix,mkx)             !  Tendency of ice specific humidity [ kg/kg/s ]
    real(r8), intent(out)   :: sten_out(mix,mkx)              !  Tendency of dry static energy [ J/kg/s ]
    real(r8), intent(out)   :: uten_out(mix,mkx)              !  Tendency of zonal wind [ m/s2 ]
    real(r8), intent(out)   :: vten_out(mix,mkx)              !  Tendency of meridional wind [ m/s2 ]
    real(r8), intent(out)   :: trten_out(mix,mkx,ncnst)       !  Tendency of tracers [ #/s, kg/kg/s ]
    real(r8), intent(out)   :: qrten_out(mix,mkx)             !  Tendency of rain water specific humidity [ kg/kg/s ]
    real(r8), intent(out)   :: qsten_out(mix,mkx)             !  Tendency of snow(mix) specific humidity [ kg/kg/s ]
    real(r8), intent(out)   :: precip_out(mix)                !  Precipitation ( rain + snow(mix) ) rate at surface [ m/s ]
    real(r8), intent(out)   :: snow_out(mix)                  !  Snow rate at surface [ m/s ]
    real(r8), intent(out)   :: evapc_out(mix,mkx)             !  Tendency of evaporation of precipitation [ kg/kg/s ]
    real(r8), intent(out)   :: slflx_out(mix,0:mkx)           !  Updraft/pen.entrainment liquid static energy flux
                                                              ! [ J/kg * kg/m2/s ]
    real(r8), intent(out)   :: qtflx_out(mix,0:mkx)           !  updraft/pen.entrainment total water flux [ kg/kg * kg/m2/s ]
    real(r8), intent(out)   :: flxprc1_out(mix,0:mkx)         ! precip(mix) (rain+snow(mix)) flux
    real(r8), intent(out)   :: flxsnow1_out(mix,0:mkx)        ! snow(mix) flux
    real(r8), intent(out)   :: cufrc_out(mix,mkx)             !  Shallow cumulus cloud fraction at the layer mid-point [ fraction ]
    real(r8), intent(out)   :: qcu_out(mix,mkx)               !  Condensate water specific humidity within cumulus updraft [ kg/kg ]
    real(r8), intent(out)   :: qlu_out(mix,mkx)               !  Liquid water specific humidity within cumulus updraft [ kg/kg ]
    real(r8), intent(out)   :: qiu_out(mix,mkx)               !  Ice specific humidity within cumulus updraft [ kg/kg ]
    real(r8), intent(out)   :: cbmf_out(mix)                  !  Cloud base mass flux [ kg/m2/s ]
    real(r8), intent(out)   :: qc_out(mix,mkx)                !  Tendency of detrained cumulus condensate
                                                              ! into the environment [ kg/kg/s ]
    real(r8), intent(out)   :: rliq_out(mix)                  !  Vertical integral of qc_out [ m/s ]
    real(r8), intent(out)   :: cnt_out(mix)                   !  Cumulus top  interface index, cnt(mix) = kpen(mix) [ no ]
    real(r8), intent(out)   :: cnb_out(mix)                   !  Cumulus base interface index, cnb(mix) = krel(mix) - 1 [ no ]

    !
    ! Internal Output Variables
    !

    real(r8)                   qtten_out(mix,mkx)             !  Tendency of qt [ kg/kg/s ]
    real(r8)                   slten_out(mix,mkx)             !  Tendency of sl [ J/kg/s ]
    real(r8)                   ufrc_out(mix,0:mkx)            !  Updraft fractional area at the interfaces [ fraction ]
    real(r8)                   uflx_out(mix,0:mkx)            !  Updraft/pen.entrainment zonal momentum flux [ m/s/m2/s ]
    real(r8)                   vflx_out(mix,0:mkx)            !  Updraft/pen.entrainment meridional momentum flux [ m/s/m2/s ]
    real(r8)                   fer_out(mix,mkx)               !  Fractional lateral entrainment rate [ 1/Pa ]
    real(r8)                   fdr_out(mix,mkx)               !  Fractional lateral detrainment rate [ 1/Pa ]
    real(r8)                   cinh_out(mix)                  !  Convective INhibition upto LFC (CIN) [ J/kg ]
    real(r8)                   trflx_out(mix,0:mkx,ncnst)     !  Updraft/pen.entrainment tracer flux [ #/m2/s, kg/kg/m2/s ]

    ! -------------------------------------------- !
    ! One-dimensional variables at each grid point !
    ! -------------------------------------------- !

    ! 1. Input variables

    real(r8)    ps0(0:mkx, mix)                                    !  Environmental pressure at the interfaces [ Pa ]
    real(r8)    zs0(0:mkx, mix)                                    !  Environmental height at the interfaces [ m ]
    real(r8)    p0(mkx, mix)                                       !  Environmental pressure at the layer mid-point [ Pa ]
    real(r8)    z0(mkx, mix)                                       !  Environmental height at the layer mid-point [ m ]
    real(r8)    dp0(mkx, mix)                                      !  Environmental layer pressure thickness [ Pa ] > 0.
    real(r8)    dpdry0(mkx, mix)                                   !  Environmental dry layer pressure thickness [ Pa ]
    real(r8)    u0(mkx, mix)                                       !  Environmental zonal wind [ m/s ]
    real(r8)    v0(mkx, mix)                                       !  Environmental meridional wind [ m/s ]
    real(r8)    tke(0:mkx, mix)                                    !  Turbulent kinetic energy at the interfaces [ m2/s2 ]
    real(r8)    cldfrct(mkx, mix)                                  !  Total cloud fraction at the previous time step [ fraction ]
    real(r8)    concldfrct(mkx, mix)                               !  Total convective cloud fraction
                                                              !  at the previous time step [ fraction ]
    real(r8)    qv0(mkx, mix)                                      !  Environmental water vapor specific humidity [ kg/kg ]
    real(r8)    ql0(mkx, mix)                                      !  Environmental liquid water specific humidity [ kg/kg ]
    real(r8)    qi0(mkx, mix)                                      !  Environmental ice specific humidity [ kg/kg ]
    real(r8)    t0(mkx, mix)                                       !  Environmental temperature [ K ]
    real(r8)    s0(mkx, mix)                                       !  Environmental dry static energy [ J/kg ]
    real(r8)    pblh(mix)                                          !  Height of PBL [ m ]
    real(r8)    cush(mix)                                          !  Convective scale height [ m ]
    real(r8)    tr0(mkx,ncnst, mix)                                !  Environmental tracers [ #, kg/kg ]

    ! 2. Environmental variables directly derived from the input variables

    real(r8)    qt0(mkx, mix)                                      !  Environmental total specific humidity [ kg/kg ]
    real(r8)    thl0(mkx, mix)                                     !  Environmental liquid potential temperature [ K ]
    real(r8)    thvl0(mkx, mix)                                    !  Environmental liquid virtual potential temperature [ K ]
    real(r8)    ssqt0(mkx, mix)                                    !  Linear internal slope
                                                              !  of environmental total specific humidity [ kg/kg/Pa ]
    real(r8)    ssthl0(mkx, mix)                                   !  Linear internal slope
                                                              ! of environmental liquid potential temperature [ K/Pa ]
    real(r8)    ssu0(mkx, mix)                                     !  Linear internal slope of environmental zonal wind [ m/s/Pa ]
    real(r8)    ssv0(mkx, mix)                                     !  Linear internal slope of environmental meridional wind [ m/s/Pa ]
    real(r8)    thv0bot(mkx, mix)                                  !  Environmental virtual potential temperature
                                                              ! at the bottom of each layer [ K ]
    real(r8)    thv0top(mkx, mix)                                  !  Environmental virtual potential temperature
                                                              ! at the top of each layer [ K ]
    real(r8)    thvl0bot(mkx, mix)                                 !  Environmental liquid virtual potential temperature
                                                              ! at the bottom of each layer [ K ]
    real(r8)    thvl0top(mkx, mix)                                 !  Environmental liquid virtual potential temperature
                                                              ! at the top of each layer [ K ]
    real(r8)    exn0(mkx, mix)                                     !  Exner function at the layer mid points [ no ]
    real(r8)    exns0(0:mkx, mix)                                  !  Exner function at the interfaces [ no ]
    real(r8)    sstr0(mkx,ncnst, mix)                              !  Linear slope of environmental tracers [ #/Pa, kg/kg/Pa ]

   ! 2-1. For preventing negative condensate at the provisional time step

    real(r8)    qv0_star(mkx, mix)                                 !  Environmental water vapor specific humidity [ kg/kg ]
    real(r8)    ql0_star(mkx, mix)                                 !  Environmental liquid water specific humidity [ kg/kg ]
    real(r8)    qi0_star(mkx, mix)                                 !  Environmental ice specific humidity [ kg/kg ]
    real(r8)    t0_star(mkx, mix)                                  !  Environmental temperature [ K ]
    real(r8)    s0_star(mkx, mix)                                  !  Environmental dry static energy [ J/kg ]

   ! 3. Variables associated with cumulus convection

    real(r8)    umf(0:mkx, mix)                                    !  Updraft mass flux at the interfaces [ kg/m2/s ]
    real(r8)    emf(0:mkx, mix)                                    !  Penetrative entrainment mass flux at the interfaces [ kg/m2/s ]
    real(r8)    qvten(mkx, mix)                                    !  Tendency of water vapor specific humidity [ kg/kg/s ]
    real(r8)    qlten(mkx, mix)                                    !  Tendency of liquid water specific humidity [ kg/kg/s ]
    real(r8)    qiten(mkx, mix)                                    !  Tendency of ice specific humidity [ kg/kg/s ]
    real(r8)    sten(mkx, mix)                                     !  Tendency of dry static energy [ J/kg ]
    real(r8)    uten(mkx, mix)                                     !  Tendency of zonal wind [ m/s2 ]
    real(r8)    vten(mkx, mix)                                     !  Tendency of meridional wind [ m/s2 ]
    real(r8)    qrten(mkx, mix)                                    !  Tendency of rain water specific humidity [ kg/kg/s ]
    real(r8)    qsten(mkx, mix)                                    !  Tendency of snow(mix) specific humidity [ kg/kg/s ]
    real(r8)    precip(mix)                                        !  Precipitation rate ( rain + snow(mix)) at the surface [ m/s ]
    real(r8)    snow(mix)                                          !  Snow rate at the surface [ m/s ]
    real(r8)    evapc(mkx, mix)                                    !  Tendency of evaporation of precipitation [ kg/kg/s ]
    real(r8)    slflx(0:mkx, mix)                                  !  Updraft/pen.entrainment liquid static energy flux
                                                              ! [ J/kg * kg/m2/s ]
    real(r8)    qtflx(0:mkx, mix)                                  !  Updraft/pen.entrainment total water flux [ kg/kg * kg/m2/s ]
    real(r8)    uflx(0:mkx, mix)                                   !  Updraft/pen.entrainment flux of zonal momentum [ m/s/m2/s ]
    real(r8)    vflx(0:mkx, mix)                                   !  Updraft/pen.entrainment flux of meridional momentum [ m/s/m2/s ]
    real(r8)    cufrc(mkx, mix)                                    !  Shallow cumulus cloud fraction at the layer mid-point [ fraction ]
    real(r8)    qcu(mkx, mix)                                      !  Condensate water specific humidity
                                                              ! within convective updraft [ kg/kg ]
    real(r8)    qlu(mkx, mix)                                      !  Liquid water specific humidity within convective updraft [ kg/kg ]
    real(r8)    qiu(mkx, mix)                                      !  Ice specific humidity within convective updraft [ kg/kg ]
    real(r8)    dwten(mkx, mix)                                    !  Detrained water tendency from cumulus updraft [ kg/kg/s ]
    real(r8)    diten(mkx, mix)                                    !  Detrained ice   tendency from cumulus updraft [ kg/kg/s ]
    real(r8)    fer(mkx, mix)                                      !  Fractional lateral entrainment rate [ 1/Pa ]
    real(r8)    fdr(mkx, mix)                                      !  Fractional lateral detrainment rate [ 1/Pa ]
    real(r8)    uf(mkx, mix)                                       !  Zonal wind at the provisional time step [ m/s ]
    real(r8)    vf(mkx, mix)                                       !  Meridional wind at the provisional time step [ m/s ]
    real(r8)    qc(mkx, mix)                                       !  Tendency due to detrained 'cloud water + cloud ice'
                                                              ! (without rain-snow(mix) contribution) [ kg/kg/s ]
    real(r8)    qc_l(mkx, mix)                                     !  Tendency due to detrained 'cloud water'
                                                              ! (without rain-snow(mix) contribution) [ kg/kg/s ]
    real(r8)    qc_i(mkx, mix)                                     !  Tendency due to detrained 'cloud ice'
                                                              ! (without rain-snow(mix) contribution) [ kg/kg/s ]
    real(r8)    qc_lm(mix)
    real(r8)    qc_im(mix)
    real(r8)    nc_lm(mix)
    real(r8)    nc_im(mix)
    real(r8)    ql_emf_kbup(mix)
    real(r8)    qi_emf_kbup(mix)
    real(r8)    nl_emf_kbup(mix)
    real(r8)    ni_emf_kbup(mix)
    real(r8)    qlten_det(mix)
    real(r8)    qiten_det(mix)
    real(r8)    rliq(mix)                                          !  Vertical integral of qc [ m/s ]
    real(r8)    cnt(mix)                                           !  Cumulus top  interface index, cnt(mix) = kpen(mix) [ no ]
    real(r8)    cnb(mix)                                           !  Cumulus base interface index, cnb(mix) = krel(mix) - 1 [ no ]
    real(r8)    qtten(mkx, mix)                                    !  Tendency of qt [ kg/kg/s ]
    real(r8)    slten(mkx, mix)                                    !  Tendency of sl [ J/kg/s ]
    real(r8)    ufrc(0:mkx, mix)                                   !  Updraft fractional area [ fraction ]
    real(r8)    trten(mkx,ncnst, mix)                              !  Tendency of tracers [ #/s, kg/kg/s ]
    real(r8)    trflx(0:mkx,ncnst, mix)                            !  Flux of tracers due to convection [ # * kg/m2/s, kg/kg * kg/m2/s ]
    real(r8)    trflx_d(0:mkx, mix)                                !  Adjustive downward flux of tracers to prevent negative tracers
    real(r8)    trflx_u(0:mkx, mix)                                !  Adjustive upward   flux of tracers to prevent negative tracers
    real(r8)    trmin(mix)                                         !  Minimum concentration of tracers allowed
    real(r8)    pdelx(mix), dum(mix)

    !----- Variables used for the calculation of condensation sink associated with compensating subsidence
    !      In the current code, this 'sink' tendency is simply set to be zero.

    real(r8)    uemf(0:mkx, mix)                                   !  Net updraft mass flux at the interface ( emf + umf ) [ kg/m2/s ]
    real(r8)    comsub(mkx, mix)                                   !  Compensating subsidence
                                                              ! at the layer mid-point ( unit of mass flux, umf ) [ kg/m2/s ]
    real(r8)    qlten_sink(mkx, mix)                               !  Liquid condensate tendency
                                                              ! by compensating subsidence/upwelling [ kg/kg/s ]
    real(r8)    qiten_sink(mkx, mix)                               !  Ice    condensate tendency
                                                              ! by compensating subsidence/upwelling [ kg/kg/s ]
    real(r8)    nlten_sink(mkx, mix)                               !  Liquid droplets # tendency
                                                              ! by compensating subsidence/upwelling [ kg/kg/s ]
    real(r8)    niten_sink(mkx, mix)                               !  Ice    droplets # tendency
                                                              ! by compensating subsidence/upwelling [ kg/kg/s ]
    real(r8)    thlten_sub(mix), qtten_sub(mix)                         !  Tendency of conservative scalars
                                                              ! by compensating subsidence/upwelling
    real(r8)    qlten_sub(mix), qiten_sub(mix)                          !  Tendency of ql0, qi0
                                                              ! by compensating subsidence/upwelling
    real(r8)    nlten_sub(mix), niten_sub(mix)                          !  Tendency of nl0, ni0
                                                              ! by compensating subsidence/upwelling
    real(r8)    thl_prog(mix), qt_prog(mix)                             !  Prognosed 'thl, qt'
                                                              ! by compensating subsidence/upwelling

    !----- Variables describing cumulus updraft

    real(r8)    wu(0:mkx, mix)                                     !  Updraft vertical velocity at the interface [ m/s ]
    real(r8)    thlu(0:mkx, mix)                                   !  Updraft liquid potential temperature at the interface [ K ]
    real(r8)    qtu(0:mkx, mix)                                    !  Updraft total specific humidity at the interface [ kg/kg ]
    real(r8)    uu(0:mkx, mix)                                     !  Updraft zonal wind at the interface [ m/s ]
    real(r8)    vu(0:mkx, mix)                                     !  Updraft meridional wind at the interface [ m/s ]
    real(r8)    thvu(0:mkx, mix)                                   !  Updraft virtual potential temperature at the interface [ m/s ]
    real(r8)    rei(mkx, mix)                                      !  Updraft fractional mixing rate with the environment [ 1/Pa ]
    real(r8)    tru(0:mkx,ncnst, mix)                              !  Updraft tracers [ #, kg/kg ]

    !----- Variables describing conservative scalars of entraining downdrafts  at the
    !      entraining interfaces, i.e., 'kbup(mix) <= k < kpen(mix)-1'. At the other interfaces,
    !      belows are simply set to equal to those of updraft for simplicity - but it
    !      does not influence numerical calculation.

    real(r8)    thlu_emf(0:mkx, mix)                               !  Penetrative downdraft liquid potential temperature
                                                              ! at entraining interfaces [ K ]
    real(r8)    qtu_emf(0:mkx, mix)                                !  Penetrative downdraft total water
                                                              ! at entraining interfaces [ kg/kg ]
    real(r8)    uu_emf(0:mkx, mix)                                 !  Penetrative downdraft zonal wind
                                                              ! at entraining interfaces [ m/s ]
    real(r8)    vu_emf(0:mkx, mix)                                 !  Penetrative downdraft meridional wind
                                                              ! at entraining interfaces [ m/s ]
    real(r8)    tru_emf(0:mkx,ncnst, mix)                          !  Penetrative Downdraft tracers
                                                              ! at entraining interfaces [ #, kg/kg ]

    !----- Variables associated with evaporations of convective 'rain' and 'snow(mix)'

    real(r8)    flxrain(0:mkx, mix)                                !  Downward rain flux at each interface [ kg/m2/s ]
    real(r8)    flxsnow(0:mkx, mix)                                !  Downward snow(mix) flux at each interface [ kg/m2/s ]
    real(r8)    ntraprd(mkx, mix)                                  !  Net production ( production - evaporation +  melting )
                                                              ! rate of rain in each layer [ kg/kg/s ]
    real(r8)    ntsnprd(mkx, mix)                                  !  Net production ( production - evaporation + freezing )
                                                              ! rate of snow(mix) in each layer [ kg/kg/s ]
    real(r8)    flxsntm(mix)                                       !  Downward snow(mix) flux
                                                              ! at the top of each layer after melting [ kg/m2/s ]
    real(r8)    snowmlt(mix)                                       !  Snow melting tendency [ kg/kg/s ]
    real(r8)    subsat(mix)                                        !  Sub-saturation ratio (1-qv/qs(mix)) [ no unit ]
    real(r8)    evprain(mix)                                       !  Evaporation rate of rain [ kg/kg/s ]
    real(r8)    evpsnow(mix)                                       !  Evaporation rate of snow(mix) [ kg/kg/s ]
    real(r8)    evplimit(mix)                                      !  Limiter of 'evprain(mix) + evpsnow(mix)' [ kg/kg/s ]
    real(r8)    evplimit_rain(mix)                                 !  Limiter of 'evprain(mix)' [ kg/kg/s ]
    real(r8)    evplimit_snow(mix)                                 !  Limiter of 'evpsnow(mix)' [ kg/kg/s ]
    real(r8)    evpint_rain(mix)                                   !  Vertically-integrated evaporative flux of rain [ kg/m2/s ]
    real(r8)    evpint_snow(mix)                                   !  Vertically-integrated evaporative flux of snow(mix) [ kg/m2/s ]
    real(r8)    kevp                                          !  Evaporative efficiency [ complex unit ]

    !----- Other internal variables

    integer     kk, mm, k, i, m, kp1, km1
    integer     iter_scaleh, iter_xc
    integer     id_check, status
    integer     klcl(mix)                                          !  Layer containing LCL of source air
    integer     kinv(mix)                                          !  Inversion layer with PBL top interface as a lower interface
    integer     krel(mix)                                          !  Release layer where buoyancy sorting mixing
                                                              ! occurs for the first time
    integer     klfc(mix)                                          !  LFC layer of cumulus source air
    integer     kbup(mix)                                          !  Top layer in which cloud buoyancy is positive at the top interface
    integer     kpen(mix)                                          !  Highest layer with positive updraft vertical velocity
                                                              ! - top layer cumulus can reach
    logical     id_exit(mix)
    logical     id_exit_ex
    logical     forcedCu(mix)                                      !  If 'true', cumulus updraft cannot overcome the buoyancy barrier
                                                              ! just above the PBL top.
    real(r8)    thlsrc(mix), qtsrc(mix), usrc(mix), vsrc(mix), thvlsrc(mix)            !  Updraft source air properties
    real(r8)    PGFc, uplus(mix), vplus(mix)
    real(r8)    trsrc(ncnst, mix), tre(ncnst, mix)
    real(r8)    plcl(mix), plfc(mix), prel(mix), wrel(mix)
    real(r8)    frc_rasn
    real(r8)    ee2(mix), ud2(mix), wtw(mix), wtwb(mix), wtwh(mix)
    real(r8)    xc(mix), xc_2(mix)
    real(r8)    cldhgt(mix), scaleh(mix), tscaleh(mix), cridis(mix), rle, rkm
    real(r8)    rkfre, sigmaw(mix), epsvarw, tkeavg(mix), dpsum(mix), dpi(mix), thvlmin(mix)
    real(r8)    thlxsat(mix), qtxsat(mix), thvxsat(mix), x_cu(mix), x_en(mix), thv_x0(mix), thv_x1(mix)
    real(r8)    thj(mix), qvj(mix), qlj(mix), qij(mix), thvj(mix), tj(mix), thv0j(mix), rho0j(mix), rhos0j(mix), qse(mix)
    real(r8)    cin(mix), cinlcl(mix)
    real(r8)    pe(mix), dpe(mix), exne(mix), thvebot(mix), thle(mix), qte(mix), ue(mix), ve(mix), thlue(mix), qtue(mix), wue(mix)
    real(r8)    mu(mix), mumin0(mix), mumin1, mumin2(mix), mulcl(mix), mulclstar(mix)
    real(r8)    cbmf(mix), wcrit(mix), winv(mix), wlcl(mix), ufrcinv(mix), ufrclcl(mix), rmaxfrac
    real(r8)    criqc, exql(mix), exqi(mix), ppen(mix)
    real(r8)    thl0top(mix), thl0bot(mix), qt0bot(mix), qt0top(mix), thvubot(mix), thvutop(mix)
    real(r8)    thlu_top(mix), qtu_top(mix), qlu_top(mix), qiu_top(mix), qlu_mid(mix), qiu_mid(mix), exntop(mix)
    real(r8)    thl0lcl(mix), qt0lcl(mix), thv0lcl(mix), thv0rel(mix), rho0inv(mix), autodet(mix)
    real(r8)    aquad(mix), bquad(mix), cquad(mix), xc1(mix), xc2(mix), excessu(mix), excess0(mix), xsat(mix), xs1(mix), xs2(mix)
    real(r8)    bogbot(mix), bogtop(mix), delbog(mix), drage(mix), expfac(mix), rbuoy, rdrag
    real(r8)    rcwp(mix), rlwp(mix), riwp(mix), qcubelow(mix), qlubelow(mix), qiubelow(mix)
    real(r8)    rainflx(mix), snowflx(mix)
    real(r8)    es(mix)
    real(r8)    qs(mix)
    real(r8)    qsat_arg(mix)
    real(r8)    xsrc(mix), xmean(mix), xtop(mix), xbot(mix), xflx(0:mkx, mix)
    real(r8)    tmp1(mix), tmp2(mix)

    !----- Some diagnostic internal output variables

    real(r8)  ufrcinvbase_out(mix)                            !  Cumulus updraft fraction at the PBL top [ fraction ]
    real(r8)  ufrclcl_out(mix)                                !  Cumulus updraft fraction at the LCL
                                                              ! ( or PBL top when LCL is below PBL top ) [ fraction ]
    real(r8)  winvbase_out(mix)                               !  Cumulus updraft velocity at the PBL top [ m/s ]
    real(r8)  wlcl_out(mix)                                   !  Cumulus updraft velocity at the LCL
                                                              ! ( or PBL top when LCL is below PBL top ) [ m/s ]
    real(r8)  plcl_out(mix)                                   !  LCL of source air [ Pa ]
    real(r8)  pinv_out(mix)                                   !  PBL top pressure [ Pa ]
    real(r8)  plfc_out(mix)                                   !  LFC of source air [ Pa ]
    real(r8)  pbup_out(mix)                                   !  Highest interface level of positive buoyancy [ Pa ]
    real(r8)  ppen_out(mix)                                   !  Highest interface evel where Cu w = 0 [ Pa ]
    real(r8)  qtsrc_out(mix)                                  !  Sourse air qt [ kg/kg ]
    real(r8)  thlsrc_out(mix)                                 !  Sourse air thl [ K ]
    real(r8)  thvlsrc_out(mix)                                !  Sourse air thvl [ K ]
    real(r8)  emfkbup_out(mix)                                !  Penetrative downward mass flux at 'kbup(mix)' interface [ kg/m2/s ]
    real(r8)  cinlclh_out(mix)                                !  Convective INhibition upto LCL (CIN) [ J/kg = m2/s2 ]
    real(r8)  tkeavg_out(mix)                                 !  Average tke over the PBL [ m2/s2 ]
    real(r8)  cbmflimit_out(mix)                              !  Cloud base mass flux limiter [ kg/m2/s ]
    real(r8)  zinv_out(mix)                                   !  PBL top height [ m ]
    real(r8)  rcwp_out(mix)                                   !  Layer mean Cumulus LWP+IWP [ kg/m2 ]
    real(r8)  rlwp_out(mix)                                   !  Layer mean Cumulus LWP [ kg/m2 ]
    real(r8)  riwp_out(mix)                                   !  Layer mean Cumulus IWP [ kg/m2 ]
    real(r8)  wu_out(mix,0:mkx)                               !  Updraft vertical velocity
                                                              ! ( defined from the release level to 'kpen(mix)-1' interface )
    real(r8)  qtu_out(mix,0:mkx)                              !  Updraft qt [ kg/kg ]
    real(r8)  thlu_out(mix,0:mkx)                             !  Updraft thl [ K ]
    real(r8)  thvu_out(mix,0:mkx)                             !  Updraft thv [ K ]
    real(r8)  uu_out(mix,0:mkx)                               !  Updraft zonal wind [ m/s ]
    real(r8)  vu_out(mix,0:mkx)                               !  Updraft meridional wind [ m/s ]
    real(r8)  qtu_emf_out(mix,0:mkx)                          !  Penetratively entrained qt [ kg/kg ]
    real(r8)  thlu_emf_out(mix,0:mkx)                         !  Penetratively entrained thl [ K ]
    real(r8)  uu_emf_out(mix,0:mkx)                           !  Penetratively entrained u [ m/s ]
    real(r8)  vu_emf_out(mix,0:mkx)                           !  Penetratively entrained v [ m/s ]
    real(r8)  uemf_out(mix,0:mkx)                             !  Net upward mass flux
                                                              ! including penetrative entrainment (umf+emf) [ kg/m2/s ]
    real(r8)  tru_out(mix,0:mkx,ncnst)                        !  Updraft tracers [ #, kg/kg ]
    real(r8)  tru_emf_out(mix,0:mkx,ncnst)                    !  Penetratively entrained tracers [ #, kg/kg ]

    real(r8)  wu_s(0:mkx, mix)                                     !  Same as above but for implicit CIN
    real(r8)  qtu_s(0:mkx, mix)
    real(r8)  thlu_s(0:mkx, mix)
    real(r8)  thvu_s(0:mkx, mix)
    real(r8)  uu_s(0:mkx, mix)
    real(r8)  vu_s(0:mkx, mix)
    real(r8)  qtu_emf_s(0:mkx, mix)
    real(r8)  thlu_emf_s(0:mkx, mix)
    real(r8)  uu_emf_s(0:mkx, mix)
    real(r8)  vu_emf_s(0:mkx, mix)
    real(r8)  uemf_s(0:mkx, mix)
    real(r8)  tru_s(0:mkx,ncnst, mix)
    real(r8)  tru_emf_s(0:mkx,ncnst, mix)

    real(r8)  dwten_out(mix,mkx)
    real(r8)  diten_out(mix,mkx)
    real(r8)  flxrain_out(mix,0:mkx)
    real(r8)  flxsnow_out(mix,0:mkx)
    real(r8)  ntraprd_out(mix,mkx)
    real(r8)  ntsnprd_out(mix,mkx)

    real(r8)  dwten_s(mkx, mix)
    real(r8)  diten_s(mkx, mix)
    real(r8)  flxrain_s(0:mkx, mix)
    real(r8)  flxsnow_s(0:mkx, mix)
    real(r8)  ntraprd_s(mkx, mix)
    real(r8)  ntsnprd_s(mkx, mix)

    real(r8)  excessu_arr_out(mix,mkx)
    real(r8)  excessu_arr(mkx, mix)
    real(r8)  excessu_arr_s(mkx, mix)
    real(r8)  excess0_arr_out(mix,mkx)
    real(r8)  excess0_arr(mkx, mix)
    real(r8)  excess0_arr_s(mkx, mix)
    real(r8)  xc_arr_out(mix,mkx)
    real(r8)  xc_arr(mkx, mix)
    real(r8)  xc_arr_s(mkx, mix)
    real(r8)  aquad_arr_out(mix,mkx)
    real(r8)  aquad_arr(mkx, mix)
    real(r8)  aquad_arr_s(mkx, mix)
    real(r8)  bquad_arr_out(mix,mkx)
    real(r8)  bquad_arr(mkx, mix)
    real(r8)  bquad_arr_s(mkx, mix)
    real(r8)  cquad_arr_out(mix,mkx)
    real(r8)  cquad_arr(mkx, mix)
    real(r8)  cquad_arr_s(mkx, mix)
    real(r8)  bogbot_arr_out(mix,mkx)
    real(r8)  bogbot_arr(mkx, mix)
    real(r8)  bogbot_arr_s(mkx, mix)
    real(r8)  bogtop_arr_out(mix,mkx)
    real(r8)  bogtop_arr(mkx, mix)
    real(r8)  bogtop_arr_s(mkx, mix)

    real(r8)  exit_UWCu(mix)
    real(r8)  exit_conden(mix)
    real(r8)  exit_klclmkx(mix)
    real(r8)  exit_klfcmkx(mix)
    real(r8)  exit_ufrc(mix)
    real(r8)  exit_wtw(mix)
    real(r8)  exit_drycore(mix)
    real(r8)  exit_wu(mix)
    real(r8)  exit_cufilter(mix)
    real(r8)  exit_kinv1(mix)
    real(r8)  exit_rei(mix)

    real(r8)  limit_shcu(mix)
    real(r8)  limit_negcon(mix)
    real(r8)  limit_ufrc(mix)
    real(r8)  limit_ppen(mix)
    real(r8)  limit_emf(mix)
    real(r8)  limit_cinlcl(mix)
    real(r8)  limit_cin(mix)
    real(r8)  limit_cbmf(mix)
    real(r8)  limit_rei(mix)
    real(r8)  ind_delcin(mix)

!!! Start of _s&_o variables

    real(r8), dimension(mix) :: ufrcinvbase_s, ufrclcl_s, winvbase_s, wlcl_s, plcl_s, pinv_s, plfc_s, &
                qtsrc_s, thlsrc_s, thvlsrc_s, emfkbup_s, cinlcl_s, pbup_s, ppen_s, cbmflimit_s, &
                tkeavg_s, zinv_s, rcwp_s, rlwp_s, riwp_s
    real(r8), dimension(mix) :: ufrcinvbase, winvbase, pinv, zinv, emfkbup, cbmflimit, rho0rel

    !----- Variables for implicit CIN computation

    real(r8), dimension(mkx, mix)         :: qv0_s  , ql0_s   , qi0_s   , s0_s    , u0_s    ,           &
                                        v0_s   , t0_s    , qt0_s   , thl0_s  , thvl0_s , qvten_s , &
                                        qlten_s, qiten_s , qrten_s , qsten_s , sten_s  , evapc_s , &
                                        uten_s , vten_s  , cufrc_s , qcu_s   , qlu_s   , qiu_s   , &
                                        fer_s  , fdr_s   , qc_s    , qtten_s , slten_s
    real(r8), dimension(0:mkx, mix)       :: umf_s  , slflx_s , qtflx_s , ufrc_s  , uflx_s , vflx_s
    real(r8), dimension(mix)                         :: cush_s , precip_s, snow_s  , cin_s   , rliq_s, cbmf_s, cnt_s, cnb_s
	real(r8), dimension(mix)                         :: cin_i,cin_f,del_CIN,ke,alpha,thlj
    real(r8), dimension(mix)                         :: cinlcl_i,cinlcl_f,del_cinlcl

    real(r8), dimension(mkx,ncnst, mix)   :: tr0_s, trten_s
    real(r8), dimension(0:mkx,ncnst, mix) :: trflx_s

    !----- Variables for temporary storages

    real(r8), dimension(mkx, mix)         :: qv0_o, ql0_o, qi0_o, t0_o, s0_o, u0_o, v0_o
    real(r8), dimension(mkx, mix)         :: qt0_o    , thl0_o   , thvl0_o   ,                         &
                                        qvten_o  , qlten_o  , qiten_o   , qrten_o   , qsten_o ,   &
                                        sten_o   , uten_o   , vten_o    , qcu_o     , qlu_o   ,   &
                                        qiu_o    , cufrc_o  , evapc_o   ,                         &
                                        thv0bot_o, thv0top_o, thvl0bot_o, thvl0top_o,             &
                                        ssthl0_o , ssqt0_o  , ssu0_o    , ssv0_o    , qc_o    ,   &
                                        qtten_o  , slten_o
    real(r8), dimension(0:mkx, mix)       :: umf_o    , slflx_o  , qtflx_o   , ufrc_o
    real(r8), dimension(mix)         :: cush_o   , precip_o , snow_o    , rliq_o, cbmf_o, cnt_o, cnb_o
    real(r8), dimension(0:mkx, mix)       :: uflx_o   , vflx_o
    real(r8), dimension(mix)                         :: tkeavg_o , thvlmin_o, qtsrc_o  , thvlsrc_o, thlsrc_o ,    &
                                        usrc_o   , vsrc_o   , plcl_o   , plfc_o   ,               &
                                        thv0lcl_o, cinlcl_o
    integer, dimension(mix)                          :: kinv_o   , klcl_o   , klfc_o

    real(r8), dimension(mkx,ncnst, mix)   :: tr0_o
    real(r8), dimension(mkx,ncnst, mix)   :: trten_o, sstr0_o
    real(r8), dimension(0:mkx,ncnst, mix) :: trflx_o
    real(r8), dimension(ncnst, mix)       :: trsrc_o

!!! End of _s&_o variables

    integer                          :: iter
	integer                          :: ixnumliq, ixnumice, ixcldliq, ixcldice

    ! ------------------ !
    !                    !
    ! Define Parameters  !
    !                    !
    ! ------------------ !

    ! ------------------------ !
    ! Iterative xc(mix) calculation !
    ! ------------------------ !

    integer , parameter              :: niter_xc = 2

    ! ----------------------------------------------------------- !
    ! Choice of 'CIN = cin(mix)' (.true.) or 'CIN = cinlcl(mix)' (.false.). !
    !                                                             !
    ! Feb 2007, Bundy: Note that use_CINcin = .false. will try to !
    ! use a variable (del_cinlcl(mix)) that is not currently set       !
    !                                                             !
    ! Sept 2012, Santos: The fact that this is still true over 5  !
    ! years later suggests that this option needs to be           !
    ! fixed or abandoned.                                         !
    ! ----------------------------------------------------------- !

    logical , parameter              :: use_CINcin = .true.

    ! --------------------------------------------------------------- !
    ! Choice of 'explicit' ( 1 ) or 'implicit' ( 2 )  CIN.            !
    !                                                                 !
    ! When choose 'CIN = cinlcl(mix)' above,  it is recommended not to use !
    ! implicit CIN, i.e., do 'NOT' choose simultaneously :            !
    !            [ 'use_CINcin=.false. & 'iter_cin=2' ]               !
    ! since 'cinlcl(mix)' will be always set to zero whenever LCL is below !
    ! the PBL top interface in the current code. So, averaging cinlcl(mix) !
    ! of two iter_cin steps is likely not so good. Except that,   all !
    ! the other combinations of  'use_CINcin'  & 'iter_cin' are OK.   !
    ! --------------------------------------------------------------- !

    integer , parameter              :: iter_cin = 2

    ! ---------------------------------------------------------------- !
    ! Choice of 'self-detrainment' by negative buoyancy in calculating !
    ! cumulus updraft mass flux at the top interface in each layer.    !
    ! ---------------------------------------------------------------- !

    logical , parameter              :: use_self_detrain = .false.

    ! --------------------------------------------------------- !
    ! Cumulus momentum flux : turn-on (.true.) or off (.false.) !
    ! --------------------------------------------------------- !

    logical , parameter              :: use_momenflx = .true.

    ! ----------------------------------------------------------------------------------------- !
    ! Penetrative Entrainment : Cumulative ( .true. , original ) or Non-Cumulative ( .false. )  !
    ! This option ( .false. ) is designed to reduce the sensitivity to the vertical resolution. !
    ! ----------------------------------------------------------------------------------------- !

    logical , parameter              :: use_cumpenent = .true.

    ! --------------------------------------------------------------------------------------------------------------- !
    ! Computation of the grid-mean condensate tendency.                                                               !
    !     use_expconten = .true.  : explcitly compute tendency by condensate detrainment and compensating subsidence  !
    !     use_expconten = .false. : use the original proportional condensate tendency equation. ( original )          !
    ! --------------------------------------------------------------------------------------------------------------- !

    logical , parameter              :: use_expconten = .true.

    ! --------------------------------------------------------------------------------------------------------------- !
    ! Treatment of reserved condensate                                                                                !
    !     use_unicondet = .true.  : detrain condensate uniformly over the environment ( original )                    !
    !     use_unicondet = .false. : detrain condensate into the pre-existing stratus                                  !
    ! --------------------------------------------------------------------------------------------------------------- !

    logical , parameter              :: use_unicondet = .false.

    ! ----------------------- !
    ! For lateral entrainment !
    ! ----------------------- !

    parameter (rle = 0.1_r8)         !  For critical stopping distance for lateral entrainment [no unit]
!   parameter (rkm = 16.0_r8)        !  Determine the amount of air that is involved in buoyancy-sorting [no unit]
    parameter (rkm = 14.0_r8)        !  Determine the amount of air that is involved in buoyancy-sorting [no unit]

    parameter (rkfre = 1.0_r8)       !  Vertical velocity variance as fraction of  tke.
    parameter (rmaxfrac = 0.10_r8)   !  Maximum allowable 'core' updraft fraction
    parameter (mumin1 = 0.906_r8)    !  Normalized CIN ('mu(mix)') corresponding to 'rmaxfrac' at the PBL top
                                     !  obtaind by inverting 'rmaxfrac = 0.5*erfc(mumin1)'.
                                     !  [rmaxfrac:mumin1]=[ 0.05:1.163, 0.075:1.018, 0.1:0.906, 0.15:0.733, 0.2:0.595, 0.25:0.477]
    parameter (rbuoy = 1.0_r8)       !  For nonhydrostatic pressure effects on updraft [no unit]
    parameter (rdrag = 1.0_r8)       !  Drag coefficient [no unit]

    parameter (epsvarw = 5.e-4_r8)   !  Variance of w at PBL top by meso-scale component [m2/s2]
    parameter (PGFc = 0.7_r8)        !  This is used for calculating vertical variations cumulus
                                     !  'u' & 'v' by horizontal PGF during upward motion [no unit]

    ! ---------------------------------------- !
    ! Bulk microphysics controlling parameters !
    ! --------------------------------------------------------------------------- !
    ! criqc    : Maximum condensate that can be hold by cumulus updraft [kg/kg]   !
    ! frc_rasn : Fraction of precipitable condensate in the expelled cloud water  !
    !            from cumulus updraft. The remaining fraction ('1-frc_rasn')  is  !
    !            'suspended condensate'.                                          !
    !                0 : all expelled condensate is 'suspended condensate'        !
    !                1 : all expelled condensate is 'precipitable condensate'     !
    ! kevp     : Evaporative efficiency                                           !
    ! noevap_krelkpen : No evaporation from 'krel(mix)' to 'kpen(mix)' layers               !
    ! --------------------------------------------------------------------------- !

    parameter ( criqc    = 0.7e-3_r8 )
    parameter ( frc_rasn = 1.0_r8    )
    parameter ( kevp     = 2.e-6_r8  )
    logical, parameter :: noevap_krelkpen = .false.

!!! end of initialization

    !------------------------!
    !                        !
    ! Start Main Calculation !
    !                        !
    !------------------------!

    call cnst_get_ind( 'NUMLIQ', ixnumliq )
    call cnst_get_ind( 'NUMICE', ixnumice )

    call cnst_get_ind( 'CLDLIQ', ixcldliq )
    call cnst_get_ind( 'CLDICE', ixcldice )




    ! ------------------------------------------------------- !
    ! Initialize output variables defined for all grid points !
    ! ------------------------------------------------------- !

    umf_out(:iend,0:mkx)         = 0.0_r8
    slflx_out(:iend,0:mkx)       = 0.0_r8
    qtflx_out(:iend,0:mkx)       = 0.0_r8
    flxprc1_out(:iend,0:mkx)     = 0.0_r8
    flxsnow1_out(:iend,0:mkx)    = 0.0_r8
    qvten_out(:iend,:mkx)        = 0.0_r8
    qlten_out(:iend,:mkx)        = 0.0_r8
    qiten_out(:iend,:mkx)        = 0.0_r8
    sten_out(:iend,:mkx)         = 0.0_r8
    uten_out(:iend,:mkx)         = 0.0_r8
    vten_out(:iend,:mkx)         = 0.0_r8
    qrten_out(:iend,:mkx)        = 0.0_r8
    qsten_out(:iend,:mkx)        = 0.0_r8
    precip_out(:iend)            = 0.0_r8
    snow_out(:iend)              = 0.0_r8
    evapc_out(:iend,:mkx)        = 0.0_r8
    cufrc_out(:iend,:mkx)        = 0.0_r8
    qcu_out(:iend,:mkx)          = 0.0_r8
    qlu_out(:iend,:mkx)          = 0.0_r8
    qiu_out(:iend,:mkx)          = 0.0_r8
    fer_out(:iend,:mkx)          = 0.0_r8
    fdr_out(:iend,:mkx)          = 0.0_r8
    cinh_out(:iend)              = -1.0_r8
    cinlclh_out(:iend)           = -1.0_r8
    cbmf_out(:iend)              = 0.0_r8
    qc_out(:iend,:mkx)           = 0.0_r8
    rliq_out(:iend)              = 0.0_r8
    cnt_out(:iend)               = real(mkx, r8)
    cnb_out(:iend)               = 0.0_r8
    qtten_out(:iend,:mkx)        = 0.0_r8
    slten_out(:iend,:mkx)        = 0.0_r8
    ufrc_out(:iend,0:mkx)        = 0.0_r8

    uflx_out(:iend,0:mkx)        = 0.0_r8
    vflx_out(:iend,0:mkx)        = 0.0_r8

    trten_out(:iend,:mkx,:ncnst) = 0.0_r8
    trflx_out(:iend,0:mkx,:ncnst)= 0.0_r8

    ufrcinvbase_out(:iend)       = 0.0_r8
    ufrclcl_out(:iend)           = 0.0_r8
    winvbase_out(:iend)          = 0.0_r8
    wlcl_out(:iend)              = 0.0_r8
    plcl_out(:iend)              = 0.0_r8
    pinv_out(:iend)              = 0.0_r8
    plfc_out(:iend)              = 0.0_r8
    pbup_out(:iend)              = 0.0_r8
    ppen_out(:iend)              = 0.0_r8
    qtsrc_out(:iend)             = 0.0_r8
    thlsrc_out(:iend)            = 0.0_r8
    thvlsrc_out(:iend)           = 0.0_r8
    emfkbup_out(:iend)           = 0.0_r8
    cbmflimit_out(:iend)         = 0.0_r8
    tkeavg_out(:iend)            = 0.0_r8
    zinv_out(:iend)              = 0.0_r8
    rcwp_out(:iend)              = 0.0_r8
    rlwp_out(:iend)              = 0.0_r8
    riwp_out(:iend)              = 0.0_r8

    wu_out(:iend,0:mkx)          = 0.0_r8
    qtu_out(:iend,0:mkx)         = 0.0_r8
    thlu_out(:iend,0:mkx)        = 0.0_r8
    thvu_out(:iend,0:mkx)        = 0.0_r8
    uu_out(:iend,0:mkx)          = 0.0_r8
    vu_out(:iend,0:mkx)          = 0.0_r8
    qtu_emf_out(:iend,0:mkx)     = 0.0_r8
    thlu_emf_out(:iend,0:mkx)    = 0.0_r8
    uu_emf_out(:iend,0:mkx)      = 0.0_r8
    vu_emf_out(:iend,0:mkx)      = 0.0_r8
    uemf_out(:iend,0:mkx)        = 0.0_r8

    tru_out(:iend,0:mkx,:ncnst)     = 0.0_r8
    tru_emf_out(:iend,0:mkx,:ncnst) = 0.0_r8

    dwten_out(:iend,:mkx)        = 0.0_r8
    diten_out(:iend,:mkx)        = 0.0_r8
    flxrain_out(:iend,0:mkx)     = 0.0_r8
    flxsnow_out(:iend,0:mkx)     = 0.0_r8
    ntraprd_out(:iend,mkx)       = 0.0_r8
    ntsnprd_out(:iend,mkx)       = 0.0_r8

    excessu_arr_out(:iend,:mkx)  = 0.0_r8
    excess0_arr_out(:iend,:mkx)  = 0.0_r8
    xc_arr_out(:iend,:mkx)       = 0.0_r8
    aquad_arr_out(:iend,:mkx)    = 0.0_r8
    bquad_arr_out(:iend,:mkx)    = 0.0_r8
    cquad_arr_out(:iend,:mkx)    = 0.0_r8
    bogbot_arr_out(:iend,:mkx)   = 0.0_r8
    bogtop_arr_out(:iend,:mkx)   = 0.0_r8

    exit_UWCu(:iend)             = 0.0_r8
    exit_conden(:iend)           = 0.0_r8
    exit_klclmkx(:iend)          = 0.0_r8
    exit_klfcmkx(:iend)          = 0.0_r8
    exit_ufrc(:iend)             = 0.0_r8
    exit_wtw(:iend)              = 0.0_r8
    exit_drycore(:iend)          = 0.0_r8
    exit_wu(:iend)               = 0.0_r8
    exit_cufilter(:iend)         = 0.0_r8
    exit_kinv1(:iend)            = 0.0_r8
    exit_rei(:iend)              = 0.0_r8

    limit_shcu(:iend)            = 0.0_r8
    limit_negcon(:iend)          = 0.0_r8
    limit_ufrc(:iend)            = 0.0_r8
    limit_ppen(:iend)            = 0.0_r8
    limit_emf(:iend)             = 0.0_r8
    limit_cinlcl(:iend)          = 0.0_r8
    limit_cin(:iend)             = 0.0_r8
    limit_cbmf(:iend)            = 0.0_r8
    limit_rei(:iend)             = 0.0_r8

    ind_delcin(:iend)            = 0.0_r8

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

      id_exit(i) = .false.
      id_exit_ex = .false.

!!!start part one
      ! -------------------------------------------- !
      ! Define 1D input variables at each grid point !
      ! -------------------------------------------- !

      ps0(0:mkx, i)       = ps0_in(i,0:mkx)
      zs0(0:mkx, i)       = zs0_in(i,0:mkx)
      p0(:mkx, i)         = p0_in(i,:mkx)
      z0(:mkx, i)         = z0_in(i,:mkx)
      dp0(:mkx, i)        = dp0_in(i,:mkx)
      dpdry0(:mkx, i)     = dpdry0_in(i,:mkx)
      u0(:mkx, i)         = u0_in(i,:mkx)
      v0(:mkx, i)         = v0_in(i,:mkx)
      qv0(:mkx, i)        = qv0_in(i,:mkx)
      ql0(:mkx, i)        = ql0_in(i,:mkx)
      qi0(:mkx, i)        = qi0_in(i,:mkx)
      t0(:mkx, i)         = t0_in(i,:mkx)
      s0(:mkx, i)         = s0_in(i,:mkx)
      tke(0:mkx, i)       = tke_in(i,0:mkx)
      cldfrct(:mkx, i)    = cldfrct_in(i,:mkx)
      concldfrct(:mkx, i) = concldfrct_in(i,:mkx)
      pblh(i)             = pblh_in(i)
      cush(i)             = cush_inout(i)
      do m = 1, ncnst
         tr0(:mkx,m, i)   = tr0_in(i,:mkx,m)
      enddo

      ! --------------------------------------------------------- !
      ! Compute other basic thermodynamic variables directly from !
      ! the input variables at each grid point                    !
      ! --------------------------------------------------------- !

      !----- 1. Compute internal environmental variables

      exn0(:mkx, i)   = (p0(:mkx, i)/p00)**rovcp
      exns0(0:mkx, i) = (ps0(0:mkx, i)/p00)**rovcp
      qt0(:mkx, i)    = (qv0(:mkx, i) + ql0(:mkx, i) + qi0(:mkx, i))
      thl0(:mkx, i)   = (t0(:mkx, i) - xlv*ql0(:mkx, i)/cp - xls*qi0(:mkx, i)/cp)/exn0(:mkx, i)
      thvl0(:mkx, i)  = (1._r8 + zvir*qt0(:mkx, i))*thl0(:mkx, i)

      !----- 2. Compute slopes of environmental variables in each layer
      !         Dimension of ssthl0(:mkx, i) is implicit.

      ssthl0(:, i)       = slope(mkx,thl0(:, i),p0(:, i))
      ssqt0(:, i)        = slope(mkx,qt0(:, i) ,p0(:, i))
      ssu0(:, i)         = slope(mkx,u0(:, i)  ,p0(:, i))
      ssv0(:, i)         = slope(mkx,v0(:, i)  ,p0(:, i))
      do m = 1, ncnst
         sstr0(:mkx,m, i) = slope(mkx,tr0(:mkx,m, i),p0(:, i))
      enddo

      !----- 3. Compute "thv0" and "thvl0(:, i)" at the top/bottom interfaces in each layer
      !         There are computed from the reconstructed thl, qt at the top/bottom.

      do k = 1, mkx

         thl0bot(i) = thl0(k, i) + ssthl0(k, i)*(ps0(k-1, i) - p0(k, i))
         qt0bot(i)  = qt0(k, i)  + ssqt0(k, i) *(ps0(k-1, i) - p0(k, i))
         call conden(ps0(k-1, i),thl0bot(i),qt0bot(i),thj(i),qvj(i),qlj(i),qij(i),qse(i),id_check)
         if( id_check .eq. 1 ) then
             exit_conden(i) = 1._r8
             id_exit(i) = .true.
             go to 333
         end if
         thv0bot(k, i)  = thj(i)*(1._r8 + zvir*qvj(i) - qlj(i) - qij(i))
         thvl0bot(k, i) = thl0bot(i)*(1._r8 + zvir*qt0bot(i))

         thl0top(i) = thl0(k, i) + ssthl0(k, i)*(ps0(k, i) - p0(k, i))
         qt0top(i)  =  qt0(k, i) + ssqt0(k, i) *(ps0(k, i) - p0(k, i))
         call conden(ps0(k, i),thl0top(i),qt0top(i),thj(i),qvj(i),qlj(i),qij(i),qse(i),id_check)
         if( id_check .eq. 1 ) then
             exit_conden(i) = 1._r8
             id_exit(i) = .true.
             go to 333
         end if
         thv0top(k, i)  = thj(i)*(1._r8 + zvir*qvj(i) - qlj(i) - qij(i))
         thvl0top(k, i) = thl0top(i)*(1._r8 + zvir*qt0top(i))

      end do

      ! ------------------------------------------------------------ !
      ! Save input and related environmental thermodynamic variables !
      ! for use at "iter_cin=2" when "del_CIN(i) >= 0"                  !
      ! ------------------------------------------------------------ !

      qv0_o(:mkx, i)          = qv0(:mkx, i)
      ql0_o(:mkx, i)          = ql0(:mkx, i)
      qi0_o(:mkx, i)          = qi0(:mkx, i)
      t0_o(:mkx, i)           = t0(:mkx, i)
      s0_o(:mkx, i)           = s0(:mkx, i)
      u0_o(:mkx, i)           = u0(:mkx, i)
      v0_o(:mkx, i)           = v0(:mkx, i)
      qt0_o(:mkx, i)          = qt0(:mkx, i)
      thl0_o(:mkx, i)         = thl0(:mkx, i)
      thvl0_o(:mkx, i)        = thvl0(:mkx, i)
      ssthl0_o(:mkx, i)       = ssthl0(:mkx, i)
      ssqt0_o(:mkx, i)        = ssqt0(:mkx, i)
      thv0bot_o(:mkx, i)      = thv0bot(:mkx, i)
      thv0top_o(:mkx, i)      = thv0top(:mkx, i)
      thvl0bot_o(:mkx, i)     = thvl0bot(:mkx, i)
      thvl0top_o(:mkx, i)     = thvl0top(:mkx, i)
      ssu0_o(:mkx, i)         = ssu0(:mkx, i)
      ssv0_o(:mkx, i)         = ssv0(:mkx, i)
      do m = 1, ncnst
         tr0_o(:mkx,m, i)     = tr0(:mkx,m, i)
         sstr0_o(:mkx,m, i)   = sstr0(:mkx,m, i)
      enddo

      ! ---------------------------------------------- !
      ! Initialize output variables at each grid point !
      ! ---------------------------------------------- !

      umf(0:mkx, i)          = 0.0_r8
      emf(0:mkx, i)          = 0.0_r8
      slflx(0:mkx, i)        = 0.0_r8
      qtflx(0:mkx, i)        = 0.0_r8
      uflx(0:mkx, i)         = 0.0_r8
      vflx(0:mkx, i)         = 0.0_r8
      qvten(:mkx, i)         = 0.0_r8
      qlten(:mkx, i)         = 0.0_r8
      qiten(:mkx, i)         = 0.0_r8
      sten(:mkx, i)          = 0.0_r8
      uten(:mkx, i)          = 0.0_r8
      vten(:mkx, i)          = 0.0_r8
      qrten(:mkx, i)         = 0.0_r8
      qsten(:mkx, i)         = 0.0_r8
      dwten(:mkx, i)         = 0.0_r8
      diten(:mkx, i)         = 0.0_r8
      precip(i)              = 0.0_r8
      snow(i)                = 0.0_r8
      evapc(:mkx, i)         = 0.0_r8
      cufrc(:mkx, i)         = 0.0_r8
      qcu(:mkx, i)           = 0.0_r8
      qlu(:mkx, i)           = 0.0_r8
      qiu(:mkx, i)           = 0.0_r8
      fer(:mkx, i)           = 0.0_r8
      fdr(:mkx, i)           = 0.0_r8
      cin(i)                 = 0.0_r8
      cbmf(i)                = 0.0_r8
      qc(:mkx, i)            = 0.0_r8
      qc_l(:mkx, i)          = 0.0_r8
      qc_i(:mkx, i)          = 0.0_r8
      rliq(i)                = 0.0_r8
      cnt(i)                 = real(mkx, r8)
      cnb(i)                 = 0.0_r8
      qtten(:mkx, i)         = 0.0_r8
      slten(:mkx, i)         = 0.0_r8
      ufrc(0:mkx, i)         = 0.0_r8

      thlu(0:mkx, i)         = 0.0_r8
      qtu(0:mkx, i)          = 0.0_r8
      uu(0:mkx, i)           = 0.0_r8
      vu(0:mkx, i)           = 0.0_r8
      wu(0:mkx, i)           = 0.0_r8
      thvu(0:mkx, i)         = 0.0_r8
      thlu_emf(0:mkx, i)     = 0.0_r8
      qtu_emf(0:mkx, i)      = 0.0_r8
      uu_emf(0:mkx, i)       = 0.0_r8
      vu_emf(0:mkx, i)       = 0.0_r8

      ufrcinvbase(i)         = 0.0_r8
      ufrclcl(i)             = 0.0_r8
      winvbase(i)            = 0.0_r8
      wlcl(i)                = 0.0_r8
      emfkbup(i)             = 0.0_r8
      cbmflimit(i)           = 0.0_r8
      excessu_arr(:mkx, i)   = 0.0_r8
      excess0_arr(:mkx, i)   = 0.0_r8
      xc_arr(:mkx, i)        = 0.0_r8
      aquad_arr(:mkx, i)     = 0.0_r8
      bquad_arr(:mkx, i)     = 0.0_r8
      cquad_arr(:mkx, i)     = 0.0_r8
      bogbot_arr(:mkx, i)    = 0.0_r8
      bogtop_arr(:mkx, i)    = 0.0_r8

      uemf(0:mkx, i)         = 0.0_r8
      comsub(:mkx, i)        = 0.0_r8
      qlten_sink(:mkx, i)    = 0.0_r8
      qiten_sink(:mkx, i)    = 0.0_r8
      nlten_sink(:mkx, i)    = 0.0_r8
      niten_sink(:mkx, i)    = 0.0_r8

      do m = 1, ncnst
         trflx(0:mkx,m, i)   = 0.0_r8
         trten(:mkx,m, i)    = 0.0_r8
         tru(0:mkx,m, i)     = 0.0_r8
         tru_emf(0:mkx,m, i) = 0.0_r8
      enddo

    !-----------------------------------------------!
    ! Below 'iter' loop is for implicit CIN closure !
    !-----------------------------------------------!

    ! ----------------------------------------------------------------------------- !
    ! It is important to note that this iterative cin(i) loop is located at the outest !
    ! shell of the code. Thus, source air properties can also be changed during the !
    ! iterative cin(i) calculation, because cumulus convection induces non-zero fluxes !
    ! even at interfaces below PBL top height through 'fluxbelowinv' subroutine.    !
    ! ----------------------------------------------------------------------------- !

    ! CODE OF goto333 MUST BE COPIED & RENUMBERED HERE !!!

    enddo ! init i loop

    do iter = 1, iter_cin

    do i = 1, iend

       if (id_exit(i)) then
          cycle
       endif

       ! ---------------------------------------------------------------------- !
       ! Cumulus scale height                                                   !
       ! In contrast to the premitive code, cumulus scale height is iteratively !
       ! calculated at each time step, and at each iterative cin(i) step.          !
       ! It is not clear whether I should locate below two lines within or  out !
       ! of the iterative cin(i) loop.                                             !
       ! ---------------------------------------------------------------------- !

       tscaleh(i) = cush(i)
       cush(i)    = -1._r8

       ! ----------------------------------------------------------------------- !
       ! Find PBL top height interface index, 'kinv(i)-1' where 'kinv(i)' is the layer !
       ! index with PBLH in it. When PBLH is exactly at interface, 'kinv(i)' is the !
       ! layer index having PBLH as a lower interface.                           !
       ! In the previous code, I set the lower limit of 'kinv(i)' by 2  in order to !
       ! be consistent with the other parts of the code. However in the modified !
       ! code, I allowed 'kinv(i)' to be 1 & if 'kinv(i) = 1', I just exit the program !
       ! without performing cumulus convection. This new approach seems to be    !
       ! more reasonable: if PBL height is within 'kinv(i)=1' layer, surface is STL !
       ! interface (bflxs <= 0) and interface just above the surface should be   !
       ! either non-turbulent (Ri>0.19) or stably turbulent (0<=Ri<0.19 but this !
       ! interface is identified as a base external interface of upperlying CL.  !
       ! Thus, when 'kinv(i)=1', PBL scheme guarantees 'bflxs <= 0'.  For this case !
       ! it is reasonable to assume that cumulus convection does not happen.     !
       ! When these is SBCL, PBL height from the PBL scheme is likely to be very !
       ! close at 'kinv(i)-1' interface, but not exactly, since 'zi' information is !
       ! changed between two model time steps. In order to ensure correct identi !
       ! fication of 'kinv(i)' for general case including SBCL, I imposed an offset !
       ! of 5 [m] in the below 'kinv(i)' finding block.                             !
       ! ----------------------------------------------------------------------- !

       do k = mkx - 1, 1, -1
          if( (pblh(i) + 5._r8 - zs0(k, i))*(pblh(i) + 5._r8 - zs0(k+1, i)) .lt. 0._r8 ) then
               kinv(i) = k + 1
               go to 15
          endif
       end do
       kinv(i) = 1
15     continue

       if( kinv(i) .le. 1 ) then
           exit_kinv1(i) = 1._r8
           id_exit(i) = .true.
           go to 333
       endif
       ! From here, it must be 'kinv(i) >= 2'.

       ! -------------------------------------------------------------------------- !
       ! Find PBL averaged tke(:, i) ('tkeavg(i)') and minimum 'thvl' ('thvlmin(i)') in the PBL !
       ! In the current code, 'tkeavg(i)' is obtained by averaging all interfacial TKE !
       ! within the PBL. However, in order to be conceptually consistent with   PBL !
       ! scheme, 'tkeavg(i)' should be calculated by considering surface buoyancy flux.!
       ! If surface buoyancy flux is positive ( bflxs >0 ), surface interfacial TKE !
       ! should be included in calculating 'tkeavg(i)', while if bflxs <= 0,   surface !
       ! interfacial TKE should not be included in calculating 'tkeavg(i)'.   I should !
       ! modify the code when 'bflxs' is available as an input of cumulus scheme.   !
       ! 'thvlmin(i)' is a minimum 'thvl' within PBL obtained by comparing top &  base !
       ! interface values of 'thvl' in each layers within the PBL.                  !
       ! -------------------------------------------------------------------------- !

       dpsum(i)    = 0._r8
       tkeavg(i)   = 0._r8
       thvlmin(i)  = 1000._r8
       do k = 0, kinv(i) - 1   ! Here, 'k' is an interfacial layer index.
          if( k .eq. 0 ) then
              dpi(i) = ps0(0, i) - p0(1, i)
          elseif( k .eq. (kinv(i)-1) ) then
              dpi(i) = p0(kinv(i)-1, i) - ps0(kinv(i)-1, i)
          else
              dpi(i) = p0(k, i) - p0(k+1, i)
          endif
          dpsum(i)  = dpsum(i)  + dpi(i)
          tkeavg(i) = tkeavg(i) + dpi(i)*tke(k, i)
          if( k .ne. 0 ) thvlmin(i) = min(thvlmin(i),min(thvl0bot(k, i),thvl0top(k, i)))
       end do
       tkeavg(i)  = tkeavg(i)/dpsum(i)

       ! ------------------------------------------------------------------ !
       ! Find characteristics of cumulus source air: qtsrc(i),thlsrc(i),usrc(i),vsrc(i) !
       ! Note that 'thlsrc(i)' was con-cocked using 'thvlsrc(i)' and 'qtsrc(i)'.     !
       ! 'qtsrc(i)' is defined as the lowest layer mid-point value;   'thlsrc(i)' !
       ! is from 'qtsrc(i)' and 'thvlmin(i)=thvlsrc(i)'; 'usrc(i)' & 'vsrc(i)' are defined !
       ! as the values just below the PBL top interface.                    !
       ! ------------------------------------------------------------------ !

       qtsrc(i)   = qt0(1, i)
       thvlsrc(i) = thvlmin(i)
       thlsrc(i)  = thvlsrc(i) / ( 1._r8 + zvir * qtsrc(i) )
       usrc(i)    = u0(kinv(i)-1, i) + ssu0(kinv(i)-1, i) * ( ps0(kinv(i)-1, i) - p0(kinv(i)-1, i) )
       vsrc(i)    = v0(kinv(i)-1, i) + ssv0(kinv(i)-1, i) * ( ps0(kinv(i)-1, i) - p0(kinv(i)-1, i) )
       do m = 1, ncnst
          trsrc(m, i) = tr0(1,m, i)
       enddo

       ! ------------------------------------------------------------------ !
       ! Find LCL of the source air and a layer index containing LCL (klcl(i)) !
       ! When the LCL is exactly at the interface, 'klcl(i)' is a layer index  !
       ! having 'plcl(i)' as the lower interface similar to the 'kinv(i)' case.   !
       ! In the previous code, I assumed that if LCL is located within the  !
       ! lowest model layer ( 1 ) or the top model layer ( mkx ), then  no  !
       ! convective adjustment is performed and just exited.   However, in  !
       ! the revised code, I relaxed the first constraint and  even though  !
       ! LCL is at the lowest model layer, I allowed cumulus convection to  !
       ! be initiated. For this case, cumulus convection should be started  !
       ! from the PBL top height, as shown in the following code.           !
       ! When source air is already saturated even at the surface, klcl(i) is  !
       ! set to 1.                                                          !
       ! ------------------------------------------------------------------ !

       plcl(i) = qsinvert(qtsrc(i),thlsrc(i),ps0(0, i))
       do k = 0, mkx
          if( ps0(k, i) .lt. plcl(i) ) then
              klcl(i) = k
              go to 25
          endif
       end do
       klcl(i) = mkx
25     continue
       klcl(i) = max(1,klcl(i))

       if( plcl(i) .lt. 30000._r8 ) then
     ! if( klcl(i) .eq. mkx ) then
           exit_klclmkx(i) = 1._r8
           id_exit(i) = .true.
           go to 333
       endif

       ! ------------------------------------------------------------- !
       ! Calculate environmental virtual potential temperature at LCL, !
       !'thv0lcl(i)' which is solely used in the 'cin(i)' calculation. Note  !
       ! that 'thv0lcl(i)' is calculated first by calculating  'thl0lcl(i)'  !
       ! and 'qt0lcl(i)' at the LCL, and performing 'conden' afterward,   !
       ! in fully consistent with the other parts of the code.         !
       ! ------------------------------------------------------------- !

       thl0lcl(i) = thl0(klcl(i), i) + ssthl0(klcl(i), i) * ( plcl(i) - p0(klcl(i), i) )
       qt0lcl(i)  = qt0(klcl(i), i)  + ssqt0(klcl(i), i)  * ( plcl(i) - p0(klcl(i), i) )
       call conden(plcl(i),thl0lcl(i),qt0lcl(i),thj(i),qvj(i),qlj(i),qij(i),qse(i),id_check)
       if( id_check .eq. 1 ) then
           exit_conden(i) = 1._r8
           id_exit(i) = .true.
           go to 333
       end if
       thv0lcl(i) = thj(i) * ( 1._r8 + zvir * qvj(i) - qlj(i) - qij(i) )

       ! ------------------------------------------------------------------------ !
       ! Compute Convective Inhibition, 'cin(i)' & 'cinlcl(i)' [J/kg]=[m2/s2] TKE unit. !
       !                                                                          !
       ! 'cin(i)' (cinlcl(i)) is computed from the PBL top interface to LFC (LCL) using !
       ! piecewisely reconstructed environmental profiles, assuming environmental !
       ! buoyancy profile within each layer ( or from LCL to upper interface in   !
       ! each layer ) is simply a linear profile. For the purpose of cin(i) (cinlcl(i)) !
       ! calculation, we simply assume that lateral entrainment does not occur in !
       ! updrafting cumulus plume, i.e., cumulus source air property is conserved.!
       ! Below explains some rules used in the calculations of cin(i) (cinlcl(i)).   In !
       ! general, both 'cin(i)' and 'cinlcl(i)' are calculated from a PBL top interface !
       ! to LCL and LFC, respectively :                                           !
       ! 1. If LCL is lower than the PBL height, cinlcl(i) = 0 and cin(i) is calculated !
       !    from PBL height to LFC.                                               !
       ! 2. If LCL is higher than PBL height,   'cinlcl(i)' is calculated by summing !
       !    both positive and negative cloud buoyancy up to LCL using 'single_cin'!
       !    From the LCL to LFC, however, only negative cloud buoyancy is counted !
       !    to calculate final 'cin(i)' upto LFC.                                    !
       ! 3. If either 'cin(i)' or 'cinlcl(i)' is negative, they are set to be zero.     !
       ! In the below code, 'klfc(i)' is the layer index containing 'LFC' similar to !
       ! 'kinv(i)' and 'klcl(i)'.                                                       !
       ! ------------------------------------------------------------------------ !

        cin(i)    = 0._r8
        cinlcl(i) = 0._r8
        plfc(i)   = 0._r8
        klfc(i)   = mkx

        ! ------------------------------------------------------------------------- !
        ! Case 1. LCL height is higher than PBL interface ( 'pLCL <= ps0(kinv(i)-1, i)' ) !
        ! ------------------------------------------------------------------------- !

        if( klcl(i) .ge. kinv(i) ) then

            do k = kinv(i), mkx - 1
               if( k .lt. klcl(i) ) then
                   thvubot(i) = thvlsrc(i)
                   thvutop(i) = thvlsrc(i)
                   cin(i)     = cin(i) + single_cin(ps0(k-1, i),thv0bot(k, i),ps0(k, i),thv0top(k, i),thvubot(i),thvutop(i))
               elseif( k .eq. klcl(i) ) then
                   !----- Bottom to LCL
                   thvubot(i) = thvlsrc(i)
                   thvutop(i) = thvlsrc(i)
                   cin(i)     = cin(i) + single_cin(ps0(k-1, i),thv0bot(k, i),plcl(i),thv0lcl(i),thvubot(i),thvutop(i))
                   if( cin(i) .lt. 0._r8 ) limit_cinlcl(i) = 1._r8
                   cinlcl(i)  = max(cin(i),0._r8)
                   cin(i)     = cinlcl(i)
                   !----- LCL to Top
                   thvubot(i) = thvlsrc(i)
                   call conden(ps0(k, i),thlsrc(i),qtsrc(i),thj(i),qvj(i),qlj(i),qij(i),qse(i),id_check)
                   if( id_check .eq. 1 ) then
                       exit_conden(i) = 1._r8
                       id_exit(i) = .true.
                       go to 333
                   end if
                   thvutop(i) = thj(i) * ( 1._r8 + zvir*qvj(i) - qlj(i) - qij(i) )
                   call getbuoy(plcl(i),thv0lcl(i),ps0(k, i),thv0top(k, i),thvubot(i),thvutop(i),plfc(i),cin(i))
                   if( plfc(i) .gt. 0._r8 ) then
                       klfc(i) = k
                       go to 35
                   end if
               else
                   thvubot(i) = thvutop(i)
                   call conden(ps0(k, i),thlsrc(i),qtsrc(i),thj(i),qvj(i),qlj(i),qij(i),qse(i),id_check)
                   if( id_check .eq. 1 ) then
                       exit_conden(i) = 1._r8
                       id_exit(i) = .true.
                       go to 333
                   end if
                   thvutop(i) = thj(i) * ( 1._r8 + zvir*qvj(i) - qlj(i) - qij(i) )
                   call getbuoy(ps0(k-1, i),thv0bot(k, i),ps0(k, i),thv0top(k, i),thvubot(i),thvutop(i),plfc(i),cin(i))
                   if( plfc(i) .gt. 0._r8 ) then
                       klfc(i) = k
                       go to 35
                   end if
               endif
            end do

       ! ----------------------------------------------------------------------- !
       ! Case 2. LCL height is lower than PBL interface ( 'pLCL > ps0(kinv(i)-1, i)' ) !
       ! ----------------------------------------------------------------------- !

       else
          cinlcl(i) = 0._r8
          do k = kinv(i), mkx - 1
             call conden(ps0(k-1, i),thlsrc(i),qtsrc(i),thj(i),qvj(i),qlj(i),qij(i),qse(i),id_check)
             if( id_check .eq. 1 ) then
                 exit_conden(i) = 1._r8
                 id_exit(i) = .true.
                 go to 333
             end if
             thvubot(i) = thj(i) * ( 1._r8 + zvir*qvj(i) - qlj(i) - qij(i) )
             call conden(ps0(k, i),thlsrc(i),qtsrc(i),thj(i),qvj(i),qlj(i),qij(i),qse(i),id_check)
             if( id_check .eq. 1 ) then
                 exit_conden(i) = 1._r8
                 id_exit(i) = .true.
                 go to 333
             end if
             thvutop(i) = thj(i) * ( 1._r8 + zvir*qvj(i) - qlj(i) - qij(i) )
             call getbuoy(ps0(k-1, i),thv0bot(k, i),ps0(k, i),thv0top(k, i),thvubot(i),thvutop(i),plfc(i),cin(i))
             if( plfc(i) .gt. 0._r8 ) then
                 klfc(i) = k
                 go to 35
             end if
          end do
       endif  ! End of CIN case selection

 35    continue
       if( cin(i) .lt. 0._r8 ) limit_cin(i) = 1._r8
       cin(i) = max(0._r8,cin(i))
       if( klfc(i) .ge. mkx ) then
           klfc(i) = mkx
         ! write(iulog,*) 'klfc(i) >= mkx'
           exit_klfcmkx(i) = 1._r8
           id_exit(i) = .true.
           go to 333
       endif

       ! ---------------------------------------------------------------------- !
       ! In order to calculate implicit 'cin(i)' (or 'cinlcl(i)'), save the initially !
       ! calculated 'cin(i)' and 'cinlcl(i)', and other related variables. These will !
       ! be restored after calculating implicit CIN.                            !
       ! ---------------------------------------------------------------------- !

       if( iter .eq. 1 ) then
           cin_i(i)       = cin(i)
           cinlcl_i(i)    = cinlcl(i)
           ke(i)          = rbuoy / ( rkfre * tkeavg(i) + epsvarw )
           kinv_o(i)      = kinv(i)
           klcl_o(i)      = klcl(i)
           klfc_o(i)      = klfc(i)
           plcl_o(i)      = plcl(i)
           plfc_o(i)      = plfc(i)
           tkeavg_o(i)    = tkeavg(i)
           thvlmin_o(i)   = thvlmin(i)
           qtsrc_o(i)     = qtsrc(i)
           thvlsrc_o(i)   = thvlsrc(i)
           thlsrc_o(i)    = thlsrc(i)
           usrc_o(i)      = usrc(i)
           vsrc_o(i)      = vsrc(i)
           thv0lcl_o(i)   = thv0lcl(i)
           do m = 1, ncnst
              trsrc_o(m, i) = trsrc(m, i)
           enddo
       endif

     ! Modification : If I impose w = max(0.1_r8, w) up to the top interface of
     !                klfc(i), I should only use cinlfc.  That is, if I want to
     !                use cinlcl(i), I should not impose w = max(0.1_r8, w).
     !                Using cinlcl(i) is equivalent to treating only 'saturated'
     !                moist convection. Note that in this sense, I should keep
     !                the functionality of both cinlfc and cinlcl(i).
     !                However, the treatment of penetrative entrainment level becomes
     !                ambiguous if I choose 'cinlcl(i)'. Thus, the best option is to use
     !                'cinlfc'.

       ! -------------------------------------------------------------------------- !
       ! Calculate implicit 'cin(i)' by averaging initial and final cins.    Note that !
       ! implicit CIN is adopted only when cumulus convection stabilized the system,!
       ! i.e., only when 'del_CIN(i) >0'. If 'del_CIN(i)<=0', just use explicit CIN. Note !
       ! also that since 'cinlcl(i)' is set to zero whenever LCL is below the PBL top, !
       ! (see above CIN calculation part), the use of 'implicit CIN=cinlcl(i)'  is not !
       ! good. Thus, when using implicit CIN, always try to only use 'implicit CIN= !
       ! cin(i)', not 'implicit CIN=cinlcl(i)'. However, both 'CIN=cin(i)' and 'CIN=cinlcl(i)'  !
       ! are good when using explicit CIN.                                          !
       ! -------------------------------------------------------------------------- !

       if( iter .ne. 1 ) then

           cin_f(i) = cin(i)
           cinlcl_f(i) = cinlcl(i)
           if( use_CINcin ) then
               del_CIN(i) = cin_f(i) - cin_i(i)
           else
               del_CIN(i) = cinlcl_f(i) - cinlcl_i(i)
           endif

           if( del_CIN(i) .gt. 0._r8 ) then

               ! -------------------------------------------------------------- !
               ! Calculate implicit 'cin(i)' and 'cinlcl(i)'. Note that when we chose !
               ! to use 'implicit CIN = cin(i)', choose 'cinlcl(i) = cinlcl_i(i)' below: !
               ! because iterative CIN only aims to obtain implicit CIN,  once  !
               ! we obtained 'implicit CIN=cin(i)', it is good to use the original !
               ! profiles information for all the other variables after that.   !
               ! Note 'cinlcl(i)' will be explicitly used in calculating  'wlcl(i)' & !
               ! 'ufrclcl(i)' after calculating 'winv(i)' & 'ufrcinv(i)'  at the PBL top !
               ! interface later, after calculating 'cbmf(i)'.                     !
               ! -------------------------------------------------------------- !

               alpha(i) = compute_alpha( del_CIN(i), ke(i) )
               cin(i)   = cin_i(i) + alpha(i) * del_CIN(i)
               if( use_CINcin ) then
                   cinlcl(i) = cinlcl_i(i)
               else
                   cinlcl(i) = cinlcl_i(i) + alpha(i) * del_cinlcl(i)
               endif

               ! ----------------------------------------------------------------- !
               ! Restore the original values from the previous 'iter_cin' step (1) !
               ! to compute correct tendencies for (n+1) time step by implicit CIN !
               ! ----------------------------------------------------------------- !

               kinv(i)      = kinv_o(i)
               klcl(i)      = klcl_o(i)
               klfc(i)      = klfc_o(i)
               plcl(i)      = plcl_o(i)
               plfc(i)      = plfc_o(i)
               tkeavg(i)    = tkeavg_o(i)
               thvlmin(i)   = thvlmin_o(i)
               qtsrc(i)     = qtsrc_o(i)
               thvlsrc(i)   = thvlsrc_o(i)
               thlsrc(i)    = thlsrc_o(i)
               usrc(i)      = usrc_o(i)
               vsrc(i)      = vsrc_o(i)
               thv0lcl(i)   = thv0lcl_o(i)
               do m = 1, ncnst
                  trsrc(m, i) = trsrc_o(m, i)
               enddo

               qv0(:mkx, i)            = qv0_o(:mkx, i)
               ql0(:mkx, i)            = ql0_o(:mkx, i)
               qi0(:mkx, i)            = qi0_o(:mkx, i)
               t0(:mkx, i)             = t0_o(:mkx, i)
               s0(:mkx, i)             = s0_o(:mkx, i)
               u0(:mkx, i)             = u0_o(:mkx, i)
               v0(:mkx, i)             = v0_o(:mkx, i)
               qt0(:mkx, i)            = qt0_o(:mkx, i)
               thl0(:mkx, i)           = thl0_o(:mkx, i)
               thvl0(:mkx, i)          = thvl0_o(:mkx, i)
               ssthl0(:mkx, i)         = ssthl0_o(:mkx, i)
               ssqt0(:mkx, i)          = ssqt0_o(:mkx, i)
               thv0bot(:mkx, i)        = thv0bot_o(:mkx, i)
               thv0top(:mkx, i)        = thv0top_o(:mkx, i)
               thvl0bot(:mkx, i)       = thvl0bot_o(:mkx, i)
               thvl0top(:mkx, i)       = thvl0top_o(:mkx, i)
               ssu0(:mkx, i)           = ssu0_o(:mkx, i)
               ssv0(:mkx, i)           = ssv0_o(:mkx, i)
               do m = 1, ncnst
                  tr0(:mkx,m, i)   = tr0_o(:mkx,m, i)
                  sstr0(:mkx,m, i) = sstr0_o(:mkx,m, i)
               enddo

               ! ------------------------------------------------------ !
               ! Initialize all fluxes, tendencies, and other variables !
               ! in association with cumulus convection.                !
               ! ------------------------------------------------------ !

               umf(0:mkx, i)          = 0.0_r8
               emf(0:mkx, i)          = 0.0_r8
               slflx(0:mkx, i)        = 0.0_r8
               qtflx(0:mkx, i)        = 0.0_r8
               uflx(0:mkx, i)         = 0.0_r8
               vflx(0:mkx, i)         = 0.0_r8
               qvten(:mkx, i)         = 0.0_r8
               qlten(:mkx, i)         = 0.0_r8
               qiten(:mkx, i)         = 0.0_r8
               sten(:mkx, i)          = 0.0_r8
               uten(:mkx, i)          = 0.0_r8
               vten(:mkx, i)          = 0.0_r8
               qrten(:mkx, i)         = 0.0_r8
               qsten(:mkx, i)         = 0.0_r8
               dwten(:mkx, i)         = 0.0_r8
               diten(:mkx, i)         = 0.0_r8
               precip(i)              = 0.0_r8
               snow(i)                = 0.0_r8
               evapc(:mkx, i)         = 0.0_r8
               cufrc(:mkx, i)         = 0.0_r8
               qcu(:mkx, i)           = 0.0_r8
               qlu(:mkx, i)           = 0.0_r8
               qiu(:mkx, i)           = 0.0_r8
               fer(:mkx, i)           = 0.0_r8
               fdr(:mkx, i)           = 0.0_r8
               qc(:mkx, i)            = 0.0_r8
               qc_l(:mkx, i)          = 0.0_r8
               qc_i(:mkx, i)          = 0.0_r8
               rliq(i)                = 0.0_r8
               cbmf(i)                = 0.0_r8
               cnt(i)                 = real(mkx, r8)
               cnb(i)                 = 0.0_r8
               qtten(:mkx, i)         = 0.0_r8
               slten(:mkx, i)         = 0.0_r8
               ufrc(0:mkx, i)         = 0.0_r8

               thlu(0:mkx, i)         = 0.0_r8
               qtu(0:mkx, i)          = 0.0_r8
               uu(0:mkx, i)           = 0.0_r8
               vu(0:mkx, i)           = 0.0_r8
               wu(0:mkx, i)           = 0.0_r8
               thvu(0:mkx, i)         = 0.0_r8
               thlu_emf(0:mkx, i)     = 0.0_r8
               qtu_emf(0:mkx, i)      = 0.0_r8
               uu_emf(0:mkx, i)       = 0.0_r8
               vu_emf(0:mkx, i)       = 0.0_r8

               do m = 1, ncnst
                  trflx(0:mkx,m, i)   = 0.0_r8
                  trten(:mkx,m, i)    = 0.0_r8
                  tru(0:mkx,m, i)     = 0.0_r8
                  tru_emf(0:mkx,m, i) = 0.0_r8
               enddo

               ! -------------------------------------------------- !
               ! Below are diagnostic output variables for detailed !
               ! analysis of cumulus scheme.                        !
               ! -------------------------------------------------- !

               ufrcinvbase(i)         = 0.0_r8
               ufrclcl(i)             = 0.0_r8
               winvbase(i)            = 0.0_r8
               wlcl(i)                = 0.0_r8
               emfkbup(i)             = 0.0_r8
               cbmflimit(i)           = 0.0_r8
               excessu_arr(:mkx, i)   = 0.0_r8
               excess0_arr(:mkx, i)   = 0.0_r8
               xc_arr(:mkx, i)        = 0.0_r8
               aquad_arr(:mkx, i)     = 0.0_r8
               bquad_arr(:mkx, i)     = 0.0_r8
               cquad_arr(:mkx, i)     = 0.0_r8
               bogbot_arr(:mkx, i)    = 0.0_r8
               bogtop_arr(:mkx, i)    = 0.0_r8

          else ! When 'del_CIN(i) < 0', use explicit CIN instead of implicit CIN.

               ! ----------------------------------------------------------- !
               ! Identifier showing whether explicit or implicit CIN is used !
               ! ----------------------------------------------------------- !

               ind_delcin(i) = 1._r8

               ! --------------------------------------------------------- !
               ! Restore original output values of "iter_cin = 1" and exit !
               ! --------------------------------------------------------- !

               umf_out(i,0:mkx)         = umf_s(0:mkx, i)
               qvten_out(i,:mkx)        = qvten_s(:mkx, i)
               qlten_out(i,:mkx)        = qlten_s(:mkx, i)
               qiten_out(i,:mkx)        = qiten_s(:mkx, i)
               sten_out(i,:mkx)         = sten_s(:mkx, i)
               uten_out(i,:mkx)         = uten_s(:mkx, i)
               vten_out(i,:mkx)         = vten_s(:mkx, i)
               qrten_out(i,:mkx)        = qrten_s(:mkx, i)
               qsten_out(i,:mkx)        = qsten_s(:mkx, i)
               precip_out(i)            = precip_s(i)
               snow_out(i)              = snow_s(i)
               evapc_out(i,:mkx)        = evapc_s(:mkx, i)
               cush_inout(i)            = cush_s(i)
               cufrc_out(i,:mkx)        = cufrc_s(:mkx, i)
               slflx_out(i,0:mkx)       = slflx_s(0:mkx, i)
               qtflx_out(i,0:mkx)       = qtflx_s(0:mkx, i)
               qcu_out(i,:mkx)          = qcu_s(:mkx, i)
               qlu_out(i,:mkx)          = qlu_s(:mkx, i)
               qiu_out(i,:mkx)          = qiu_s(:mkx, i)
               cbmf_out(i)              = cbmf_s(i)
               qc_out(i,:mkx)           = qc_s(:mkx, i)
               rliq_out(i)              = rliq_s(i)
               cnt_out(i)               = cnt_s(i)
               cnb_out(i)               = cnb_s(i)
               do m = 1, ncnst
                  trten_out(i,:mkx,m)   = trten_s(:mkx,m, i)
               enddo

               ! ------------------------------------------------------------------------------ !
               ! Below are diagnostic output variables for detailed analysis of cumulus scheme. !
               ! The order of vertical index is reversed for this internal diagnostic output.   !
               ! ------------------------------------------------------------------------------ !

               fer_out(i,mkx:1:-1)      = fer_s(:mkx, i)
               fdr_out(i,mkx:1:-1)      = fdr_s(:mkx, i)
               cinh_out(i)              = cin_s(i)
               cinlclh_out(i)           = cinlcl_s(i)
               qtten_out(i,mkx:1:-1)    = qtten_s(:mkx, i)
               slten_out(i,mkx:1:-1)    = slten_s(:mkx, i)
               ufrc_out(i,mkx:0:-1)     = ufrc_s(0:mkx, i)
               uflx_out(i,mkx:0:-1)     = uflx_s(0:mkx, i)
               vflx_out(i,mkx:0:-1)     = vflx_s(0:mkx, i)

               ufrcinvbase_out(i)       = ufrcinvbase_s(i)
               ufrclcl_out(i)           = ufrclcl_s(i)
               winvbase_out(i)          = winvbase_s(i)
               wlcl_out(i)              = wlcl_s(i)
               plcl_out(i)              = plcl_s(i)
               pinv_out(i)              = pinv_s(i)
               plfc_out(i)              = plfc_s(i)
               pbup_out(i)              = pbup_s(i)
               ppen_out(i)              = ppen_s(i)
               qtsrc_out(i)             = qtsrc_s(i)
               thlsrc_out(i)            = thlsrc_s(i)
               thvlsrc_out(i)           = thvlsrc_s(i)
               emfkbup_out(i)           = emfkbup_s(i)
               cbmflimit_out(i)         = cbmflimit_s(i)
               tkeavg_out(i)            = tkeavg_s(i)
               zinv_out(i)              = zinv_s(i)
               rcwp_out(i)              = rcwp_s(i)
               rlwp_out(i)              = rlwp_s(i)
               riwp_out(i)              = riwp_s(i)

               wu_out(i,mkx:0:-1)       = wu_s(0:mkx, i)
               qtu_out(i,mkx:0:-1)      = qtu_s(0:mkx, i)
               thlu_out(i,mkx:0:-1)     = thlu_s(0:mkx, i)
               thvu_out(i,mkx:0:-1)     = thvu_s(0:mkx, i)
               uu_out(i,mkx:0:-1)       = uu_s(0:mkx, i)
               vu_out(i,mkx:0:-1)       = vu_s(0:mkx, i)
               qtu_emf_out(i,mkx:0:-1)  = qtu_emf_s(0:mkx, i)
               thlu_emf_out(i,mkx:0:-1) = thlu_emf_s(0:mkx, i)
               uu_emf_out(i,mkx:0:-1)   = uu_emf_s(0:mkx, i)
               vu_emf_out(i,mkx:0:-1)   = vu_emf_s(0:mkx, i)
               uemf_out(i,mkx:0:-1)     = uemf_s(0:mkx, i)

               dwten_out(i,mkx:1:-1)    = dwten_s(:mkx, i)
               diten_out(i,mkx:1:-1)    = diten_s(:mkx, i)
               flxrain_out(i,mkx:0:-1)  = flxrain_s(0:mkx, i)
               flxsnow_out(i,mkx:0:-1)  = flxsnow_s(0:mkx, i)
               ntraprd_out(i,mkx:1:-1)  = ntraprd_s(:mkx, i)
               ntsnprd_out(i,mkx:1:-1)  = ntsnprd_s(:mkx, i)

               excessu_arr_out(i,mkx:1:-1)  = excessu_arr_s(:mkx, i)
               excess0_arr_out(i,mkx:1:-1)  = excess0_arr_s(:mkx, i)
               xc_arr_out(i,mkx:1:-1)       = xc_arr_s(:mkx, i)
               aquad_arr_out(i,mkx:1:-1)    = aquad_arr_s(:mkx, i)
               bquad_arr_out(i,mkx:1:-1)    = bquad_arr_s(:mkx, i)
               cquad_arr_out(i,mkx:1:-1)    = cquad_arr_s(:mkx, i)
               bogbot_arr_out(i,mkx:1:-1)   = bogbot_arr_s(:mkx, i)
               bogtop_arr_out(i,mkx:1:-1)   = bogtop_arr_s(:mkx, i)

               do m = 1, ncnst
                  trflx_out(i,mkx:0:-1,m)   = trflx_s(0:mkx,m, i)
                  tru_out(i,mkx:0:-1,m)     = tru_s(0:mkx,m, i)
                  tru_emf_out(i,mkx:0:-1,m) = tru_emf_s(0:mkx,m, i)
               enddo

               id_exit(i) = .false.
               go to 333

          endif

       endif

       ! ------------------------------------------------------------------ !
       ! Define a release level, 'prel(i)' and release layer, 'krel(i)'.          !
       ! 'prel(i)' is the lowest level from which buoyancy sorting occurs, and !
       ! 'krel(i)' is the layer index containing 'prel(i)' in it, similar to  the !
       ! previous definitions of 'kinv(i)', 'klcl(i)', and 'klfc(i)'.    In order to !
       ! ensure that only PBL scheme works within the PBL,  if LCL is below !
       ! PBL top height, then 'krel(i) = kinv(i)', while if LCL is above  PBL top !
       ! height, then 'krel(i) = klcl(i)'.   Note however that regardless of  the !
       ! definition of 'krel(i)', cumulus convection induces fluxes within PBL !
       ! through 'fluxbelowinv'.  We can make cumulus convection start from !
       ! any level, even within the PBL by appropriately defining 'krel(i)'  & !
       ! 'prel(i)' here. Then it must be accompanied by appropriate definition !
       ! of source air properties, CIN, and re-setting of 'fluxbelowinv', & !
       ! many other stuffs.                                                 !
       ! Note that even when 'prel(i)' is located above the PBL top height, we !
       ! still have cumulus convection between PBL top height and 'prel(i)':   !
       ! we simply assume that no lateral mixing occurs in this range.      !
       ! ------------------------------------------------------------------ !

       if( klcl(i) .lt. kinv(i) ) then
           krel(i)    = kinv(i)
           prel(i)    = ps0(krel(i)-1, i)
           thv0rel(i) = thv0bot(krel(i), i)
       else
           krel(i)    = klcl(i)
           prel(i)    = plcl(i)
           thv0rel(i) = thv0lcl(i)
       endif

       ! --------------------------------------------------------------------------- !
       ! Calculate cumulus base mass flux ('cbmf(i)'), fractional area ('ufrcinv(i)'), and !
       ! and mean vertical velocity (winv(i)) of cumulus updraft at PBL top interface.  !
       ! Also, calculate updraft fractional area (ufrclcl(i)) and vertical velocity  at !
       ! the LCL (wlcl(i)). When LCL is below PBLH, cinlcl(i) = 0 and 'ufrclcl(i) = ufrcinv(i)', !
       ! and 'wlcl(i) = winv(i).                                                           !
       ! Only updrafts strong enough to overcome CIN can rise over PBL top interface.!
       ! Thus,  in order to calculate cumulus mass flux at PBL top interface, 'cbmf(i)',!
       ! we need to know 'CIN' ( the strength of potential energy barrier ) and      !
       ! 'sigmaw(i)' ( a standard deviation of updraft vertical velocity at the PBL top !
       ! interface, a measure of turbulentce strength in the PBL ).   Naturally, the !
       ! ratio of these two variables, 'mu(i)' - normalized CIN by TKE- is key variable !
       ! controlling 'cbmf(i)'.  If 'mu(i)' becomes large, only small fraction of updrafts !
       ! with very strong TKE can rise over the PBL - both 'cbmf(i)' and 'ufrc(:, i)' becomes !
       ! small, but 'winv(i)' becomes large ( this can be easily understood by PDF of w !
       ! at PBL top ).  If 'mu(i)' becomes small, lots of updraft can rise over the PBL !
       ! top - both 'cbmf(i)' and 'ufrc(:, i)' becomes large, but 'winv(i)' becomes small. Thus, !
       ! all of the key variables associated with cumulus convection  at the PBL top !
       ! - 'cbmf(i)', 'ufrc(:, i)', 'winv(i)' where 'cbmf(i) = rho*ufrc(:, i)*winv(i)' - are a unique functi !
       ! ons of 'mu(i)', normalized CIN. Although these are uniquely determined by 'mu(i)',!
       ! we usually impose two comstraints on 'cbmf(i)' and 'ufrc(:, i)': (1) because we will !
       ! simply assume that subsidence warming and drying of 'kinv(i)-1' layer in assoc !
       ! iation with 'cbmf(i)' at PBL top interface is confined only in 'kinv(i)-1' layer, !
       ! cbmf(i) must not be larger than the mass within the 'kinv(i)-1' layer. Otherwise, !
       ! instability will occur due to the breaking of stability con. If we consider !
       ! semi-Lagrangian vertical advection scheme and explicitly consider the exten !
       ! t of vertical movement of each layer in association with cumulus mass flux, !
       ! we don't need to impose this constraint. However,  using a  semi-Lagrangian !
       ! scheme is a future research subject. Note that this constraint should be ap !
       ! plied for all interfaces above PBL top as well as PBL top interface.   As a !
       ! result, this 'cbmf(i)' constraint impose a 'lower' limit on mu(i) - 'mumin0(i)'. (2) !
       ! in order for mass flux parameterization - rho*(w'a')= M*(a_c-a_e) - to   be !
       ! valid, cumulus updraft fractional area should be much smaller than 1.    In !
       ! current code, we impose 'rmaxfrac = 0.1 ~ 0.2'   through the whole vertical !
       ! layers where cumulus convection occurs. At the PBL top interface,  the same !
       ! constraint is made by imposing another lower 'lower' limit on mu(i), 'mumin1'. !
       ! After that, also limit 'ufrclcl(i)' to be smaller than 'rmaxfrac' by 'mumin2(i)'. !
       ! --------------------------------------------------------------------------- !

       ! --------------------------------------------------------------------------- !
       ! Calculate normalized CIN, 'mu(i)' satisfying all the three constraints imposed !
       ! on 'cbmf(i)'('mumin0(i)'), 'ufrc(:, i)' at the PBL top - 'ufrcinv(i)' - ( by 'mumin1' from !
       ! a parameter sentence), and 'ufrc(:, i)' at the LCL - 'ufrclcl(i)' ( by 'mumin2(i)').    !
       ! Note that 'cbmf(i)' does not change between PBL top and LCL  because we assume !
       ! that buoyancy sorting does not occur when cumulus updraft is unsaturated.   !
       ! --------------------------------------------------------------------------- !

       if( use_CINcin ) then
           wcrit(i) = sqrt( 2._r8 * cin(i) * rbuoy )
       else
           wcrit(i) = sqrt( 2._r8 * cinlcl(i) * rbuoy )
       endif
       sigmaw(i) = sqrt( rkfre * tkeavg(i) + epsvarw )
       mu(i) = wcrit(i)/sigmaw(i)/1.4142_r8
       if( mu(i) .ge. 3._r8 ) then
         ! write(iulog,*) 'mu(i) >= 3'
           id_exit(i) = .true.
           go to 333
       endif
       rho0inv(i) = ps0(kinv(i)-1, i)/(r*thv0top(kinv(i)-1, i)*exns0(kinv(i)-1, i))
       cbmf(i) = (rho0inv(i)*sigmaw(i)/2.5066_r8)*exp(-mu(i)**2)
       ! 1. 'cbmf(i)' constraint
       cbmflimit(i) = 0.9_r8*dp0(kinv(i)-1, i)/g/dt
       mumin0(i) = 0._r8
       if( cbmf(i) .gt. cbmflimit(i) ) mumin0(i) = sqrt(-log(2.5066_r8*cbmflimit(i)/rho0inv(i)/sigmaw(i)))
       ! 2. 'ufrcinv(i)' constraint
       mu(i) = max(max(mu(i),mumin0(i)),mumin1)
       ! 3. 'ufrclcl(i)' constraint
       mulcl(i) = sqrt(2._r8*cinlcl(i)*rbuoy)/1.4142_r8/sigmaw(i)
       mulclstar(i) = sqrt(max(0._r8,2._r8*(exp(-mu(i)**2)/2.5066_r8)**2*(1._r8/erfc(mu(i))**2-0.25_r8/rmaxfrac**2)))
       if( mulcl(i) .gt. 1.e-8_r8 .and. mulcl(i) .gt. mulclstar(i) ) then
           mumin2(i) = compute_mumin2(mulcl(i),rmaxfrac,mu(i))
           if( mu(i) .gt. mumin2(i) ) then
               write(iulog,*) 'Critical error in mu(i) calculation in UW_ShCu'
               call endrun
           endif
           mu(i) = max(mu(i),mumin2(i))
           if( mu(i) .eq. mumin2(i) ) limit_ufrc(i) = 1._r8
       endif
       if( mu(i) .eq. mumin0(i) ) limit_cbmf(i) = 1._r8
       if( mu(i) .eq. mumin1 ) limit_ufrc(i) = 1._r8

       ! ------------------------------------------------------------------- !
       ! Calculate final ['cbmf(i)','ufrcinv(i)','winv(i)'] at the PBL top interface. !
       ! Note that final 'cbmf(i)' here is obtained in such that 'ufrcinv(i)' and  !
       ! 'ufrclcl(i)' are smaller than ufrcmax with no instability.             !
       ! ------------------------------------------------------------------- !

       cbmf(i) = (rho0inv(i)*sigmaw(i)/2.5066_r8)*exp(-mu(i)**2)
       winv(i) = sigmaw(i)*(2._r8/2.5066_r8)*exp(-mu(i)**2)/erfc(mu(i))
       ufrcinv(i) = cbmf(i)/winv(i)/rho0inv(i)

       ! ------------------------------------------------------------------- !
       ! Calculate ['ufrclcl(i)','wlcl(i)'] at the LCL. When LCL is below PBL top, !
       ! it automatically becomes 'ufrclcl(i) = ufrcinv(i)' & 'wlcl(i) = winv(i)', since !
       ! it was already set to 'cinlcl(i)=0' if LCL is below PBL top interface. !
       ! Note 'cbmf(i)' at the PBL top is the same as 'cbmf(i)' at the LCL.  Note  !
       ! also that final 'cbmf(i)' here is obtained in such that 'ufrcinv(i)' and  !
       ! 'ufrclcl(i)' are smaller than ufrcmax and there is no instability.     !
       ! By construction, it must be 'wlcl(i) > 0' but for assurance, I checked !
       ! this again in the below block. If 'ufrclcl(i) < 0.1%', just exit.      !
       ! ------------------------------------------------------------------- !

       wtw(i) = winv(i) * winv(i) - 2._r8 * cinlcl(i) * rbuoy
       if( wtw(i) .le. 0._r8 ) then
         ! write(iulog,*) 'wlcl(i) < 0 at the LCL'
           exit_wtw(i) = 1._r8
           id_exit(i) = .true.
           go to 333
       endif
       wlcl(i) = sqrt(wtw(i))
       ufrclcl(i) = cbmf(i)/wlcl(i)/rho0inv(i)
       wrel(i) = wlcl(i)
       if( ufrclcl(i) .le. 0.0001_r8 ) then
         ! write(iulog,*) 'ufrclcl(i) <= 0.0001'
           exit_ufrc(i) = 1._r8
           id_exit(i) = .true.
           go to 333
       endif
       ufrc(krel(i)-1, i) = ufrclcl(i)

       ! ----------------------------------------------------------------------- !
       ! Below is just diagnostic output for detailed analysis of cumulus scheme !
       ! ----------------------------------------------------------------------- !

       ufrcinvbase(i)        = ufrcinv(i)
       winvbase(i)           = winv(i)
       umf(kinv(i)-1:krel(i)-1, i) = cbmf(i)
       wu(kinv(i)-1:krel(i)-1, i)  = winv(i)

       ! -------------------------------------------------------------------------- !
       ! Define updraft properties at the level where buoyancy sorting starts to be !
       ! happening, i.e., by definition, at 'prel(i)' level within the release layer.  !
       ! Because no lateral entrainment occurs upto 'prel(i)', conservative scalars of !
       ! cumulus updraft at release level is same as those of source air.  However, !
       ! horizontal momentums of source air are modified by horizontal PGF forcings !
       ! from PBL top interface to 'prel(i)'.  For this case, we should add additional !
       ! horizontal momentum from PBL top interface to 'prel(i)' as will be done below !
       ! to 'usrc(i)' and 'vsrc(i)'. Note that below cumulus updraft properties - umf(:, i), wu(:, i),!
       ! thlu(:, i), qtu(:, i), thvu(:, i), uu(:, i), vu(:, i) - are defined all interfaces not at the layer mid- !
       ! point. From the index notation of cumulus scheme, wu(k, i) is the cumulus up- !
       ! draft vertical velocity at the top interface of k layer.                   !
       ! Diabatic horizontal momentum forcing should be treated as a kind of 'body' !
       ! forcing without actual mass exchange between convective updraft and        !
       ! environment, but still taking horizontal momentum from the environment to  !
       ! the convective updrafts. Thus, diabatic convective momentum transport      !
       ! vertically redistributes environmental horizontal momentum.                !
       ! -------------------------------------------------------------------------- !

       emf(krel(i)-1, i)  = 0._r8
       umf(krel(i)-1, i)  = cbmf(i)
       wu(krel(i)-1, i)   = wrel(i)
       thlu(krel(i)-1, i) = thlsrc(i)
       qtu(krel(i)-1, i)  = qtsrc(i)
       call conden(prel(i),thlsrc(i),qtsrc(i),thj(i),qvj(i),qlj(i),qij(i),qse(i),id_check)
       if( id_check .eq. 1 ) then
           exit_conden(i) = 1._r8
           id_exit(i) = .true.
           go to 333
       endif
       thvu(krel(i)-1, i) = thj(i) * ( 1._r8 + zvir*qvj(i) - qlj(i) - qij(i) )

       uplus(i) = 0._r8
       vplus(i) = 0._r8
       if( krel(i) .eq. kinv(i) ) then
           uplus(i) = PGFc * ssu0(kinv(i), i) * ( prel(i) - ps0(kinv(i)-1, i) )
           vplus(i) = PGFc * ssv0(kinv(i), i) * ( prel(i) - ps0(kinv(i)-1, i) )
       else
           do k = kinv(i), max(krel(i)-1,kinv(i))
              uplus(i) = uplus(i) + PGFc * ssu0(k, i) * ( ps0(k, i) - ps0(k-1, i) )
              vplus(i) = vplus(i) + PGFc * ssv0(k, i) * ( ps0(k, i) - ps0(k-1, i) )
           end do
           uplus(i) = uplus(i) + PGFc * ssu0(krel(i), i) * ( prel(i) - ps0(krel(i)-1, i) )
           vplus(i) = vplus(i) + PGFc * ssv0(krel(i), i) * ( prel(i) - ps0(krel(i)-1, i) )
       end if
       uu(krel(i)-1, i) = usrc(i) + uplus(i)
       vu(krel(i)-1, i) = vsrc(i) + vplus(i)

       do m = 1, ncnst
          tru(krel(i)-1,m, i)  = trsrc(m, i)
       enddo

       ! -------------------------------------------------------------------------- !
       ! Define environmental properties at the level where buoyancy sorting occurs !
       ! ('pe(i)', normally, layer midpoint except in the 'krel(i)' layer). In the 'krel(i)' !
       ! layer where buoyancy sorting starts to occur, however, 'pe(i)' is defined     !
       ! differently because LCL is regarded as lower interface for mixing purpose. !
       ! -------------------------------------------------------------------------- !

       pe(i)      = 0.5_r8 * ( prel(i) + ps0(krel(i), i) )
       dpe(i)     = prel(i) - ps0(krel(i), i)
       exne(i)    = exnf(pe(i))
       thvebot(i) = thv0rel(i)
       thle(i)    = thl0(krel(i), i) + ssthl0(krel(i), i) * ( pe(i) - p0(krel(i), i) )
       qte(i)     = qt0(krel(i), i)  + ssqt0(krel(i), i)  * ( pe(i) - p0(krel(i), i) )
       ue(i)      = u0(krel(i), i)   + ssu0(krel(i), i)   * ( pe(i) - p0(krel(i), i) )
       ve(i)      = v0(krel(i), i)   + ssv0(krel(i), i)   * ( pe(i) - p0(krel(i), i) )
       do m = 1, ncnst
          tre(m, i) = tr0(krel(i),m, i)  + sstr0(krel(i),m, i) * ( pe(i) - p0(krel(i), i) )
       enddo

       !-------------------------!
       ! Buoyancy-Sorting Mixing !
       !-------------------------!------------------------------------------------ !
       !                                                                           !
       !  In order to complete buoyancy-sorting mixing at layer mid-point, and so  !
       !  calculate 'updraft mass flux, updraft w velocity, conservative scalars'  !
       !  at the upper interface of each layer, we need following 3 information.   !
       !                                                                           !
       !  1. Pressure where mixing occurs ('pe(i)'), and temperature at 'pe(i)' which is !
       !     necessary to calculate various thermodynamic coefficients at pe(i). This !
       !     temperature is obtained by undiluted cumulus properties lifted to pe(i). !
       !  2. Undiluted updraft properties at pe(i) - conservative scalar and vertical !
       !     velocity -which are assumed to be the same as the properties at lower !
       !     interface only for calculation of fractional lateral entrainment  and !
       !     detrainment rate ( fer(k, i) and fdr(k, i) [Pa-1] ), respectively.    Final !
       !     values of cumulus conservative scalars and w at the top interface are !
       !     calculated afterward after obtaining fer(k, i) & fdr(k, i).                 !
       !  3. Environmental properties at pe(i).                                       !
       ! ------------------------------------------------------------------------- !

       ! ------------------------------------------------------------------------ !
       ! Define cumulus scale height.                                             !
       ! Cumulus scale height is defined as the maximum height cumulus can reach. !
       ! In case of premitive code, cumulus scale height ('cush(i)')  at the current !
       ! time step was assumed to be the same as 'cush(i)' of previous time step.    !
       ! However, I directly calculated cush(i) at each time step using an iterative !
       ! method. Note that within the cumulus scheme, 'cush(i)' information is  used !
       ! only at two places during buoyancy-sorting process:                      !
       ! (1) Even negatively buoyancy mixtures with strong vertical velocity      !
       !     enough to rise up to 'rle*scaleh(i)' (rle = 0.1) from pe(i) are entrained  !
       !     into cumulus updraft,                                                !
       ! (2) The amount of mass that is involved in buoyancy-sorting mixing       !
       !      process at pe(i) is rei(k, i) = rkm/scaleh(i)/rho*g [Pa-1]                   !
       ! In terms of (1), I think critical stopping distance might be replaced by !
       ! layer thickness. In future, we will use rei(k, i) = (0.5*rkm/z0(k, i)/rho/g).  !
       ! In the premitive code,  'scaleh(i)' was largely responsible for the jumping !
       ! variation of precipitation amount.                                       !
       ! ------------------------------------------------------------------------ !

       scaleh(i) = tscaleh(i)
       if( tscaleh(i) .lt. 0.0_r8 ) scaleh(i) = 1000._r8
!!!end part one
     ! Save time : Set iter_scaleh = 1. This will automatically use 'cush(i)' from the previous time step
     !             at the first implicit iteration. At the second implicit iteration, it will use
     !             the updated 'cush(i)' by the first implicit cin(i). So, this updating has an effect of
     !             doing one iteration for cush(i) calculation, which is good.
     !             So, only this setting of 'iter_scaleh = 1' is sufficient-enough to save computation time.
     ! OK
!!!end part one
       do iter_scaleh = 1, 3

       ! ---------------------------------------------------------------- !
       ! Initialization of 'kbup(i)' and 'kpen(i)'                              !
       ! ---------------------------------------------------------------- !
       ! 'kbup(i)' is the top-most layer in which cloud buoyancy is positive !
       ! both at the top and bottom interface of the layer. 'kpen(i)' is the !
       ! layer upto which cumulus panetrates ,i.e., cumulus w at the base !
       ! interface is positive, but becomes negative at the top interface.!
       ! Here, we initialize 'kbup(i)' and 'kpen(i)'. These initializations are !
       ! not trivial but important, expecially   in calculating turbulent !
       ! fluxes without confliction among several physics as explained in !
       ! detail in the part of turbulent fluxes calculation later.   Note !
       ! that regardless of whether 'kbup(i)' and 'kpen(i)' are updated or  not !
       ! during updraft motion,  penetrative entrainments are dumped down !
       ! across the top interface of 'kbup(i)' later.      More specifically,!
       ! penetrative entrainment heat and moisture fluxes are  calculated !
       ! from the top interface of 'kbup(i)' layer  to the base interface of !
       ! 'kpen(i)' layer. Because of this, initialization of 'kbup(i)' & 'kpen(i)' !
       ! influence the convection system when there are not updated.  The !
       ! below initialization of 'kbup(i) = krel(i)' assures  that  penetrative !
       ! entrainment fluxes always occur at interfaces above the PBL  top !
       ! interfaces (i.e., only at interfaces k >=kinv(i) ), which seems  to !
       ! be attractable considering that the most correct fluxes  at  the !
       ! PBL top interface can be ontained from the 'fluxbelowinv'  using !
       ! reconstructed PBL height.                                        !
       ! The 'kbup(i) = krel(i)'(after going through the whole buoyancy sorting !
       ! proces during updraft motion) implies that cumulus updraft  from !
       ! the PBL top interface can not reach to the LFC,so that 'kbup(i)' is !
       ! not updated during upward. This means that cumulus updraft   did !
       ! not fully overcome the buoyancy barrier above just the PBL top.  !
       ! If 'kpen(i)' is not updated either ( i.e., cumulus cannot rise over !
       ! the top interface of release layer),penetrative entrainment will !
       ! not happen at any interfaces.  If cumulus updraft can rise above !
       ! the release layer but cannot fully overcome the buoyancy barrier !
       ! just above PBL top interface, penetratve entrainment   occurs at !
       ! several above interfaces, including the top interface of release !
       ! layer. In the latter case, warming and drying tendencies will be !
       ! be initiated in 'krel(i)' layer. Note current choice of 'kbup(i)=krel(i)' !
       ! is completely compatible with other flux physics without  double !
       ! or miss counting turbulent fluxes at any interface. However, the !
       ! alternative choice of 'kbup(i)=krel(i)-1' also has itw own advantage - !
       ! when cumulus updraft cannot overcome buoyancy barrier just above !
       ! PBL top, entrainment warming and drying are concentrated in  the !
       ! 'kinv(i)-1' layer instead of 'kinv(i)' layer for this case. This might !
       ! seems to be more dynamically reasonable, but I will choose the   !
       ! 'kbup(i) = krel(i)' choice since it is more compatible  with the other !
       ! parts of the code, expecially, when we chose ' use_emf=.false. ' !
       ! as explained in detail in turbulent flux calculation part.       !
       ! ---------------------------------------------------------------- !

       kbup(i)    = krel(i)
       kpen(i)    = krel(i)

       ! ------------------------------------------------------------ !
       ! Since 'wtw(i)' is continuously updated during vertical motion,  !
       ! I need below initialization command within this 'iter_scaleh'!
       ! do loop. Similarily, I need initializations of environmental !
       ! properties at 'krel(i)' layer as below.                         !
       ! ------------------------------------------------------------ !

       wtw(i)     = wlcl(i) * wlcl(i)
       pe(i)      = 0.5_r8 * ( prel(i) + ps0(krel(i), i) )
       dpe(i)     = prel(i) - ps0(krel(i), i)
       exne(i)    = exnf(pe(i))
       thvebot(i) = thv0rel(i)
       thle(i)    = thl0(krel(i), i) + ssthl0(krel(i), i) * ( pe(i) - p0(krel(i), i) )
       qte(i)     = qt0(krel(i), i)  + ssqt0(krel(i), i)  * ( pe(i) - p0(krel(i), i) )
       ue(i)      = u0(krel(i), i)   + ssu0(krel(i), i)   * ( pe(i) - p0(krel(i), i) )
       ve(i)      = v0(krel(i), i)   + ssv0(krel(i), i)   * ( pe(i) - p0(krel(i), i) )
       do m = 1, ncnst
          tre(m, i) = tr0(krel(i),m, i)  + sstr0(krel(i),m, i)  * ( pe(i) - p0(krel(i), i) )
       enddo

       ! ----------------------------------------------------------------------- !
       ! Cumulus rises upward from 'prel(i)' ( or base interface of  'krel(i)' layer ) !
       ! until updraft vertical velocity becomes zero.                           !
       ! Buoyancy sorting is performed via two stages. (1) Using cumulus updraft !
       ! properties at the base interface of each layer,perform buoyancy sorting !
       ! at the layer mid-point, 'pe(i)',  and update cumulus properties at the top !
       ! interface, and then  (2) by averaging updated cumulus properties at the !
       ! top interface and cumulus properties at the base interface,   calculate !
       ! cumulus updraft properties at pe(i) that will be used  in buoyancy sorting !
       ! mixing - thlue(i), qtue(i) and, wue(i).  Using this averaged properties, perform !
       ! buoyancy sorting again at pe(i), and re-calculate fer(k, i) and fdr(k, i). Using !
       ! this recalculated fer(k, i) and fdr(k, i),  finally calculate cumulus updraft !
       ! properties at the top interface - thlu(:, i), qtu(:, i), thvu(:, i), uu(:, i), vu(:, i). In the below,!
       ! 'iter_xc = 1' performs the first stage, while 'iter_xc= 2' performs the !
       ! second stage. We can increase the number of iterations, 'nter_xc'.as we !
       ! want, but a sample test indicated that about 3 - 5 iterations  produced !
       ! satisfactory converent solution. Finally, identify 'kbup(i)' and 'kpen(i)'.   !
       ! ----------------------------------------------------------------------- !

       do k = krel(i), mkx - 1 ! Here, 'k' is a layer index.

          km1 = k - 1

          thlue(i) = thlu(km1, i)
          qtue(i)  = qtu(km1, i)
          wue(i)   = wu(km1, i)
          wtwb(i)  = wtw(i)

       do iter_xc = 1, niter_xc

          wtw(i) = wu(km1, i) * wu(km1, i)

          ! ---------------------------------------------------------------- !
          ! Calculate environmental and cumulus saturation 'excess' at 'pe(i)'. !
          ! Note that in order to calculate saturation excess, we should use !
          ! liquid water temperature instead of temperature  as the argument !
          ! of "qsat". But note normal argument of "qsat" is temperature.    !
          ! ---------------------------------------------------------------- !

          call conden(pe(i),thle(i),qte(i),thj(i),qvj(i),qlj(i),qij(i),qse(i),id_check)
          if( id_check .eq. 1 ) then
              exit_conden(i) = 1._r8
              id_exit(i) = .true.
              go to 333
          end if
          thv0j(i)    = thj(i) * ( 1._r8 + zvir*qvj(i) - qlj(i) - qij(i) )
          rho0j(i)    = pe(i) / ( r * thv0j(i) * exne(i) )
          qsat_arg(i) = thle(i)*exne(i)
          call t_startf('qsat_o')
!#define QSAT_WV
#ifdef QSAT_WV
          call qsat_wv(qsat_arg(i), pe(i), es(i), qs(i))
#else
          call qsat_o(qsat_arg(i), pe(i), es(i), qs(i))
#endif
          call t_stopf('qsat_o')
          excess0(i)  = qte(i) - qs(i)

          call conden(pe(i),thlue(i),qtue(i),thj(i),qvj(i),qlj(i),qij(i),qse(i),id_check)
          if( id_check .eq. 1 ) then
              exit_conden(i) = 1._r8
              id_exit(i) = .true.
              go to 333
          end if
          ! ----------------------------------------------------------------- !
          ! Detrain excessive condensate larger than 'criqc' from the cumulus !
          ! updraft before performing buoyancy sorting. All I should to do is !
          ! to update 'thlue(i)' &  'que' here. Below modification is completely !
          ! compatible with the other part of the code since 'thule' & 'qtue(i)' !
          ! are used only for buoyancy sorting. I found that as long as I use !
          ! 'niter_xc >= 2',  detraining excessive condensate before buoyancy !
          ! sorting has negligible influence on the buoyancy sorting results. !
          ! ----------------------------------------------------------------- !
          if( (qlj(i) + qij(i)) .gt. criqc ) then
               exql(i)  = ( ( qlj(i) + qij(i) ) - criqc ) * qlj(i) / ( qlj(i) + qij(i) )
               exqi(i)  = ( ( qlj(i) + qij(i) ) - criqc ) * qij(i) / ( qlj(i) + qij(i) )
               qtue(i)  = qtue(i) - exql(i) - exqi(i)
               thlue(i) = thlue(i) + (xlv/cp/exne(i))*exql(i) + (xls/cp/exne(i))*exqi(i)
          endif
          call conden(pe(i),thlue(i),qtue(i),thj(i),qvj(i),qlj(i),qij(i),qse(i),id_check)
          if( id_check .eq. 1 ) then
              exit_conden(i) = 1._r8
              id_exit(i) = .true.
              go to 333
          end if
          thvj(i)     = thj(i) * ( 1._r8 + zvir * qvj(i) - qlj(i) - qij(i) )
          tj(i)       = thj(i) * exne(i) ! This 'tj(i)' is used for computing thermo. coeffs. below
          qsat_arg(i) = thlue(i)*exne(i)
          call t_startf('qsat_o')
#ifdef QSAT_WV
          call qsat_wv(qsat_arg(i), pe(i), es(i), qs(i))
#else
          call qsat_o(qsat_arg(i), pe(i), es(i), qs(i))
#endif
          call t_stopf('qsat_o')
          excessu(i)  = qtue(i) - qs(i)

          ! ------------------------------------------------------------------- !
          ! Calculate critical mixing fraction, 'xc(i)'. Mixture with mixing ratio !
          ! smaller than 'xc(i)' will be entrained into cumulus updraft.  Both the !
          ! saturated updrafts with 'positive buoyancy' or 'negative buoyancy + !
          ! strong vertical velocity enough to rise certain threshold distance' !
          ! are kept into the updraft in the below program. If the core updraft !
          ! is unsaturated, we can set 'xc(i) = 0' and let the cumulus  convection !
          ! still works or we may exit.                                         !
          ! Current below code does not entrain unsaturated mixture. However it !
          ! should be modified such that it also entrain unsaturated mixture.   !
          ! ------------------------------------------------------------------- !

          ! ----------------------------------------------------------------- !
          ! cridis(i) : Critical stopping distance for buoyancy sorting purpose. !
          !          scaleh(i) is only used here.                                !
          ! ----------------------------------------------------------------- !

            cridis(i) = rle*scaleh(i)                 ! Original code
          ! cridis(i) = 1._r8*(zs0(k, i) - zs0(k-1, i))  ! New code

          ! ---------------- !
          ! Buoyancy Sorting !
          ! ---------------- !

          ! ----------------------------------------------------------------- !
          ! Case 1 : When both cumulus and env. are unsaturated or saturated. !
          ! ----------------------------------------------------------------- !

          if( ( excessu(i) .le. 0._r8 .and. excess0(i) .le. 0._r8 ) .or. ( excessu(i) .ge. 0._r8 .and. excess0(i) .ge. 0._r8 ) ) then
                xc(i) = min(1._r8,max(0._r8,1._r8-2._r8*rbuoy*g*cridis(i)/wue(i)**2._r8*(1._r8-thvj(i)/thv0j(i))))
              ! Below 3 lines are diagnostic output not influencing
              ! numerical calculations.
                aquad(i) = 0._r8
                bquad(i) = 0._r8
                cquad(i) = 0._r8
          else
          ! -------------------------------------------------- !
          ! Case 2 : When either cumulus or env. is saturated. !
          ! -------------------------------------------------- !
              xsat(i)    = excessu(i) / ( excessu(i) - excess0(i) );
              thlxsat(i) = thlue(i) + xsat(i) * ( thle(i) - thlue(i) );
              qtxsat(i)  = qtue(i)  + xsat(i) * ( qte(i) - qtue(i) );
              call conden(pe(i),thlxsat(i),qtxsat(i),thj(i),qvj(i),qlj(i),qij(i),qse(i),id_check)
              if( id_check .eq. 1 ) then
                  exit_conden(i) = 1._r8
                  id_exit(i) = .true.
                  go to 333
              end if
              thvxsat(i) = thj(i) * ( 1._r8 + zvir * qvj(i) - qlj(i) - qij(i) )
              ! -------------------------------------------------- !
              ! kk=1 : Cumulus Segment, kk=2 : Environment Segment !
              ! -------------------------------------------------- !
              do kk = 1, 2
                   if( kk .eq. 1 ) then
                       thv_x0(i) = thvj(i);
                       thv_x1(i) = ( 1._r8 - 1._r8/xsat(i) ) * thvj(i) + ( 1._r8/xsat(i) ) * thvxsat(i);
                   else
                       thv_x1(i) = thv0j(i);
                       thv_x0(i) = ( xsat(i) / ( xsat(i) - 1._r8 ) ) * thv0j(i) + ( 1._r8/( 1._r8 - xsat(i) ) ) * thvxsat(i);
                   endif
                   aquad(i) =  wue(i)**2;
                   bquad(i) =  2._r8*rbuoy*g*cridis(i)*(thv_x1(i) - thv_x0(i))/thv0j(i) - 2._r8*wue(i)**2;
                   cquad(i) =  2._r8*rbuoy*g*cridis(i)*(thv_x0(i) -  thv0j(i))/thv0j(i) +       wue(i)**2;
                   if( kk .eq. 1 ) then
                       if( ( bquad(i)**2-4._r8*aquad(i)*cquad(i) ) .ge. 0._r8 ) then
                             call roots(aquad(i),bquad(i),cquad(i),xs1(i),xs2(i),status)
                             x_cu(i) = min(1._r8,max(0._r8,min(xsat(i),min(xs1(i),xs2(i)))))
                       else
                             x_cu(i) = xsat(i);
                       endif
                   else
                       if( ( bquad(i)**2-4._r8*aquad(i)*cquad(i)) .ge. 0._r8 ) then
                             call roots(aquad(i),bquad(i),cquad(i),xs1(i),xs2(i),status)
                             x_en(i) = min(1._r8,max(0._r8,max(xsat(i),min(xs1(i),xs2(i)))))
                       else
                             x_en(i) = 1._r8;
                       endif
                   endif
              enddo
              if( x_cu(i) .eq. xsat(i) ) then
                  xc(i) = max(x_cu(i), x_en(i));
              else
                  xc(i) = x_cu(i);
              endif
          endif

          ! ------------------------------------------------------------------------ !
          ! Compute fractional lateral entrainment & detrainment rate in each layers.!
          ! The unit of rei(k, i), fer(k, i), and fdr(k, i) is [Pa-1].  Alternative choice of !
          ! 'rei(k, i)' is also shown below, where coefficient 0.5 was from approximate !
          ! tuning against the BOMEX case.                                           !
          ! In order to prevent the onset of instability in association with cumulus !
          ! induced subsidence advection, cumulus mass flux at the top interface  in !
          ! any layer should be smaller than ( 90% of ) total mass within that layer.!
          ! I imposed limits on 'rei(k, i)' as below,  in such that stability condition !
          ! is always satisfied.                                                     !
          ! Below limiter of 'rei(k, i)' becomes negative for some cases, causing error.!
          ! So, for the time being, I came back to the original limiter.             !
          ! ------------------------------------------------------------------------ !
          ee2(i)    = xc(i)**2
          ud2(i)    = 1._r8 - 2._r8*xc(i) + xc(i)**2
        ! rei(k, i) = ( rkm / scaleh(i) / g / rho0j(i) )        ! Default.
          rei(k, i) = ( 0.5_r8 * rkm / z0(k, i) / g /rho0j(i) ) ! Alternative.
          if( xc(i) .gt. 0.5_r8 ) rei(k, i) = min(rei(k, i),0.9_r8*log(dp0(k, i)/g/dt/umf(km1, i) + 1._r8)/dpe(i)/(2._r8*xc(i)-1._r8))
          fer(k, i) = rei(k, i) * ee2(i)
          fdr(k, i) = rei(k, i) * ud2(i)

          ! ------------------------------------------------------------------------------ !
          ! Iteration Start due to 'maxufrc' constraint [ ****************************** ] !
          ! ------------------------------------------------------------------------------ !

          ! -------------------------------------------------------------------------- !
          ! Calculate cumulus updraft mass flux and penetrative entrainment mass flux. !
          ! Note that  non-zero penetrative entrainment mass flux will be asigned only !
          ! to interfaces from the top interface of 'kbup(i)' layer to the base interface !
          ! of 'kpen(i)' layer as will be shown later.                                    !
          ! -------------------------------------------------------------------------- !

          umf(k, i) = umf(km1, i) * exp( dpe(i) * ( fer(k, i) - fdr(k, i) ) )
          emf(k, i) = 0._r8

          ! --------------------------------------------------------- !
          ! Compute cumulus updraft properties at the top interface.  !
          ! Also use Tayler expansion in order to treat limiting case !
          ! --------------------------------------------------------- !

          if( fer(k, i)*dpe(i) .lt. 1.e-4_r8 ) then
              thlu(k, i) = thlu(km1, i) + ( thle(i) + ssthl0(k, i) * dpe(i) / 2._r8 - thlu(km1, i) ) * fer(k, i) * dpe(i)
              qtu(k, i)  =  qtu(km1, i) + ( qte(i)  +  ssqt0(k, i) * dpe(i) / 2._r8 -  qtu(km1, i) ) * fer(k, i) * dpe(i)
              uu(k, i)   =   uu(km1, i) + ( ue(i)   +   ssu0(k, i) * dpe(i) / 2._r8 -   uu(km1, i) ) * fer(k, i) * dpe(i) - PGFc * ssu0(k, i) * dpe(i)
              vu(k, i)   =   vu(km1, i) + ( ve(i)   +   ssv0(k, i) * dpe(i) / 2._r8 -   vu(km1, i) ) * fer(k, i) * dpe(i) - PGFc * ssv0(k, i) * dpe(i)
              do m = 1, ncnst
                 tru(k,m, i)  =  tru(km1,m, i) + ( tre(m, i)  + sstr0(k,m, i) * dpe(i) / 2._r8  -  tru(km1,m, i) ) * fer(k, i) * dpe(i)
              enddo
          else
              thlu(k, i) = ( thle(i) + ssthl0(k, i) / fer(k, i) - ssthl0(k, i) * dpe(i) / 2._r8 ) -          &
                        ( thle(i) + ssthl0(k, i) * dpe(i) / 2._r8 - thlu(km1, i) + ssthl0(k, i) / fer(k, i) ) * exp(-fer(k, i) * dpe(i))
              qtu(k, i)  = ( qte(i)  +  ssqt0(k, i) / fer(k, i) -  ssqt0(k, i) * dpe(i) / 2._r8 ) -          &
                        ( qte(i)  +  ssqt0(k, i) * dpe(i) / 2._r8 -  qtu(km1, i) +  ssqt0(k, i) / fer(k, i) ) * exp(-fer(k, i) * dpe(i))
              uu(k, i) =   ( ue(i) + ( 1._r8 - PGFc ) * ssu0(k, i) / fer(k, i) - ssu0(k, i) * dpe(i) / 2._r8 ) - &
                        ( ue(i) +     ssu0(k, i) * dpe(i) / 2._r8 -   uu(km1, i) + ( 1._r8 - PGFc ) * ssu0(k, i) / fer(k, i) ) * exp(-fer(k, i) * dpe(i))
              vu(k, i) =   ( ve(i) + ( 1._r8 - PGFc ) * ssv0(k, i) / fer(k, i) - ssv0(k, i) * dpe(i) / 2._r8 ) - &
                        ( ve(i) +     ssv0(k, i) * dpe(i) / 2._r8 -   vu(km1, i) + ( 1._r8 - PGFc ) * ssv0(k, i) / fer(k, i) ) * exp(-fer(k, i) * dpe(i))
              do m = 1, ncnst
                 tru(k,m, i)  = ( tre(m, i)  + sstr0(k,m, i) / fer(k, i) - sstr0(k,m, i) * dpe(i) / 2._r8 ) - &
                             ( tre(m, i)  + sstr0(k,m, i) * dpe(i) / 2._r8 - tru(km1,m, i) + sstr0(k,m, i) / fer(k, i) ) * exp(-fer(k, i) * dpe(i))
              enddo
          end if

          !------------------------------------------------------------------- !
          ! Expel some of cloud water and ice from cumulus  updraft at the top !
          ! interface.  Note that this is not 'detrainment' term  but a 'sink' !
          ! term of cumulus updraft qt ( or one part of 'source' term of  mean !
          ! environmental qt ). At this stage, as the most simplest choice, if !
          ! condensate amount within cumulus updraft is larger than a critical !
          ! value, 'criqc', expels the surplus condensate from cumulus updraft !
          ! to the environment. A certain fraction ( e.g., 'frc_sus' ) of this !
          ! expelled condesnate will be in a form that can be suspended in the !
          ! layer k where it was formed, while the other fraction, '1-frc_sus' !
          ! will be in a form of precipitatble (e.g.,can potentially fall down !
          ! across the base interface of layer k ). In turn we should describe !
          ! subsequent falling of precipitable condensate ('1-frc_sus') across !
          ! the base interface of the layer k, &  evaporation of precipitating !
          ! water in the below layer k-1 and associated evaporative cooling of !
          ! the later, k-1, and falling of 'non-evaporated precipitating water !
          ! ( which was initially formed in layer k ) and a newly-formed preci !
          ! pitable water in the layer, k-1', across the base interface of the !
          ! lower layer k-1.  Cloud microphysics should correctly describe all !
          ! of these process.  In a near future, I should significantly modify !
          ! this cloud microphysics, including precipitation-induced downdraft !
          ! also.                                                              !
          ! ------------------------------------------------------------------ !

          call conden(ps0(k, i),thlu(k, i),qtu(k, i),thj(i),qvj(i),qlj(i),qij(i),qse(i),id_check)
          if( id_check .eq. 1 ) then
              exit_conden(i) = 1._r8
              id_exit(i) = .true.
              go to 333
          end if
          if( (qlj(i) + qij(i)) .gt. criqc ) then
               exql(i)    = ( ( qlj(i) + qij(i) ) - criqc ) * qlj(i) / ( qlj(i) + qij(i) )
               exqi(i)    = ( ( qlj(i) + qij(i) ) - criqc ) * qij(i) / ( qlj(i) + qij(i) )
               ! ---------------------------------------------------------------- !
               ! It is very important to re-update 'qtu(:, i)' and 'thlu(:, i)'  at the upper !
               ! interface after expelling condensate from cumulus updraft at the !
               ! top interface of the layer. As mentioned above, this is a 'sink' !
               ! of cumulus qt (or equivalently, a 'source' of environmentasl qt),!
               ! not a regular convective'detrainment'.                           !
               ! ---------------------------------------------------------------- !
               qtu(k, i)  = qtu(k, i) - exql(i) - exqi(i)
               thlu(k, i) = thlu(k, i) + (xlv/cp/exns0(k, i))*exql(i) + (xls/cp/exns0(k, i))*exqi(i)
               ! ---------------------------------------------------------------- !
               ! Expelled cloud condensate into the environment from the updraft. !
               ! After all the calculation later, 'dwten(:, i)' and 'diten(:, i)' will have a !
               ! unit of [ kg/kg/s ], because it is a tendency of qt. Restoration !
               ! of 'dwten(:, i)' and 'diten(:, i)' to this correct unit through  multiplying !
               ! 'umf(k, i)*g/dp0(k, i)' will be performed later after finally updating !
               ! 'umf(:, i)' using a 'rmaxfrac' constraint near the end of this updraft !
               ! buoyancy sorting loop.                                           !
               ! ---------------------------------------------------------------- !
               dwten(k, i) = exql(i)
               diten(k, i) = exqi(i)
          else
               dwten(k, i) = 0._r8
               diten(k, i) = 0._r8
          endif
          ! ----------------------------------------------------------------- !
          ! Update 'thvu(k, i)' after detraining condensate from cumulus updraft.!
          ! ----------------------------------------------------------------- !
          call conden(ps0(k, i),thlu(k, i),qtu(k, i),thj(i),qvj(i),qlj(i),qij(i),qse(i),id_check)
          if( id_check .eq. 1 ) then
              exit_conden(i) = 1._r8
              id_exit(i) = .true.
              go to 333
          end if
          thvu(k, i) = thj(i) * ( 1._r8 + zvir * qvj(i) - qlj(i) - qij(i) )

          ! ----------------------------------------------------------- !
          ! Calculate updraft vertical velocity at the upper interface. !
          ! In order to calculate 'wtw(i)' at the upper interface, we use  !
          ! 'wtw(i)' at the lower interface. Note  'wtw(i)'  is continuously  !
          ! updated as cumulus updraft rises.                           !
          ! ----------------------------------------------------------- !

          bogbot(i) = rbuoy * ( thvu(km1, i) / thvebot(i)  - 1._r8 ) ! Cloud buoyancy at base interface
          bogtop(i) = rbuoy * ( thvu(k, i) / thv0top(k, i) - 1._r8 ) ! Cloud buoyancy at top  interface

          delbog(i) = bogtop(i) - bogbot(i)
          drage(i)  = fer(k, i) * ( 1._r8 + rdrag )
          expfac(i) = exp(-2._r8*drage(i)*dpe(i))

          wtwb(i) = wtw(i)
          if( drage(i)*dpe(i) .gt. 1.e-3_r8 ) then
              wtw(i) = wtw(i)*expfac(i) + (delbog(i) + (1._r8-expfac(i))*(bogbot(i) + delbog(i)/(-2._r8*drage(i)*dpe(i))))/(rho0j(i)*drage(i))
          else
              wtw(i) = wtw(i) + dpe(i) * ( bogbot(i) + bogtop(i) ) / rho0j(i)
          endif

        ! Force the plume rise at least to klfc(i) of the undiluted plume.
        ! Because even the below is not complete, I decided not to include this.

        ! if( k .le. klfc(i) ) then
        !     wtw(i) = max( 1.e-2_r8, wtw(i) )
        ! endif

          ! -------------------------------------------------------------- !
          ! Repeat 'iter_xc' iteration loop until 'iter_xc = niter_xc'.    !
          ! Also treat the case even when wtw(i) < 0 at the 'kpen(i)' interface. !
          ! -------------------------------------------------------------- !

          if( wtw(i) .gt. 0._r8 ) then
              thlue(i) = 0.5_r8 * ( thlu(km1, i) + thlu(k, i) )
              qtue(i)  = 0.5_r8 * ( qtu(km1, i)  +  qtu(k, i) )
              wue(i)   = 0.5_r8 *   sqrt( max( wtwb(i) + wtw(i), 0._r8 ) )
          else
              go to 111
          endif

       enddo ! End of 'iter_xc' loop

   111 continue

          ! --------------------------------------------------------------------------- !
          ! Add the contribution of self-detrainment  to vertical variations of cumulus !
          ! updraft mass flux. The reason why we are trying to include self-detrainment !
          ! is as follows.  In current scheme,  vertical variation of updraft mass flux !
          ! is not fully consistent with the vertical variation of updraft vertical w.  !
          ! For example, within a given layer, let's assume that  cumulus w is positive !
          ! at the base interface, while negative at the top interface. This means that !
          ! cumulus updraft cannot reach to the top interface of the layer. However,    !
          ! cumulus updraft mass flux at the top interface is not zero according to the !
          ! vertical tendency equation of cumulus mass flux.   Ideally, cumulus updraft !
          ! mass flux at the top interface should be zero for this case. In order to    !
          ! assures that cumulus updraft mass flux goes to zero when cumulus updraft    !
          ! vertical velocity goes to zero, we are imposing self-detrainment term as    !
          ! below by considering layer-mean cloud buoyancy and cumulus updraft vertical !
          ! velocity square at the top interface. Use of auto-detrainment term will  be !
          ! determined by setting 'use_self_detrain=.true.' in the parameter sentence.  !
          ! --------------------------------------------------------------------------- !

          if( use_self_detrain ) then
              autodet(i) = min( 0.5_r8*g*(bogbot(i)+bogtop(i))/(max(wtw(i),0._r8)+1.e-4_r8), 0._r8 )
              umf(k, i)  = umf(k, i) * exp( 0.637_r8*(dpe(i)/rho0j(i)/g) * autodet(i) )
          end if
          if( umf(k, i) .eq. 0._r8 ) wtw(i) = -1._r8

          ! -------------------------------------- !
          ! Below block is just a dignostic output !
          ! -------------------------------------- !

          excessu_arr(k, i) = excessu(i)
          excess0_arr(k, i) = excess0(i)
          xc_arr(k, i)      = xc(i)
          aquad_arr(k, i)   = aquad(i)
          bquad_arr(k, i)   = bquad(i)
          cquad_arr(K, i)   = cquad(i)
          bogbot_arr(k, i)  = bogbot(i)
          bogtop_arr(k, i)  = bogtop(i)

          ! ------------------------------------------------------------------- !
          ! 'kbup(i)' is the upper most layer in which cloud buoyancy  is positive !
          ! both at the base and top interface.  'kpen(i)' is the upper most layer !
          ! up to cumulus can reach. Usually, 'kpen(i)' is located higher than the !
          ! 'kbup(i)'. Note we initialized these by 'kbup(i) = krel(i)' & 'kpen(i) = krel(i)'. !
          ! As explained before, it is possible that only 'kpen(i)' is updated,    !
          ! while 'kbup(i)' keeps its initialization value. For this case, current !
          ! scheme will simply turns-off penetrative entrainment fluxes and use !
          ! normal buoyancy-sorting fluxes for 'kbup(i) <= k <= kpen(i)-1' interfaces,!
          ! in order to describe shallow continental cumulus convection.        !
          ! ------------------------------------------------------------------- !

        ! if( bogbot(i) .gt. 0._r8 .and. bogtop(i) .gt. 0._r8 ) then
        ! if( bogtop(i) .gt. 0._r8 ) then
          if( bogtop(i) .gt. 0._r8 .and. wtw(i) .gt. 0._r8 ) then
              kbup(i) = k
          end if

          if( wtw(i) .le. 0._r8 ) then
              kpen(i) = k
              go to 45
          end if

          wu(k, i) = sqrt(wtw(i))
          if( wu(k, i) .gt. 100._r8 ) then
              exit_wu(i) = 1._r8
              id_exit(i) = .true.
              go to 333
          endif

          ! ---------------------------------------------------------------------------- !
          ! Iteration end due to 'rmaxfrac' constraint [ ***************************** ] !
          ! ---------------------------------------------------------------------------- !

          ! ---------------------------------------------------------------------- !
          ! Calculate updraft fractional area at the upper interface and set upper !
          ! limit to 'ufrc(:, i)' by 'rmaxfrac'. In order to keep the consistency  among !
          ! ['ufrc(:, i)','umf(:, i)','wu(:, i) (or wtw(i))'], if ufrc(:, i) is limited by 'rmaxfrac', either !
          ! 'umf(:, i)' or 'wu(:, i)' should be changed. Although both 'umf(:, i)' and 'wu(:, i) (wtw(i))' at !
          ! the current upper interface are used for updating 'umf(:, i)' & 'wu(:, i)'  at the !
          ! next upper interface, 'umf(:, i)' is a passive variable not influencing  the !
          ! buoyancy sorting process in contrast to 'wtw(i)'. This is a reason why we !
          ! adjusted 'umf(:, i)' instead of 'wtw(i)'. In turn we updated 'fdr(:, i)' here instead !
          ! of 'fer(:, i)',  which guarantees  that all previously updated thermodynamic !
          ! variables at the upper interface before applying 'rmaxfrac' constraint !
          ! are already internally consistent,  even though 'ufrc(:, i)'  is  limited by !
          ! 'rmaxfrac'. Thus, we don't need to go through interation loop again.If !
          ! If we update 'fer(:, i)' however, we should go through above iteration loop. !
          ! ---------------------------------------------------------------------- !

          rhos0j(i)  = ps0(k, i) / ( r * 0.5_r8 * ( thv0bot(k+1, i) + thv0top(k, i) ) * exns0(k, i) )
          ufrc(k, i) = umf(k, i) / ( rhos0j(i) * wu(k, i) )
          if( ufrc(k, i) .gt. rmaxfrac ) then
              limit_ufrc(i) = 1._r8
              ufrc(k, i) = rmaxfrac
              umf(k, i)  = rmaxfrac * rhos0j(i) * wu(k, i)
              fdr(k, i)  = fer(k, i) - log( umf(k, i) / umf(km1, i) ) / dpe(i)
          endif

          ! ------------------------------------------------------------ !
          ! Update environmental properties for at the mid-point of next !
          ! upper layer for use in buoyancy sorting.                     !
          ! ------------------------------------------------------------ !

          pe(i)      = p0(k+1, i)
          dpe(i)     = dp0(k+1, i)
          exne(i)    = exn0(k+1, i)
          thvebot(i) = thv0bot(k+1, i)
          thle(i)    = thl0(k+1, i)
          qte(i)     = qt0(k+1, i)
          ue(i)      = u0(k+1, i)
          ve(i)      = v0(k+1, i)
          do m = 1, ncnst
             tre(m, i)  = tr0(k+1,m, i)
          enddo

       end do   ! End of cumulus updraft loop from the 'krel(i)' layer to 'kpen(i)' layer.

       ! ------------------------------------------------------------------------------- !
       ! Up to this point, we finished all of buoyancy sorting processes from the 'krel(i)' !
       ! layer to 'kpen(i)' layer: at the top interface of individual layers, we calculated !
       ! updraft and penetrative mass fluxes [ umf(k, i) & emf(k, i) = 0 ], updraft fractional !
       ! area [ ufrc(k, i) ],  updraft vertical velocity [ wu(k, i) ],  updraft  thermodynamic !
       ! variables [thlu(k, i),qtu(k, i),uu(k, i),vu(k, i),thvu(k, i)]. In the layer,we also calculated !
       ! fractional entrainment-detrainment rate [ fer(k, i), fdr(k, i) ], and detrainment ten !
       ! dency of water and ice from cumulus updraft [ dwten(k, i), diten(k, i) ]. In addition,!
       ! we updated and identified 'krel(i)' and 'kpen(i)' layer index, if any.  In the 'kpen(i)' !
       ! layer, we calculated everything mentioned above except the 'wu(k, i)' and 'ufrc(k, i)'!
       ! since a real value of updraft vertical velocity is not defined at the kpen(i)  top !
       ! interface (note 'ufrc(:, i)' at the top interface of layer is calculated from 'umf(k, i)'!
       ! and 'wu(k, i)'). As mentioned before, special treatment is required when 'kbup(i)' is !
       ! not updated and so 'kbup(i) = krel(i)'.                                               !
       ! ------------------------------------------------------------------------------- !

       ! ------------------------------------------------------------------------------ !
       ! During the 'iter_scaleh' iteration loop, non-physical ( with non-zero values ) !
       ! values can remain in the variable arrays above (also 'including' in case of wu(:, i) !
       ! and ufrc(:, i) at the top interface) the 'kpen(i)' layer. This can happen when the kpen(i) !
       ! layer index identified from the 'iter_scaleh = 1' iteration loop is located at !
       ! above the kpen(i) layer index identified from   'iter_scaleh = 3' iteration loop. !
       ! Thus, in the following calculations, we should only use the values in each     !
       ! variables only up to finally identified 'kpen(i)' layer & 'kpen(i)' interface except !
       ! 'wu(:, i)' and 'ufrc(:, i)' at the top interface of 'kpen(i)' layer.    Note that in order to !
       ! prevent any problems due to these non-physical values, I re-initialized    the !
       ! values of [ umf(kpen(i):mkx, i), emf(kpen(i):mkx, i), dwten(kpen(i)+1:mkx, i), diten(kpen(i)+1:mkx, i),!
       ! fer(kpen(i):mkx, i), fdr(kpen(i)+1:mkx, i), ufrc(kpen(i):mkx, i) ] to be zero after 'iter_scaleh'!
       ! do loop.                                                                       !
       ! ------------------------------------------------------------------------------ !

 45    continue

       ! ------------------------------------------------------------------------------ !
       ! Calculate 'ppen( < 0 , i)', updarft penetrative distance from the lower interface !
       ! of 'kpen(i)' layer. Note that bogbot(i) & bogtop(i) at the 'kpen(i)' layer either when fer(:, i) !
       ! is zero or non-zero was already calculated above.                              !
       ! It seems that below qudarature solving formula is valid only when bogbot(i) < 0.  !
       ! Below solving equation is clearly wrong ! I should revise this !               !
       ! ------------------------------------------------------------------------------ !

       if( drage(i) .eq. 0._r8 ) then
           aquad(i) =  ( bogtop(i) - bogbot(i) ) / ( ps0(kpen(i), i) - ps0(kpen(i)-1, i) )
           bquad(i) =  2._r8 * bogbot(i)
           cquad(i) = -wu(kpen(i)-1, i)**2 * rho0j(i)
           call roots(aquad(i),bquad(i),cquad(i),xc1(i),xc2(i),status)
           if( status .eq. 0 ) then
               if( xc1(i) .le. 0._r8 .and. xc2(i) .le. 0._r8 ) then
                   ppen(i) = max( xc1(i), xc2(i) )
                   ppen(i) = min( 0._r8,max( -dp0(kpen(i), i), ppen(i) ) )
               elseif( xc1(i) .gt. 0._r8 .and. xc2(i) .gt. 0._r8 ) then
                   ppen(i) = -dp0(kpen(i), i)
                   write(iulog,*) 'Warning : UW-Cumulus penetrates upto kpen(i) interface'
               else
                   ppen(i) = min( xc1(i), xc2(i) )
                   ppen(i) = min( 0._r8,max( -dp0(kpen(i), i), ppen(i) ) )
               endif
           else
               ppen(i) = -dp0(kpen(i), i)
               write(iulog,*) 'Warning : UW-Cumulus penetrates upto kpen(i) interface'
           endif
       else
           ppen(i) = compute_ppen(wtwb(i),drage(i),bogbot(i),bogtop(i),rho0j(i),dp0(kpen(i), i))
       endif
       if( ppen(i) .eq. -dp0(kpen(i), i) .or. ppen(i) .eq. 0._r8 ) limit_ppen(i) = 1._r8

       ! -------------------------------------------------------------------- !
       ! Re-calculate the amount of expelled condensate from cloud updraft    !
       ! at the cumulus top. This is necessary for refined calculations of    !
       ! bulk cloud microphysics at the cumulus top. Note that ppen(i) < 0._r8   !
       ! In the below, I explicitly calculate 'thlu_top(i)' & 'qtu_top(i)' by       !
       ! using non-zero 'fer(kpen(i), i)'.                                          !
       ! -------------------------------------------------------------------- !

       if( fer(kpen(i), i)*(-ppen(i)) .lt. 1.e-4_r8 ) then
           thlu_top(i) = thlu(kpen(i)-1, i) + ( thl0(kpen(i), i) + ssthl0(kpen(i), i) * (-ppen(i)) / 2._r8 - thlu(kpen(i)-1, i) ) * fer(kpen(i), i) * (-ppen(i))
           qtu_top(i)  =  qtu(kpen(i)-1, i) + (  qt0(kpen(i), i) +  ssqt0(kpen(i), i) * (-ppen(i)) / 2._r8  - qtu(kpen(i)-1, i) ) * fer(kpen(i), i) * (-ppen(i))
       else
           thlu_top(i) = ( thl0(kpen(i), i) + ssthl0(kpen(i), i) / fer(kpen(i), i) - ssthl0(kpen(i), i) * (-ppen(i)) / 2._r8 ) - &
                      ( thl0(kpen(i), i) + ssthl0(kpen(i), i) * (-ppen(i)) / 2._r8 - thlu(kpen(i)-1, i) + ssthl0(kpen(i), i) / fer(kpen(i), i) ) &
                      * exp(-fer(kpen(i), i) * (-ppen(i)))
           qtu_top(i)  = ( qt0(kpen(i), i)  +  ssqt0(kpen(i), i) / fer(kpen(i), i) -  ssqt0(kpen(i), i) * (-ppen(i)) / 2._r8 ) - &
                      ( qt0(kpen(i), i)  +  ssqt0(kpen(i), i) * (-ppen(i)) / 2._r8 -  qtu(kpen(i)-1, i) +  ssqt0(kpen(i), i) / fer(kpen(i), i) ) &
                      * exp(-fer(kpen(i), i) * (-ppen(i)))
       end if

       call conden(ps0(kpen(i)-1, i)+ppen(i),thlu_top(i),qtu_top(i),thj(i),qvj(i),qlj(i),qij(i),qse(i),id_check)
       if( id_check .eq. 1 ) then
           exit_conden(i) = 1._r8
           id_exit(i) = .true.
           go to 333
       end if
       exntop(i) = ((ps0(kpen(i)-1, i)+ppen(i))/p00)**rovcp
       if( (qlj(i) + qij(i)) .gt. criqc ) then
            dwten(kpen(i), i) = ( ( qlj(i) + qij(i) ) - criqc ) * qlj(i) / ( qlj(i) + qij(i) )
            diten(kpen(i), i) = ( ( qlj(i) + qij(i) ) - criqc ) * qij(i) / ( qlj(i) + qij(i) )
            qtu_top(i)  = qtu_top(i) - dwten(kpen(i), i) - diten(kpen(i), i)
            thlu_top(i) = thlu_top(i) + (xlv/cp/exntop(i))*dwten(kpen(i), i) + (xls/cp/exntop(i))*diten(kpen(i), i)
       else
            dwten(kpen(i), i) = 0._r8
            diten(kpen(i), i) = 0._r8
       endif

       ! ----------------------------------------------------------------------- !
       ! Calculate cumulus scale height as the top height that cumulus can reach.!
       ! ----------------------------------------------------------------------- !

       rhos0j(i) = ps0(kpen(i)-1, i)/(r*0.5_r8*(thv0bot(kpen(i), i)+thv0top(kpen(i)-1, i))*exns0(kpen(i)-1, i))
       cush(i)   = zs0(kpen(i)-1, i) - ppen(i)/rhos0j(i)/g
       scaleh(i) = cush(i)

    end do   ! End of 'iter_scaleh' loop.

       ! -------------------------------------------------------------------- !
       ! The 'forcedCu(i)' is logical identifier saying whether cumulus updraft  !
       ! overcome the buoyancy barrier just above the PBL top. If it is true, !
       ! cumulus did not overcome the barrier -  this is a shallow convection !
       ! with negative cloud buoyancy, mimicking  shallow continental cumulus !
       ! convection. Depending on 'forcedCu(i)' parameter, treatment of heat  &  !
       ! moisture fluxes at the entraining interfaces, 'kbup(i) <= k < kpen(i) - 1' !
       ! will be set up in a different ways, as will be shown later.          !
       ! -------------------------------------------------------------------- !

       if( kbup(i) .eq. krel(i) ) then
           forcedCu(i) = .true.
           limit_shcu(i) = 1._r8
       else
           forcedCu(i) = .false.
           limit_shcu(i) = 0._r8
       endif

       ! ------------------------------------------------------------------ !
       ! Filtering of unerasonable cumulus adjustment here.  This is a very !
       ! important process which should be done cautiously. Various ways of !
       ! filtering are possible depending on cases mainly using the indices !
       ! of key layers - 'klcl(i)','kinv(i)','krel(i)','klfc(i)','kbup(i)','kpen(i)'. At this !
       ! stage, the followings are all possible : 'kinv(i) >= 2', 'klcl(i) >= 1', !
       ! 'krel(i) >= kinv(i)', 'kbup(i) >= krel(i)', 'kpen(i) >= krel(i)'. I must design this !
       ! filtering very cautiously, in such that none of  realistic cumulus !
       ! convection is arbitrarily turned-off. Potentially, I might turn-off!
       ! cumulus convection if layer-mean 'ql > 0' in the 'kinv(i)-1' layer,in !
       ! order to suppress cumulus convection growing, based at the Sc top. !
       ! This is one of potential future modifications. Note that ppen(i) < 0. !
       ! ------------------------------------------------------------------ !

       cldhgt(i) = ps0(kpen(i)-1, i) + ppen(i)
       if( forcedCu(i) ) then
           ! write(iulog,*) 'forcedCu(i) - did not overcome initial buoyancy barrier'
           exit_cufilter(i) = 1._r8
           id_exit(i) = .true.
           go to 333
       end if
       ! Limit 'additional shallow cumulus' for DYCOMS simulation.
       ! if( cldhgt(i).ge.88000._r8 ) then
       !     id_exit(i) = .true.
       !     go to 333
       ! end if

       ! ------------------------------------------------------------------------------ !
       ! Re-initializing some key variables above the 'kpen(i)' layer in order to suppress !
       ! the influence of non-physical values above 'kpen(i)', in association with the use !
       ! of 'iter_scaleh' loop. Note that umf(:, i), emf(:, i),  ufrc(:, i) are defined at the interfaces !
       ! (0:mkx), while 'dwten(:, i)','diten(:, i)', 'fer(:, i)', 'fdr(:, i)' are defined at layer mid-points.  !
       ! Initialization of 'fer(:, i)' and 'fdr(:, i)' is for correct writing purpose of diagnostic !
       ! output. Note that we set umf(kpen(i), i)=emf(kpen(i), i)=ufrc(kpen(i), i)=0, in consistent  with !
       ! wtw(i) < 0  at the top interface of 'kpen(i)' layer. However, we still have non-zero !
       ! expelled cloud condensate in the 'kpen(i)' layer.                                 !
       ! ------------------------------------------------------------------------------ !

       umf(kpen(i):mkx, i)     = 0._r8
       emf(kpen(i):mkx, i)     = 0._r8
       ufrc(kpen(i):mkx, i)    = 0._r8
       dwten(kpen(i)+1:mkx, i) = 0._r8
       diten(kpen(i)+1:mkx, i) = 0._r8
       fer(kpen(i)+1:mkx, i)   = 0._r8
       fdr(kpen(i)+1:mkx, i)   = 0._r8

       ! ------------------------------------------------------------------------ !
       ! Calculate downward penetrative entrainment mass flux, 'emf(k, i) < 0',  and !
       ! thermodynamic properties of penetratively entrained airs at   entraining !
       ! interfaces. emf(k, i) is defined from the top interface of the  layer  kbup(i) !
       ! to the bottom interface of the layer 'kpen(i)'. Note even when  kbup(i) = krel(i),!
       ! i.e.,even when 'kbup(i)' was not updated in the above buoyancy  sorting  do !
       ! loop (i.e., 'kbup(i)' remains as the initialization value),   below do loop !
       ! of penetrative entrainment flux can be performed without  any conceptual !
       ! or logical problems, because we have already computed all  the variables !
       ! necessary for performing below penetrative entrainment block.            !
       ! In the below 'do' loop, 'k' is an interface index at which non-zero 'emf(:, i)'!
       ! (penetrative entrainment mass flux) is calculated. Since cumulus updraft !
       ! is negatively buoyant in the layers between the top interface of 'kbup(i)'  !
       ! layer (interface index, kbup(i)) and the top interface of 'kpen(i)' layer, the !
       ! fractional lateral entrainment, fer(k, i) within these layers will be close !
       ! to zero - so it is likely that only strong lateral detrainment occurs in !
       ! thses layers. Under this situation,we can easily calculate the amount of !
       ! detrainment cumulus air into these negatively buoyanct layers by  simply !
       ! comparing cumulus updraft mass fluxes between the base and top interface !
       ! of each layer: emf(k, i) = emf(k-1, i)*exp(-fdr(k, i)*dp0(k, i))                     !
       !                       ~ emf(k-1, i)*(1-rei(k, i)*dp0(k, i))                       !
       !                emf(k-1, i)-emf(k, i) ~ emf(k-1, i)*rei(k, i)*dp0(k, i)                  !
       ! Current code assumes that about 'rpen~10' times of these detrained  mass !
       ! are penetratively re-entrained down into the 'k-1' interface. And all of !
       ! these detrained masses are finally dumped down into the top interface of !
       ! 'kbup(i)' layer. Thus, the amount of penetratively entrained air across the !
       ! top interface of 'kbup(i)' layer with 'rpen~10' becomes too large.          !
       ! Note that this penetrative entrainment part can be completely turned-off !
       ! and we can simply use normal buoyancy-sorting involved turbulent  fluxes !
       ! by modifying 'penetrative entrainment fluxes' part below.                !
       ! ------------------------------------------------------------------------ !

       ! -----------------------------------------------------------------------!
       ! Calculate entrainment mass flux and conservative scalars of entraining !
       ! free air at interfaces of 'kbup(i) <= k < kpen(i) - 1'                       !
       ! ---------------------------------------------------------------------- !

       do k = 0, mkx
          thlu_emf(k, i) = thlu(k, i)
          qtu_emf(k, i)  = qtu(k, i)
          uu_emf(k, i)   = uu(k, i)
          vu_emf(k, i)   = vu(k, i)
          do m = 1, ncnst
             tru_emf(k,m, i)  = tru(k,m, i)
          enddo
       end do

       do k = kpen(i) - 1, kbup(i), -1  ! Here, 'k' is an interface index at which
                                  ! penetrative entrainment fluxes are calculated.

          rhos0j(i) = ps0(k, i) / ( r * 0.5_r8 * ( thv0bot(k+1, i) + thv0top(k, i) ) * exns0(k, i) )

          if( k .eq. kpen(i) - 1 ) then

             ! ------------------------------------------------------------------------ !
             ! Note that 'ppen(i)' has already been calculated in the above 'iter_scaleh'  !
             ! loop assuming zero lateral entrainmentin the layer 'kpen(i)'.               !
             ! ------------------------------------------------------------------------ !

             ! -------------------------------------------------------------------- !
             ! Calculate returning mass flux, emf(:, i) ( < 0 )                           !
             ! Current penetrative entrainment rate with 'rpen~10' is too large and !
             ! future refinement is necessary including the definition of 'thl','qt'!
             ! of penetratively entrained air.  Penetratively entrained airs across !
             ! the 'kpen(i)-1' interface is assumed to have the properties of the base !
             ! interface of 'kpen(i)' layer. Note that 'emf(:, i) ~ - umf(:, i)/ufrc(:, i) = - w * rho'. !
             ! Thus, below limit sets an upper limit of |emf(:, i)| to be ~ 10cm/s, which !
             ! is very loose constraint. Here, I used more restricted constraint on !
             ! the limit of emf(:, i), assuming 'emf(:, i)' cannot exceed a net mass within the !
             ! layer above the interface. Similar to the case of warming and drying !
             ! due to cumulus updraft induced compensating subsidence,  penetrative !
             ! entrainment induces compensating upwelling -     in order to prevent !
             ! numerical instability in association with compensating upwelling, we !
             ! should similarily limit the amount of penetrative entrainment at the !
             ! interface by the amount of masses within the layer just above the    !
             ! penetratively entraining interface.                                  !
             ! -------------------------------------------------------------------- !

             if( ( umf(k, i)*ppen(i)*rei(kpen(i), i)*rpen ) .lt. -0.1_r8*rhos0j(i) )         limit_emf(i) = 1._r8
             if( ( umf(k, i)*ppen(i)*rei(kpen(i), i)*rpen ) .lt. -0.9_r8*dp0(kpen(i), i)/g/dt ) limit_emf(i) = 1._r8

             emf(k, i) = max( max( umf(k, i)*ppen(i)*rei(kpen(i), i)*rpen, -0.1_r8*rhos0j(i)), -0.9_r8*dp0(kpen(i), i)/g/dt)
             thlu_emf(k, i) = thl0(kpen(i), i) + ssthl0(kpen(i), i) * ( ps0(k, i) - p0(kpen(i), i) )
             qtu_emf(k, i)  = qt0(kpen(i), i)  + ssqt0(kpen(i), i)  * ( ps0(k, i) - p0(kpen(i), i) )
             uu_emf(k, i)   = u0(kpen(i), i)   + ssu0(kpen(i), i)   * ( ps0(k, i) - p0(kpen(i), i) )
             vu_emf(k, i)   = v0(kpen(i), i)   + ssv0(kpen(i), i)   * ( ps0(k, i) - p0(kpen(i), i) )
             do m = 1, ncnst
                tru_emf(k,m, i)  = tr0(kpen(i),m, i)  + sstr0(kpen(i),m, i)  * ( ps0(k, i) - p0(kpen(i), i) )
             enddo

          else ! if(k.lt.kpen(i)-1).

             ! --------------------------------------------------------------------------- !
             ! Note we are coming down from the higher interfaces to the lower interfaces. !
             ! Also note that 'emf(:, i) < 0'. So, below operation is a summing not subtracting. !
             ! In order to ensure numerical stability, I imposed a modified correct limit  !
             ! of '-0.9*dp0(k+1, i)/g/dt' on emf(k, i).                                          !
             ! --------------------------------------------------------------------------- !

             if( use_cumpenent ) then  ! Original Cumulative Penetrative Entrainment

                 if( ( emf(k+1, i)-umf(k, i)*dp0(k+1, i)*rei(k+1, i)*rpen ) .lt. -0.1_r8*rhos0j(i) )        limit_emf(i) = 1
                 if( ( emf(k+1, i)-umf(k, i)*dp0(k+1, i)*rei(k+1, i)*rpen ) .lt. -0.9_r8*dp0(k+1, i)/g/dt ) limit_emf(i) = 1
                 emf(k, i) = max(max(emf(k+1, i)-umf(k, i)*dp0(k+1, i)*rei(k+1, i)*rpen, -0.1_r8*rhos0j(i)), -0.9_r8*dp0(k+1, i)/g/dt )
                 if( abs(emf(k, i)) .gt. abs(emf(k+1, i)) ) then
                     thlu_emf(k, i) = ( thlu_emf(k+1, i) * emf(k+1, i) + thl0(k+1, i) * ( emf(k, i) - emf(k+1, i) ) ) / emf(k, i)
                     qtu_emf(k, i)  = ( qtu_emf(k+1, i)  * emf(k+1, i) + qt0(k+1, i)  * ( emf(k, i) - emf(k+1, i) ) ) / emf(k, i)
                     uu_emf(k, i)   = ( uu_emf(k+1, i)   * emf(k+1, i) + u0(k+1, i)   * ( emf(k, i) - emf(k+1, i) ) ) / emf(k, i)
                     vu_emf(k, i)   = ( vu_emf(k+1, i)   * emf(k+1, i) + v0(k+1, i)   * ( emf(k, i) - emf(k+1, i) ) ) / emf(k, i)
                     do m = 1, ncnst
                        tru_emf(k,m, i)  = ( tru_emf(k+1,m, i)  * emf(k+1, i) + tr0(k+1,m, i)  * ( emf(k, i) - emf(k+1, i) ) ) / emf(k, i)
                     enddo
                 else
                     thlu_emf(k, i) = thl0(k+1, i)
                     qtu_emf(k, i)  =  qt0(k+1, i)
                     uu_emf(k, i)   =   u0(k+1, i)
                     vu_emf(k, i)   =   v0(k+1, i)
                     do m = 1, ncnst
                        tru_emf(k,m, i)  =  tr0(k+1,m, i)
                     enddo
                 endif

             else ! Alternative Non-Cumulative Penetrative Entrainment

                 if( ( -umf(k, i)*dp0(k+1, i)*rei(k+1, i)*rpen ) .lt. -0.1_r8*rhos0j(i) )        limit_emf(i) = 1
                 if( ( -umf(k, i)*dp0(k+1, i)*rei(k+1, i)*rpen ) .lt. -0.9_r8*dp0(k+1, i)/g/dt ) limit_emf(i) = 1
                 emf(k, i) = max(max(-umf(k, i)*dp0(k+1, i)*rei(k+1, i)*rpen, -0.1_r8*rhos0j(i)), -0.9_r8*dp0(k+1, i)/g/dt )
                 thlu_emf(k, i) = thl0(k+1, i)
                 qtu_emf(k, i)  =  qt0(k+1, i)
                 uu_emf(k, i)   =   u0(k+1, i)
                 vu_emf(k, i)   =   v0(k+1, i)
                 do m = 1, ncnst
                    tru_emf(k,m, i)  =  tr0(k+1,m, i)
                 enddo

             endif

          endif

          ! ---------------------------------------------------------------------------- !
          ! In this GCM modeling framework,  all what we should do is to calculate  heat !
          ! and moisture fluxes at the given geometrically-fixed height interfaces -  we !
          ! don't need to worry about movement of material height surface in association !
          ! with compensating subsidence or unwelling, in contrast to the bulk modeling. !
          ! In this geometrically fixed height coordinate system, heat and moisture flux !
          ! at the geometrically fixed height handle everything - a movement of material !
          ! surface is implicitly treated automatically. Note that in terms of turbulent !
          ! heat and moisture fluxes at model interfaces, both the cumulus updraft  mass !
          ! flux and penetratively entraining mass flux play the same role -both of them !
          ! warms and dries the 'kbup(i)' layer, cools and moistens the 'kpen(i)' layer,   and !
          ! cools and moistens any intervening layers between 'kbup(i)' and 'kpen(i)' layers.  !
          ! It is important to note these identical roles on turbulent heat and moisture !
          ! fluxes of 'umf(:, i)' and 'emf(:, i)'.                                                   !
          ! When 'kbup(i)' is a stratocumulus-topped PBL top interface,  increase of 'rpen' !
          ! is likely to strongly diffuse stratocumulus top interface,  resulting in the !
          ! reduction of cloud fraction. In this sense, the 'kbup(i)' interface has a  very !
          ! important meaning and role : across the 'kbup(i)' interface, strong penetrative !
          ! entrainment occurs, thus any sharp gradient properties across that interface !
          ! are easily diffused through strong mass exchange. Thus, an initialization of !
          ! 'kbup(i)' (and also 'kpen(i)') should be done very cautiously as mentioned before. !
          ! In order to prevent this stron diffusion for the shallow cumulus convection  !
          ! based at the Sc top, it seems to be good to initialize 'kbup(i) = krel(i)', rather !
          ! that 'kbup(i) = krel(i)-1'.                                                        !
          ! ---------------------------------------------------------------------------- !

       end do

       !------------------------------------------------------------------ !
       !                                                                   !
       ! Compute turbulent heat, moisture, momentum flux at all interfaces !
       !                                                                   !
       !------------------------------------------------------------------ !
       ! It is very important to note that in calculating turbulent fluxes !
       ! below, we must not double count turbulent flux at any interefaces.!
       ! In the below, turbulent fluxes at the interfaces (interface index !
       ! k) are calculated by the following 4 blocks in consecutive order: !
       !                                                                   !
       ! (1) " 0 <= k <= kinv(i) - 1 "  : PBL fluxes.                         !
       !     From 'fluxbelowinv' using reconstructed PBL height. Currently,!
       !     the reconstructed PBLs are independently calculated for  each !
       !     individual conservative scalar variables ( qt, thl, u, v ) in !
       !     each 'fluxbelowinv',  instead of being uniquely calculated by !
       !     using thvl. Turbulent flux at the surface is assumed to be 0. !
       ! (2) " kinv(i) <= k <= krel(i) - 1 " : Non-buoyancy sorting fluxes       !
       !     Assuming cumulus mass flux  and cumulus updraft thermodynamic !
       !     properties (except u, v which are modified by the PGFc during !
       !     upward motion) are conserved during a updraft motion from the !
       !     PBL top interface to the release level. If these layers don't !
       !     exist (e,g, when 'krel(i) = kinv(i)'), then  current routine do not !
       !     perform this routine automatically. So I don't need to modify !
       !     anything.                                                     !
       ! (3) " krel(i) <= k <= kbup(i) - 1 " : Buoyancy sorting fluxes           !
       !     From laterally entraining-detraining buoyancy sorting plumes. !
       ! (4) " kbup(i) <= k < kpen(i)-1 " : Penetrative entrainment fluxes       !
       !     From penetratively entraining plumes,                         !
       !                                                                   !
       ! In case of normal situation, turbulent interfaces  in each groups !
       ! are mutually independent of each other. Thus double flux counting !
       ! or ambiguous flux counting requiring the choice among the above 4 !
       ! groups do not occur normally. However, in case that cumulus plume !
       ! could not completely overcome the buoyancy barrier just above the !
       ! PBL top interface and so 'kbup(i) = krel(i)' (.forcedCu(i)=.true.) ( here, !
       ! it can be either 'kpen(i) = krel(i)' as the initialization, or ' kpen(i) > !
       ! krel(i)' if cumulus updraft just penetrated over the top of  release !
       ! layer ). If this happens, we should be very careful in organizing !
       ! the sequence of the 4 calculation routines above -  note that the !
       ! routine located at the later has the higher priority.  Additional !
       ! feature I must consider is that when 'kbup(i) = kinv(i) - 1' (this is a !
       ! combined situation of 'kbup(i)=krel(i)-1' & 'krel(i) = kinv(i)' when I  chose !
       ! 'kbup(i)=krel(i)-1' instead of current choice of 'kbup(i)=krel(i)'), a strong !
       ! penetrative entrainment fluxes exists at the PBL top interface, & !
       ! all of these fluxes are concentrated (deposited) within the layer !
       ! just below PBL top interface (i.e., 'kinv(i)-1' layer). On the other !
       ! hand, in case of 'fluxbelowinv', only the compensating subsidence !
       ! effect is concentrated in the 'kinv(i)-1' layer and 'pure' turbulent !
       ! heat and moisture fluxes ( 'pure' means the fluxes not associated !
       ! with compensating subsidence) are linearly distributed throughout !
       ! the whole PBL. Thus different choice of the above flux groups can !
       ! produce very different results. Output variable should be written !
       ! consistently to the choice of computation sequences.              !
       ! When the case of 'kbup(i) = krel(-1, i)' happens,another way to dealing !
       ! with this case is to simply ' exit ' the whole cumulus convection !
       ! calculation without performing any cumulus convection.     We can !
       ! choose this approach by specifying a condition in the  'Filtering !
       ! of unreasonable cumulus adjustment' just after 'iter_scaleh'. But !
       ! this seems not to be a good choice (although this choice was used !
       ! previous code ), since it might arbitrary damped-out  the shallow !
       ! cumulus convection over the continent land, where shallow cumulus !
       ! convection tends to be negatively buoyant.                        !
       ! ----------------------------------------------------------------- !

       ! --------------------------------------------------- !
       ! 1. PBL fluxes :  0 <= k <= kinv(i) - 1                 !
       !    All the information necessary to reconstruct PBL !
       !    height are passed to 'fluxbelowinv'.             !
       ! --------------------------------------------------- !

       xsrc(i)  = qtsrc(i)
       xmean(i) = qt0(kinv(i), i)
       xtop(i)  = qt0(kinv(i)+1, i) + ssqt0(kinv(i)+1, i) * ( ps0(kinv(i), i)   - p0(kinv(i)+1, i) )
       xbot(i)  = qt0(kinv(i)-1, i) + ssqt0(kinv(i)-1, i) * ( ps0(kinv(i)-1, i) - p0(kinv(i)-1, i) )
       call fluxbelowinv( cbmf(i), ps0(0:mkx, i), mkx, kinv(i), dt, xsrc(i), xmean(i), xtop(i), xbot(i), xflx(:, i) )
       qtflx(0:kinv(i)-1, i) = xflx(0:kinv(i)-1, i)

       xsrc(i)  = thlsrc(i)
       xmean(i) = thl0(kinv(i), i)
       xtop(i)  = thl0(kinv(i)+1, i) + ssthl0(kinv(i)+1, i) * ( ps0(kinv(i), i)   - p0(kinv(i)+1, i) )
       xbot(i)  = thl0(kinv(i)-1, i) + ssthl0(kinv(i)-1, i) * ( ps0(kinv(i)-1, i) - p0(kinv(i)-1, i) )
       call fluxbelowinv( cbmf(i), ps0(0:mkx, i), mkx, kinv(i), dt, xsrc(i), xmean(i), xtop(i), xbot(i), xflx(:, i) )
       slflx(0:kinv(i)-1, i) = cp * exns0(0:kinv(i)-1, i) * xflx(0:kinv(i)-1, i)

       xsrc(i)  = usrc(i)
       xmean(i) = u0(kinv(i), i)
       xtop(i)  = u0(kinv(i)+1, i) + ssu0(kinv(i)+1, i) * ( ps0(kinv(i), i)   - p0(kinv(i)+1, i) )
       xbot(i)  = u0(kinv(i)-1, i) + ssu0(kinv(i)-1, i) * ( ps0(kinv(i)-1, i) - p0(kinv(i)-1, i) )
       call fluxbelowinv( cbmf(i), ps0(0:mkx, i), mkx, kinv(i), dt, xsrc(i), xmean(i), xtop(i), xbot(i), xflx(:, i) )
       uflx(0:kinv(i)-1, i) = xflx(0:kinv(i)-1, i)

       xsrc(i)  = vsrc(i)
       xmean(i) = v0(kinv(i), i)
       xtop(i)  = v0(kinv(i)+1, i) + ssv0(kinv(i)+1, i) * ( ps0(kinv(i), i)   - p0(kinv(i)+1, i) )
       xbot(i)  = v0(kinv(i)-1, i) + ssv0(kinv(i)-1, i) * ( ps0(kinv(i)-1, i) - p0(kinv(i)-1, i) )
       call fluxbelowinv( cbmf(i), ps0(0:mkx, i), mkx, kinv(i), dt, xsrc(i), xmean(i), xtop(i), xbot(i), xflx(:, i) )
       vflx(0:kinv(i)-1, i) = xflx(0:kinv(i)-1, i)

       do m = 1, ncnst
          xsrc(i)  = trsrc(m, i)
          xmean(i) = tr0(kinv(i),m, i)
          xtop(i)  = tr0(kinv(i)+1,m, i) + sstr0(kinv(i)+1,m, i) * ( ps0(kinv(i), i)   - p0(kinv(i)+1, i) )
          xbot(i)  = tr0(kinv(i)-1,m, i) + sstr0(kinv(i)-1,m, i) * ( ps0(kinv(i)-1, i) - p0(kinv(i)-1, i) )
          call fluxbelowinv( cbmf(i), ps0(0:mkx, i), mkx, kinv(i), dt, xsrc(i), xmean(i), xtop(i), xbot(i), xflx(:, i) )
          trflx(0:kinv(i)-1,m, i) = xflx(0:kinv(i)-1, i)
       enddo

       ! -------------------------------------------------------------- !
       ! 2. Non-buoyancy sorting fluxes : kinv(i) <= k <= krel(i) - 1         !
       !    Note that when 'krel(i) = kinv(i)', below block is never executed !
       !    as in a desirable, expected way ( but I must check  if this !
       !    is the case ). The non-buoyancy sorting fluxes are computed !
       !    only when 'krel(i) > kinv(i)'.                                    !
       ! -------------------------------------------------------------- !

       uplus(i) = 0._r8
       vplus(i) = 0._r8
       do k = kinv(i), krel(i) - 1
          kp1 = k + 1
          qtflx(k, i) = cbmf(i) * ( qtsrc(i)  - (  qt0(kp1, i) +  ssqt0(kp1, i) * ( ps0(k, i) - p0(kp1, i) ) ) )
          slflx(k, i) = cbmf(i) * ( thlsrc(i) - ( thl0(kp1, i) + ssthl0(kp1, i) * ( ps0(k, i) - p0(kp1, i) ) ) ) * cp * exns0(k, i)
          uplus(i)    = uplus(i) + PGFc * ssu0(k, i) * ( ps0(k, i) - ps0(k-1, i) )
          vplus(i)    = vplus(i) + PGFc * ssv0(k, i) * ( ps0(k, i) - ps0(k-1, i) )
          uflx(k, i)  = cbmf(i) * ( usrc(i) + uplus(i) -  (  u0(kp1, i)  +   ssu0(kp1, i) * ( ps0(k, i) - p0(kp1, i) ) ) )
          vflx(k, i)  = cbmf(i) * ( vsrc(i) + vplus(i) -  (  v0(kp1, i)  +   ssv0(kp1, i) * ( ps0(k, i) - p0(kp1, i) ) ) )
          do m = 1, ncnst
             trflx(k,m, i) = cbmf(i) * ( trsrc(m, i)  - (  tr0(kp1,m, i) +  sstr0(kp1,m, i) * ( ps0(k, i) - p0(kp1, i) ) ) )
          enddo
       end do

       ! ------------------------------------------------------------------------ !
       ! 3. Buoyancy sorting fluxes : krel(i) <= k <= kbup(i) - 1                       !
       !    In case that 'kbup(i) = krel(i) - 1 ' ( or even in case 'kbup(i) = krel(i)' ),    !
       !    buoyancy sorting fluxes are not calculated, which is consistent,      !
       !    desirable feature.                                                    !
       ! ------------------------------------------------------------------------ !

       do k = krel(i), kbup(i) - 1
          kp1 = k + 1
          slflx(k, i) = cp * exns0(k, i) * umf(k, i) * ( thlu(k, i) - ( thl0(kp1, i) + ssthl0(kp1, i) * ( ps0(k, i) - p0(kp1, i) ) ) )
          qtflx(k, i) = umf(k, i) * ( qtu(k, i) - ( qt0(kp1, i) + ssqt0(kp1, i) * ( ps0(k, i) - p0(kp1, i) ) ) )
          uflx(k, i)  = umf(k, i) * ( uu(k, i) - ( u0(kp1, i) + ssu0(kp1, i) * ( ps0(k, i) - p0(kp1, i) ) ) )
          vflx(k, i)  = umf(k, i) * ( vu(k, i) - ( v0(kp1, i) + ssv0(kp1, i) * ( ps0(k, i) - p0(kp1, i) ) ) )
          do m = 1, ncnst
             trflx(k,m, i) = umf(k, i) * ( tru(k,m, i) - ( tr0(kp1,m, i) + sstr0(kp1,m, i) * ( ps0(k, i) - p0(kp1, i) ) ) )
          enddo
       end do

       ! ------------------------------------------------------------------------- !
       ! 4. Penetrative entrainment fluxes : kbup(i) <= k <= kpen(i) - 1                 !
       !    The only confliction that can happen is when 'kbup(i) = kinv(i)-1'. For this !
       !    case, turbulent flux at kinv(i)-1 is calculated  both from 'fluxbelowinv' !
       !    and here as penetrative entrainment fluxes.  Since penetrative flux is !
       !    calculated later, flux at 'kinv(i) - 1 ' will be that of penetrative flux.!
       !    However, turbulent flux calculated at 'kinv(i) - 1' from penetrative entr.!
       !    is less attractable,  since more reasonable turbulent flux at 'kinv(i)-1' !
       !    should be obtained from 'fluxbelowinv', by considering  re-constructed !
       !    inversion base height. This conflicting problem can be solved if we can!
       !    initialize 'kbup(i) = krel(i)', instead of kbup(i) = krel(i) - 1. This choice seems!
       !    to be more reasonable since it is not conflicted with 'fluxbelowinv' in!
       !    calculating fluxes at 'kinv(i) - 1' ( for this case, flux at 'kinv(i)-1' is  !
       !    always from 'fluxbelowinv' ), and flux at 'krel(i)-1' is calculated from  !
       !    the non-buoyancy sorting flux without being competed with penetrative  !
       !    entrainment fluxes. Even when we use normal cumulus flux instead of    !
       !    penetrative entrainment fluxes at 'kbup(i) <= k <= kpen(i)-1' interfaces,    !
       !    the initialization of kbup(i)=krel(i) perfectly works without any conceptual !
       !    confliction. Thus it seems to be much better to choose 'kbup(i) = krel(i)'   !
       !    initialization of 'kbup(i)', which is current choice.                     !
       !    Note that below formula uses conventional updraft cumulus fluxes for   !
       !    shallow cumulus which did not overcome the first buoyancy barrier above!
       !    PBL top while uses penetrative entrainment fluxes for the other cases  !
       !    'kbup(i) <= k <= kpen(i)-1' interfaces. Depending on cases, however, I can   !
       !    selelct different choice.                                              !
       ! ------------------------------------------------------------------------------------------------------------------ !
       !   if( forcedCu(i) ) then                                                                                              !
       !       slflx(k, i) = cp * exns0(k, i) * umf(k, i) * ( thlu(k, i) - ( thl0(kp1, i) + ssthl0(kp1, i) * ( ps0(k, i) - p0(kp1, i) ) ) )         !
       !       qtflx(k, i) =                 umf(k, i) * (  qtu(k, i) - (  qt0(kp1, i) +  ssqt0(kp1, i) * ( ps0(k, i) - p0(kp1, i) ) ) )         !
       !       uflx(k, i)  =                 umf(k, i) * (   uu(k, i) - (   u0(kp1, i) +   ssu0(kp1, i) * ( ps0(k, i) - p0(kp1, i) ) ) )         !
       !       vflx(k, i)  =                 umf(k, i) * (   vu(k, i) - (   v0(kp1, i) +   ssv0(kp1, i) * ( ps0(k, i) - p0(kp1, i) ) ) )         !
       !       do m = 1, ncnst                                                                                              !
       !          trflx(k,m, i) = umf(k, i) * ( tru(k,m, i) - ( tr0(kp1,m, i) + sstr0(kp1,m, i) * ( ps0(k, i) - p0(kp1, i) ) ) )                 !
       !       enddo                                                                                                        !
       !   else                                                                                                             !
       !       slflx(k, i) = cp * exns0(k, i) * emf(k, i) * ( thlu_emf(k, i) - ( thl0(k, i) + ssthl0(k, i) * ( ps0(k, i) - p0(k, i) ) ) )           !
       !       qtflx(k, i) =                 emf(k, i) * (  qtu_emf(k, i) - (  qt0(k, i) +  ssqt0(k, i) * ( ps0(k, i) - p0(k, i) ) ) )           !
       !       uflx(k, i)  =                 emf(k, i) * (   uu_emf(k, i) - (   u0(k, i) +   ssu0(k, i) * ( ps0(k, i) - p0(k, i) ) ) )           !
       !       vflx(k, i)  =                 emf(k, i) * (   vu_emf(k, i) - (   v0(k, i) +   ssv0(k, i) * ( ps0(k, i) - p0(k, i) ) ) )           !
       !       do m = 1, ncnst                                                                                              !
       !          trflx(k,m, i) = emf(k, i) * ( tru_emf(k,m, i) - ( tr0(k,m, i) + sstr0(k,m, i) * ( ps0(k, i) - p0(k, i) ) ) )                   !
       !       enddo                                                                                                        !
       !   endif                                                                                                            !
       !                                                                                                                    !
       !   if( use_uppenent ) then ! Combined Updraft + Penetrative Entrainment Flux                                        !
       !       slflx(k, i) = cp * exns0(k, i) * umf(k, i) * ( thlu(k, i)     - ( thl0(kp1, i) + ssthl0(kp1, i) * ( ps0(k, i) - p0(kp1, i) ) ) ) + & !
       !                  cp * exns0(k, i) * emf(k, i) * ( thlu_emf(k, i) - (   thl0(k, i) +   ssthl0(k, i) * ( ps0(k, i) - p0(k, i) ) ) )       !
       !       qtflx(k, i) =                 umf(k, i) * (  qtu(k, i)     - (  qt0(kp1, i) +  ssqt0(kp1, i) * ( ps0(k, i) - p0(kp1, i) ) ) ) + & !
       !                                  emf(k, i) * (  qtu_emf(k, i) - (    qt0(k, i) +    ssqt0(k, i) * ( ps0(k, i) - p0(k, i) ) ) )       !
       !       uflx(k, i)  =                 umf(k, i) * (   uu(k, i)     - (   u0(kp1, i) +   ssu0(kp1, i) * ( ps0(k, i) - p0(kp1, i) ) ) ) + & !
       !                                  emf(k, i) * (   uu_emf(k, i) - (     u0(k, i) +     ssu0(k, i) * ( ps0(k, i) - p0(k, i) ) ) )       !
       !       vflx(k, i)  =                 umf(k, i) * (   vu(k, i)     - (   v0(kp1, i) +   ssv0(kp1, i) * ( ps0(k, i) - p0(kp1, i) ) ) ) + & !
       !                                  emf(k, i) * (   vu_emf(k, i) - (     v0(k, i) +     ssv0(k, i) * ( ps0(k, i) - p0(k, i) ) ) )       !
       !       do m = 1, ncnst                                                                                              !
       !          trflx(k,m, i) = umf(k, i) * ( tru(k,m, i) - ( tr0(kp1,m, i) + sstr0(kp1,m, i) * ( ps0(k, i) - p0(kp1, i) ) ) ) + &             !
       !                       emf(k, i) * ( tru_emf(k,m, i) - ( tr0(k,m, i) + sstr0(k,m, i) * ( ps0(k, i) - p0(k, i) ) ) )                   !
       !       enddo                                                                                                        !
       ! ------------------------------------------------------------------------------------------------------------------ !

       do k = kbup(i), kpen(i) - 1
          kp1 = k + 1
          slflx(k, i) = cp * exns0(k, i) * emf(k, i) * ( thlu_emf(k, i) - ( thl0(k, i) + ssthl0(k, i) * ( ps0(k, i) - p0(k, i) ) ) )
          qtflx(k, i) =                 emf(k, i) * (  qtu_emf(k, i) - (  qt0(k, i) +  ssqt0(k, i) * ( ps0(k, i) - p0(k, i) ) ) )
          uflx(k, i)  =                 emf(k, i) * (   uu_emf(k, i) - (   u0(k, i) +   ssu0(k, i) * ( ps0(k, i) - p0(k, i) ) ) )
          vflx(k, i)  =                 emf(k, i) * (   vu_emf(k, i) - (   v0(k, i) +   ssv0(k, i) * ( ps0(k, i) - p0(k, i) ) ) )
          do m = 1, ncnst
             trflx(k,m, i) = emf(k, i) * ( tru_emf(k,m, i) - ( tr0(k,m, i) + sstr0(k,m, i) * ( ps0(k, i) - p0(k, i) ) ) )
          enddo
       end do

       ! ------------------------------------------- !
       ! Turn-off cumulus momentum flux as an option !
       ! ------------------------------------------- !

       if( .not. use_momenflx ) then
           uflx(0:mkx, i) = 0._r8
           vflx(0:mkx, i) = 0._r8
       endif

       ! -------------------------------------------------------- !
       ! Condensate tendency by compensating subsidence/upwelling !
       ! -------------------------------------------------------- !

       uemf(0:mkx, i)         = 0._r8
       do k = 0, kinv(i) - 2  ! Assume linear updraft mass flux within the PBL.
          uemf(k, i) = cbmf(i) * ( ps0(0, i) - ps0(k, i) ) / ( ps0(0, i) - ps0(kinv(i)-1, i) )
       end do
       uemf(kinv(i)-1:krel(i)-1, i) = cbmf(i)
       uemf(krel(i):kbup(i)-1, i)   = umf(krel(i):kbup(i)-1, i)
       uemf(kbup(i):kpen(i)-1, i)   = emf(kbup(i):kpen(i)-1, i) ! Only use penetrative entrainment flux consistently.

       comsub(1:mkx, i) = 0._r8
       do k = 1, kpen(i)
          comsub(k, i)  = 0.5_r8 * ( uemf(k, i) + uemf(k-1, i) )
       end do

       do k = 1, kpen(i)
          if( comsub(k, i) .ge. 0._r8 ) then
              if( k .eq. mkx ) then
                  thlten_sub(i) = 0._r8
                  qtten_sub(i)  = 0._r8
                  qlten_sub(i)  = 0._r8
                  qiten_sub(i)  = 0._r8
                  nlten_sub(i)  = 0._r8
                  niten_sub(i)  = 0._r8
              else
                  thlten_sub(i) = g * comsub(k, i) * ( thl0(k+1, i) - thl0(k, i) ) / ( p0(k, i) - p0(k+1, i) )
                  qtten_sub(i)  = g * comsub(k, i) * (  qt0(k+1, i) -  qt0(k, i) ) / ( p0(k, i) - p0(k+1, i) )
                  qlten_sub(i)  = g * comsub(k, i) * (  ql0(k+1, i) -  ql0(k, i) ) / ( p0(k, i) - p0(k+1, i) )
                  qiten_sub(i)  = g * comsub(k, i) * (  qi0(k+1, i) -  qi0(k, i) ) / ( p0(k, i) - p0(k+1, i) )
                  nlten_sub(i)  = g * comsub(k, i) * (  tr0(k+1,ixnumliq, i) -  tr0(k,ixnumliq, i) ) / ( p0(k, i) - p0(k+1, i) )
                  niten_sub(i)  = g * comsub(k, i) * (  tr0(k+1,ixnumice, i) -  tr0(k,ixnumice, i) ) / ( p0(k, i) - p0(k+1, i) )
              endif
          else
              if( k .eq. 1 ) then
                  thlten_sub(i) = 0._r8
                  qtten_sub(i)  = 0._r8
                  qlten_sub(i)  = 0._r8
                  qiten_sub(i)  = 0._r8
                  nlten_sub(i)  = 0._r8
                  niten_sub(i)  = 0._r8
              else
                  thlten_sub(i) = g * comsub(k, i) * ( thl0(k, i) - thl0(k-1, i) ) / ( p0(k-1, i) - p0(k, i) )
                  qtten_sub(i)  = g * comsub(k, i) * (  qt0(k, i) -  qt0(k-1, i) ) / ( p0(k-1, i) - p0(k, i) )
                  qlten_sub(i)  = g * comsub(k, i) * (  ql0(k, i) -  ql0(k-1, i) ) / ( p0(k-1, i) - p0(k, i) )
                  qiten_sub(i)  = g * comsub(k, i) * (  qi0(k, i) -  qi0(k-1, i) ) / ( p0(k-1, i) - p0(k, i) )
                  nlten_sub(i)  = g * comsub(k, i) * (  tr0(k,ixnumliq, i) -  tr0(k-1,ixnumliq, i) ) / ( p0(k-1, i) - p0(k, i) )
                  niten_sub(i)  = g * comsub(k, i) * (  tr0(k,ixnumice, i) -  tr0(k-1,ixnumice, i) ) / ( p0(k-1, i) - p0(k, i) )
              endif
          endif
          thl_prog(i) = thl0(k, i) + thlten_sub(i) * dt
          qt_prog(i)  = max( qt0(k, i) + qtten_sub(i) * dt, 1.e-12_r8 )
          call conden(p0(k, i),thl_prog(i),qt_prog(i),thj(i),qvj(i),qlj(i),qij(i),qse(i),id_check)
          if( id_check .eq. 1 ) then
              id_exit(i) = .true.
              go to 333
          endif
        ! qlten_sink(k, i) = ( qlj(i) - ql0(k, i) ) / dt
        ! qiten_sink(k, i) = ( qij(i) - qi0(k, i) ) / dt
          qlten_sink(k, i) = max( qlten_sub(i), - ql0(k, i) / dt ) ! For consistency with prognostic macrophysics scheme
          qiten_sink(k, i) = max( qiten_sub(i), - qi0(k, i) / dt ) ! For consistency with prognostic macrophysics scheme
          nlten_sink(k, i) = max( nlten_sub(i), - tr0(k,ixnumliq, i) / dt )
          niten_sink(k, i) = max( niten_sub(i), - tr0(k,ixnumice, i) / dt )
       end do

       ! --------------------------------------------- !
       !                                               !
       ! Calculate convective tendencies at each layer !
       !                                               !
       ! --------------------------------------------- !

       ! ----------------- !
       ! Momentum tendency !
       ! ----------------- !

       do k = 1, kpen(i)
          km1 = k - 1
          uten(k, i) = ( uflx(km1, i) - uflx(k, i) ) * g / dp0(k, i)
          vten(k, i) = ( vflx(km1, i) - vflx(k, i) ) * g / dp0(k, i)
          uf(k, i)   = u0(k, i) + uten(k, i) * dt
          vf(k, i)   = v0(k, i) + vten(k, i) * dt
        ! do m = 1, ncnst
        !    trten(k,m, i) = ( trflx(km1,m, i) - trflx(k,m, i) ) * g / dp0(k, i)
        !  ! Limit trten(k,m, i) such that negative value is not developed.
        !  ! This limitation does not conserve grid-mean tracers and future
        !  ! refinement is required for tracer-conserving treatment.
        !    trten(k,m, i) = max(trten(k,m, i),-tr0(k,m, i)/dt)
        ! enddo
       end do

       ! ----------------------------------------------------------------- !
       ! Tendencies of thermodynamic variables.                            !
       ! This part requires a careful treatment of bulk cloud microphysics.!
       ! Relocations of 'precipitable condensates' either into the surface !
       ! or into the tendency of 'krel(i)' layer will be performed just after !
       ! finishing the below 'do-loop'.                                    !
       ! ----------------------------------------------------------------- !

       rliq(i)    = 0._r8
       rainflx(i) = 0._r8
       snowflx(i) = 0._r8

       do k = 1, kpen(i)

          km1 = k - 1

          ! ------------------------------------------------------------------------------ !
          ! Compute 'slten(:, i)', 'qtten(:, i)', 'qvten(:, i)', 'qlten(:, i)', 'qiten(:, i)', and 'sten(:, i)'                !
          !                                                                                !
          ! Key assumptions made in this 'cumulus scheme' are :                            !
          ! 1. Cumulus updraft expels condensate into the environment at the top interface !
          !    of each layer. Note that in addition to this expel process ('source' term), !
          !    cumulus updraft can modify layer mean condensate through normal detrainment !
          !    forcing or compensating subsidence.                                         !
          ! 2. Expelled water can be either 'sustaining' or 'precipitating' condensate. By !
          !    definition, 'suataining condensate' will remain in the layer where it was   !
          !    formed, while 'precipitating condensate' will fall across the base of the   !
          !    layer where it was formed.                                                  !
          ! 3. All precipitating condensates are assumed to fall into the release layer or !
          !    ground as soon as it was formed without being evaporated during the falling !
          !    process down to the desinated layer ( either release layer of surface ).    !
          ! ------------------------------------------------------------------------------ !

          ! ------------------------------------------------------------------------- !
          ! 'dwten(k, i)','diten(k, i)' : Production rate of condensate  within the layer k !
          !      [ kg/kg/s ]        by the expels of condensate from cumulus updraft. !
          ! It is important to note that in terms of moisture tendency equation, this !
          ! is a 'source' term of enviromental 'qt'.  More importantly,  these source !
          ! are already counted in the turbulent heat and moisture fluxes we computed !
          ! until now, assuming all the expelled condensate remain in the layer where !
          ! it was formed. Thus, in calculation of 'qtten(:, i)' and 'slten(:, i)' below, we MUST !
          ! NOT add or subtract these terms explicitly in order not to double or miss !
          ! count, unless some expelled condensates fall down out of the layer.  Note !
          ! this falling-down process ( i.e., precipitation process ) and  associated !
          ! 'qtten(:, i)' and 'slten(:, i)' and production of surface precipitation flux  will be !
          ! treated later in 'zm_conv_evap' in 'convect_shallow_tend' subroutine.     !
          ! In below, we are converting expelled cloud condensate into correct unit.  !
          ! I found that below use of '0.5 * (umf(k-1, i) + umf(k, i))' causes conservation !
          ! errors at some columns in global simulation. So, I returned to originals. !
          ! This will cause no precipitation flux at 'kpen(i)' layer since umf(kpen(i), i)=0.  !
          ! ------------------------------------------------------------------------- !

          dwten(k, i) = dwten(k, i) * 0.5_r8 * ( umf(k-1, i) + umf(k, i) ) * g / dp0(k, i) ! [ kg/kg/s ]
          diten(k, i) = diten(k, i) * 0.5_r8 * ( umf(k-1, i) + umf(k, i) ) * g / dp0(k, i) ! [ kg/kg/s ]

          ! dwten(k, i) = dwten(k, i) * umf(k, i) * g / dp0(k, i) ! [ kg/kg/s ]
          ! diten(k, i) = diten(k, i) * umf(k, i) * g / dp0(k, i) ! [ kg/kg/s ]

          ! --------------------------------------------------------------------------- !
          ! 'qrten(k, i)','qsten(k, i)' : Production rate of rain and snow(i) within the layer k !
          !     [ kg/kg/s ]         by cumulus expels of condensates to the environment.!
          ! This will be falled-out of the layer where it was formed and will be dumped !
          ! dumped into the release layer assuming that there is no evaporative cooling !
          ! while precipitable condensate moves to the relaes level. This is reasonable !
          ! assumtion if cumulus is purely vertical and so the path along which precita !
          ! ble condensate falls is fully saturared. This 're-allocation' process of    !
          ! precipitable condensate into the release layer is fully described in this   !
          ! convection scheme. After that, the dumped water into the release layer will !
          ! falling down across the base of release layer ( or LCL, if  exact treatment !
          ! is required ) and will be allowed to be evaporated in layers below  release !
          ! layer, and finally non-zero surface precipitation flux will be calculated.  !
          ! This latter process will be separately treated 'zm_conv_evap' routine.      !
          ! --------------------------------------------------------------------------- !

          qrten(k, i) = frc_rasn * dwten(k, i)
          qsten(k, i) = frc_rasn * diten(k, i)

          ! ----------------------------------------------------------------------- !
          ! 'rainflx(i)','snowflx(i)' : Cumulative rain and snow(i) flux integrated from the !
          !     [ kg/m2/s ]       release leyer to the 'kpen(i)' layer. Note that even !
          ! though wtw(kpen(i), i) < 0 (and umf(kpen(i), i) = 0) at the top interface of 'kpen(i)' !
          ! layer, 'dwten(kpen(i), i)' and diten(kpen(i), i)  were calculated after calculating !
          ! explicit cloud top height. Thus below calculation of precipitation flux !
          ! is correct. Note that  precipitating condensates are formed only in the !
          ! layers from 'krel(i)' to 'kpen(i)', including the two layers.                 !
          ! ----------------------------------------------------------------------- !

          rainflx(i) = rainflx(i) + qrten(k, i) * dp0(k, i) / g
          snowflx(i) = snowflx(i) + qsten(k, i) * dp0(k, i) / g

          ! ------------------------------------------------------------------------ !
          ! 'slten(k, i)','qtten(k, i)'                                                    !
          !  Note that 'slflx(k, i)' and 'qtflx(k, i)' we have calculated already included !
          !  all the contributions of (1) expels of condensate (dwten(k, i), diten(k, i)), !
          !  (2) mass detrainment ( delta * umf(:, i) * ( qtu(:, i) - qt ) ), & (3) compensating !
          !  subsidence ( M * dqt / dz ). Thus 'slflx(k, i)' and 'qtflx(k, i)' we computed !
          !  is a hybrid turbulent flux containing one part of 'source' term - expel !
          !  of condensate. In order to calculate 'slten(:, i)' and 'qtten(:, i)', we should add !
          !  additional 'source' term, if any. If the expelled condensate falls down !
          !  across the base of the layer, it will be another sink (negative source) !
          !  term.  Note also that we included frictional heating terms in the below !
          !  calculation of 'slten(:, i)'.                                                 !
          ! ------------------------------------------------------------------------ !

          slten(k, i) = ( slflx(km1, i) - slflx(k, i) ) * g / dp0(k, i)
          if( k .eq. 1 ) then
              slten(k, i) = slten(k, i) - g / 4._r8 / dp0(k, i) * (                            &
                                    uflx(k, i)*(uf(k+1, i) - uf(k, i) + u0(k+1, i) - u0(k, i)) +     &
                                    vflx(k, i)*(vf(k+1, i) - vf(k, i) + v0(k+1, i) - v0(k, i)))
          elseif( k .ge. 2 .and. k .le. kpen(i)-1 ) then
              slten(k, i) = slten(k, i) - g / 4._r8 / dp0(k, i) * (                            &
                                    uflx(k, i)*(uf(k+1, i) - uf(k, i) + u0(k+1, i) - u0(k, i)) +     &
                                    uflx(k-1, i)*(uf(k, i) - uf(k-1, i) + u0(k, i) - u0(k-1, i)) +   &
                                    vflx(k, i)*(vf(k+1, i) - vf(k, i) + v0(k+1, i) - v0(k, i)) +     &
                                    vflx(k-1, i)*(vf(k, i) - vf(k-1, i) + v0(k, i) - v0(k-1, i)))
          elseif( k .eq. kpen(i) ) then
              slten(k, i) = slten(k, i) - g / 4._r8 / dp0(k, i) * (                            &
                                    uflx(k-1, i)*(uf(k, i) - uf(k-1, i) + u0(k, i) - u0(k-1, i)) +   &
                                    vflx(k-1, i)*(vf(k, i) - vf(k-1, i) + v0(k, i) - v0(k-1, i)))
          endif
          qtten(k, i) = ( qtflx(km1, i) - qtflx(k, i) ) * g / dp0(k, i)

          ! ---------------------------------------------------------------------------- !
          ! Compute condensate tendency, including reserved condensate                   !
          ! We assume that eventual detachment and detrainment occurs in kbup(i) layer  due !
          ! to downdraft buoyancy sorting. In the layer above the kbup(i), only penetrative !
          ! entrainment exists. Penetrative entrained air is assumed not to contain any  !
          ! condensate.                                                                  !
          ! ---------------------------------------------------------------------------- !

          ! Compute in-cumulus condensate at the layer mid-point.

          if( k .lt. krel(i) .or. k .gt. kpen(i) ) then
              qlu_mid(i) = 0._r8
              qiu_mid(i) = 0._r8
              qlj(i)     = 0._r8
              qij(i)     = 0._r8
          elseif( k .eq. krel(i) ) then
              call conden(prel(i),thlu(krel(i)-1, i),qtu(krel(i)-1, i),thj(i),qvj(i),qlj(i),qij(i),qse(i),id_check)
              if( id_check .eq. 1 ) then
                  exit_conden(i) = 1._r8
                  id_exit(i) = .true.
                  go to 333
              endif
              qlubelow(i) = qlj(i)
              qiubelow(i) = qij(i)
              call conden(ps0(k, i),thlu(k, i),qtu(k, i),thj(i),qvj(i),qlj(i),qij(i),qse(i),id_check)
              if( id_check .eq. 1 ) then
                  exit_conden(i) = 1._r8
                  id_exit(i) = .true.
                  go to 333
              end if
              qlu_mid(i) = 0.5_r8 * ( qlubelow(i) + qlj(i) ) * ( prel(i) - ps0(k, i) )/( ps0(k-1, i) - ps0(k, i) )
              qiu_mid(i) = 0.5_r8 * ( qiubelow(i) + qij(i) ) * ( prel(i) - ps0(k, i) )/( ps0(k-1, i) - ps0(k, i) )
          elseif( k .eq. kpen(i) ) then
              call conden(ps0(k-1, i)+ppen(i),thlu_top(i),qtu_top(i),thj(i),qvj(i),qlj(i),qij(i),qse(i),id_check)
              if( id_check .eq. 1 ) then
                  exit_conden(i) = 1._r8
                  id_exit(i) = .true.
                  go to 333
              end if
              qlu_mid(i) = 0.5_r8 * ( qlubelow(i) + qlj(i) ) * ( -ppen(i) )        /( ps0(k-1, i) - ps0(k, i) )
              qiu_mid(i) = 0.5_r8 * ( qiubelow(i) + qij(i) ) * ( -ppen(i) )        /( ps0(k-1, i) - ps0(k, i) )
              qlu_top(i) = qlj(i)
              qiu_top(i) = qij(i)
          else
              call conden(ps0(k, i),thlu(k, i),qtu(k, i),thj(i),qvj(i),qlj(i),qij(i),qse(i),id_check)
              if( id_check .eq. 1 ) then
                  exit_conden(i) = 1._r8
                  id_exit(i) = .true.
                  go to 333
              end if
              qlu_mid(i) = 0.5_r8 * ( qlubelow(i) + qlj(i) )
              qiu_mid(i) = 0.5_r8 * ( qiubelow(i) + qij(i) )
          endif
          qlubelow(i) = qlj(i)
          qiubelow(i) = qij(i)

          ! 1. Sustained Precipitation

          qc_l(k, i) = ( 1._r8 - frc_rasn ) * dwten(k, i) ! [ kg/kg/s ]
          qc_i(k, i) = ( 1._r8 - frc_rasn ) * diten(k, i) ! [ kg/kg/s ]

          ! 2. Detrained Condensate

          if( k .le. kbup(i) ) then
              qc_l(k, i) = qc_l(k, i) + g * 0.5_r8 * ( umf(k-1, i) + umf(k, i) ) * fdr(k, i) * qlu_mid(i) ! [ kg/kg/s ]
              qc_i(k, i) = qc_i(k, i) + g * 0.5_r8 * ( umf(k-1, i) + umf(k, i) ) * fdr(k, i) * qiu_mid(i) ! [ kg/kg/s ]
              qc_lm(i)   =         - g * 0.5_r8 * ( umf(k-1, i) + umf(k, i) ) * fdr(k, i) * ql0(k, i)
              qc_im(i)   =         - g * 0.5_r8 * ( umf(k-1, i) + umf(k, i) ) * fdr(k, i) * qi0(k, i)
            ! Below 'nc_lm(i)', 'nc_im(i)' should be used only when frc_rasn = 1.
              nc_lm(i)   =         - g * 0.5_r8 * ( umf(k-1, i) + umf(k, i) ) * fdr(k, i) * tr0(k,ixnumliq, i)
              nc_im(i)   =         - g * 0.5_r8 * ( umf(k-1, i) + umf(k, i) ) * fdr(k, i) * tr0(k,ixnumice, i)
          else
              qc_lm(i)   = 0._r8
              qc_im(i)   = 0._r8
              nc_lm(i)   = 0._r8
              nc_im(i)   = 0._r8
          endif

          ! 3. Detached Updraft

          if( k .eq. kbup(i) ) then
              qc_l(k, i) = qc_l(k, i) + g * umf(k, i) * qlj(i)     / ( ps0(k-1, i) - ps0(k, i) ) ! [ kg/kg/s ]
              qc_i(k, i) = qc_i(k, i) + g * umf(k, i) * qij(i)     / ( ps0(k-1, i) - ps0(k, i) ) ! [ kg/kg/s ]
              qc_lm(i)   = qc_lm(i)   - g * umf(k, i) * ql0(k, i)  / ( ps0(k-1, i) - ps0(k, i) ) ! [ kg/kg/s ]
              qc_im(i)   = qc_im(i)   - g * umf(k, i) * qi0(k, i)  / ( ps0(k-1, i) - ps0(k, i) ) ! [ kg/kg/s ]
              nc_lm(i)   = nc_lm(i)   - g * umf(k, i) * tr0(k,ixnumliq, i)  / ( ps0(k-1, i) - ps0(k, i) ) ! [ kg/kg/s ]
              nc_im(i)   = nc_im(i)   - g * umf(k, i) * tr0(k,ixnumice, i)  / ( ps0(k-1, i) - ps0(k, i) ) ! [ kg/kg/s ]
          endif

          ! 4. Cumulative Penetrative entrainment detrained in the 'kbup(i)' layer
          !    Explicitly compute the properties detrained penetrative entrained airs in k = kbup(i) layer.

          if( k .eq. kbup(i) ) then
              call conden(p0(k, i),thlu_emf(k, i),qtu_emf(k, i),thj(i),qvj(i),ql_emf_kbup(i),qi_emf_kbup(i),qse(i),id_check)
              if( id_check .eq. 1 ) then
                  id_exit(i) = .true.
                  go to 333
              endif
              if( ql_emf_kbup(i) .gt. 0._r8 ) then
                  nl_emf_kbup(i) = tru_emf(k,ixnumliq, i)
              else
                  nl_emf_kbup(i) = 0._r8
              endif
              if( qi_emf_kbup(i) .gt. 0._r8 ) then
                  ni_emf_kbup(i) = tru_emf(k,ixnumice, i)
              else
                  ni_emf_kbup(i) = 0._r8
              endif
              qc_lm(i)   = qc_lm(i)   - g * emf(k, i) * ( ql_emf_kbup(i) - ql0(k, i) ) / ( ps0(k-1, i) - ps0(k, i) ) ! [ kg/kg/s ]
              qc_im(i)   = qc_im(i)   - g * emf(k, i) * ( qi_emf_kbup(i) - qi0(k, i) ) / ( ps0(k-1, i) - ps0(k, i) ) ! [ kg/kg/s ]
              nc_lm(i)   = nc_lm(i)   - g * emf(k, i) * ( nl_emf_kbup(i) - tr0(k,ixnumliq, i) ) / ( ps0(k-1, i) - ps0(k, i) ) ! [ kg/kg/s ]
              nc_im(i)   = nc_im(i)   - g * emf(k, i) * ( ni_emf_kbup(i) - tr0(k,ixnumice, i) ) / ( ps0(k-1, i) - ps0(k, i) ) ! [ kg/kg/s ]
          endif

          qlten_det(i)   = qc_l(k, i) + qc_lm(i)
          qiten_det(i)   = qc_i(k, i) + qc_im(i)

          ! --------------------------------------------------------------------------------- !
          ! 'qlten(k, i)','qiten(k, i)','qvten(k, i)','sten(k, i)'                                        !
          ! Note that falling of precipitation will be treated later.                         !
          ! The prevension of negative 'qv,ql,qi' will be treated later in positive_moisture. !
          ! --------------------------------------------------------------------------------- !

          if( use_expconten ) then
              if( use_unicondet ) then
                  qc_l(k, i) = 0._r8
                  qc_i(k, i) = 0._r8
                  qlten(k, i) = frc_rasn * dwten(k, i) + qlten_sink(k, i) + qlten_det(i)
                  qiten(k, i) = frc_rasn * diten(k, i) + qiten_sink(k, i) + qiten_det(i)
              else
                  qlten(k, i) = qc_l(k, i) + frc_rasn * dwten(k, i) + ( max( 0._r8, ql0(k, i) + ( qc_lm(i) + qlten_sink(k, i) ) * dt ) - ql0(k, i) ) / dt
                  qiten(k, i) = qc_i(k, i) + frc_rasn * diten(k, i) + ( max( 0._r8, qi0(k, i) + ( qc_im(i) + qiten_sink(k, i) ) * dt ) - qi0(k, i) ) / dt
                  trten(k,ixnumliq, i) = max( nc_lm(i) + nlten_sink(k, i), - tr0(k,ixnumliq, i) / dt )
                  trten(k,ixnumice, i) = max( nc_im(i) + niten_sink(k, i), - tr0(k,ixnumice, i) / dt )
              endif
          else
              if( use_unicondet ) then
                  qc_l(k, i) = 0._r8
                  qc_i(k, i) = 0._r8
              endif
              qlten(k, i) = dwten(k, i) + ( qtten(k, i) - dwten(k, i) - diten(k, i) ) * ( ql0(k, i) / qt0(k, i) )
              qiten(k, i) = diten(k, i) + ( qtten(k, i) - dwten(k, i) - diten(k, i) ) * ( qi0(k, i) / qt0(k, i) )
          endif

          qvten(k, i) = qtten(k, i) - qlten(k, i) - qiten(k, i)
          sten(k, i)  = slten(k, i) + xlv * qlten(k, i) + xls * qiten(k, i)

          ! -------------------------------------------------------------------------- !
          ! 'rliq(i)' : Verticall-integrated 'suspended cloud condensate'                 !
          !  [m/s]   This is so called 'reserved liquid water'  in other subroutines   !
          ! of CAM3, since the contribution of this term should not be included into   !
          ! the tendency of each layer or surface flux (precip(i))  within this cumulus   !
          ! scheme. The adding of this term to the layer tendency will be done inthe   !
          ! 'stratiform_tend', just after performing sediment process there.           !
          ! The main problem of these rather going-back-and-forth and stupid-seeming   !
          ! approach is that the sediment process of suspendened condensate will not   !
          ! be treated at all in the 'stratiform_tend'.                                !
          ! Note that 'precip(i)' [m/s] is vertically-integrated total 'rain+snow(i)' formed !
          ! from the cumulus updraft. Important : in the below, 1000 is rhoh2o ( water !
          ! density ) [ kg/m^3 ] used for unit conversion from [ kg/m^2/s ] to [ m/s ] !
          ! for use in stratiform.F90.                                                 !
          ! -------------------------------------------------------------------------- !

          qc(k, i)  =  qc_l(k, i) +  qc_i(k, i)
          rliq(i)   =  rliq(i)    + qc(k, i) * dp0(k, i) / g / 1000._r8    ! [ m/s ]

       end do

          precip(i)  =  rainflx(i) + snowflx(i)                       ! [ kg/m2/s ]
          snow(i)    =  snowflx(i)                                 ! [ kg/m2/s ]

       ! ---------------------------------------------------------------- !
       ! Now treats the 'evaporation' and 'melting' of rain ( qrten(:, i) ) and !
       ! snow(i) ( qsten(:, i) ) during falling process. Below algorithms are from !
       ! 'zm_conv_evap' but with some modification, which allows separate !
       ! treatment of 'rain' and 'snow(i)' condensates. Note that I included !
       ! the evaporation dynamics into the convection scheme for complete !
       ! development of cumulus scheme especially in association with the !
       ! implicit CIN closure. In compatible with this internal treatment !
       ! of evaporation, I should modify 'convect_shallow',  in such that !
       ! 'zm_conv_evap' is not performed when I choose UW PBL-Cu schemes. !
       ! ---------------------------------------------------------------- !

       evpint_rain(i)    = 0._r8
       evpint_snow(i)    = 0._r8
       flxrain(0:mkx, i) = 0._r8
       flxsnow(0:mkx, i) = 0._r8
       ntraprd(:mkx, i)  = 0._r8
       ntsnprd(:mkx, i)  = 0._r8

       do k = mkx, 1, -1  ! 'k' is a layer index : 'mkx'('1') is the top ('bottom') layer

          ! ----------------------------------------------------------------------------- !
          ! flxsntm(i) [kg/m2/s] : Downward snow(i) flux at the top of each layer after melting.!
          ! snowmlt(i) [kg/kg/s] : Snow melting tendency.                                    !
          ! Below allows melting of snow(i) when it goes down into the warm layer below.     !
          ! ----------------------------------------------------------------------------- !

          if( t0(k, i) .gt. 273.16_r8 ) then
              snowmlt(i) = max( 0._r8, flxsnow(k, i) * g / dp0(k, i) )
          else
              snowmlt(i) = 0._r8
          endif

          ! ----------------------------------------------------------------- !
          ! Evaporation rate of 'rain' and 'snow(i)' in the layer k, [ kg/kg/s ] !
          ! where 'rain' and 'snow(i)' are coming down from the upper layers.    !
          ! I used the same evaporative efficiency both for 'rain' and 'snow(i)'.!
          ! Note that evaporation is not allowed in the layers 'k >= krel(i)' by !
          ! assuming that inside of cumulus cloud, across which precipitation !
          ! is falling down, is fully saturated.                              !
          ! The asumptions in association with the 'evplimit_rain(snow(i), i)' are  !
          !   1. Do not allow evaporation to supersate the layer              !
          !   2. Do not evaporate more than the flux falling into the layer   !
          !   3. Total evaporation cannot exceed the input total surface flux !
          ! ----------------------------------------------------------------- !

          call t_startf('qsat_o')
#ifdef QSAT_WV
          call qsat_wv(t0(k, i), p0(k, i), es(i), qs(i))
#else
          call qsat_o(t0(k, i), p0(k, i), es(i), qs(i))
#endif
          call t_stopf('qsat_o')
          subsat(i) = max( ( 1._r8 - qv0(k, i)/qs(i) ), 0._r8 )
          if( noevap_krelkpen ) then
              if( k .ge. krel(i) ) subsat(i) = 0._r8
          endif

          evprain(i)  = kevp * subsat(i) * sqrt(flxrain(k, i)+snowmlt(i)*dp0(k, i)/g)
          evpsnow(i)  = kevp * subsat(i) * sqrt(max(flxsnow(k, i)-snowmlt(i)*dp0(k, i)/g,0._r8))

          evplimit(i) = max( 0._r8, ( qw0_in(i,k) - qv0(k, i) ) / dt )

          evplimit_rain(i) = min( evplimit(i),      ( flxrain(k, i) + snowmlt(i) * dp0(k, i) / g ) * g / dp0(k, i) )
          evplimit_rain(i) = min( evplimit_rain(i), ( rainflx(i) - evpint_rain(i) ) * g / dp0(k, i) )
          evprain(i) = max(0._r8,min( evplimit_rain(i), evprain(i) ))

          evplimit_snow(i) = min( evplimit(i),   max( flxsnow(k, i) - snowmlt(i) * dp0(k, i) / g , 0._r8 ) * g / dp0(k, i) )
          evplimit_snow(i) = min( evplimit_snow(i), ( snowflx(i) - evpint_snow(i) ) * g / dp0(k, i) )
          evpsnow(i) = max(0._r8,min( evplimit_snow(i), evpsnow(i) ))

          if( ( evprain(i) + evpsnow(i) ) .gt. evplimit(i) ) then
                tmp1(i) = evprain(i) * evplimit(i) / ( evprain(i) + evpsnow(i) )
                tmp2(i) = evpsnow(i) * evplimit(i) / ( evprain(i) + evpsnow(i) )
                evprain(i) = tmp1(i)
                evpsnow(i) = tmp2(i)
          endif

          evapc(k, i) = evprain(i) + evpsnow(i)

          ! ------------------------------------------------------------- !
          ! Vertically-integrated evaporative fluxes of 'rain' and 'snow(i)' !
          ! ------------------------------------------------------------- !

          evpint_rain(i) = evpint_rain(i) + evprain(i) * dp0(k, i) / g
          evpint_snow(i) = evpint_snow(i) + evpsnow(i) * dp0(k, i) / g

          ! -------------------------------------------------------------- !
          ! Net 'rain' and 'snow(i)' production rate in the layer [ kg/kg/s ] !
          ! -------------------------------------------------------------- !

          ntraprd(k, i) = qrten(k, i) - evprain(i) + snowmlt(i)
          ntsnprd(k, i) = qsten(k, i) - evpsnow(i) - snowmlt(i)

          ! -------------------------------------------------------------------------------- !
          ! Downward fluxes of 'rain' and 'snow(i)' fluxes at the base of the layer [ kg/m2/s ] !
          ! Note that layer index increases with height.                                     !
          ! -------------------------------------------------------------------------------- !

          flxrain(k-1, i) = flxrain(k, i) + ntraprd(k, i) * dp0(k, i) / g
          flxsnow(k-1, i) = flxsnow(k, i) + ntsnprd(k, i) * dp0(k, i) / g
          flxrain(k-1, i) = max( flxrain(k-1, i), 0._r8 )
          if( flxrain(k-1, i) .eq. 0._r8 ) ntraprd(k, i) = -flxrain(k, i) * g / dp0(k, i)
          flxsnow(k-1, i) = max( flxsnow(k-1, i), 0._r8 )
          if( flxsnow(k-1, i) .eq. 0._r8 ) ntsnprd(k, i) = -flxsnow(k, i) * g / dp0(k, i)

          ! ---------------------------------- !
          ! Calculate thermodynamic tendencies !
          ! --------------------------------------------------------------------------- !
          ! Note that equivalently, we can write tendency formula of 'sten(:, i)' and 'slten(:, i)' !
          ! by 'sten(k, i)  = sten(k, i) - xlv*evprain(i)  - xls*evpsnow(i) - (xls-xlv)*snowmlt(i)' &  !
          !    'slten(k, i) = sten(k, i) - xlv*qlten(k, i) - xls*qiten(k, i)'.                      !
          ! The above formula is equivalent to the below formula. However below formula !
          ! is preferred since we have already imposed explicit constraint on 'ntraprd(:, i)' !
          ! and 'ntsnprd(:, i)' in case that flxrain(k-1, i) < 0 & flxsnow(k-1, i) < 0._r8          !
          ! Note : In future, I can elborate the limiting of 'qlten(:, i)','qvten(:, i)','qiten(:, i)'    !
          !        such that that energy and moisture conservation error is completely  !
          !        suppressed.                                                          !
          ! Re-storation to the positive condensate will be performed later below       !
          ! --------------------------------------------------------------------------- !

          qlten(k, i) = qlten(k, i) - qrten(k, i)
          qiten(k, i) = qiten(k, i) - qsten(k, i)
          qvten(k, i) = qvten(k, i) + evprain(i)  + evpsnow(i)
          qtten(k, i) = qlten(k, i) + qiten(k, i) + qvten(k, i)
          if( ( qv0(k, i) + qvten(k, i)*dt ) .lt. qmin(1) .or. &
              ( ql0(k, i) + qlten(k, i)*dt ) .lt. qmin(ixcldliq) .or. &
              ( qi0(k, i) + qiten(k, i)*dt ) .lt. qmin(ixcldice) ) then
               limit_negcon(i) = 1._r8
          end if
          sten(k, i)  = sten(k, i) - xlv*evprain(i)  - xls*evpsnow(i) - (xls-xlv)*snowmlt(i)
          slten(k, i) = sten(k, i) - xlv*qlten(k, i) - xls*qiten(k, i)

        !  slten(k, i) = slten(k, i) + xlv * ntraprd(k, i) + xls * ntsnprd(k, i)
        !  sten(k, i)  = slten(k, i) + xlv * qlten(k, i)   + xls * qiten(k, i)

       end do

       ! ------------------------------------------------------------- !
       ! Calculate final surface flux of precipitation, rain, and snow(i) !
       ! Convert unit to [m/s] for use in 'check_energy_chng'.         !
       ! ------------------------------------------------------------- !

       precip(i)  = ( flxrain(0, i) + flxsnow(0, i) ) / 1000._r8
       snow(i)    =   flxsnow(0, i) / 1000._r8

       ! --------------------------------------------------------------------------- !
       ! Until now, all the calculations are done completely in this shallow cumulus !
       ! scheme. If you want to use this cumulus scheme other than CAM3, then do not !
       ! perform below block. However, for compatible use with the other subroutines !
       ! in CAM3, I should subtract the effect of 'qc(k, i)' ('rliq(i)') from the tendency !
       ! equation in each layer, since this effect will be separately added later in !
       ! in 'stratiform_tend' just after performing sediment process there. In order !
       ! to be consistent with 'stratiform_tend', just subtract qc(k, i)  from tendency !
       ! equation of each layer, but do not add it to the 'precip(i)'. Apprently,  this !
       ! will violate energy and moisture conservations.    However, when performing !
       ! conservation check in 'tphysbc.F90' just after 'convect_shallow_tend',   we !
       ! will add 'qc(k, i)' ( rliq(i) ) to the surface flux term just for the purpose  of !
       ! passing the energy-moisture conservation check. Explicit adding-back of 'qc(:, i)'!
       ! to the individual layer tendency equation will be done in 'stratiform_tend' !
       ! after performing sediment process there. Simply speaking, in 'tphysbc' just !
       ! after 'convect_shallow_tend', we will dump 'rliq(i)' into surface as a  'rain' !
       ! in order to satisfy energy and moisture conservation, and  in the following !
       ! 'stratiform_tend', we will restore it back to 'qlten(k, i)' ( 'ice' will go to !
       ! 'water' there) from surface precipitation. This is a funny but conceptually !
       ! entertaining procedure. One concern I have for this complex process is that !
       ! output-writed stratiform precipitation amount will be underestimated due to !
       ! arbitrary subtracting of 'rliq(i)' in stratiform_tend, where                   !
       ! ' prec_str = prec_sed + prec_pcw - rliq(i)' and 'rliq(i)' is not real but fake.   !
       ! However, as shown in 'srfxfer.F90', large scale precipitation amount (PRECL)!
       ! that is writed-output is corrected written since in 'srfxfer.F90',  PRECL = !
       ! 'prec_sed + prec_pcw', without including 'rliq(i)'. So current code is correct.!
       ! Note also in 'srfxfer.F90', convective precipitation amount is 'PRECC =     !
       ! prec_zmc(i) + prec_cmf(i)' which is also correct.                           !
       ! --------------------------------------------------------------------------- !

       do k = 1, kpen(i)
          qtten(k, i) = qtten(k, i) - qc(k, i)
          qlten(k, i) = qlten(k, i) - qc_l(k, i)
          qiten(k, i) = qiten(k, i) - qc_i(k, i)
          slten(k, i) = slten(k, i) + ( xlv * qc_l(k, i) + xls * qc_i(k, i) )
          ! ---------------------------------------------------------------------- !
          ! Since all reserved condensates will be treated  as liquid water in the !
          ! 'check_energy_chng' & 'stratiform_tend' without an explicit conversion !
          ! algorithm, I should consider explicitly the energy conversions between !
          ! 'ice' and 'liquid' - i.e., I should convert 'ice' to 'liquid'  and the !
          ! necessary energy for this conversion should be subtracted from 'sten(:, i)'. !
          ! Without this conversion here, energy conservation error come out. Note !
          ! that there should be no change of 'qvten(k, i)'.                          !
          ! ---------------------------------------------------------------------- !
          sten(k, i)  = sten(k, i)  - ( xls - xlv ) * qc_i(k, i)
       end do

       ! --------------------------------------------------------------- !
       ! Prevent the onset-of negative condensate at the next time step  !
       ! Potentially, this block can be moved just in front of the above !
       ! block.                                                          !
       ! --------------------------------------------------------------- !

       ! Modification : I should check whether this 'positive_moisture_single' routine is
       !                consistent with the one used in UW PBL and cloud macrophysics schemes.
       ! Modification : Below may overestimate resulting 'ql, qi' if we use the new 'qc_l(:, i)', 'qc_i(:, i)'
       !                in combination with the original computation of qlten(:, i), qiten(:, i). However,
       !                if we use new 'qlten(:, i),qiten(:, i)', there is no problem.

        qv0_star(:mkx, i) = qv0(:mkx, i) + qvten(:mkx, i) * dt
        ql0_star(:mkx, i) = ql0(:mkx, i) + qlten(:mkx, i) * dt
        qi0_star(:mkx, i) = qi0(:mkx, i) + qiten(:mkx, i) * dt
        s0_star(:mkx, i)  =  s0(:mkx, i) +  sten(:mkx, i) * dt
        call positive_moisture_single( xlv, xls, mkx, dt, qmin(1), qmin(ixcldliq), qmin(ixcldice), &
             dp0(:, i), qv0_star(:, i), ql0_star(:, i), qi0_star(:, i), s0_star(:, i), qvten(:, i), qlten(:, i), qiten(:, i), sten(:, i) )
        qtten(:mkx, i)    = qvten(:mkx, i) + qlten(:mkx, i) + qiten(:mkx, i)
        slten(:mkx, i)    = sten(:mkx, i)  - xlv * qlten(:mkx, i) - xls * qiten(:mkx, i)

       ! --------------------- !
       ! Tendencies of tracers !
       ! --------------------- !

       do m = 4, ncnst

       if( m .ne. ixnumliq .and. m .ne. ixnumice ) then

          trmin(i) = qmin(m)
          trflx_d(0:mkx, i) = 0._r8
          trflx_u(0:mkx, i) = 0._r8
          do k = 1, mkx-1
             if( cnst_get_type_byind(m) .eq. 'wet' ) then
                 pdelx(i) = dp0(k, i)
             else
                 pdelx(i) = dpdry0(k, i)
             endif
             km1 = k - 1
             dum(i) = ( tr0(k,m, i) - trmin(i) ) *  pdelx(i) / g / dt + trflx(km1,m, i) - trflx(k,m, i) + trflx_d(km1, i)
             trflx_d(k, i) = min( 0._r8, dum(i) )
          enddo
          do k = mkx, 2, -1
             if( cnst_get_type_byind(m) .eq. 'wet' ) then
                 pdelx(i) = dp0(k, i)
             else
                 pdelx(i) = dpdry0(k, i)
             endif
             km1 = k - 1
             dum(i) = ( tr0(k,m, i) - trmin(i) ) * pdelx(i) / g / dt + trflx(km1,m, i) - trflx(k,m, i) + &
                                                           trflx_d(km1, i) - trflx_d(k, i) - trflx_u(k, i)
             trflx_u(km1, i) = max( 0._r8, -dum(i) )
          enddo
          do k = 1, mkx
             if( cnst_get_type_byind(m) .eq. 'wet' ) then
                 pdelx(i) = dp0(k, i)
             else
                 pdelx(i) = dpdry0(k, i)
             endif
             km1 = k - 1
           ! Check : I should re-check whether '_u', '_d' are correctly ordered in
           !         the below tendency computation.
             trten(k,m, i) = ( trflx(km1,m, i) - trflx(k,m, i) + &
                            trflx_d(km1, i) - trflx_d(k, i) + &
                            trflx_u(km1, i) - trflx_u(k, i) ) * g / pdelx(i)
          enddo

       endif

       enddo

       ! ---------------------------------------------------------------- !
       ! Cumpute default diagnostic outputs                               !
       ! Note that since 'qtu(krel(i)-1:kpen(i)-1, i)' & 'thlu(krel(i)-1:kpen(i)-1, i)' has !
       ! been adjusted after detraining cloud condensate into environment !
       ! during cumulus updraft motion,  below calculations will  exactly !
       ! reproduce in-cloud properties as shown in the output analysis.   !
       ! ---------------------------------------------------------------- !

       call conden(prel(i),thlu(krel(i)-1, i),qtu(krel(i)-1, i),thj(i),qvj(i),qlj(i),qij(i),qse(i),id_check)
       if( id_check .eq. 1 ) then
           exit_conden(i) = 1._r8
           id_exit(i) = .true.
           go to 333
       end if
       qcubelow(i) = qlj(i) + qij(i)
       qlubelow(i) = qlj(i)
       qiubelow(i) = qij(i)
       rcwp(i)     = 0._r8
       rlwp(i)     = 0._r8
       riwp(i)     = 0._r8

       ! --------------------------------------------------------------------- !
       ! In the below calculations, I explicitly considered cloud base ( LCL ) !
       ! and cloud top height ( ps0(kpen(i)-1, i) + ppen(i) )                           !
       ! ----------------------------------------------------------------------!
       do k = krel(i), kpen(i) ! This is a layer index
          ! ------------------------------------------------------------------ !
          ! Calculate cumulus condensate at the upper interface of each layer. !
          ! Note 'ppen(i) < 0' and at 'k=kpen(i)' layer, I used 'thlu_top(i)'&'qtu_top(i)' !
          ! which explicitly considered zero or non-zero 'fer(kpen(i), i)'.          !
          ! ------------------------------------------------------------------ !
          if( k .eq. kpen(i) ) then
              call conden(ps0(k-1, i)+ppen(i),thlu_top(i),qtu_top(i),thj(i),qvj(i),qlj(i),qij(i),qse(i),id_check)
          else
              call conden(ps0(k, i),thlu(k, i),qtu(k, i),thj(i),qvj(i),qlj(i),qij(i),qse(i),id_check)
          endif
          if( id_check .eq. 1 ) then
              exit_conden(i) = 1._r8
              id_exit(i) = .true.
              go to 333
          end if
          ! ---------------------------------------------------------------- !
          ! Calculate in-cloud mean LWC ( qlu(k, i) ), IWC ( qiu(k, i) ),  & layer !
          ! mean cumulus fraction ( cufrc(k, i) ),  vertically-integrated layer !
          ! mean LWP and IWP. Expel some of in-cloud condensate at the upper !
          ! interface if it is largr than criqc. Note cumulus cloud fraction !
          ! is assumed to be twice of core updraft fractional area. Thus LWP !
          ! and IWP will be twice of actual value coming from our scheme.    !
          ! ---------------------------------------------------------------- !
          qcu(k, i)   = 0.5_r8 * ( qcubelow(i) + qlj(i) + qij(i) )
          qlu(k, i)   = 0.5_r8 * ( qlubelow(i) + qlj(i) )
          qiu(k, i)   = 0.5_r8 * ( qiubelow(i) + qij(i) )
          cufrc(k, i) = ( ufrc(k-1, i) + ufrc(k, i) )
          if( k .eq. krel(i) ) then
              cufrc(k, i) = ( ufrclcl(i) + ufrc(k, i) )*( prel(i) - ps0(k, i) )/( ps0(k-1, i) - ps0(k, i) )
          else if( k .eq. kpen(i) ) then
              cufrc(k, i) = ( ufrc(k-1, i) + 0._r8 )*( -ppen(i) )        /( ps0(k-1, i) - ps0(k, i) )
              if( (qlj(i) + qij(i)) .gt. criqc ) then
                   qcu(k, i) = 0.5_r8 * ( qcubelow(i) + criqc )
                   qlu(k, i) = 0.5_r8 * ( qlubelow(i) + criqc * qlj(i) / ( qlj(i) + qij(i) ) )
                   qiu(k, i) = 0.5_r8 * ( qiubelow(i) + criqc * qij(i) / ( qlj(i) + qij(i) ) )
              endif
          endif
          rcwp(i) = rcwp(i) + ( qlu(k, i) + qiu(k, i) ) * ( ps0(k-1, i) - ps0(k, i) ) / g * cufrc(k, i)
          rlwp(i) = rlwp(i) +   qlu(k, i)            * ( ps0(k-1, i) - ps0(k, i) ) / g * cufrc(k, i)
          riwp(i) = riwp(i) +   qiu(k, i)            * ( ps0(k-1, i) - ps0(k, i) ) / g * cufrc(k, i)
          qcubelow(i) = qlj(i) + qij(i)
          qlubelow(i) = qlj(i)
          qiubelow(i) = qij(i)
       end do
       ! ------------------------------------ !
       ! Cloud top and base interface indices !
       ! ------------------------------------ !
       cnt(i) = real( kpen(i), r8 )
       cnb(i) = real( krel(i) - 1, r8 )

       ! ------------------------------------------------------------------------- !
       ! End of formal calculation. Below blocks are for implicit CIN calculations !
       ! with re-initialization and save variables at iter_cin = 1._r8             !
       ! ------------------------------------------------------------------------- !

       ! --------------------------------------------------------------- !
       ! Adjust the original input profiles for implicit CIN calculation !
       ! --------------------------------------------------------------- !

       if( iter .ne. iter_cin ) then

          ! ------------------------------------------------------------------- !
          ! Save the output from "iter_cin = 1"                                 !
          ! These output will be writed-out if "iter_cin = 1" was not performed !
          ! for some reasons.                                                   !
          ! ------------------------------------------------------------------- !

          qv0_s(:mkx, i)           = qv0(:mkx, i) + qvten(:mkx, i) * dt
          ql0_s(:mkx, i)           = ql0(:mkx, i) + qlten(:mkx, i) * dt
          qi0_s(:mkx, i)           = qi0(:mkx, i) + qiten(:mkx, i) * dt
          s0_s(:mkx, i)            = s0(:mkx, i)  +  sten(:mkx, i) * dt
          u0_s(:mkx, i)            = u0(:mkx, i)  +  uten(:mkx, i) * dt
          v0_s(:mkx, i)            = v0(:mkx, i)  +  vten(:mkx, i) * dt
          qt0_s(:mkx, i)           = qv0_s(:mkx, i) + ql0_s(:mkx, i) + qi0_s(:mkx, i)
          t0_s(:mkx, i)            = t0(:mkx, i)  +  sten(:mkx, i) * dt / cp
          do m = 1, ncnst
             tr0_s(:mkx,m, i)      = tr0(:mkx,m, i) + trten(:mkx,m, i) * dt
          enddo

          umf_s(0:mkx, i)          = umf(0:mkx, i)
          qvten_s(:mkx, i)         = qvten(:mkx, i)
          qlten_s(:mkx, i)         = qlten(:mkx, i)
          qiten_s(:mkx, i)         = qiten(:mkx, i)
          sten_s(:mkx, i)          = sten(:mkx, i)
          uten_s(:mkx, i)          = uten(:mkx, i)
          vten_s(:mkx, i)          = vten(:mkx, i)
          qrten_s(:mkx, i)         = qrten(:mkx, i)
          qsten_s(:mkx, i)         = qsten(:mkx, i)
          precip_s(i)              = precip(i)
          snow_s(i)                = snow(i)
          evapc_s(:mkx, i)         = evapc(:mkx, i)
          cush_s(i)                = cush(i)
          cufrc_s(:mkx, i)         = cufrc(:mkx, i)
          slflx_s(0:mkx, i)        = slflx(0:mkx, i)
          qtflx_s(0:mkx, i)        = qtflx(0:mkx, i)
          qcu_s(:mkx, i)           = qcu(:mkx, i)
          qlu_s(:mkx, i)           = qlu(:mkx, i)
          qiu_s(:mkx, i)           = qiu(:mkx, i)
          fer_s(:mkx, i)           = fer(:mkx, i)
          fdr_s(:mkx, i)           = fdr(:mkx, i)
          cin_s(i)                 = cin(i)
          cinlcl_s(i)              = cinlcl(i)
          cbmf_s(i)                = cbmf(i)
          rliq_s(i)                = rliq(i)
          qc_s(:mkx, i)            = qc(:mkx, i)
          cnt_s(i)                 = cnt(i)
          cnb_s(i)                 = cnb(i)
          qtten_s(:mkx, i)         = qtten(:mkx, i)
          slten_s(:mkx, i)         = slten(:mkx, i)
          ufrc_s(0:mkx, i)         = ufrc(0:mkx, i)

          uflx_s(0:mkx, i)         = uflx(0:mkx, i)
          vflx_s(0:mkx, i)         = vflx(0:mkx, i)

          ufrcinvbase_s(i)         = ufrcinvbase(i)
          ufrclcl_s(i)             = ufrclcl(i)
          winvbase_s(i)            = winvbase(i)
          wlcl_s(i)                = wlcl(i)
          plcl_s(i)                = plcl(i)
          pinv_s(i)                = ps0(kinv(i)-1, i)
          plfc_s(i)                = plfc(i)
          pbup_s(i)                = ps0(kbup(i), i)
          ppen_s(i)                = ps0(kpen(i)-1, i) + ppen(i)
          qtsrc_s(i)               = qtsrc(i)
          thlsrc_s(i)              = thlsrc(i)
          thvlsrc_s(i)             = thvlsrc(i)
          emfkbup_s(i)             = emf(kbup(i), i)
          cbmflimit_s(i)           = cbmflimit(i)
          tkeavg_s(i)              = tkeavg(i)
          zinv_s(i)                = zs0(kinv(i)-1, i)
          rcwp_s(i)                = rcwp(i)
          rlwp_s(i)                = rlwp(i)
          riwp_s(i)                = riwp(i)

          wu_s(0:mkx, i)           = wu(0:mkx, i)
          qtu_s(0:mkx, i)          = qtu(0:mkx, i)
          thlu_s(0:mkx, i)         = thlu(0:mkx, i)
          thvu_s(0:mkx, i)         = thvu(0:mkx, i)
          uu_s(0:mkx, i)           = uu(0:mkx, i)
          vu_s(0:mkx, i)           = vu(0:mkx, i)
          qtu_emf_s(0:mkx, i)      = qtu_emf(0:mkx, i)
          thlu_emf_s(0:mkx, i)     = thlu_emf(0:mkx, i)
          uu_emf_s(0:mkx, i)       = uu_emf(0:mkx, i)
          vu_emf_s(0:mkx, i)       = vu_emf(0:mkx, i)
          uemf_s(0:mkx, i)         = uemf(0:mkx, i)

          dwten_s(:mkx, i)         = dwten(:mkx, i)
          diten_s(:mkx, i)         = diten(:mkx, i)
          flxrain_s(0:mkx, i)      = flxrain(0:mkx, i)
          flxsnow_s(0:mkx, i)      = flxsnow(0:mkx, i)
          ntraprd_s(:mkx, i)       = ntraprd(:mkx, i)
          ntsnprd_s(:mkx, i)       = ntsnprd(:mkx, i)

          excessu_arr_s(:mkx, i)   = excessu_arr(:mkx, i)
          excess0_arr_s(:mkx, i)   = excess0_arr(:mkx, i)
          xc_arr_s(:mkx, i)        = xc_arr(:mkx, i)
          aquad_arr_s(:mkx, i)     = aquad_arr(:mkx, i)
          bquad_arr_s(:mkx, i)     = bquad_arr(:mkx, i)
          cquad_arr_s(:mkx, i)     = cquad_arr(:mkx, i)
          bogbot_arr_s(:mkx, i)    = bogbot_arr(:mkx, i)
          bogtop_arr_s(:mkx, i)    = bogtop_arr(:mkx, i)

          do m = 1, ncnst
             trten_s(:mkx,m, i)    = trten(:mkx,m, i)
             trflx_s(0:mkx,m, i)   = trflx(0:mkx,m, i)
             tru_s(0:mkx,m, i)     = tru(0:mkx,m, i)
             tru_emf_s(0:mkx,m, i) = tru_emf(0:mkx,m, i)
          enddo

          ! ----------------------------------------------------------------------------- !
          ! Recalculate environmental variables for new cin(i) calculation at "iter_cin = 2" !
          ! using the updated state variables. Perform only for variables necessary  for  !
          ! the new cin(i) calculation.                                                      !
          ! ----------------------------------------------------------------------------- !

          qv0(:mkx, i)   = qv0_s(:mkx, i)
          ql0(:mkx, i)   = ql0_s(:mkx, i)
          qi0(:mkx, i)   = qi0_s(:mkx, i)
          s0(:mkx, i)    = s0_s(:mkx, i)
          t0(:mkx, i)    = t0_s(:mkx, i)

          qt0(:mkx, i)   = (qv0(:mkx, i) + ql0(:mkx, i) + qi0(:mkx, i))
          thl0(:mkx, i)  = (t0(:mkx, i) - xlv*ql0(:mkx, i)/cp - xls*qi0(:mkx, i)/cp)/exn0(:mkx, i)
          thvl0(:mkx, i) = (1._r8 + zvir*qt0(:mkx, i))*thl0(:mkx, i)

          ssthl0(:, i)      = slope(mkx,thl0(:, i),p0(:, i)) ! Dimension of ssthl0(:mkx, i) is implicit
          ssqt0(:, i)       = slope(mkx,qt0(:, i) ,p0(:, i))
          ssu0(:, i)        = slope(mkx,u0(:, i)  ,p0(:, i))
          ssv0(:, i)        = slope(mkx,v0(:, i)  ,p0(:, i))
          do m = 1, ncnst
             sstr0(:mkx,m, i) = slope(mkx,tr0(:mkx,m, i),p0(:, i))
          enddo

          do k = 1, mkx

             thl0bot(i) = thl0(k, i) + ssthl0(k, i) * ( ps0(k-1, i) - p0(k, i) )
             qt0bot(i)  = qt0(k, i)  + ssqt0(k, i)  * ( ps0(k-1, i) - p0(k, i) )
             call conden(ps0(k-1, i),thl0bot(i),qt0bot(i),thj(i),qvj(i),qlj(i),qij(i),qse(i),id_check)
             if( id_check .eq. 1 ) then
                 exit_conden(i) = 1._r8
                 id_exit(i) = .true.
                 go to 333
             end if
             thv0bot(k, i)  = thj(i) * ( 1._r8 + zvir*qvj(i) - qlj(i) - qij(i) )
             thvl0bot(k, i) = thl0bot(i) * ( 1._r8 + zvir*qt0bot(i) )

             thl0top(i) = thl0(k, i) + ssthl0(k, i) * ( ps0(k, i) - p0(k, i) )
             qt0top(i)  =  qt0(k, i) + ssqt0(k, i)  * ( ps0(k, i) - p0(k, i) )
             call conden(ps0(k, i),thl0top(i),qt0top(i),thj(i),qvj(i),qlj(i),qij(i),qse(i),id_check)
             if( id_check .eq. 1 ) then
                 exit_conden(i) = 1._r8
                 id_exit(i) = .true.
                 go to 333
             end if
             thv0top(k, i)  = thj(i) * ( 1._r8 + zvir*qvj(i) - qlj(i) - qij(i) )
             thvl0top(k, i) = thl0top(i) * ( 1._r8 + zvir*qt0top(i) )

          end do

       endif               ! End of 'if(iter .ne. iter_cin)' if sentence.

     ! end do                ! End of implicit CIN loop (cin_iter)

     if (iter == iter_cin) then

     ! ----------------------- !
     ! Update Output Variables !
     ! ----------------------- !

     umf_out(i,0:mkx)             = umf(0:mkx, i)
     slflx_out(i,0:mkx)           = slflx(0:mkx, i)
     qtflx_out(i,0:mkx)           = qtflx(0:mkx, i)
!the indices are not reversed, these variables go into compute_mcshallow_inv, this is why they are called "flxprc1" and "flxsnow1".
     flxprc1_out(i,0:mkx)         = flxrain(0:mkx, i) + flxsnow(0:mkx, i)
     flxsnow1_out(i,0:mkx)        = flxsnow(0:mkx, i)
     qvten_out(i,:mkx)            = qvten(:mkx, i)
     qlten_out(i,:mkx)            = qlten(:mkx, i)
     qiten_out(i,:mkx)            = qiten(:mkx, i)
     sten_out(i,:mkx)             = sten(:mkx, i)
     uten_out(i,:mkx)             = uten(:mkx, i)
     vten_out(i,:mkx)             = vten(:mkx, i)
     qrten_out(i,:mkx)            = qrten(:mkx, i)
     qsten_out(i,:mkx)            = qsten(:mkx, i)
     precip_out(i)                = precip(i)
     snow_out(i)                  = snow(i)
     evapc_out(i,:mkx)            = evapc(:mkx, i)
     cufrc_out(i,:mkx)            = cufrc(:mkx, i)
     qcu_out(i,:mkx)              = qcu(:mkx, i)
     qlu_out(i,:mkx)              = qlu(:mkx, i)
     qiu_out(i,:mkx)              = qiu(:mkx, i)
     cush_inout(i)                = cush(i)
     cbmf_out(i)                  = cbmf(i)
     rliq_out(i)                  = rliq(i)
     qc_out(i,:mkx)               = qc(:mkx, i)
     cnt_out(i)                   = cnt(i)
     cnb_out(i)                   = cnb(i)

     do m = 1, ncnst
        trten_out(i,:mkx,m)       = trten(:mkx,m, i)
     enddo

     ! ------------------------------------------------- !
     ! Below are specific diagnostic output for detailed !
     ! analysis of cumulus scheme                        !
     ! ------------------------------------------------- !

     fer_out(i,mkx:1:-1)          = fer(:mkx, i)
     fdr_out(i,mkx:1:-1)          = fdr(:mkx, i)
     cinh_out(i)                  = cin(i)
     cinlclh_out(i)               = cinlcl(i)
     qtten_out(i,mkx:1:-1)        = qtten(:mkx, i)
     slten_out(i,mkx:1:-1)        = slten(:mkx, i)
     ufrc_out(i,mkx:0:-1)         = ufrc(0:mkx, i)
     uflx_out(i,mkx:0:-1)         = uflx(0:mkx, i)
     vflx_out(i,mkx:0:-1)         = vflx(0:mkx, i)

     ufrcinvbase_out(i)           = ufrcinvbase(i)
     ufrclcl_out(i)               = ufrclcl(i)
     winvbase_out(i)              = winvbase(i)
     wlcl_out(i)                  = wlcl(i)
     plcl_out(i)                  = plcl(i)
     pinv_out(i)                  = ps0(kinv(i)-1, i)
     plfc_out(i)                  = plfc(i)
     pbup_out(i)                  = ps0(kbup(i), i)
     ppen_out(i)                  = ps0(kpen(i)-1, i) + ppen(i)
     qtsrc_out(i)                 = qtsrc(i)
     thlsrc_out(i)                = thlsrc(i)
     thvlsrc_out(i)               = thvlsrc(i)
     emfkbup_out(i)               = emf(kbup(i), i)
     cbmflimit_out(i)             = cbmflimit(i)
     tkeavg_out(i)                = tkeavg(i)
     zinv_out(i)                  = zs0(kinv(i)-1, i)
     rcwp_out(i)                  = rcwp(i)
     rlwp_out(i)                  = rlwp(i)
     riwp_out(i)                  = riwp(i)

     wu_out(i,mkx:0:-1)           = wu(0:mkx, i)
     qtu_out(i,mkx:0:-1)          = qtu(0:mkx, i)
     thlu_out(i,mkx:0:-1)         = thlu(0:mkx, i)
     thvu_out(i,mkx:0:-1)         = thvu(0:mkx, i)
     uu_out(i,mkx:0:-1)           = uu(0:mkx, i)
     vu_out(i,mkx:0:-1)           = vu(0:mkx, i)
     qtu_emf_out(i,mkx:0:-1)      = qtu_emf(0:mkx, i)
     thlu_emf_out(i,mkx:0:-1)     = thlu_emf(0:mkx, i)
     uu_emf_out(i,mkx:0:-1)       = uu_emf(0:mkx, i)
     vu_emf_out(i,mkx:0:-1)       = vu_emf(0:mkx, i)
     uemf_out(i,mkx:0:-1)         = uemf(0:mkx, i)

     dwten_out(i,mkx:1:-1)        = dwten(:mkx, i)
     diten_out(i,mkx:1:-1)        = diten(:mkx, i)
     flxrain_out(i,mkx:0:-1)      = flxrain(0:mkx, i)
     flxsnow_out(i,mkx:0:-1)      = flxsnow(0:mkx, i)
     ntraprd_out(i,mkx:1:-1)      = ntraprd(:mkx, i)
     ntsnprd_out(i,mkx:1:-1)      = ntsnprd(:mkx, i)

     excessu_arr_out(i,mkx:1:-1)  = excessu_arr(:mkx, i)
     excess0_arr_out(i,mkx:1:-1)  = excess0_arr(:mkx, i)
     xc_arr_out(i,mkx:1:-1)       = xc_arr(:mkx, i)
     aquad_arr_out(i,mkx:1:-1)    = aquad_arr(:mkx, i)
     bquad_arr_out(i,mkx:1:-1)    = bquad_arr(:mkx, i)
     cquad_arr_out(i,mkx:1:-1)    = cquad_arr(:mkx, i)
     bogbot_arr_out(i,mkx:1:-1)   = bogbot_arr(:mkx, i)
     bogtop_arr_out(i,mkx:1:-1)   = bogtop_arr(:mkx, i)

     do m = 1, ncnst
        trflx_out(i,mkx:0:-1,m)   = trflx(0:mkx,m, i)
        tru_out(i,mkx:0:-1,m)     = tru(0:mkx,m, i)
        tru_emf_out(i,mkx:0:-1,m) = tru_emf(0:mkx,m, i)
     enddo

     endif ! (iter == iter_cin)

 333 if(id_exit(i)) then ! Exit without cumulus convection

     exit_UWCu(i) = 1._r8

     ! --------------------------------------------------------------------- !
     ! Initialize output variables when cumulus convection was not performed.!
     ! --------------------------------------------------------------------- !

     umf_out(i,0:mkx)             = 0._r8
     slflx_out(i,0:mkx)           = 0._r8
     qtflx_out(i,0:mkx)           = 0._r8
     qvten_out(i,:mkx)            = 0._r8
     qlten_out(i,:mkx)            = 0._r8
     qiten_out(i,:mkx)            = 0._r8
     sten_out(i,:mkx)             = 0._r8
     uten_out(i,:mkx)             = 0._r8
     vten_out(i,:mkx)             = 0._r8
     qrten_out(i,:mkx)            = 0._r8
     qsten_out(i,:mkx)            = 0._r8
     precip_out(i)                = 0._r8
     snow_out(i)                  = 0._r8
     evapc_out(i,:mkx)            = 0._r8
     cufrc_out(i,:mkx)            = 0._r8
     qcu_out(i,:mkx)              = 0._r8
     qlu_out(i,:mkx)              = 0._r8
     qiu_out(i,:mkx)              = 0._r8
     cush_inout(i)                = -1._r8
     cbmf_out(i)                  = 0._r8
     rliq_out(i)                  = 0._r8
     qc_out(i,:mkx)               = 0._r8
     cnt_out(i)                   = 1._r8
     cnb_out(i)                   = real(mkx, r8)

     fer_out(i,mkx:1:-1)          = 0._r8
     fdr_out(i,mkx:1:-1)          = 0._r8
     cinh_out(i)                  = -1._r8
     cinlclh_out(i)               = -1._r8
     qtten_out(i,mkx:1:-1)        = 0._r8
     slten_out(i,mkx:1:-1)        = 0._r8
     ufrc_out(i,mkx:0:-1)         = 0._r8
     uflx_out(i,mkx:0:-1)         = 0._r8
     vflx_out(i,mkx:0:-1)         = 0._r8

     ufrcinvbase_out(i)           = 0._r8
     ufrclcl_out(i)               = 0._r8
     winvbase_out(i)              = 0._r8
     wlcl_out(i)                  = 0._r8
     plcl_out(i)                  = 0._r8
     pinv_out(i)                  = 0._r8
     plfc_out(i)                  = 0._r8
     pbup_out(i)                  = 0._r8
     ppen_out(i)                  = 0._r8
     qtsrc_out(i)                 = 0._r8
     thlsrc_out(i)                = 0._r8
     thvlsrc_out(i)               = 0._r8
     emfkbup_out(i)               = 0._r8
     cbmflimit_out(i)             = 0._r8
     tkeavg_out(i)                = 0._r8
     zinv_out(i)                  = 0._r8
     rcwp_out(i)                  = 0._r8
     rlwp_out(i)                  = 0._r8
     riwp_out(i)                  = 0._r8

     wu_out(i,mkx:0:-1)           = 0._r8
     qtu_out(i,mkx:0:-1)          = 0._r8
     thlu_out(i,mkx:0:-1)         = 0._r8
     thvu_out(i,mkx:0:-1)         = 0._r8
     uu_out(i,mkx:0:-1)           = 0._r8
     vu_out(i,mkx:0:-1)           = 0._r8
     qtu_emf_out(i,mkx:0:-1)      = 0._r8
     thlu_emf_out(i,mkx:0:-1)     = 0._r8
     uu_emf_out(i,mkx:0:-1)       = 0._r8
     vu_emf_out(i,mkx:0:-1)       = 0._r8
     uemf_out(i,mkx:0:-1)         = 0._r8

     dwten_out(i,mkx:1:-1)        = 0._r8
     diten_out(i,mkx:1:-1)        = 0._r8
     flxrain_out(i,mkx:0:-1)      = 0._r8
     flxsnow_out(i,mkx:0:-1)      = 0._r8
     ntraprd_out(i,mkx:1:-1)      = 0._r8
     ntsnprd_out(i,mkx:1:-1)      = 0._r8

     excessu_arr_out(i,mkx:1:-1)  = 0._r8
     excess0_arr_out(i,mkx:1:-1)  = 0._r8
     xc_arr_out(i,mkx:1:-1)       = 0._r8
     aquad_arr_out(i,mkx:1:-1)    = 0._r8
     bquad_arr_out(i,mkx:1:-1)    = 0._r8
     cquad_arr_out(i,mkx:1:-1)    = 0._r8
     bogbot_arr_out(i,mkx:1:-1)   = 0._r8
     bogtop_arr_out(i,mkx:1:-1)   = 0._r8

     do m = 1, ncnst
        trten_out(i,:mkx,m)       = 0._r8
        trflx_out(i,mkx:0:-1,m)   = 0._r8
        tru_out(i,mkx:0:-1,m)     = 0._r8
        tru_emf_out(i,mkx:0:-1,m) = 0._r8
     enddo

     end if

     end do                  ! end of big i loop for each column.

     end do                ! End of implicit CIN loop (iter_cin)

     ! ---------------------------------------- !
     ! Writing main diagnostic output variables !
     ! ---------------------------------------- !

     call outfld( 'qtflx_Cu'        , qtflx_out(:,mkx:0:-1),    mix,    lchnk )
     call outfld( 'slflx_Cu'        , slflx_out(:,mkx:0:-1),    mix,    lchnk )
     call outfld( 'uflx_Cu'         , uflx_out,                 mix,    lchnk )
     call outfld( 'vflx_Cu'         , vflx_out,                 mix,    lchnk )

     call outfld( 'qtten_Cu'        , qtten_out,                mix,    lchnk )
     call outfld( 'slten_Cu'        , slten_out,                mix,    lchnk )
     call outfld( 'uten_Cu'         , uten_out(:,mkx:1:-1),     mix,    lchnk )
     call outfld( 'vten_Cu'         , vten_out(:,mkx:1:-1),     mix,    lchnk )
     call outfld( 'qvten_Cu'        , qvten_out(:,mkx:1:-1),    mix,    lchnk )
     call outfld( 'qlten_Cu'        , qlten_out(:,mkx:1:-1),    mix,    lchnk )
     call outfld( 'qiten_Cu'        , qiten_out(:,mkx:1:-1),    mix,    lchnk )

     call outfld( 'cbmf_Cu'         , cbmf_out,                 mix,    lchnk )
     call outfld( 'ufrcinvbase_Cu'  , ufrcinvbase_out,          mix,    lchnk )
     call outfld( 'ufrclcl_Cu'      , ufrclcl_out,              mix,    lchnk )
     call outfld( 'winvbase_Cu'     , winvbase_out,             mix,    lchnk )
     call outfld( 'wlcl_Cu'         , wlcl_out,                 mix,    lchnk )
     call outfld( 'plcl_Cu'         , plcl_out,                 mix,    lchnk )
     call outfld( 'pinv_Cu'         , pinv_out,                 mix,    lchnk )
     call outfld( 'plfc_Cu'         , plfc_out,                 mix,    lchnk )
     call outfld( 'pbup_Cu'         , pbup_out,                 mix,    lchnk )
     call outfld( 'ppen_Cu'         , ppen_out,                 mix,    lchnk )
     call outfld( 'qtsrc_Cu'        , qtsrc_out,                mix,    lchnk )
     call outfld( 'thlsrc_Cu'       , thlsrc_out,               mix,    lchnk )
     call outfld( 'thvlsrc_Cu'      , thvlsrc_out,              mix,    lchnk )
     call outfld( 'emfkbup_Cu'      , emfkbup_out,              mix,    lchnk )
     call outfld( 'cin_Cu'          , cinh_out,                 mix,    lchnk )
     call outfld( 'cinlcl_Cu'       , cinlclh_out,              mix,    lchnk )
     call outfld( 'cbmflimit_Cu'    , cbmflimit_out,            mix,    lchnk )
     call outfld( 'tkeavg_Cu'       , tkeavg_out,               mix,    lchnk )
     call outfld( 'zinv_Cu'         , zinv_out,                 mix,    lchnk )
     call outfld( 'rcwp_Cu'         , rcwp_out,                 mix,    lchnk )
     call outfld( 'rlwp_Cu'         , rlwp_out,                 mix,    lchnk )
     call outfld( 'riwp_Cu'         , riwp_out,                 mix,    lchnk )
     call outfld( 'tophgt_Cu'       , cush_inout,               mix,    lchnk )

     call outfld( 'wu_Cu'           , wu_out,                   mix,    lchnk )
     call outfld( 'ufrc_Cu'         , ufrc_out,                 mix,    lchnk )
     call outfld( 'qtu_Cu'          , qtu_out,                  mix,    lchnk )
     call outfld( 'thlu_Cu'         , thlu_out,                 mix,    lchnk )
     call outfld( 'thvu_Cu'         , thvu_out,                 mix,    lchnk )
     call outfld( 'uu_Cu'           , uu_out,                   mix,    lchnk )
     call outfld( 'vu_Cu'           , vu_out,                   mix,    lchnk )
     call outfld( 'qtu_emf_Cu'      , qtu_emf_out,              mix,    lchnk )
     call outfld( 'thlu_emf_Cu'     , thlu_emf_out,             mix,    lchnk )
     call outfld( 'uu_emf_Cu'       , uu_emf_out,               mix,    lchnk )
     call outfld( 'vu_emf_Cu'       , vu_emf_out,               mix,    lchnk )
     call outfld( 'umf_Cu'          , umf_out(:,mkx:0:-1),      mix,    lchnk )
     call outfld( 'uemf_Cu'         , uemf_out,                 mix,    lchnk )
     call outfld( 'qcu_Cu'          , qcu_out(:,mkx:1:-1),      mix,    lchnk )
     call outfld( 'qlu_Cu'          , qlu_out(:,mkx:1:-1),      mix,    lchnk )
     call outfld( 'qiu_Cu'          , qiu_out(:,mkx:1:-1),      mix,    lchnk )
     call outfld( 'cufrc_Cu'        , cufrc_out(:,mkx:1:-1),    mix,    lchnk )
     call outfld( 'fer_Cu'          , fer_out,                  mix,    lchnk )
     call outfld( 'fdr_Cu'          , fdr_out,                  mix,    lchnk )

     call outfld( 'dwten_Cu'        , dwten_out,                mix,    lchnk )
     call outfld( 'diten_Cu'        , diten_out,                mix,    lchnk )
     call outfld( 'qrten_Cu'        , qrten_out(:,mkx:1:-1),    mix,    lchnk )
     call outfld( 'qsten_Cu'        , qsten_out(:,mkx:1:-1),    mix,    lchnk )
     call outfld( 'flxrain_Cu'      , flxrain_out,              mix,    lchnk )
     call outfld( 'flxsnow_Cu'      , flxsnow_out,              mix,    lchnk )
     call outfld( 'ntraprd_Cu'      , ntraprd_out,              mix,    lchnk )
     call outfld( 'ntsnprd_Cu'      , ntsnprd_out,              mix,    lchnk )

     call outfld( 'excessu_Cu'      , excessu_arr_out,          mix,    lchnk )
     call outfld( 'excess0_Cu'      , excess0_arr_out,          mix,    lchnk )
     call outfld( 'xc_Cu'           , xc_arr_out,               mix,    lchnk )
     call outfld( 'aquad_Cu'        , aquad_arr_out,            mix,    lchnk )
     call outfld( 'bquad_Cu'        , bquad_arr_out,            mix,    lchnk )
     call outfld( 'cquad_Cu'        , cquad_arr_out,            mix,    lchnk )
     call outfld( 'bogbot_Cu'       , bogbot_arr_out,           mix,    lchnk )
     call outfld( 'bogtop_Cu'       , bogtop_arr_out,           mix,    lchnk )

     call outfld( 'exit_UWCu_Cu'    , exit_UWCu,                mix,    lchnk )
     call outfld( 'exit_conden_Cu'  , exit_conden,              mix,    lchnk )
     call outfld( 'exit_klclmkx_Cu' , exit_klclmkx,             mix,    lchnk )
     call outfld( 'exit_klfcmkx_Cu' , exit_klfcmkx,             mix,    lchnk )
     call outfld( 'exit_ufrc_Cu'    , exit_ufrc,                mix,    lchnk )
     call outfld( 'exit_wtw_Cu'     , exit_wtw,                 mix,    lchnk )
     call outfld( 'exit_drycore_Cu' , exit_drycore,             mix,    lchnk )
     call outfld( 'exit_wu_Cu'      , exit_wu,                  mix,    lchnk )
     call outfld( 'exit_cufilter_Cu', exit_cufilter,            mix,    lchnk )
     call outfld( 'exit_kinv1_Cu'   , exit_kinv1,               mix,    lchnk )
     call outfld( 'exit_rei_Cu'     , exit_rei,                 mix,    lchnk )

     call outfld( 'limit_shcu_Cu'   , limit_shcu,               mix,    lchnk )
     call outfld( 'limit_negcon_Cu' , limit_negcon,             mix,    lchnk )
     call outfld( 'limit_ufrc_Cu'   , limit_ufrc,               mix,    lchnk )
     call outfld( 'limit_ppen_Cu'   , limit_ppen,               mix,    lchnk )
     call outfld( 'limit_emf_Cu'    , limit_emf,                mix,    lchnk )
     call outfld( 'limit_cinlcl_Cu' , limit_cinlcl,             mix,    lchnk )
     call outfld( 'limit_cin_Cu'    , limit_cin,                mix,    lchnk )
     call outfld( 'limit_cbmf_Cu'   , limit_cbmf,               mix,    lchnk )
     call outfld( 'limit_rei_Cu'    , limit_rei,                mix,    lchnk )
     call outfld( 'ind_delcin_Cu'   , ind_delcin,               mix,    lchnk )

    return

  end subroutine compute_uwshcu

!!! end of uwshcu

  ! ------------------------------ !
  !                                !
  ! Beginning of subroutine blocks !
  !                                !
  ! ------------------------------ !

  subroutine getbuoy(pbot,thv0bot,ptop,thv0top,thvubot,thvutop,plfc,cin)
  ! ----------------------------------------------------------- !
  ! Subroutine to calculate integrated CIN [ J/kg = m2/s2 ] and !
  ! 'cinlcl, plfc' if any. Assume 'thv' is linear in each layer !
  ! both for cumulus and environment. Note that this subroutine !
  ! only include positive CIN in calculation - if there are any !
  ! negative CIN, it is assumed to be zero.    This is slightly !
  ! different from 'single_cin' below, where both positive  and !
  ! negative CIN are included.                                  !
  ! ----------------------------------------------------------- !
    real(r8) pbot,thv0bot,ptop,thv0top,thvubot,thvutop,plfc,cin,frc

    if( thvubot .gt. thv0bot .and. thvutop .gt. thv0top ) then
        plfc = pbot
        return
    elseif( thvubot .le. thv0bot .and. thvutop .le. thv0top ) then
        cin  = cin - ( (thvubot/thv0bot - 1._r8) + (thvutop/thv0top - 1._r8)) * (pbot - ptop) /        &
                     ( pbot/(r*thv0bot*exnf(pbot)) + ptop/(r*thv0top*exnf(ptop)) )
    elseif( thvubot .gt. thv0bot .and. thvutop .le. thv0top ) then
        frc  = ( thvutop/thv0top - 1._r8 ) / ( (thvutop/thv0top - 1._r8) - (thvubot/thv0bot - 1._r8) )
        cin  = cin - ( thvutop/thv0top - 1._r8 ) * ( (ptop + frc*(pbot - ptop)) - ptop ) /             &
                     ( pbot/(r*thv0bot*exnf(pbot)) + ptop/(r*thv0top*exnf(ptop)) )
    else
        frc  = ( thvubot/thv0bot - 1._r8 ) / ( (thvubot/thv0bot - 1._r8) - (thvutop/thv0top - 1._r8) )
        plfc = pbot - frc * ( pbot - ptop )
        cin  = cin - ( thvubot/thv0bot - 1._r8)*(pbot - plfc)/                                         &
                     ( pbot/(r*thv0bot*exnf(pbot)) + ptop/(r*thv0top * exnf(ptop)))
    endif

    return
  end subroutine getbuoy

  function single_cin(pbot,thv0bot,ptop,thv0top,thvubot,thvutop)
  ! ------------------------------------------------------- !
  ! Function to calculate a single layer CIN by summing all !
  ! positive and negative CIN.                              !
  ! ------------------------------------------------------- !
    real(r8) :: single_cin
    real(r8)    pbot,thv0bot,ptop,thv0top,thvubot,thvutop

    single_cin = ( (1._r8 - thvubot/thv0bot) + (1._r8 - thvutop/thv0top)) * ( pbot - ptop ) / &
                 ( pbot/(r*thv0bot*exnf(pbot)) + ptop/(r*thv0top*exnf(ptop)) )
    return
  end function single_cin


  subroutine conden(p,thl,qt,th,qv,ql,qi,rvls,id_check)
  ! --------------------------------------------------------------------- !
  ! Calculate thermodynamic properties from a given set of ( p, thl, qt ) !
  ! --------------------------------------------------------------------- !
    implicit none
    real(r8), intent(in)  :: p
    real(r8), intent(in)  :: thl
    real(r8), intent(in)  :: qt
    real(r8), intent(out) :: th
    real(r8), intent(out) :: qv
    real(r8), intent(out) :: ql
    real(r8), intent(out) :: qi
    real(r8), intent(out) :: rvls
    integer , intent(out) :: id_check
    real(r8)              :: tc,temps,t
    real(r8)              :: leff, nu, qc
    integer               :: iteration
    real(r8)              :: es              ! Saturation vapor pressure
    real(r8)              :: qs              ! Saturation spec. humidity


    tc   = thl*exnf(p)
  ! Modification : In order to be compatible with the dlf treatment in stratiform.F90,
  !                we may use ( 268.15, 238.15 ) with 30K ramping instead of 20 K,
  !                in computing ice fraction below.
  !                Note that 'cldfrc_fice' uses ( 243.15, 263.15 ) with 20K ramping for stratus.
    nu   = max(min((268._r8 - tc)/20._r8,1.0_r8),0.0_r8)  ! Fraction of ice in the condensate.
    leff = (1._r8 - nu)*xlv + nu*xls                      ! This is an estimate that hopefully speeds convergence

    ! --------------------------------------------------------------------------- !
    ! Below "temps" and "rvls" are just initial guesses for iteration loop below. !
    ! Note that the output "temps" from the below iteration loop is "temperature" !
    ! NOT "liquid temperature".                                                   !
    ! --------------------------------------------------------------------------- !

    temps  = tc
    call t_startf('qsat_o')
#ifdef QSAT_WV
    call qsat_wv(temps, p, es, qs)
#else
    call qsat_o(temps, p, es, qs)
#endif
    call t_stopf('qsat_o')
    rvls   = qs

    if( qs .ge. qt ) then
        id_check = 0
        qv = qt
        qc = 0._r8
        ql = 0._r8
        qi = 0._r8
        th = tc/exnf(p)
    else
        do iteration = 1, 10
           temps  = temps + ( (tc-temps)*cp/leff + qt - rvls )/( cp/leff + ep2*leff*rvls/r/temps/temps )
           call t_startf('qsat_o')
#ifdef QSAT_WV
           call qsat_wv(temps, p, es, qs)
#else
           call qsat_o(temps, p, es, qs)
#endif
           call t_stopf('qsat_o')
           rvls   = qs
        end do
        qc = max(qt - qs,0._r8)
        qv = qt - qc
        ql = qc*(1._r8 - nu)
        qi = nu*qc
        th = temps/exnf(p)
        if( abs((temps-(leff/cp)*qc)-tc) .ge. 1._r8 ) then
            id_check = 1
        else
            id_check = 0
        end if
    end if

    return
  end subroutine conden

  subroutine roots(a,b,c,r1,r2,status)
  ! --------------------------------------------------------- !
  ! Subroutine to solve the second order polynomial equation. !
  ! I should check this subroutine later.                     !
  ! --------------------------------------------------------- !
    real(r8), intent(in)  :: a
    real(r8), intent(in)  :: b
    real(r8), intent(in)  :: c
    real(r8), intent(out) :: r1
    real(r8), intent(out) :: r2
    integer , intent(out) :: status
    real(r8)              :: q

    status = 0

    if( a .eq. 0._r8 ) then                            ! Form b*x + c = 0
        if( b .eq. 0._r8 ) then                        ! Failure: c = 0
            status = 1
        else                                           ! b*x + c = 0
            r1 = -c/b
        endif
        r2 = r1
    else
        if( b .eq. 0._r8 ) then                        ! Form a*x**2 + c = 0
            if( a*c .gt. 0._r8 ) then                  ! Failure: x**2 = -c/a < 0
                status = 2
            else                                       ! x**2 = -c/a
                r1 = sqrt(-c/a)
            endif
            r2 = -r1
       else                                            ! Form a*x**2 + b*x + c = 0
            if( (b**2 - 4._r8*a*c) .lt. 0._r8 ) then   ! Failure, no real roots
                 status = 3
            else
                 q  = -0.5_r8*(b + sign(1.0_r8,b)*sqrt(b**2 - 4._r8*a*c))
                 r1 =  q/a
                 r2 =  c/q
            endif
       endif
    endif

    return
  end subroutine roots

  function slope(mkx,field,p0)
  ! ------------------------------------------------------------------ !
  ! Function performing profile reconstruction of conservative scalars !
  ! in each layer. This is identical to profile reconstruction used in !
  ! UW-PBL scheme but from bottom to top layer here.     At the lowest !
  ! layer near to surface, slope is defined using the two lowest layer !
  ! mid-point values. I checked this subroutine and it is correct.     !
  ! ------------------------------------------------------------------ !
    integer,  intent(in) :: mkx
    real(r8)             :: slope(mkx)
    real(r8), intent(in) :: field(mkx)
    real(r8), intent(in) :: p0(mkx)

    real(r8)             :: below
    real(r8)             :: above
    integer              :: k

    below = ( field(2) - field(1) ) / ( p0(2) - p0(1) )
    do k = 2, mkx
       above = ( field(k) - field(k-1) ) / ( p0(k) - p0(k-1) )
       if( above .gt. 0._r8 ) then
           slope(k-1) = max(0._r8,min(above,below))
       else
           slope(k-1) = min(0._r8,max(above,below))
       end if
       below = above
    end do
    slope(mkx) = slope(mkx-1)

    return
  end function slope

  function qsinvert(qt,thl,psfc)
  ! ----------------------------------------------------------------- !
  ! Function calculating saturation pressure ps (or pLCL) from qt and !
  ! thl ( liquid potential temperature,  NOT liquid virtual potential !
  ! temperature) by inverting Bolton formula. I should check later if !
  ! current use of 'leff' instead of 'xlv' here is reasonable or not. !
  ! ----------------------------------------------------------------- !
    real(r8)          :: qsinvert
    real(r8)             qt, thl, psfc
    real(r8)             ps, Pis, Ts, err, dlnqsdT, dTdPis
    real(r8)             dPisdps, dlnqsdps, derrdps, dps
    real(r8)             Ti, rhi, TLCL, PiLCL, psmin, dpsmax
    integer              i
    real(r8)          :: es                     ! saturation vapor pressure
    real(r8)          :: qs                     ! saturation spec. humidity
    real(r8)          :: gam                    ! (L/cp)*dqs/dT
    real(r8)          :: leff, nu

    psmin  = 100._r8*100._r8 ! Default saturation pressure [Pa] if iteration does not converge
    dpsmax = 1._r8           ! Tolerance [Pa] for convergence of iteration

    ! ------------------------------------ !
    ! Calculate best initial guess of pLCL !
    ! ------------------------------------ !

    Ti       =  thl*(psfc/p00)**rovcp
    call t_startf('qsat_o')
#ifdef QSAT_WV
    call qsat_wv(Ti, psfc, es, qs)
#else
    call qsat_o(Ti, psfc, es, qs)
#endif
    call t_stopf('qsat_o')

    rhi      =  qt/qs
    if( rhi .le. 0.01_r8 ) then
      !  write(iulog,*) 'Source air is too dry and pLCL is set to psmin in uwshcu.F90'
        qsinvert = psmin
        return
    end if
    TLCL     =  55._r8 + 1._r8/(1._r8/(Ti-55._r8)-log(rhi)/2840._r8); ! Bolton's formula. MWR.1980.Eq.(22)
    PiLCL    =  TLCL/thl
    ps       =  p00*(PiLCL)**(1._r8/rovcp)

    do i = 1, 10
       Pis      =  (ps/p00)**rovcp
       Ts       =  thl*Pis
       call t_startf('qsat_o')
#ifdef QSAT_WV
       call qsat_wv_gam(Ts, ps, es, qs, gam=gam)
#else
       call qsat_o(Ts, ps, es, qs, gam=gam)
#endif
       call t_stopf('qsat_o')
       err      =  qt - qs
       nu       =  max(min((268._r8 - Ts)/20._r8,1.0_r8),0.0_r8)
       leff     =  (1._r8 - nu)*xlv + nu*xls
       dlnqsdT  =  gam*(cp/leff)/qs
       dTdPis   =  thl
       dPisdps  =  rovcp*Pis/ps
       dlnqsdps = -1._r8/(ps - (1._r8 - ep2)*es)
       derrdps  = -qs*(dlnqsdT * dTdPis * dPisdps + dlnqsdps)
       dps      = -err/derrdps
       ps       =  ps + dps
       if( ps .lt. 0._r8 ) then
           write(iulog,*) 'pLCL iteration is negative and set to psmin in uwshcu.F90', qt, thl, psfc
           qsinvert = psmin
           return
       end if
       if( abs(dps) .le. dpsmax ) then
           qsinvert = ps
           return
       end if
    end do
    write(iulog,*) 'pLCL does not converge and is set to psmin in uwshcu.F90', qt, thl, psfc
    qsinvert = psmin
    return
  end function qsinvert

  real(r8) function compute_alpha(del_CIN,ke)
  ! ------------------------------------------------ !
  ! Subroutine to compute proportionality factor for !
  ! implicit CIN calculation.                        !
  ! ------------------------------------------------ !
    real(r8) :: del_CIN, ke
    real(r8) :: x0, x1

    integer  :: iteration

    x0 = 0._r8
    do iteration = 1, 10
       x1 = x0 - (exp(-x0*ke*del_CIN) - x0)/(-ke*del_CIN*exp(-x0*ke*del_CIN) - 1._r8)
       x0 = x1
    end do
    compute_alpha = x0

    return

  end function compute_alpha

  real(r8) function compute_mumin2(mulcl,rmaxfrac,mulow)
  ! --------------------------------------------------------- !
  ! Subroutine to compute critical 'mu' (normalized CIN) such !
  ! that updraft fraction at the LCL is equal to 'rmaxfrac'.  !
  ! --------------------------------------------------------- !
    real(r8) :: mulcl, rmaxfrac, mulow
    real(r8) :: x0, x1, ex, ef, exf, f, fs
    integer  :: iteration

    x0 = mulow
    do iteration = 1, 10
       ex = exp(-x0**2)
       ef = erfc(x0)
       ! if(x0.ge.3._r8) then
       !    compute_mumin2 = 3._r8
       !    goto 20
       ! endif
       exf = ex/ef
       f  = 0.5_r8*exf**2 - 0.5_r8*(ex/2._r8/rmaxfrac)**2 - (mulcl*2.5066_r8/2._r8)**2
       fs = (2._r8*exf**2)*(exf/sqrt(3.141592_r8)-x0) + (0.5_r8*x0*ex**2)/(rmaxfrac**2)
       x1 = x0 - f/fs
       x0 = x1
    end do
    compute_mumin2 = x0

 20 return

  end function compute_mumin2

  real(r8) function compute_ppen(wtwb,D,bogbot,bogtop,rho0j,dpen)
  ! ----------------------------------------------------------- !
  ! Subroutine to compute critical 'ppen[Pa]<0' ( pressure dis. !
  ! from 'ps0(kpen-1)' to the cumulus top where cumulus updraft !
  ! vertical velocity is exactly zero ) by considering exact    !
  ! non-zero fer(kpen).                                         !
  ! ----------------------------------------------------------- !
    real(r8) :: wtwb, D, bogbot, bogtop, rho0j, dpen
    real(r8) :: x0, x1, f, fs, SB, s00
    integer  :: iteration

    ! Buoyancy slope
      SB = ( bogtop - bogbot ) / dpen
    ! Sign of slope, 'f' at x = 0
    ! If 's00>0', 'w' increases with height.
      s00 = bogbot / rho0j - D * wtwb

    if( D*dpen .lt. 1.e-8_r8 ) then
        if( s00 .ge. 0._r8 ) then
            x0 = dpen
        else
            x0 = max(0._r8,min(dpen,-0.5_r8*wtwb/s00))
        endif
    else
        if( s00 .ge. 0._r8 ) then
            x0 = dpen
        else
            x0 = 0._r8
        endif
        do iteration = 1, 5
           f  = exp(-2._r8*D*x0)*(wtwb-(bogbot-SB/(2._r8*D))/(D*rho0j)) + &
                                 (SB*x0+bogbot-SB/(2._r8*D))/(D*rho0j)
           fs = -2._r8*D*exp(-2._r8*D*x0)*(wtwb-(bogbot-SB/(2._r8*D))/(D*rho0j)) + &
                                 (SB)/(D*rho0j)
           if( fs .ge. 0._r8 ) then
		fs = max(fs, 1.e-10_r8)
           else
           	fs = min(fs,-1.e-10_r8)
           endif
           x1 = x0 - f/fs
           x0 = x1
      end do

    endif

    compute_ppen = -max(0._r8,min(dpen,x0))

  end function compute_ppen

  subroutine fluxbelowinv(cbmf,ps0,mkx,kinv,dt,xsrc,xmean,xtopin,xbotin,xflx)
  ! ------------------------------------------------------------------------- !
  ! Subroutine to calculate turbulent fluxes at and below 'kinv-1' interfaces.!
  ! Check in the main program such that input 'cbmf' should not be zero.      !
  ! If the reconstructed inversion height does not go down below the 'kinv-1' !
  ! interface, then turbulent flux at 'kinv-1' interface  is simply a product !
  ! of 'cmbf' and 'qtsrc-xbot' where 'xbot' is the value at the top interface !
  ! of 'kinv-1' layer. This flux is linearly interpolated down to the surface !
  ! assuming turbulent fluxes at surface are zero. If reconstructed inversion !
  ! height goes down below the 'kinv-1' interface, subsidence warming &drying !
  ! measured by 'xtop-xbot', where  'xtop' is the value at the base interface !
  ! of 'kinv+1' layer, is added ONLY to the 'kinv-1' layer, using appropriate !
  ! mass weighting ( rpinv and rcbmf, or rr = rpinv / rcbmf ) between current !
  ! and next provisional time step. Also impose a limiter to enforce outliers !
  ! of thermodynamic variables in 'kinv' layer  to come back to normal values !
  ! at the next step.                                                         !
  ! ------------------------------------------------------------------------- !
    integer,  intent(in)                     :: mkx, kinv
    real(r8), intent(in)                     :: cbmf, dt, xsrc, xmean, xtopin, xbotin
    real(r8), intent(in),  dimension(0:mkx)  :: ps0
    real(r8), intent(out), dimension(0:mkx)  :: xflx
    integer k
    real(r8) rcbmf, rpeff, dp, rr, pinv_eff, xtop, xbot, pinv, xtop_ori, xbot_ori

    xflx(0:mkx) = 0._r8
    dp = ps0(kinv-1) - ps0(kinv)
    xbot = xbotin
    xtop = xtopin

    ! -------------------------------------- !
    ! Compute reconstructed inversion height !
    ! -------------------------------------- !
    xtop_ori = xtop
    xbot_ori = xbot
    rcbmf = ( cbmf * g * dt ) / dp                  ! Can be larger than 1 : 'OK'

    if( xbot .ge. xtop ) then
        rpeff = ( xmean - xtop ) / max(  1.e-20_r8, xbot - xtop )
    else
        rpeff = ( xmean - xtop ) / min( -1.e-20_r8, xbot - xtop )
    endif

    rpeff = min( max(0._r8,rpeff), 1._r8 )          ! As of this, 0<= rpeff <= 1
    if( rpeff .eq. 0._r8 .or. rpeff .eq. 1._r8 ) then
        xbot = xmean
        xtop = xmean
    endif
    ! Below two commented-out lines are the old code replacing the above 'if' block.
    ! if(rpeff.eq.1) xbot = xmean
    ! if(rpeff.eq.0) xtop = xmean
    rr       = rpeff / rcbmf
    pinv     = ps0(kinv-1) - rpeff * dp             ! "pinv" before detraining mass
    pinv_eff = ps0(kinv-1) + ( rcbmf - rpeff ) * dp ! Effective "pinv" after detraining mass
    ! ----------------------------------------------------------------------- !
    ! Compute turbulent fluxes.                                               !
    ! Below two cases exactly converges at 'kinv-1' interface when rr = 1._r8 !
    ! ----------------------------------------------------------------------- !
    do k = 0, kinv - 1
       xflx(k) = cbmf * ( xsrc - xbot ) * ( ps0(0) - ps0(k) ) / ( ps0(0) - pinv )
    end do
    if( rr .le. 1._r8 ) then
        xflx(kinv-1) =  xflx(kinv-1) - ( 1._r8 - rr ) * cbmf * ( xtop_ori - xbot_ori )
    endif

    return
  end subroutine fluxbelowinv

  subroutine positive_moisture_single( xlv, xls, mkx, dt, qvmin, qlmin, qimin, dp, qv, ql, qi, s, qvten, qlten, qiten, sten )
  ! ------------------------------------------------------------------------------- !
  ! If any 'ql < qlmin, qi < qimin, qv < qvmin' are developed in any layer,         !
  ! force them to be larger than minimum value by (1) condensating water vapor      !
  ! into liquid or ice, and (2) by transporting water vapor from the very lower     !
  ! layer. '2._r8' is multiplied to the minimum values for safety.                  !
  ! Update final state variables and tendencies associated with this correction.    !
  ! If any condensation happens, update (s,t) too.                                  !
  ! Note that (qv,ql,qi,s) are final state variables after applying corresponding   !
  ! input tendencies and corrective tendencies                                      !
  ! ------------------------------------------------------------------------------- !
    implicit none
    integer,  intent(in)     :: mkx
    real(r8), intent(in)     :: xlv, xls
    real(r8), intent(in)     :: dt, qvmin, qlmin, qimin
    real(r8), intent(in)     :: dp(mkx)
    real(r8), intent(inout)  :: qv(mkx), ql(mkx), qi(mkx), s(mkx)
    real(r8), intent(inout)  :: qvten(mkx), qlten(mkx), qiten(mkx), sten(mkx)
    integer   k
    real(r8)  dql, dqi, dqv, sum, aa, dum

    do k = mkx, 1, -1        ! From the top to the 1st (lowest) layer from the surface
       dql = max(0._r8,1._r8*qlmin-ql(k))
       dqi = max(0._r8,1._r8*qimin-qi(k))
       qlten(k) = qlten(k) +  dql/dt
       qiten(k) = qiten(k) +  dqi/dt
       qvten(k) = qvten(k) - (dql+dqi)/dt
       sten(k)  = sten(k)  + xlv * (dql/dt) + xls * (dqi/dt)
       ql(k)    = ql(k) +  dql
       qi(k)    = qi(k) +  dqi
       qv(k)    = qv(k) -  dql - dqi
       s(k)     = s(k)  +  xlv * dql + xls * dqi
       dqv      = max(0._r8,1._r8*qvmin-qv(k))
       qvten(k) = qvten(k) + dqv/dt
       qv(k)    = qv(k)   + dqv
       if( k .ne. 1 ) then
           qv(k-1)    = qv(k-1)    - dqv*dp(k)/dp(k-1)
           qvten(k-1) = qvten(k-1) - dqv*dp(k)/dp(k-1)/dt
       endif
       qv(k) = max(qv(k),qvmin)
       ql(k) = max(ql(k),qlmin)
       qi(k) = max(qi(k),qimin)
    end do
    ! Extra moisture used to satisfy 'qv(i,1)=qvmin' is proportionally
    ! extracted from all the layers that has 'qv > 2*qvmin'. This fully
    ! preserves column moisture.
    if( dqv .gt. 1.e-20_r8 ) then
        sum = 0._r8
        do k = 1, mkx
           if( qv(k) .gt. 2._r8*qvmin ) sum = sum + qv(k)*dp(k)
        enddo
        aa = dqv*dp(1)/max(1.e-20_r8,sum)
        if( aa .lt. 0.5_r8 ) then
            do k = 1, mkx
               if( qv(k) .gt. 2._r8*qvmin ) then
                   dum      = aa*qv(k)
                   qv(k)    = qv(k) - dum
                   qvten(k) = qvten(k) - dum/dt
               endif
            enddo
        else
            write(iulog,*) 'Full positive_moisture is impossible in uwshcu'
        endif
    endif

    return
  end subroutine positive_moisture_single

  ! ------------------------ !
  !                          !
  ! End of subroutine blocks !
  !                          !
  ! ------------------------ !

  end module uwshcu
