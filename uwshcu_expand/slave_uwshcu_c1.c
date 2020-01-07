
#include <slave.h>
#include "dma_macros.h"

typedef struct {
  double *ptr_param_uw1;
  int iend;
} param_uw_t;

typedef struct {
    double *zs0
    double *ps0
    double *p0
    double *tke
    double *thvl0bot
    double *thvl0top
    double *qt0
    double *u0
    double *ssu0
    double *v0
    double *ssv0
    double *trsrc
    double *thl0
    double *ssthl0
    double *ssqt0
    double *thv0bot
    double *thv0top
    double *trsrc_o
    double *qv0
    double *qv0_o
    double *ql0
    double *ql0_o
    double *qi0
    double *qi0_o
    double *t0
    double *t0_o
    double *s0
    double *s0_o
    double *u0_o
    double *v0_o
    double *qt0_o
    double *thl0_o
    double *thvl0
    double *thvl0_o
    double *ssthl0_o
    double *ssqt0_o
    double *thv0bot_o
    double *thv0top_o
    double *thvl0bot_o
    double *thvl0top_o
    double *ssu0_o
    double *ssv0_o
    double *umf
    double *emf
    double *slflx
    double *qtflx
    double *uflx
    double *vflx
    double *qvten
    double *qlten
    double *qiten
    double *sten
    double *uten
    double *vten
    double *qrten
    double *qsten
    double *dwten
    double *diten
    double *evapc
    double *cufrc
    double *qcu
    double *qlu
    double *qiu
    double *fer
    double *fdr
    double *qc
    double *qc_l
    double *qc_i
    double *qtten
    double *slten
    double *ufrc
    double *thlu
    double *qtu
    double *uu
    double *vu
    double *wu
    double *thvu
    double *thlu_emf
    double *qtu_emf
    double *uu_emf
    double *vu_emf
    double *excessu_arr
    double *excess0_arr
    double *xc_arr
    double *aquad_arr
    double *bquad_arr
    double *cquad_arr
    double *bogbot_arr
    double *bogtop_arr
    double *umf_out
    double *umf_s
    double *qvten_out
    double *qvten_s
    double *qlten_out
    double *qlten_s
    double *qiten_out
    double *qiten_s
    double *sten_out
    double *sten_s
    double *uten_out
    double *uten_s
    double *vten_out
    double *vten_s
    double *qrten_out
    double *qrten_s
    double *qsten_out
    double *qsten_s
    double *evapc_out
    double *evapc_s
    double *cufrc_out
    double *cufrc_s
    double *slflx_out
    double *slflx_s
    double *qtflx_out
    double *qtflx_s
    double *qcu_out
    double *qcu_s
    double *qlu_out
    double *qlu_s
    double *qiu_out
    double *qiu_s
    double *qc_out
    double *qc_s
    double *fer_out
    double *fer_s
    double *fdr_out
    double *fdr_s
    double *qtten_out
    double *qtten_s
    double *slten_out
    double *slten_s
    double *ufrc_out
    double *ufrc_s
    double *uflx_out
    double *uflx_s
    double *vflx_out
    double *vflx_s
    double *wu_out
    double *wu_s
    double *qtu_out
    double *qtu_s
    double *thlu_out
    double *thlu_s
    double *thvu_out
    double *thvu_s
    double *uu_out
    double *uu_s
    double *vu_out
    double *vu_s
    double *qtu_emf_out
    double *qtu_emf_s
    double *thlu_emf_out
    double *thlu_emf_s
    double *uu_emf_out
    double *uu_emf_s
    double *vu_emf_out
    double *vu_emf_s
    double *uemf_out
    double *uemf_s
    double *dwten_out
    double *dwten_s
    double *diten_out
    double *diten_s
    double *flxrain_out
    double *flxrain_s
    double *flxsnow_out
    double *flxsnow_s
    double *ntraprd_out
    double *ntraprd_s
    double *ntsnprd_out
    double *ntsnprd_s
    double *excessu_arr_out
    double *excessu_arr_s
    double *excess0_arr_out
    double *excess0_arr_s
    double *xc_arr_out
    double *xc_arr_s
    double *aquad_arr_out
    double *aquad_arr_s
    double *bquad_arr_out
    double *bquad_arr_s
    double *cquad_arr_out
    double *cquad_arr_s
    double *bogbot_arr_out
    double *bogbot_arr_s
    double *bogtop_arr_out
    double *bogtop_arr_s
    double *exns0
    double *dp0
    double *write
    double *tre
    double *rei
    double *exn0
    double *xflx
    double *uemf
    double *comsub
    double *nlten_sink
    double *niten_sink
    double *uf
    double *vf
    double *qlten_sink
    double *qiten_sink
    double *flxrain
    double *flxsnow
    double *ntraprd
    double *ntsnprd
    double *qw0_in
    double *qv0_star
    double *ql0_star
    double *qi0_star
    double *s0_star
    double *trflx_d
    double *trflx_u
    double *dpdry0
    double *qv0_s
    double *ql0_s
    double *qi0_s
    double *s0_s
    double *u0_s
    double *v0_s
    double *qt0_s
    double *t0_s
    double *flxprc1_out
    double *flxsnow1_out
    double *tr0
    double *tr0_o
    double *sstr0
    double *sstr0_o
    double *trflx
    double *trten
    double *tru
    double *tru_emf
    double *trten_out
    double *trten_s
    double *trflx_out
    double *trflx_s
    double *tru_out
    double *tru_s
    double *tru_emf_out
    double *tru_emf_s
    double *tr0_s
    double tscaleh
    double cush
    double pblh
    double exit_kinv1
    double id_exit
    double dpsum
    double tkeavg
    double thvlmin
    double dpi
    double qtsrc
    double thvlsrc
    double thlsrc
    double usrc
    double vsrc
    double plcl
    double exit_klclmkx
    double thl0lcl
    double qt0lcl
    double thj
    double qvj
    double qlj
    double qij
    double qse
    double exit_conden
    double thv0lcl
    double cin
    double cinlcl
    double plfc
    double thvubot
    double thvutop
    double limit_cinlcl
    double limit_cin
    double exit_klfcmkx
    double cin_i
    double cinlcl_i
    double ke
    double kinv_o
    double klcl_o
    double klfc_o
    double plcl_o
    double plfc_o
    double tkeavg_o
    double thvlmin_o
    double qtsrc_o
    double thvlsrc_o
    double thlsrc_o
    double usrc_o
    double vsrc_o
    double thv0lcl_o
    double cin_f
    double cinlcl_f
    double del_CIN
    double alpha
    double del_cinlcl
    double precip
    double snow
    double rliq
    double cbmf
    double cnt
    double cnb
    double ufrcinvbase
    double ufrclcl
    double winvbase
    double wlcl
    double emfkbup
    double cbmflimit
    double ind_delcin
    double precip_out
    double precip_s
    double snow_out
    double snow_s
    double cush_inout
    double cush_s
    double cbmf_out
    double cbmf_s
    double rliq_out
    double rliq_s
    double cnt_out
    double cnt_s
    double cnb_out
    double cnb_s
    double cinh_out
    double cin_s
    double cinlclh_out
    double cinlcl_s
    double ufrcinvbase_out
    double ufrcinvbase_s
    double ufrclcl_out
    double ufrclcl_s
    double winvbase_out
    double winvbase_s
    double wlcl_out
    double wlcl_s
    double plcl_out
    double plcl_s
    double pinv_out
    double pinv_s
    double plfc_out
    double plfc_s
    double pbup_out
    double pbup_s
    double ppen_out
    double ppen_s
    double qtsrc_out
    double qtsrc_s
    double thlsrc_out
    double thlsrc_s
    double thvlsrc_out
    double thvlsrc_s
    double emfkbup_out
    double emfkbup_s
    double cbmflimit_out
    double cbmflimit_s
    double tkeavg_out
    double tkeavg_s
    double zinv_out
    double zinv_s
    double rcwp_out
    double rcwp_s
    double rlwp_out
    double rlwp_s
    double riwp_out
    double riwp_s
    double id_exit_ex
    double prel
    double thv0rel
    double wcrit
    double sigmaw
    double mu
    double rho0inv
    double mumin0
    double mulcl
    double mulclstar
    double erfc
    double mumin2
    double limit_ufrc
    double limit_cbmf
    double winv
    double ufrcinv
    double wtw
    double exit_wtw
    double wrel
    double exit_ufrc
    double uplus
    double vplus
    double pe
    double dpe
    double exne
    double thvebot
    double thle
    double qte
    double ue
    double ve
    double scaleh
    double thlue
    double qtue
    double wue
    double wtwb
    double thv0j
    double rho0j
    double qsat_arg
    double es
    double weight
    double estbl
    double qs
    double excess0
    double exql
    double exqi
    double thvj
    double excessu
    double xc
    double cridis
    double aquad
    double bquad
    double cquad
    double xsat
    double thlxsat
    double qtxsat
    double thvxsat
    double thv_x0
    double thv_x1
    double xs1
    double xs2
    double x_cu
    double x_en
    double ee2
    double ud2
    double delbog
    double bogtop
    double bogbot
    double drage
    double expfac
    double autodet
    double exit_wu
    double rhos0j
    double xc1
    double xc2
    double ppen
    double limit_ppen
    double thlu_top
    double qtu_top
    double exntop
    double forcedCu
    double limit_shcu
    double cldhgt
    double exit_cufilter
    double limit_emf
    double xsrc
    double xmean
    double xtop
    double xbot
    double thlten_sub
    double qtten_sub
    double qlten_sub
    double qiten_sub
    double nlten_sub
    double niten_sub
    double thl_prog
    double qt_prog
    double rainflx
    double snowflx
    double qlu_mid
    double qiu_mid
    double qlubelow
    double qiubelow
    double qlu_top
    double qiu_top
    double qc_lm
    double qc_im
    double nc_lm
    double nc_im
    double ql_emf_kbup
    double qi_emf_kbup
    double nl_emf_kbup
    double ni_emf_kbup
    double qlten_det
    double qiten_det
    double evpint_rain
    double evpint_snow
    double snowmlt
    double subsat
    double evprain
    double evpsnow
    double evplimit
    double evplimit_rain
    double evplimit_snow
    double tmp1
    double tmp2
    double limit_negcon
    double trmin
    double pdelx
    double dum
    double pdel
    double qcubelow
    double rcwp
    double rlwp
    double riwp
    double thl0bot
    double qt0bot
    double thl0top
    double qt0top
    int kinv
    int klcl
    int klfc
    int krel
    int kbup
} param_t;

void slave_uwsch_c1_(param_uw_t *param_uw_s)
{
  volatile int id = athread_get_id(-1);
  volatile int unsigned long get_reply, put_reply;
  dma_init();
  param_uw_t param_uw_d;
  param_t param_d;
  pe_get(param_uw_s, &param_uw_d, sizeof(param_uw_t));
  dma_syn();
  int limit = param_uw_d.iend;
  if(id < limit) {
    pe_get(param_uw_d.ptr_param_uw1 + id, &param_d, sizeof(param_t));
    dma_syn();
    double tscaleh = param_d.tscaleh;
    double cush = param_d.cush;
    double pblh = param_d.pblh;
    double exit_kinv1 = param_d.exit_kinv1;
    double id_exit = param_d.id_exit;
    double dpsum = param_d.dpsum;
    double tkeavg = param_d.tkeavg;
    double thvlmin = param_d.thvlmin;
    double dpi = param_d.dpi;
    double qtsrc = param_d.qtsrc;
    double thvlsrc = param_d.thvlsrc;
    double thlsrc = param_d.thlsrc;
    double usrc = param_d.usrc;
    double vsrc = param_d.vsrc;
    double plcl = param_d.plcl;
    double exit_klclmkx = param_d.exit_klclmkx;
    double thl0lcl = param_d.thl0lcl;
    double qt0lcl = param_d.qt0lcl;
    double thj = param_d.thj;
    double qvj = param_d.qvj;
    double qlj = param_d.qlj;
    double qij = param_d.qij;
    double qse = param_d.qse;
    double exit_conden = param_d.exit_conden;
    double thv0lcl = param_d.thv0lcl;
    double cin = param_d.cin;
    double cinlcl = param_d.cinlcl;
    double plfc = param_d.plfc;
    double thvubot = param_d.thvubot;
    double thvutop = param_d.thvutop;
    double limit_cinlcl = param_d.limit_cinlcl;
    double limit_cin = param_d.limit_cin;
    double exit_klfcmkx = param_d.exit_klfcmkx;
    double cin_i = param_d.cin_i;
    double cinlcl_i = param_d.cinlcl_i;
    double ke = param_d.ke;
    double kinv_o = param_d.kinv_o;
    double klcl_o = param_d.klcl_o;
    double klfc_o = param_d.klfc_o;
    double plcl_o = param_d.plcl_o;
    double plfc_o = param_d.plfc_o;
    double tkeavg_o = param_d.tkeavg_o;
    double thvlmin_o = param_d.thvlmin_o;
    double qtsrc_o = param_d.qtsrc_o;
    double thvlsrc_o = param_d.thvlsrc_o;
    double thlsrc_o = param_d.thlsrc_o;
    double usrc_o = param_d.usrc_o;
    double vsrc_o = param_d.vsrc_o;
    double thv0lcl_o = param_d.thv0lcl_o;
    double cin_f = param_d.cin_f;
    double cinlcl_f = param_d.cinlcl_f;
    double del_CIN = param_d.del_CIN;
    double alpha = param_d.alpha;
    double del_cinlcl = param_d.del_cinlcl;
    double precip = param_d.precip;
    double snow = param_d.snow;
    double rliq = param_d.rliq;
    double cbmf = param_d.cbmf;
    double cnt = param_d.cnt;
    double cnb = param_d.cnb;
    double ufrcinvbase = param_d.ufrcinvbase;
    double ufrclcl = param_d.ufrclcl;
    double winvbase = param_d.winvbase;
    double wlcl = param_d.wlcl;
    double emfkbup = param_d.emfkbup;
    double cbmflimit = param_d.cbmflimit;
    double ind_delcin = param_d.ind_delcin;
    double precip_out = param_d.precip_out;
    double precip_s = param_d.precip_s;
    double snow_out = param_d.snow_out;
    double snow_s = param_d.snow_s;
    double cush_inout = param_d.cush_inout;
    double cush_s = param_d.cush_s;
    double cbmf_out = param_d.cbmf_out;
    double cbmf_s = param_d.cbmf_s;
    double rliq_out = param_d.rliq_out;
    double rliq_s = param_d.rliq_s;
    double cnt_out = param_d.cnt_out;
    double cnt_s = param_d.cnt_s;
    double cnb_out = param_d.cnb_out;
    double cnb_s = param_d.cnb_s;
    double cinh_out = param_d.cinh_out;
    double cin_s = param_d.cin_s;
    double cinlclh_out = param_d.cinlclh_out;
    double cinlcl_s = param_d.cinlcl_s;
    double ufrcinvbase_out = param_d.ufrcinvbase_out;
    double ufrcinvbase_s = param_d.ufrcinvbase_s;
    double ufrclcl_out = param_d.ufrclcl_out;
    double ufrclcl_s = param_d.ufrclcl_s;
    double winvbase_out = param_d.winvbase_out;
    double winvbase_s = param_d.winvbase_s;
    double wlcl_out = param_d.wlcl_out;
    double wlcl_s = param_d.wlcl_s;
    double plcl_out = param_d.plcl_out;
    double plcl_s = param_d.plcl_s;
    double pinv_out = param_d.pinv_out;
    double pinv_s = param_d.pinv_s;
    double plfc_out = param_d.plfc_out;
    double plfc_s = param_d.plfc_s;
    double pbup_out = param_d.pbup_out;
    double pbup_s = param_d.pbup_s;
    double ppen_out = param_d.ppen_out;
    double ppen_s = param_d.ppen_s;
    double qtsrc_out = param_d.qtsrc_out;
    double qtsrc_s = param_d.qtsrc_s;
    double thlsrc_out = param_d.thlsrc_out;
    double thlsrc_s = param_d.thlsrc_s;
    double thvlsrc_out = param_d.thvlsrc_out;
    double thvlsrc_s = param_d.thvlsrc_s;
    double emfkbup_out = param_d.emfkbup_out;
    double emfkbup_s = param_d.emfkbup_s;
    double cbmflimit_out = param_d.cbmflimit_out;
    double cbmflimit_s = param_d.cbmflimit_s;
    double tkeavg_out = param_d.tkeavg_out;
    double tkeavg_s = param_d.tkeavg_s;
    double zinv_out = param_d.zinv_out;
    double zinv_s = param_d.zinv_s;
    double rcwp_out = param_d.rcwp_out;
    double rcwp_s = param_d.rcwp_s;
    double rlwp_out = param_d.rlwp_out;
    double rlwp_s = param_d.rlwp_s;
    double riwp_out = param_d.riwp_out;
    double riwp_s = param_d.riwp_s;
    double id_exit_ex = param_d.id_exit_ex;
    double prel = param_d.prel;
    double thv0rel = param_d.thv0rel;
    double wcrit = param_d.wcrit;
    double sigmaw = param_d.sigmaw;
    double mu = param_d.mu;
    double rho0inv = param_d.rho0inv;
    double mumin0 = param_d.mumin0;
    double mulcl = param_d.mulcl;
    double mulclstar = param_d.mulclstar;
    double erfc = param_d.erfc;
    double mumin2 = param_d.mumin2;
    double limit_ufrc = param_d.limit_ufrc;
    double limit_cbmf = param_d.limit_cbmf;
    double winv = param_d.winv;
    double ufrcinv = param_d.ufrcinv;
    double wtw = param_d.wtw;
    double exit_wtw = param_d.exit_wtw;
    double wrel = param_d.wrel;
    double exit_ufrc = param_d.exit_ufrc;
    double uplus = param_d.uplus;
    double vplus = param_d.vplus;
    double pe = param_d.pe;
    double dpe = param_d.dpe;
    double exne = param_d.exne;
    double thvebot = param_d.thvebot;
    double thle = param_d.thle;
    double qte = param_d.qte;
    double ue = param_d.ue;
    double ve = param_d.ve;
    double scaleh = param_d.scaleh;
    double thlue = param_d.thlue;
    double qtue = param_d.qtue;
    double wue = param_d.wue;
    double wtwb = param_d.wtwb;
    double thv0j = param_d.thv0j;
    double rho0j = param_d.rho0j;
    double qsat_arg = param_d.qsat_arg;
    double es = param_d.es;
    double weight = param_d.weight;
    double estbl = param_d.estbl;
    double qs = param_d.qs;
    double excess0 = param_d.excess0;
    double exql = param_d.exql;
    double exqi = param_d.exqi;
    double thvj = param_d.thvj;
    double excessu = param_d.excessu;
    double xc = param_d.xc;
    double cridis = param_d.cridis;
    double aquad = param_d.aquad;
    double bquad = param_d.bquad;
    double cquad = param_d.cquad;
    double xsat = param_d.xsat;
    double thlxsat = param_d.thlxsat;
    double qtxsat = param_d.qtxsat;
    double thvxsat = param_d.thvxsat;
    double thv_x0 = param_d.thv_x0;
    double thv_x1 = param_d.thv_x1;
    double xs1 = param_d.xs1;
    double xs2 = param_d.xs2;
    double x_cu = param_d.x_cu;
    double x_en = param_d.x_en;
    double ee2 = param_d.ee2;
    double ud2 = param_d.ud2;
    double delbog = param_d.delbog;
    double bogtop = param_d.bogtop;
    double bogbot = param_d.bogbot;
    double drage = param_d.drage;
    double expfac = param_d.expfac;
    double autodet = param_d.autodet;
    double exit_wu = param_d.exit_wu;
    double rhos0j = param_d.rhos0j;
    double xc1 = param_d.xc1;
    double xc2 = param_d.xc2;
    double ppen = param_d.ppen;
    double limit_ppen = param_d.limit_ppen;
    double thlu_top = param_d.thlu_top;
    double qtu_top = param_d.qtu_top;
    double exntop = param_d.exntop;
    double forcedCu = param_d.forcedCu;
    double limit_shcu = param_d.limit_shcu;
    double cldhgt = param_d.cldhgt;
    double exit_cufilter = param_d.exit_cufilter;
    double limit_emf = param_d.limit_emf;
    double xsrc = param_d.xsrc;
    double xmean = param_d.xmean;
    double xtop = param_d.xtop;
    double xbot = param_d.xbot;
    double thlten_sub = param_d.thlten_sub;
    double qtten_sub = param_d.qtten_sub;
    double qlten_sub = param_d.qlten_sub;
    double qiten_sub = param_d.qiten_sub;
    double nlten_sub = param_d.nlten_sub;
    double niten_sub = param_d.niten_sub;
    double thl_prog = param_d.thl_prog;
    double qt_prog = param_d.qt_prog;
    double rainflx = param_d.rainflx;
    double snowflx = param_d.snowflx;
    double qlu_mid = param_d.qlu_mid;
    double qiu_mid = param_d.qiu_mid;
    double qlubelow = param_d.qlubelow;
    double qiubelow = param_d.qiubelow;
    double qlu_top = param_d.qlu_top;
    double qiu_top = param_d.qiu_top;
    double qc_lm = param_d.qc_lm;
    double qc_im = param_d.qc_im;
    double nc_lm = param_d.nc_lm;
    double nc_im = param_d.nc_im;
    double ql_emf_kbup = param_d.ql_emf_kbup;
    double qi_emf_kbup = param_d.qi_emf_kbup;
    double nl_emf_kbup = param_d.nl_emf_kbup;
    double ni_emf_kbup = param_d.ni_emf_kbup;
    double qlten_det = param_d.qlten_det;
    double qiten_det = param_d.qiten_det;
    double evpint_rain = param_d.evpint_rain;
    double evpint_snow = param_d.evpint_snow;
    double snowmlt = param_d.snowmlt;
    double subsat = param_d.subsat;
    double evprain = param_d.evprain;
    double evpsnow = param_d.evpsnow;
    double evplimit = param_d.evplimit;
    double evplimit_rain = param_d.evplimit_rain;
    double evplimit_snow = param_d.evplimit_snow;
    double tmp1 = param_d.tmp1;
    double tmp2 = param_d.tmp2;
    double limit_negcon = param_d.limit_negcon;
    double trmin = param_d.trmin;
    double pdelx = param_d.pdelx;
    double dum = param_d.dum;
    double pdel = param_d.pdel;
    double qcubelow = param_d.qcubelow;
    double rcwp = param_d.rcwp;
    double rlwp = param_d.rlwp;
    double riwp = param_d.riwp;
    double thl0bot = param_d.thl0bot;
    double qt0bot = param_d.qt0bot;
    double thl0top = param_d.thl0top;
    double qt0top = param_d.qt0top;
    int kinv = param_d.kinv;
    int klcl = param_d.klcl;
    int klfc = param_d.klfc;
    int krel = param_d.krel;
    int kbup = param_d.kbup;
    double zs0[mkx + 1];
    double ps0[mkx + 1];
    double p0[mkx];
    double tke[mkx + 1];
    double thvl0bot[mkx];
    double thvl0top[mkx];
    double qt0[mkx];
    double u0[mkx];
    double ssu0[mkx];
    double v0[mkx];
    double ssv0[mkx];
    double trsrc[mkx];
    double thl0[mkx];
    double ssthl0[mkx];
    double ssqt0[mkx];
    double thv0bot[mkx];
    double thv0top[mkx];
    double trsrc_o[mkx];
    double qv0[mkx];
    double qv0_o[mkx];
    double ql0[mkx];
    double ql0_o[mkx];
    double qi0[mkx];
    double qi0_o[mkx];
    double t0[mkx];
    double t0_o[mkx];
    double s0[mkx];
    double s0_o[mkx];
    double u0_o[mkx];
    double v0_o[mkx];
    double qt0_o[mkx];
    double thl0_o[mkx];
    double thvl0[mkx];
    double thvl0_o[mkx];
    double ssthl0_o[mkx];
    double ssqt0_o[mkx];
    double thv0bot_o[mkx];
    double thv0top_o[mkx];
    double thvl0bot_o[mkx];
    double thvl0top_o[mkx];
    double ssu0_o[mkx];
    double ssv0_o[mkx];
    double umf[mkx + 1];
    double emf[mkx + 1];
    double slflx[mkx + 1];
    double qtflx[mkx + 1];
    double uflx[mkx + 1];
    double vflx[mkx + 1];
    double qvten[mkx];
    double qlten[mkx];
    double qiten[mkx];
    double sten[mkx];
    double uten[mkx];
    double vten[mkx];
    double qrten[mkx];
    double qsten[mkx];
    double dwten[mkx];
    double diten[mkx];
    double evapc[mkx];
    double cufrc[mkx];
    double qcu[mkx];
    double qlu[mkx];
    double qiu[mkx];
    double fer[mkx];
    double fdr[mkx];
    double qc[mkx];
    double qc_l[mkx];
    double qc_i[mkx];
    double qtten[mkx];
    double slten[mkx];
    double ufrc[mkx + 1];
    double thlu[mkx + 1];
    double qtu[mkx + 1];
    double uu[mkx + 1];
    double vu[mkx + 1];
    double wu[mkx + 1];
    double thvu[mkx + 1];
    double thlu_emf[mkx + 1];
    double qtu_emf[mkx + 1];
    double uu_emf[mkx + 1];
    double vu_emf[mkx + 1];
    double excessu_arr[mkx];
    double excess0_arr[mkx];
    double xc_arr[mkx];
    double aquad_arr[mkx];
    double bquad_arr[mkx];
    double cquad_arr[mkx];
    double bogbot_arr[mkx];
    double bogtop_arr[mkx];
    double umf_out[mkx + 1];
    double umf_s[mkx];
    double qvten_out[mkx];
    double qvten_s[mkx];
    double qlten_out[mkx];
    double qlten_s[mkx];
    double qiten_out[mkx];
    double qiten_s[mkx];
    double sten_out[mkx];
    double sten_s[mkx];
    double uten_out[mkx];
    double uten_s[mkx];
    double vten_out[mkx];
    double vten_s[mkx];
    double qrten_out[mkx];
    double qrten_s[mkx];
    double qsten_out[mkx];
    double qsten_s[mkx];
    double evapc_out[mkx];
    double evapc_s[mkx];
    double cufrc_out[mkx];
    double cufrc_s[mkx];
    double slflx_out[mkx + 1];
    double slflx_s[mkx];
    double qtflx_out[mkx + 1];
    double qtflx_s[mkx];
    double qcu_out[mkx];
    double qcu_s[mkx];
    double qlu_out[mkx];
    double qlu_s[mkx];
    double qiu_out[mkx];
    double qiu_s[mkx];
    double qc_out[mkx];
    double qc_s[mkx];
    double fer_out[mkx];
    double fer_s[mkx];
    double fdr_out[mkx];
    double fdr_s[mkx];
    double qtten_out[mkx];
    double qtten_s[mkx];
    double slten_out[mkx];
    double slten_s[mkx];
    double ufrc_out[mkx + 1];
    double ufrc_s[mkx];
    double uflx_out[mkx + 1];
    double uflx_s[mkx];
    double vflx_out[mkx + 1];
    double vflx_s[mkx];
    double wu_out[mkx + 1];
    double wu_s[mkx];
    double qtu_out[mkx + 1];
    double qtu_s[mkx];
    double thlu_out[mkx + 1];
    double thlu_s[mkx];
    double thvu_out[mkx];
    double thvu_s[mkx];
    double uu_out[mkx];
    double uu_s[mkx];
    double vu_out[mkx];
    double vu_s[mkx];
    double qtu_emf_out[mkx];
    double qtu_emf_s[mkx];
    double thlu_emf_out[mkx];
    double thlu_emf_s[mkx];
    double uu_emf_out[mkx];
    double uu_emf_s[mkx];
    double vu_emf_out[mkx];
    double vu_emf_s[mkx];
    double uemf_out[mkx];
    double uemf_s[mkx];
    double dwten_out[mkx];
    double dwten_s[mkx];
    double diten_out[mkx];
    double diten_s[mkx];
    double flxrain_out[mkx];
    double flxrain_s[mkx];
    double flxsnow_out[mkx];
    double flxsnow_s[mkx];
    double ntraprd_out[mkx];
    double ntraprd_s[mkx];
    double ntsnprd_out[mkx];
    double ntsnprd_s[mkx];
    double excessu_arr_out[mkx];
    double excessu_arr_s[mkx];
    double excess0_arr_out[mkx];
    double excess0_arr_s[mkx];
    double xc_arr_out[mkx];
    double xc_arr_s[mkx];
    double aquad_arr_out[mkx];
    double aquad_arr_s[mkx];
    double bquad_arr_out[mkx];
    double bquad_arr_s[mkx];
    double cquad_arr_out[mkx];
    double cquad_arr_s[mkx];
    double bogbot_arr_out[mkx];
    double bogbot_arr_s[mkx];
    double bogtop_arr_out[mkx];
    double bogtop_arr_s[mkx];
    double exns0[mkx + 1];
    double dp0[mkx];
    double write[mkx];
    double tre[mkx];
    double rei[mkx];
    double exn0[mkx];
    double xflx[mkx + 1];
    double uemf[mkx + 1];
    double comsub[mkx];
    double nlten_sink[mkx];
    double niten_sink[mkx];
    double uf[mkx];
    double vf[mkx];
    double qlten_sink[mkx];
    double qiten_sink[mkx];
    double flxrain[mkx + 1];
    double flxsnow[mkx];
    double ntraprd[mkx];
    double ntsnprd[mkx];
    double qw0_in[mkx];
    double qv0_star[mkx];
    double ql0_star[mkx];
    double qi0_star[mkx];
    double s0_star[mkx];
    double trflx_d[mkx];
    double trflx_u[mkx + 1];
    double dpdry0[mkx];
    double qv0_s[mkx];
    double ql0_s[mkx];
    double qi0_s[mkx];
    double s0_s[mkx];
    double u0_s[mkx];
    double v0_s[mkx];
    double qt0_s[mkx];
    double t0_s[mkx];
    double flxprc1_out[mkx + 1];
    double flxsnow1_out[mkx + 1];
    double tr0[ncnst][mkx];
    double tr0_o[ncnst][mkx];
    double sstr0[ncnst][mkx];
    double sstr0_o[ncnst][mkx];
    double trflx[ncnst][mkx];
    double trten[ncnst][mkx];
    double tru[ncnst][mkx];
    double tru_emf[ncnst][mkx];
    double trten_out[ncnst][mkx];
    double trten_s[ncnst][mkx];
    double trflx_out[ncnst][mkx];
    double trflx_s[ncnst][mkx];
    double tru_out[ncnst][mkx];
    double tru_s[ncnst][mkx];
    double tru_emf_out[ncnst][mkx];
    double tru_emf_s[ncnst][mkx];
    double tr0_s[ncnst][mkx];
    pe_get(param_d.zs0, zs0, sizeof(double)*(mkx + 1));
    pe_get(param_d.ps0, ps0, sizeof(double)*(mkx + 1));
    pe_get(param_d.p0, p0, sizeof(double)*mkx);
    pe_get(param_d.tke, tke, sizeof(double)*(mkx + 1));
    pe_get(param_d.thvl0bot, thvl0bot, sizeof(double)*mkx);
    pe_get(param_d.thvl0top, thvl0top, sizeof(double)*mkx);
    pe_get(param_d.qt0, qt0, sizeof(double)*mkx);
    pe_get(param_d.u0, u0, sizeof(double)*mkx);
    pe_get(param_d.ssu0, ssu0, sizeof(double)*mkx);
    pe_get(param_d.v0, v0, sizeof(double)*mkx);
    pe_get(param_d.ssv0, ssv0, sizeof(double)*mkx);
    pe_get(param_d.trsrc, trsrc, sizeof(double)*mkx);
    pe_get(param_d.thl0, thl0, sizeof(double)*mkx);
    pe_get(param_d.ssthl0, ssthl0, sizeof(double)*mkx);
    pe_get(param_d.ssqt0, ssqt0, sizeof(double)*mkx);
    pe_get(param_d.thv0bot, thv0bot, sizeof(double)*mkx);
    pe_get(param_d.thv0top, thv0top, sizeof(double)*mkx);
    pe_get(param_d.trsrc_o, trsrc_o, sizeof(double)*mkx);
    pe_get(param_d.qv0, qv0, sizeof(double)*mkx);
    pe_get(param_d.qv0_o, qv0_o, sizeof(double)*mkx);
    pe_get(param_d.ql0, ql0, sizeof(double)*mkx);
    pe_get(param_d.ql0_o, ql0_o, sizeof(double)*mkx);
    pe_get(param_d.qi0, qi0, sizeof(double)*mkx);
    pe_get(param_d.qi0_o, qi0_o, sizeof(double)*mkx);
    pe_get(param_d.t0, t0, sizeof(double)*mkx);
    pe_get(param_d.t0_o, t0_o, sizeof(double)*mkx);
    pe_get(param_d.s0, s0, sizeof(double)*mkx);
    pe_get(param_d.s0_o, s0_o, sizeof(double)*mkx);
    pe_get(param_d.u0_o, u0_o, sizeof(double)*mkx);
    pe_get(param_d.v0_o, v0_o, sizeof(double)*mkx);
    pe_get(param_d.qt0_o, qt0_o, sizeof(double)*mkx);
    pe_get(param_d.thl0_o, thl0_o, sizeof(double)*mkx);
    pe_get(param_d.thvl0, thvl0, sizeof(double)*mkx);
    pe_get(param_d.thvl0_o, thvl0_o, sizeof(double)*mkx);
    pe_get(param_d.ssthl0_o, ssthl0_o, sizeof(double)*mkx);
    pe_get(param_d.ssqt0_o, ssqt0_o, sizeof(double)*mkx);
    pe_get(param_d.thv0bot_o, thv0bot_o, sizeof(double)*mkx);
    pe_get(param_d.thv0top_o, thv0top_o, sizeof(double)*mkx);
    pe_get(param_d.thvl0bot_o, thvl0bot_o, sizeof(double)*mkx);
    pe_get(param_d.thvl0top_o, thvl0top_o, sizeof(double)*mkx);
    pe_get(param_d.ssu0_o, ssu0_o, sizeof(double)*mkx);
    pe_get(param_d.ssv0_o, ssv0_o, sizeof(double)*mkx);
    pe_get(param_d.umf, umf, sizeof(double)*(mkx + 1));
    pe_get(param_d.emf, emf, sizeof(double)*(mkx + 1));
    pe_get(param_d.slflx, slflx, sizeof(double)*(mkx + 1));
    pe_get(param_d.qtflx, qtflx, sizeof(double)*(mkx + 1));
    pe_get(param_d.uflx, uflx, sizeof(double)*(mkx + 1));
    pe_get(param_d.vflx, vflx, sizeof(double)*(mkx + 1));
    pe_get(param_d.qvten, qvten, sizeof(double)*mkx);
    pe_get(param_d.qlten, qlten, sizeof(double)*mkx);
    pe_get(param_d.qiten, qiten, sizeof(double)*mkx);
    pe_get(param_d.sten, sten, sizeof(double)*mkx);
    pe_get(param_d.uten, uten, sizeof(double)*mkx);
    pe_get(param_d.vten, vten, sizeof(double)*mkx);
    pe_get(param_d.qrten, qrten, sizeof(double)*mkx);
    pe_get(param_d.qsten, qsten, sizeof(double)*mkx);
    pe_get(param_d.dwten, dwten, sizeof(double)*mkx);
    pe_get(param_d.diten, diten, sizeof(double)*mkx);
    pe_get(param_d.evapc, evapc, sizeof(double)*mkx);
    pe_get(param_d.cufrc, cufrc, sizeof(double)*mkx);
    pe_get(param_d.qcu, qcu, sizeof(double)*mkx);
    pe_get(param_d.qlu, qlu, sizeof(double)*mkx);
    pe_get(param_d.qiu, qiu, sizeof(double)*mkx);
    pe_get(param_d.fer, fer, sizeof(double)*mkx);
    pe_get(param_d.fdr, fdr, sizeof(double)*mkx);
    pe_get(param_d.qc, qc, sizeof(double)*mkx);
    pe_get(param_d.qc_l, qc_l, sizeof(double)*mkx);
    pe_get(param_d.qc_i, qc_i, sizeof(double)*mkx);
    pe_get(param_d.qtten, qtten, sizeof(double)*mkx);
    pe_get(param_d.slten, slten, sizeof(double)*mkx);
    pe_get(param_d.ufrc, ufrc, sizeof(double)*(mkx + 1));
    pe_get(param_d.thlu, thlu, sizeof(double)*(mkx + 1));
    pe_get(param_d.qtu, qtu, sizeof(double)*(mkx + 1));
    pe_get(param_d.uu, uu, sizeof(double)*(mkx + 1));
    pe_get(param_d.vu, vu, sizeof(double)*(mkx + 1));
    pe_get(param_d.wu, wu, sizeof(double)*(mkx + 1));
    pe_get(param_d.thvu, thvu, sizeof(double)*(mkx + 1));
    pe_get(param_d.thlu_emf, thlu_emf, sizeof(double)*(mkx + 1));
    pe_get(param_d.qtu_emf, qtu_emf, sizeof(double)*(mkx + 1));
    pe_get(param_d.uu_emf, uu_emf, sizeof(double)*(mkx + 1));
    pe_get(param_d.vu_emf, vu_emf, sizeof(double)*(mkx + 1));
    pe_get(param_d.excessu_arr, excessu_arr, sizeof(double)*mkx);
    pe_get(param_d.excess0_arr, excess0_arr, sizeof(double)*mkx);
    pe_get(param_d.xc_arr, xc_arr, sizeof(double)*mkx);
    pe_get(param_d.aquad_arr, aquad_arr, sizeof(double)*mkx);
    pe_get(param_d.bquad_arr, bquad_arr, sizeof(double)*mkx);
    pe_get(param_d.cquad_arr, cquad_arr, sizeof(double)*mkx);
    pe_get(param_d.bogbot_arr, bogbot_arr, sizeof(double)*mkx);
    pe_get(param_d.bogtop_arr, bogtop_arr, sizeof(double)*mkx);
    pe_get(param_d.umf_out, umf_out, sizeof(double)*(mkx + 1));
    pe_get(param_d.umf_s, umf_s, sizeof(double)*mkx);
    pe_get(param_d.qvten_out, qvten_out, sizeof(double)*mkx);
    pe_get(param_d.qvten_s, qvten_s, sizeof(double)*mkx);
    pe_get(param_d.qlten_out, qlten_out, sizeof(double)*mkx);
    pe_get(param_d.qlten_s, qlten_s, sizeof(double)*mkx);
    pe_get(param_d.qiten_out, qiten_out, sizeof(double)*mkx);
    pe_get(param_d.qiten_s, qiten_s, sizeof(double)*mkx);
    pe_get(param_d.sten_out, sten_out, sizeof(double)*mkx);
    pe_get(param_d.sten_s, sten_s, sizeof(double)*mkx);
    pe_get(param_d.uten_out, uten_out, sizeof(double)*mkx);
    pe_get(param_d.uten_s, uten_s, sizeof(double)*mkx);
    pe_get(param_d.vten_out, vten_out, sizeof(double)*mkx);
    pe_get(param_d.vten_s, vten_s, sizeof(double)*mkx);
    pe_get(param_d.qrten_out, qrten_out, sizeof(double)*mkx);
    pe_get(param_d.qrten_s, qrten_s, sizeof(double)*mkx);
    pe_get(param_d.qsten_out, qsten_out, sizeof(double)*mkx);
    pe_get(param_d.qsten_s, qsten_s, sizeof(double)*mkx);
    pe_get(param_d.evapc_out, evapc_out, sizeof(double)*mkx);
    pe_get(param_d.evapc_s, evapc_s, sizeof(double)*mkx);
    pe_get(param_d.cufrc_out, cufrc_out, sizeof(double)*mkx);
    pe_get(param_d.cufrc_s, cufrc_s, sizeof(double)*mkx);
    pe_get(param_d.slflx_out, slflx_out, sizeof(double)*(mkx + 1));
    pe_get(param_d.slflx_s, slflx_s, sizeof(double)*mkx);
    pe_get(param_d.qtflx_out, qtflx_out, sizeof(double)*(mkx + 1));
    pe_get(param_d.qtflx_s, qtflx_s, sizeof(double)*mkx);
    pe_get(param_d.qcu_out, qcu_out, sizeof(double)*mkx);
    pe_get(param_d.qcu_s, qcu_s, sizeof(double)*mkx);
    pe_get(param_d.qlu_out, qlu_out, sizeof(double)*mkx);
    pe_get(param_d.qlu_s, qlu_s, sizeof(double)*mkx);
    pe_get(param_d.qiu_out, qiu_out, sizeof(double)*mkx);
    pe_get(param_d.qiu_s, qiu_s, sizeof(double)*mkx);
    pe_get(param_d.qc_out, qc_out, sizeof(double)*mkx);
    pe_get(param_d.qc_s, qc_s, sizeof(double)*mkx);
    pe_get(param_d.fer_out, fer_out, sizeof(double)*mkx);
    pe_get(param_d.fer_s, fer_s, sizeof(double)*mkx);
    pe_get(param_d.fdr_out, fdr_out, sizeof(double)*mkx);
    pe_get(param_d.fdr_s, fdr_s, sizeof(double)*mkx);
    pe_get(param_d.qtten_out, qtten_out, sizeof(double)*mkx);
    pe_get(param_d.qtten_s, qtten_s, sizeof(double)*mkx);
    pe_get(param_d.slten_out, slten_out, sizeof(double)*mkx);
    pe_get(param_d.slten_s, slten_s, sizeof(double)*mkx);
    pe_get(param_d.ufrc_out, ufrc_out, sizeof(double)*(mkx + 1));
    pe_get(param_d.ufrc_s, ufrc_s, sizeof(double)*mkx);
    pe_get(param_d.uflx_out, uflx_out, sizeof(double)*(mkx + 1));
    pe_get(param_d.uflx_s, uflx_s, sizeof(double)*mkx);
    pe_get(param_d.vflx_out, vflx_out, sizeof(double)*(mkx + 1));
    pe_get(param_d.vflx_s, vflx_s, sizeof(double)*mkx);
    pe_get(param_d.wu_out, wu_out, sizeof(double)*(mkx + 1));
    pe_get(param_d.wu_s, wu_s, sizeof(double)*mkx);
    pe_get(param_d.qtu_out, qtu_out, sizeof(double)*(mkx + 1));
    pe_get(param_d.qtu_s, qtu_s, sizeof(double)*mkx);
    pe_get(param_d.thlu_out, thlu_out, sizeof(double)*(mkx + 1));
    pe_get(param_d.thlu_s, thlu_s, sizeof(double)*mkx);
    pe_get(param_d.thvu_out, thvu_out, sizeof(double)*mkx);
    pe_get(param_d.thvu_s, thvu_s, sizeof(double)*mkx);
    pe_get(param_d.uu_out, uu_out, sizeof(double)*mkx);
    pe_get(param_d.uu_s, uu_s, sizeof(double)*mkx);
    pe_get(param_d.vu_out, vu_out, sizeof(double)*mkx);
    pe_get(param_d.vu_s, vu_s, sizeof(double)*mkx);
    pe_get(param_d.qtu_emf_out, qtu_emf_out, sizeof(double)*mkx);
    pe_get(param_d.qtu_emf_s, qtu_emf_s, sizeof(double)*mkx);
    pe_get(param_d.thlu_emf_out, thlu_emf_out, sizeof(double)*mkx);
    pe_get(param_d.thlu_emf_s, thlu_emf_s, sizeof(double)*mkx);
    pe_get(param_d.uu_emf_out, uu_emf_out, sizeof(double)*mkx);
    pe_get(param_d.uu_emf_s, uu_emf_s, sizeof(double)*mkx);
    pe_get(param_d.vu_emf_out, vu_emf_out, sizeof(double)*mkx);
    pe_get(param_d.vu_emf_s, vu_emf_s, sizeof(double)*mkx);
    pe_get(param_d.uemf_out, uemf_out, sizeof(double)*mkx);
    pe_get(param_d.uemf_s, uemf_s, sizeof(double)*mkx);
    pe_get(param_d.dwten_out, dwten_out, sizeof(double)*mkx);
    pe_get(param_d.dwten_s, dwten_s, sizeof(double)*mkx);
    pe_get(param_d.diten_out, diten_out, sizeof(double)*mkx);
    pe_get(param_d.diten_s, diten_s, sizeof(double)*mkx);
    pe_get(param_d.flxrain_out, flxrain_out, sizeof(double)*mkx);
    pe_get(param_d.flxrain_s, flxrain_s, sizeof(double)*mkx);
    pe_get(param_d.flxsnow_out, flxsnow_out, sizeof(double)*mkx);
    pe_get(param_d.flxsnow_s, flxsnow_s, sizeof(double)*mkx);
    pe_get(param_d.ntraprd_out, ntraprd_out, sizeof(double)*mkx);
    pe_get(param_d.ntraprd_s, ntraprd_s, sizeof(double)*mkx);
    pe_get(param_d.ntsnprd_out, ntsnprd_out, sizeof(double)*mkx);
    pe_get(param_d.ntsnprd_s, ntsnprd_s, sizeof(double)*mkx);
    pe_get(param_d.excessu_arr_out, excessu_arr_out, sizeof(double)*mkx);
    pe_get(param_d.excessu_arr_s, excessu_arr_s, sizeof(double)*mkx);
    pe_get(param_d.excess0_arr_out, excess0_arr_out, sizeof(double)*mkx);
    pe_get(param_d.excess0_arr_s, excess0_arr_s, sizeof(double)*mkx);
    pe_get(param_d.xc_arr_out, xc_arr_out, sizeof(double)*mkx);
    pe_get(param_d.xc_arr_s, xc_arr_s, sizeof(double)*mkx);
    pe_get(param_d.aquad_arr_out, aquad_arr_out, sizeof(double)*mkx);
    pe_get(param_d.aquad_arr_s, aquad_arr_s, sizeof(double)*mkx);
    pe_get(param_d.bquad_arr_out, bquad_arr_out, sizeof(double)*mkx);
    pe_get(param_d.bquad_arr_s, bquad_arr_s, sizeof(double)*mkx);
    pe_get(param_d.cquad_arr_out, cquad_arr_out, sizeof(double)*mkx);
    pe_get(param_d.cquad_arr_s, cquad_arr_s, sizeof(double)*mkx);
    pe_get(param_d.bogbot_arr_out, bogbot_arr_out, sizeof(double)*mkx);
    pe_get(param_d.bogbot_arr_s, bogbot_arr_s, sizeof(double)*mkx);
    pe_get(param_d.bogtop_arr_out, bogtop_arr_out, sizeof(double)*mkx);
    pe_get(param_d.bogtop_arr_s, bogtop_arr_s, sizeof(double)*mkx);
    pe_get(param_d.exns0, exns0, sizeof(double)*(mkx + 1));
    pe_get(param_d.dp0, dp0, sizeof(double)*mkx);
    pe_get(param_d.write, write, sizeof(double)*mkx);
    pe_get(param_d.tre, tre, sizeof(double)*mkx);
    pe_get(param_d.rei, rei, sizeof(double)*mkx);
    pe_get(param_d.exn0, exn0, sizeof(double)*mkx);
    pe_get(param_d.xflx, xflx, sizeof(double)*(mkx + 1));
    pe_get(param_d.uemf, uemf, sizeof(double)*(mkx + 1));
    pe_get(param_d.comsub, comsub, sizeof(double)*mkx);
    pe_get(param_d.nlten_sink, nlten_sink, sizeof(double)*mkx);
    pe_get(param_d.niten_sink, niten_sink, sizeof(double)*mkx);
    pe_get(param_d.uf, uf, sizeof(double)*mkx);
    pe_get(param_d.vf, vf, sizeof(double)*mkx);
    pe_get(param_d.qlten_sink, qlten_sink, sizeof(double)*mkx);
    pe_get(param_d.qiten_sink, qiten_sink, sizeof(double)*mkx);
    pe_get(param_d.flxrain, flxrain, sizeof(double)*(mkx + 1));
    pe_get(param_d.flxsnow, flxsnow, sizeof(double)*mkx);
    pe_get(param_d.ntraprd, ntraprd, sizeof(double)*mkx);
    pe_get(param_d.ntsnprd, ntsnprd, sizeof(double)*mkx);
    pe_get(param_d.qw0_in, qw0_in, sizeof(double)*mkx);
    pe_get(param_d.qv0_star, qv0_star, sizeof(double)*mkx);
    pe_get(param_d.ql0_star, ql0_star, sizeof(double)*mkx);
    pe_get(param_d.qi0_star, qi0_star, sizeof(double)*mkx);
    pe_get(param_d.s0_star, s0_star, sizeof(double)*mkx);
    pe_get(param_d.trflx_d, trflx_d, sizeof(double)*mkx);
    pe_get(param_d.trflx_u, trflx_u, sizeof(double)*(mkx + 1));
    pe_get(param_d.dpdry0, dpdry0, sizeof(double)*mkx);
    pe_get(param_d.qv0_s, qv0_s, sizeof(double)*mkx);
    pe_get(param_d.ql0_s, ql0_s, sizeof(double)*mkx);
    pe_get(param_d.qi0_s, qi0_s, sizeof(double)*mkx);
    pe_get(param_d.s0_s, s0_s, sizeof(double)*mkx);
    pe_get(param_d.u0_s, u0_s, sizeof(double)*mkx);
    pe_get(param_d.v0_s, v0_s, sizeof(double)*mkx);
    pe_get(param_d.qt0_s, qt0_s, sizeof(double)*mkx);
    pe_get(param_d.t0_s, t0_s, sizeof(double)*mkx);
    pe_get(param_d.flxprc1_out, flxprc1_out, sizeof(double)*(mkx + 1));
    pe_get(param_d.flxsnow1_out, flxsnow1_out, sizeof(double)*(mkx + 1));
    pe_get(param_d.tr0, tr0, sizeof(double)*mkx*ncnst);
    pe_get(param_d.tr0_o, tr0_o, sizeof(double)*mkx*ncnst);
    pe_get(param_d.sstr0, sstr0, sizeof(double)*mkx*ncnst);
    pe_get(param_d.sstr0_o, sstr0_o, sizeof(double)*mkx*ncnst);
    pe_get(param_d.trflx, trflx, sizeof(double)*mkx*ncnst);
    pe_get(param_d.trten, trten, sizeof(double)*mkx*ncnst);
    pe_get(param_d.tru, tru, sizeof(double)*mkx*ncnst);
    pe_get(param_d.tru_emf, tru_emf, sizeof(double)*mkx*ncnst);
    pe_get(param_d.trten_out, trten_out, sizeof(double)*mkx*ncnst);
    pe_get(param_d.trten_s, trten_s, sizeof(double)*mkx*ncnst);
    pe_get(param_d.trflx_out, trflx_out, sizeof(double)*mkx*ncnst);
    pe_get(param_d.trflx_s, trflx_s, sizeof(double)*mkx*ncnst);
    pe_get(param_d.tru_out, tru_out, sizeof(double)*mkx*ncnst);
    pe_get(param_d.tru_s, tru_s, sizeof(double)*mkx*ncnst);
    pe_get(param_d.tru_emf_out, tru_emf_out, sizeof(double)*mkx*ncnst);
    pe_get(param_d.tru_emf_s, tru_emf_s, sizeof(double)*mkx*ncnst);
    pe_get(param_d.tr0_s, tr0_s, sizeof(double)*mkx*ncnst);
    dma_syn();
    uw_compute_loop1_(&tscaleh, &cush, &pblh, &kinv, &exit_kinv1, &id_exit, &dpsum, &tkeavg, &thvlmin, &dpi, &qtsrc, &thvlsrc, &thlsrc, &usrc, &vsrc, &plcl, &klcl, &exit_klclmkx, &thl0lcl, &qt0lcl, &thj, &qvj, &qlj, &qij, &qse, &exit_conden, &thv0lcl, &cin, &cinlcl, &plfc, &klfc, &thvubot, &thvutop, &limit_cinlcl, &limit_cin, &exit_klfcmkx, &cin_i, &cinlcl_i, &ke, &kinv_o, &klcl_o, &klfc_o, &plcl_o, &plfc_o, &tkeavg_o, &thvlmin_o, &qtsrc_o, &thvlsrc_o, &thlsrc_o, &usrc_o, &vsrc_o, &thv0lcl_o, &cin_f, &cinlcl_f, &del_CIN, &alpha, &del_cinlcl, &precip, &snow, &rliq, &cbmf, &cnt, &cnb, &ufrcinvbase, &ufrclcl, &winvbase, &wlcl, &emfkbup, &cbmflimit, &ind_delcin, &precip_out, &precip_s, &snow_out, &snow_s, &cush_inout, &cush_s, &cbmf_out, &cbmf_s, &rliq_out, &rliq_s, &cnt_out, &cnt_s, &cnb_out, &cnb_s, &cinh_out, &cin_s, &cinlclh_out, &cinlcl_s, &ufrcinvbase_out, &ufrcinvbase_s, &ufrclcl_out, &ufrclcl_s, &winvbase_out, &winvbase_s, &wlcl_out, &wlcl_s, &plcl_out, &plcl_s, &pinv_out, &pinv_s, &plfc_out, &plfc_s, &pbup_out, &pbup_s, &ppen_out, &ppen_s, &qtsrc_out, &qtsrc_s, &thlsrc_out, &thlsrc_s, &thvlsrc_out, &thvlsrc_s, &emfkbup_out, &emfkbup_s, &cbmflimit_out, &cbmflimit_s, &tkeavg_out, &tkeavg_s, &zinv_out, &zinv_s, &rcwp_out, &rcwp_s, &rlwp_out, &rlwp_s, &riwp_out, &riwp_s, &id_exit_ex, &krel, &prel, &thv0rel, &wcrit, &sigmaw, &mu, &rho0inv, &mumin0, &mulcl, &mulclstar, &erfc, &mumin2, &limit_ufrc, &limit_cbmf, &winv, &ufrcinv, &wtw, &exit_wtw, &wrel, &exit_ufrc, &uplus, &vplus, &pe, &dpe, &exne, &thvebot, &thle, &qte, &ue, &ve, &scaleh, &kbup, &thlue, &qtue, &wue, &wtwb, &thv0j, &rho0j, &qsat_arg, &es, &weight, &estbl, &qs, &excess0, &exql, &exqi, &thvj, &excessu, &xc, &cridis, &aquad, &bquad, &cquad, &xsat, &thlxsat, &qtxsat, &thvxsat, &thv_x0, &thv_x1, &xs1, &xs2, &x_cu, &x_en, &ee2, &ud2, &delbog, &bogtop, &bogbot, &drage, &expfac, &autodet, &exit_wu, &rhos0j, &xc1, &xc2, &ppen, &limit_ppen, &thlu_top, &qtu_top, &exntop, &forcedCu, &limit_shcu, &cldhgt, &exit_cufilter, &limit_emf, &xsrc, &xmean, &xtop, &xbot, &thlten_sub, &qtten_sub, &qlten_sub, &qiten_sub, &nlten_sub, &niten_sub, &thl_prog, &qt_prog, &rainflx, &snowflx, &qlu_mid, &qiu_mid, &qlubelow, &qiubelow, &qlu_top, &qiu_top, &qc_lm, &qc_im, &nc_lm, &nc_im, &ql_emf_kbup, &qi_emf_kbup, &nl_emf_kbup, &ni_emf_kbup, &qlten_det, &qiten_det, &evpint_rain, &evpint_snow, &snowmlt, &subsat, &evprain, &evpsnow, &evplimit, &evplimit_rain, &evplimit_snow, &tmp1, &tmp2, &limit_negcon, &trmin, &pdelx, &dum, &pdel, &qcubelow, &rcwp, &rlwp, &riwp, &thl0bot, &qt0bot, &thl0top, &qt0top, zs0, ps0, p0, tke, thvl0bot, thvl0top, qt0, u0, ssu0, v0, ssv0, trsrc, thl0, ssthl0, ssqt0, thv0bot, thv0top, trsrc_o, qv0, qv0_o, ql0, ql0_o, qi0, qi0_o, t0, t0_o, s0, s0_o, u0_o, v0_o, qt0_o, thl0_o, thvl0, thvl0_o, ssthl0_o, ssqt0_o, thv0bot_o, thv0top_o, thvl0bot_o, thvl0top_o, ssu0_o, ssv0_o, umf, emf, slflx, qtflx, uflx, vflx, qvten, qlten, qiten, sten, uten, vten, qrten, qsten, dwten, diten, evapc, cufrc, qcu, qlu, qiu, fer, fdr, qc, qc_l, qc_i, qtten, slten, ufrc, thlu, qtu, uu, vu, wu, thvu, thlu_emf, qtu_emf, uu_emf, vu_emf, excessu_arr, excess0_arr, xc_arr, aquad_arr, bquad_arr, cquad_arr, bogbot_arr, bogtop_arr, umf_out, umf_s, qvten_out, qvten_s, qlten_out, qlten_s, qiten_out, qiten_s, sten_out, sten_s, uten_out, uten_s, vten_out, vten_s, qrten_out, qrten_s, qsten_out, qsten_s, evapc_out, evapc_s, cufrc_out, cufrc_s, slflx_out, slflx_s, qtflx_out, qtflx_s, qcu_out, qcu_s, qlu_out, qlu_s, qiu_out, qiu_s, qc_out, qc_s, fer_out, fer_s, fdr_out, fdr_s, qtten_out, qtten_s, slten_out, slten_s, ufrc_out, ufrc_s, uflx_out, uflx_s, vflx_out, vflx_s, wu_out, wu_s, qtu_out, qtu_s, thlu_out, thlu_s, thvu_out, thvu_s, uu_out, uu_s, vu_out, vu_s, qtu_emf_out, qtu_emf_s, thlu_emf_out, thlu_emf_s, uu_emf_out, uu_emf_s, vu_emf_out, vu_emf_s, uemf_out, uemf_s, dwten_out, dwten_s, diten_out, diten_s, flxrain_out, flxrain_s, flxsnow_out, flxsnow_s, ntraprd_out, ntraprd_s, ntsnprd_out, ntsnprd_s, excessu_arr_out, excessu_arr_s, excess0_arr_out, excess0_arr_s, xc_arr_out, xc_arr_s, aquad_arr_out, aquad_arr_s, bquad_arr_out, bquad_arr_s, cquad_arr_out, cquad_arr_s, bogbot_arr_out, bogbot_arr_s, bogtop_arr_out, bogtop_arr_s, exns0, dp0, write, tre, rei, exn0, xflx, uemf, comsub, nlten_sink, niten_sink, uf, vf, qlten_sink, qiten_sink, flxrain, flxsnow, ntraprd, ntsnprd, qw0_in, qv0_star, ql0_star, qi0_star, s0_star, trflx_d, trflx_u, dpdry0, qv0_s, ql0_s, qi0_s, s0_s, u0_s, v0_s, qt0_s, t0_s, flxprc1_out, flxsnow1_out, tr0, tr0_o, sstr0, sstr0_o, trflx, trten, tru, tru_emf, trten_out, trten_s, trflx_out, trflx_s, tru_out, tru_s, tru_emf_out, tru_emf_s, tr0_s);

  }
}
