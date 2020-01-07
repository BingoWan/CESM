# Script: For Variables Substitution & Loop Exchange for uwshcu.F90 of CESM
# Author: Wan Wubing & Xie Jiayu

import os
import sys
import re
import random

# 数组声明替换
def array_expand00(m):
	return m.group(0).replace(")", ", mix)")

# 单变量声明替换
def array_expand01(m):
	name  = re.sub("\W", "", m.group(0)[1:])
	other = m.group(0)[1:].replace(name, "")
	return m.group(0)[0] + name + "(mix)" + other

# 数组引用替换
def array_expand10(m):
	return m.group(0).replace(")", ", i)")

# 数组参数替换
def array_expand11(m):
	name  = re.sub("\W", "", m.group(0)[1:])
	other = m.group(0)[1:].replace(name, "")
	return m.group(0)[0] + name + "(:, i)" + other

# 单变量使用替换
def array_expand12(m):
	name  = re.sub("\W", "", m.group(0)[1:])
	other = m.group(0)[1:].replace(name, "")
	return m.group(0)[0] + name + "(i)" + other

# _o数组声明替换
def array_expand20(m):
	return m.group(0).replace(")", ", mix)")

# _o变量声明替换
def array_expand21(m):
	#str = m.group(0) + '~'
	name  = re.sub("[ :]", "", m.group(0)[1:])
	other = m.group(0)[1:].replace(name, "")
	return m.group(0)[0] + name + ", dimension(mix)" + other

# 交换循环i和iter次序
def loop_exchange1(m):
	return "! CODE OF goto333 MUST BE COPIED & RENUMBERED HERE !!!\n\n    " \
			+ "enddo ! init i loop\n\n    " + m.group(0) \
			+ "\n\n    do i = 1, iend\n\n       " \
			+ "if (id_exit) then\n          cycle\n       endif"

def loop_exchange20(m):
	return "! " + m.group(0) + "\n\n     if (iter == iter_cin) then"

def loop_exchange21(m):
	return "     endif ! (iter == iter_cin)\n\n" + m.group(0)

def loop_exchange22(m):
	return m.group(0) + "\n\n     end do                ! End of implicit CIN loop (iter_cin)"


# 字符串匹配策略
def write_code(value_list, nrplc_list, aflag_list):
    flag1 = False
    flag2 = False
    flag3 = False
    with open("uwshcu_bac.F90", "w") as fw:
        with open("uwshcu_o.F90", "r") as fr:
            files = fr.readlines()
            for line in files:
                #交换循环顺序
                line = re.sub("do iter = 1, iter_cin", loop_exchange1, line)
                line = re.sub("end do                ! End of implicit CIN loop \(cin_iter\)", \
                              loop_exchange20, line)
                line = re.sub(" 333 if\(id_exit\) then ! Exit without cumulus convection", \
                              loop_exchange21, line)
                line = re.sub("end do                  ! end of big i loop for each column.", \
                              loop_exchange22, line)

                #变量扩展
                if re.match("!!! start of uwshcu", line) != None:
                    flag1 = True
                if re.match("!!! end of uwshcu", line) != None:
                    flag1 = False
                if re.match("!!! start of initialization", line) != None:
                    flag2 = True
                if re.match("!!! end of initialization", line) != None:
                    flag2 = False
                if re.match("!!! Start of _s&_o variables", line) != None:
                    flag3 = True
                if re.match("!!! End of _s&_o variables", line) != None:
                    flag3 = False
                if (flag1 == True and flag2 == True):
                    for index in range(len(value_list)):
                        value = value_list[index]
                        nrplc_list_temp = nrplc_list[index] + len(re.compile("\W%s\(.+?\)"%value).findall(line))
                        if (nrplc_list_temp > nrplc_list[index]):
                            aflag_list[index] = 1
                        nrplc_list[index] = nrplc_list_temp
                        line  = re.sub("\W%s\(.+?\)"%value, array_expand00, line)
                if (flag1 == True and flag2 == True and flag3 == False):
                    for index in range(len(value_list)):
                        value = value_list[index]
                        if (aflag_list[index] == 0):
                            nrplc_list[index] = nrplc_list[index] + len(re.compile("\W%s[^A-Za-z0-9_\(]"%value).findall(line))
                            line  = re.sub("\W%s[^A-Za-z0-9_\(]"%value, array_expand01, line)
                if (flag1 == True and flag2 == False):
                    for index in range(len(value_list)):
                        value = value_list[index]
                        nrplc_list[index] = nrplc_list[index] + len(re.compile("\W%s\(.+?\)"%value).findall(line))
                        line  = re.sub("\W%s\(.+?\)"%value, array_expand10, line)
                if (flag1 == True and flag2 == False):
                    for index in range(len(value_list)):
                        value = value_list[index]
                        if (aflag_list[index] == 1):
                            nrplc_list[index] = nrplc_list[index] + len(re.compile("\W%s[^A-Za-z0-9_\(]"%value).findall(line))
                            line  = re.sub("\W%s[^A-Za-z0-9_\(]"%value, array_expand11, line)
                if (flag1 == True and flag2 == False):
                    for index in range(len(value_list)):
                        value = value_list[index]
                        if (aflag_list[index] == 0):
                            nrplc_list[index] = nrplc_list[index] + len(re.compile("\W%s[^A-Za-z0-9_\(]"%value).findall(line))
                            line  = re.sub("\W%s[^A-Za-z0-9_\(]"%value, array_expand12, line)
                if (flag3 == True):
                    line  = re.sub("\Wdimension\([^i]+?\)", array_expand20, line)
                    line  = re.sub("\Wreal\(r8\) +::", array_expand21, line)
                    line  = re.sub("\Winteger +::", array_expand21, line)
                fw.write(line)
    for index in range(len(value_list)):
        print("%15s"%value_list[index], " 被替换 ", "%4d"%nrplc_list[index], " 次")
    print("--->finish")

# 主程序
def main(argv):
    value_list = [
				  # Input variables (All Arrays)
				  "ps0", "zs0", "p0", "z0", "dp0", "dpdry0", "u0", "v0", "tke", "cldfrct", "concldfrct", \
				  "qv0", "ql0", "qi0", "t0", "s0", "pblh", "cush", "tr0", \

				  # Environmental variables directly derived from the input variables (All Arrays)
				  "qt0", "thl0", "thvl0", "ssqt0", "ssthl0", "ssu0", "ssv0", "thv0bot", "thv0top", "thvl0bot", "thvl0top", "exn0", "exns0", "sstr0", \

				  # Variables associated with cumulus convection (Mainly Arrays)
				  "qv0_star", "ql0_star", "qi0_star", "t0_star", "s0_star", "umf", "emf", "qvten", "qlten", \
				  "qiten", "sten", "uten", "vten", "qrten", "qsten", "precip", "snow", "evapc", "slflx", \
				  "qtflx", "uflx", "vflx", "cufrc", "qcu", "qlu", "qiu", "dwten", "diten", "fer", "fdr", "uf", \
				  "vf", "qc", "qc_l", "qc_i", "qc_lm", "qc_im", "nc_lm", "nc_im", "ql_emf_kbup", \
				  "qi_emf_kbup", "nl_emf_kbup", "ni_emf_kbup", "qlten_det", "qiten_det", "rliq", "cnt", "cnb", \
				  "qtten", "slten", "ufrc", "trten", "trflx", "trflx_d", "trflx_u", "trmin", "pdelx", "dum", \

				  # Variables used for the calculation of condensation sink (Arrays Except _sub & _prog)
				  "uemf", "comsub", "qlten_sink", "qiten_sink", "nlten_sink", "niten_sink", "thlten_sub", \
				  "qtten_sub", "qlten_sub", "qiten_sub", "nlten_sub", "niten_sub", "thl_prog", "qt_prog", \

				  # Variables describing cumulus updraft (All Arrays)
				  "wu", "thlu", "qtu", "uu", "vu", "thvu", "rei", "tru", \

				  # Variables describing conservative scalars of entraining downdrafts (All Arrays)
				  "thlu_emf", "qtu_emf", "uu_emf", "vu_emf", "tru_emf", \

				  # Variables associated with evaporations of convective 'rain' and 'snow'
				  # (flxrain, flxsnow, ntraprd and ntsnprd are Arrays)
				  "flxrain", "flxsnow", "ntraprd", "ntsnprd", "flxsntm", "snowmlt", "subsat", "evprain", \
				  "evpsnow", "evplimit", "evplimit_rain", "evplimit_snow", "evpint_rain", "evpint_snow", \

				  # Other internal variables (All Single Vars Except trsrc)
				  "id_exit", "forcedCu", "klcl", "kinv", "krel", "klfc", "kbup", "kpen", "thlsrc", "qtsrc", \
				  "usrc", "vsrc", "thvlsrc", "uplus", "vplus", "trsrc", "tre", "plcl", "plfc", "prel", "wrel", \
				  "ee2", "ud2", "wtw", "wtwb", "wtwh", "xc", "xc_2", "cldhgt", "scaleh", "tscaleh", "cridis", \
				  "sigmaw", "tkeavg", "dpsum", "dpi", "thvlmin", "thlxsat", "qtxsat", "thvxsat", "x_cu", \
				  "x_en", "thv_x0", "thv_x1", "thj", "qvj", "qlj", "qij", "thvj", "tj", "thv0j", "rho0j", \
				  "rhos0j", "qse", "cin", "cinlcl", "pe", "dpe", "exne", "thvebot", "thle", "qte", "ue", "ve", \
				  "thlue", "qtue", "wue", "mu", "mumin0", "mumin2", "mulcl", "mulclstar", "cbmf", "wcrit", \
				  "winv", "wlcl", "ufrcinv", "ufrclcl", "exql", "exqi", "ppen", "thl0top", "thl0bot", \
				  "qt0bot", "qt0top", "thvubot", "thvutop", "thlu_top", "qtu_top", "qlu_top", "qiu_top", \
				  "qlu_mid", "qiu_mid", "exntop", "thl0lcl", "qt0lcl", "thv0lcl", "thv0rel", "rho0inv", \
				  "autodet", "aquad", "bquad", "cquad", "xc1", "xc2", "excessu", "excess0", "xsat", "xs1", \
				  "xs2", "bogbot", "bogtop", "delbog", "drage", "expfac", "rcwp", "rlwp", "riwp", "qcubelow", \
				  "qlubelow", "qiubelow", "rainflx", "snowflx", "es", "qs", "qsat_arg", "xsrc", "xmean", \
				  "xtop", "xbot", "xflx", "tmp1", "tmp2", \

				  # Variables for implicit CIN computation (Mainly Arrays)
				  "ufrcinvbase_s", "ufrclcl_s", "winvbase_s", "wlcl_s", "plcl_s", "pinv_s", "plfc_s", \
				  "qtsrc_s", "thlsrc_s", "thvlsrc_s", "emfkbup_s", "cinlcl_s", "pbup_s", "ppen_s", \
				  "cbmflimit_s", "tkeavg_s", "zinv_s", "rcwp_s", "rlwp_s", "riwp_s", "ufrcinvbase", \
				  "winvbase", "pinv", "zinv", "emfkbup", "cbmflimit", "rho0rel", "qv0_s", "ql0_s", "qi0_s", \
				  "s0_s", "u0_s", "v0_s", "t0_s", "qt0_s", "thl0_s", "thvl0_s", "qvten_s", "qlten_s", \
				  "qiten_s", "qrten_s", "qsten_s", "sten_s", "evapc_s", "uten_s", "vten_s", "cufrc_s", \
				  "qcu_s", "qlu_s", "qiu_s", "fer_s", "fdr_s", "qc_s", "qtten_s", "slten_s", "umf_s", \
				  "slflx_s", "qtflx_s", "ufrc_s", "uflx_s", "vflx_s", "cush_s", "precip_s", "snow_s", \
				  "cin_s", "rliq_s", "cbmf_s", "cnt_s", "cnb_s", "tr0_s", "trten_s", "trflx_s", \
				  "cin_i", "cin_f", "del_CIN", "ke", "alpha", "thlj", "cinlcl_i", "cinlcl_f", "del_cinlcl", \
				  "wu_s", "qtu_s", "thlu_s", "thvu_s", "uu_s", "vu_s", "qtu_emf_s", "thlu_emf_s", "uu_emf_s", \
				  "vu_emf_s", "uemf_s", "tru_s", "tru_emf_s", "dwten_s", "diten_s", "flxrain_s", "flxsnow_s", \
				  "ntraprd_s", "ntsnprd_s", "excessu_arr", "excessu_arr_s", "excess0_arr", "excess0_arr_s", \
				  "xc_arr", "xc_arr_s", "aquad_arr", "aquad_arr_s", "bquad_arr", "bquad_arr_s", "cquad_arr", \
				  "cquad_arr_s", "bogbot_arr", "bogbot_arr_s", "bogtop_arr", "bogtop_arr_s", \

				  # Variables for temporary storages (Mainly Arrays)
				  "qv0_o", "ql0_o", "qi0_o", "t0_o", "s0_o", "u0_o", "v0_o", "qt0_o", "thl0_o", "thvl0_o", \
				  "qvten_o", "qlten_o", "qiten_o", "qrten_o", "qsten_o", "sten_o", "uten_o", "vten_o", \
				  "qcu_o", "qlu_o", "qiu_o", "cufrc_o", "evapc_o", "thv0bot_o", "thv0top_o", "thvl0bot_o", \
				  "thvl0top_o", "ssthl0_o", "ssqt0_o", "ssu0_o", "ssv0_o", "qc_o", "qtten_o", "slten_o", \
				  "umf_o", "slflx_o", "qtflx_o", "ufrc_o", "uflx_o", "vflx_o", "tkeavg_o", "thvlmin_o", \
				  "qtsrc_o", "thvlsrc_o", "thlsrc_o", "usrc_o", "vsrc_o", "plcl_o", "plfc_o", "thv0lcl_o", \
				  "cinlcl_o", "kinv_o", "klcl_o", "klfc_o", "tr0_o", "trten_o", "sstr0_o", "trflx_o", \
				  "trsrc_o", \
				  ]
    nrplc_list = [0] * len(value_list)
    aflag_list = [0] * len(value_list)
    write_code(value_list, nrplc_list, aflag_list)

if __name__ == "__main__":
   sys.exit(main(sys.argv))
