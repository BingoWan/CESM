#! /home/ruo/anaconda3/bin/python
# Script: For Variables Substitution & Loop Exchange for uwshcu.F90 of CESM
# Author: Wan Wubing & Xie Jiayu

import os
import sys
import re
import random


flag1 = False
flag2 = False
flag3 = False
flag4 = False
flag5 = False
flag6 = False
flag7 = False
flag8 = False

value_list = [
    # Input variables (All Arrays)
    "ps0", "zs0", "p0", "z0", "dp0", "dpdry0", "u0", "v0", "tke", "cldfrct", "concldfrct", \
    "qv0", "ql0", "qi0", "t0", "s0", "pblh", "cush", "tr0", "t_tmp", "weight",\

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
    "id_exit", "id_exit_ex", "forcedCu", "klcl", "kinv", "krel", "klfc", "kbup", "kpen", "thlsrc", "qtsrc", \
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

def exchange_cycle(line):
    line = re.sub("do iter = 1, iter_cin", loop_exchange1, line)
    line = re.sub("end do                ! End of implicit CIN loop \(cin_iter\)", \
                  loop_exchange20, line)
    line = re.sub("333 if\(id_exit\) then ! Exit without cumulus convection", \
                  loop_exchange21, line)
    line = re.sub("end do                  ! end of big i loop for each column.", \
                              loop_exchange22, line)
    return line


def var_expand(line):
    global flag1, flag2, flag3
    global value_list, nrplc_list, aflag_list
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
    return line

extra_lines = '''
    end do

    do i = 1, iend
      if (id_exit .or. id_exit_ex) then
        cycle
      endif
'''
def split_cycle(line, lines):
    global extra_lines, flag4, flag5, flag6, flag_ex
    if re.search("! CODE OF goto333 MUST BE COPIED & RENUMBERED HERE !!!", line) != None:
        str = "333   continue\n"
        lines.append(str)
        lines.append("!!!end part init\n")
        flag4 = True
    if re.match("!!!end part one", line) != None:
        flag4 = False
        flag5 = True
        str = "3331   continue\n"
        lines.append(str)
        for line_ex in extra_lines.split("\n"):
            lines.append(line_ex+"\n")
    if re.match("!!!end part two", line) != None:
        flag5 = False
        flag6 = True
        str = "3332   continue\n"
        lines.append(str)
        for line_ex in extra_lines.split("\n"):
            lines.append(line_ex+"\n")
    if re.match("!!!end part three\n", line) != None:
        flag6 = False
    if (flag4 == True):
        if re.search("id_exit = .false.", line) != None:
            lines.append("               id_exit_ex = .true.\n")
        line = re.sub('go to 333', 'go to 3331', line)
    if (flag5 == True):
        line = re.sub('go to 333', 'go to 3332', line)
    if (flag6 == True):
        line = re.sub('go to 333', 'go to 3333', line)
    line = re.sub("333 if\(id_exit", "3333 if(id_exit", line)
    return (line, lines)

key_words = ["conden", "if", "else", "kpen", "max", "min", "sqrt", "int"        \
             , "aint", "exp", "roots", "log", "slope", "real", "getbuoy"        \
             , "single_cin", "fluxbelowinv", "positive_moisture_single"         \
             , "uwshcu_readnl", "findsp_vc", "cnst_get_type_byind"              \
             , "cnst_get_ind", "phys_decomp", "addfld", "outfld", "qmin"        \
             , "qsinvert", "compute_alpha", "compute_mumin2", "compute_ppen"    \
             , "exnf", "qsat_o", "abs", "elseif", "status"]

tail_line = {}   # dic for line number of words tail
lines_subset = []
lines_s = 0
def match_value(line, count_line):
    global tail_line
    stack = []
    pos = []
    brackets = {')':'('}
    count = 0
    for char in line:
        if char in brackets.values():
            stack.append(char)
            pos.append(count)
        elif char in brackets.keys():
            try:
                if brackets[char] == stack.pop():
                    pos_e = pos.pop()
                    for i in range(pos_e - 1, 0, -1):
                        if not (line[i] == "_" or line[i].isdigit() or line[i].isalpha()):
                            break
                    word = line[i + 1 : pos_e]
                    if (word not in key_words) and len(word) != 0:
                        tail_line[word] = count_line
            except IndexError:
                continue
        count += 1

def is_skip_comments(line):
    flag = True
    for char in line:
        if char == " ":
            continue
        elif char == "!":
            flag = False
            break
    return flag

def tail_value(line, count_line):
    global  flag8, lines_s
    flag = is_skip_comments(line)
    if re.match("!!!start count", line) != None:
        flag8 = True
        lines_s = count_line
    if re.match("!!!end count", line) != None:
        flag8 = False
    if (flag == True and flag8 == True):
        match_value(line, count_line)
    if(flag8 == True):
        lines_subset.append(line)


vlist = []
sum_val = 0
value_1 = []
value_2 = []
value_3 = []
def vol(ki, word = None):
    dot_c = 0
    for s in ki:
        if "," == s:
            dot_c += 1
    if (dot_c == 0 and word != None):
        value_1.append(word)
    if (dot_c == 1 and word != None):
        value_2.append(word)
    if (dot_c == 2 and word != None):
        value_3.append(word)
    if (dot_c == 0):
        return 8
    if (dot_c == 1):
        return 8*30
    if (dot_c == 2):
        return 8*25*30

def count_vol(line, count_line):
    global key_words, vlist, sum_val, tail_line
    #print(line)
    stack = []
    pos = []
    brackets = {')':'('}
    count = 0
    word_ki = {}

    for char in line:
        if char in brackets.values():
            stack.append(char)
            pos.append(count)
        elif char in brackets.keys():
            try:
                if brackets[char] == stack.pop():
                    pos_e = pos.pop()
                    for i in range(pos_e - 1, 0, -1):
                        if not (line[i] == "_" or line[i].isdigit() or line[i].isalpha()):
                            break
                    word = line[i + 1 : pos_e]
                    ki = line[pos_e : count + 1]
                    if (word not in key_words) and (word != ""):
                        word_ki[word] = ki

            except IndexError:
                continue
        count += 1
    for word in word_ki.keys():
        #print(word)
        if (word not in vlist):
            vlist.append(word)
            sum_val += vol(word_ki[word], word)
        if (tail_line[word] == count_line):
            sum_val -= vol(word_ki[word])
    #print(line)
    print("line: %d, volume: %d B" %(count_line, sum_val))

vlist_s = []
sum_val_s = 0
def sum_vol(line, count_line):
    global key_words, vlist_s, sum_val_s, tail_line
    #print(line)
    stack = []
    pos = []
    brackets = {')':'('}
    count = 0
    word_ki = {}

    for char in line:
        if char in brackets.values():
            stack.append(char)
            pos.append(count)
        elif char in brackets.keys():
            try:
                if brackets[char] == stack.pop():
                    pos_e = pos.pop()
                    for i in range(pos_e - 1, 0, -1):
                        if not (line[i] == "_" or line[i].isdigit() or line[i].isalpha()):
                            break
                    word = line[i + 1 : pos_e]
                    ki = line[pos_e : count + 1]
                    if (word not in key_words) and (word != ""):
                        word_ki[word] = ki

            except IndexError:
                continue
        count += 1
    for word in word_ki.keys():
        #print(word)
        if (word not in vlist_s):
            vlist_s.append(word)
            sum_val_s += vol(word_ki[word])
    #print(line)
    print("line: %d, volume: %d B" %(count_line, sum_val_s))

def count_value():
    global tail_line, lines_s, lines_subset, value_1, value_2, value_3
    count = lines_s
    for word in value_1:
        print("v1:", word)
    for word in value_2:
        print("v2:", word)
    for word in value_3:
        print("v3:", word)
    #for key, value in tail_line.items():
    #    print('{key}:{value}'.format(key = key, value = value))
    #for line in lines_subset:
    #    print(count, line)
    #    count += 1

# 字符串匹配策略
def write_code():
    global value_list, nrplc_list, aflag_list
    with open("uwshcu_t.F90", "w") as fw:
        with open("uwshcu_o.F90", "r") as fr:
            files = fr.readlines()
            count_line = 1
            for line in files:
                #交换循环顺序
                lines = []
                line = exchange_cycle(line)
                # split cycle into three part, and this will bring extra codes which also need to expand
                (line, lines) = split_cycle(line, lines)
                lines.append(line)
                #变量扩展
                for line in lines:
                    line = var_expand(line)
                    fw.write(line)
                    #trans_array(line, count_line)
                    tail_value(line, count_line)
                    count_line += 1
    print("--->finish")

def line_volume():
    count = lines_s
    print("************************* max ******************************")
    for line in lines_subset:
        if(is_skip_comments(line)):
            count_vol(line, count)
        count += 1

    count = lines_s
    print("************************* sum ******************************")
    for line in lines_subset:
        if(is_skip_comments(line)):
            sum_vol(line, count)
        count += 1
inter_words = ["iter_scaleh", "iter_xc", "id_check", "status", "klcl", "kinv"   \
              , "krel", "klfc", "kbup", "kpen", "iter", "ixnumliq", "ixnumice"  \
              , "ixcldli1", "ixcldice"]
# without out value
value_k_0 = ["ps0_in", "zs0_in", "tke_in", "umf_out", "slflx_out", "qtflx_out"  \
             , "flxprc1_out", "flxsnow1_out", "ufrc_out", "uflx_out", "vflx_out"\
             , "trflx_out", "ps0", "zs0", "tke", "exns0", "umf", "emf", "slflx" \
             , "qtflx", "uflx", "vflx", "ufrc", "trflx", "trlfx_d", "trflx_u"   \
             , "uemf", "wu", "thlu", "qtu", "uu", "vu", "thvu", "tru"           \
             , "thlu_emf", "qtu_emf", "uu_emf", "vu_emf", "tru_emf", "flxrain"  \
             , "xflx", "wu_out", "qtu_out", "thlu_out"]
key_logical = ["id_exit", "id_exit_ex", "forceCu"]

def cons_struct():
    global tail_line, lines_s, lines_subset, value_1, value_2, value_3
    count = lines_s
    struct_s = '''
#define SW_UWSHCU_01
#ifdef SW_UWSHCU_01
    external :: slave_uwshcu_c1
    type param_uw_com_t1
'''
    struct_e = '''
    end type param_uw_com_t1
    type(param_uw_com_t1) param_uw_com_s1(iend)
'''
    struct_s1 = '''
    type param_uw_t1
      integer*8 :: ptr_param_uw1
      integer :: iend
    end type param_uw_t1
    type(param_uw_t1) param_uw_s1
'''
    struct_init1 = '''
    param_uw_s1%ptr_param_uw1 = loc(param_uw_com_s1)
    param_uw_s1%iend = iend
'''
    with open("struct.F90", "w") as fw:
        fw.write(struct_s)
        for word in value_2:
            fw.write("      integer*8 :: %s\n"%word)
        for word in value_3:
            fw.write("      integer*8 :: %s\n"%word)
        for word in value_1:
            if (word not in inter_words) and (word not in key_logical):
                fw.write("      real(r8) :: %s\n"%word)
        for word in value_1:
            if (word in inter_words) or (word in key_logical):
                fw.write("      integer :: %s\n"%word)
        fw.write("      integer :: mkx\n")
        fw.write(struct_e)
        fw.write(struct_s1)
        fw.write("#endif")
        fw.write(struct_init1)
        fw.write("    do i = 1, iend\n")
        for word in value_2:
            fw.write("      param_uw_com_s1(i)%s%s = loc(%s)\n"%("%", word      \
                      , word))
        for word in value_3:
            fw.write("      param_uw_com_s1(i)%s%s = loc(%s)\n"%("%", word      \
                      , word))
        for word in value_1:
            if word in key_logical:
                fw.write('      if (%s(i) .eq. .true.) then\n'%word)
                fw.write('        param_uw_com_s1(i)%s%s = 1\n'%("%", word))
                fw.write('      else\n')
                fw.write('        param_uw_com_s1(i)%s%s = 0\n'%("%", word))
                fw.write('      endif\n')
            else:
                fw.write('      param_uw_com_s1(i)%s%s = %s(i)\n'%("%", word    \
                          , word))
        fw.write("      param_uw_com_s1(i)%mkx = mkx\n")
        fw.write("    end do\n")

    inc = '''
#include <slave.h>
#include "dma_macros.h"
'''
    struct_s1 = '''
typedef struct {
  double *ptr_param_uw1;
  int iend;
} param_uw_t;
'''
    struct_s = '''
typedef struct {
'''
    struct_e = '''} param_t;
'''
    func_s = '''
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
'''
    func_e = '''
  }
}
'''
    with open("slave_uwshcu_c1.c", "w") as fw:
        fw.write(inc)
        fw.write(struct_s1)
        fw.write(struct_s)
        for word in value_2:
            fw.write("    double *%s\n"%word)
        for word in value_3:
            fw.write("    double *%s\n"%word)
        for word in value_1:
            if word not in inter_words:
                fw.write("    double %s\n"%word)
        for word in value_1:
            if word in inter_words:
                fw.write("    int %s\n"%word)
        fw.write(struct_e)
        fw.write(func_s)
        for word in value_1:
            if word not in inter_words:
                fw.write("    double %s = param_d.%s;\n"%(word, word))
        for word in value_1:
            if word in inter_words:
                fw.write("    int %s = param_d.%s;\n"%(word, word))
        for word in value_2:
            if word in value_k_0:
                fw.write("    double %s[mkx + 1];\n"%word)
            else:
                fw.write("    double %s[mkx];\n"%word)
        for word in value_3:
            fw.write("    double %s[ncnst][mkx];\n"%word)
        for word in value_2:
            if word in value_k_0:
                fw.write("    pe_get(param_d.%s, %s, sizeof(double)*(mkx + 1));\n"%(word, word))
            else:
                fw.write("    pe_get(param_d.%s, %s, sizeof(double)*mkx);\n"%(word, word))
        for word in value_3:
            fw.write("    pe_get(param_d.%s, %s, sizeof(double)*mkx*ncnst);\n"%(word, word))
        fw.write("    dma_syn();\n")
        call_f_s = "    uw_compute_loop1_("
        call_f_e = ");\n"
        fw.write(call_f_s)
        for word in value_1:
            fw.write("&%s, "%word)
        for word in value_2:
            fw.write("%s, "%word)
        for word in value_3[:len(value_3) - 1]:
            fw.write("%s, "%word)
        fw.write(value_3[len(value_3) - 1])
        fw.write(call_f_e)
        fw.write(func_e)

def trans_array(line, count_line):
    global key_words, vlist, sum_val, tail_line
    #print(line)
    stack = []
    pos = []
    brackets = {')':'('}
    count = 0
    word_ki = {}

    for char in line:
        if char in brackets.values():
            stack.append(char)
            pos.append(count)
        elif char in brackets.keys():
            try:
                if brackets[char] == stack.pop():
                    pos_e = pos.pop()
                    for i in range(pos_e - 1, 0, -1):
                        if not (line[i] == "_" or line[i].isdigit() or line[i].isalpha()):
                            break
                    word = line[i + 1 : pos_e]
                    ki = line[pos_e + 1 : count]
                    if (word not in key_words) and (word != ""):
                        word_ki[word] = ki
            except IndexError:
                continue
        count += 1
    word_list = []
    lis = 1
    ris = 1
    if len(word_ki) == 2:
        for word in word_ki.keys():
            word_list.append(word)
        #print(line, word_ki[word_list[0]], word_ki[word_list[1]])
        li = word_ki[word_list[0]].split(",")[0]
        if ":" in li:
            lis = li.split(":")[0]
            if lis == "":
                lis = 1
        else:
            lis = li
        ri = word_ki[word_list[0]].split(",")[1]
        if ":" in ri:
            ris = ri.split(":")[0]
            if ris == "":
                ris = 1
        print(line, lis, ris)



# 主程序
def main(argv):
    write_code()
    line_volume()
    count_value()
    cons_struct()

if __name__ == "__main__":
   sys.exit(main(sys.argv))
