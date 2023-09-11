import string
import os
import re
import json
from decimal import *
import pandas as pd
import openpyxl 
from joblib import Parallel, delayed
import multiprocessing
import xlwings as xw
import time
import chardet
import interval

# coding=utf-8

#



pd.set_option('display.max_rows', 500)  
pd.set_option('display.max_columns', 500)  
pd.set_option('display.width', 1000)  


with open("config_peak.json", "r") as jn:
    config_peak = json.loads(jn.read())
    

with open("config_control.json", "r") as jn_control:
    config_control = json.loads(jn_control.read())
    
folders = ['infile', 'outfile','lastfile','error_log','seg_file','lastfile/ALL','lastfile/HP','lastfile/MP','lastfile/LP','lastfile/xDMA','ucsc_bed','ucsc_bed/bedGf','ucsc_bed/region','file_mer','peak_type']

for folder in folders:
    if not os.path.isdir(folder):
        os.mkdir(folder)

infile="in_outfile\\"
outfile="outfile\\"
lastfile="lastfile\\"
hpfile = "HP\\"
mpfile = "MP\\"
lpfile = "LP\\"
all_file = "ALL\\"
ssr_all_file = r"ssr_all\\"
xDMA = "xDMA\\"
xDMAfile = "lastfile\\xDMA\\"
error_log = "error_log\\"
seg_file ="seg_file\\"
bed_outfile = "ucsc_bed\\"
ALL_mer = 'ALL_mer\\'
mer_file = 'ssr_all\\ALL_mer\\'
land_bgf = 'ucsc_bed\\bedGf\\'
ucsc_regfile = 'ucsc_bed\\region\\'
file_mers = 'file_mer\\'
excel_subpf = "lastfile\\excel_subp\\"
peak_typef = "peak_type\\"
in_outfile = r"in_outfile\\"




flag_replace = '1'

seg_list = os.listdir(seg_file)  

global data_seg

if len(seg_list) != 0:
    data_seg = pd.read_excel(seg_file + seg_list[0], engine="openpyxl")  
    old_name = 'G38'
    new_name = 'GRCh38'
else:
    old_name = 'HS'
    new_name = 'CHM13'
    
flag_name = 0
peak_name_new = ''

with open("config.json", "r") as fp:
    config = json.loads(fp.read())

print("\n")
print("Processing》》》》》")
print("\n")



def get_encoding(file):
    with open(file,'rb') as f:
        file_format = chardet.detect(f.read())['encoding']
        #print(file_format)
        return file_format
        
        
def csv_peak(file_name):
    LP, MP, sHP =0, 0, 0
    name_i = file_name
    outname = name_i+"_MP+"+".csv"
    inname = name_i+"_all"+".csv"
    flag_ins = 0
    flag_insl = 0
    LMDP_zoom = interval.Interval(90,150,upper_closed=False)
    MMDP_zoom = interval.Interval(150,300,upper_closed=False)
    HMDP_zoom = interval.Interval(300,1001,upper_closed=False)
    with open('in_outfile\\'+outname,"w") as pf:
        pf.write('Peak_name,Peak_type,Type,position,label,ratio_mono,ratio_di,ratio_tri,'\
        'ratio_tetra,ratio_penta,ratio_hexa,ratio_total,ratio_len_mono,ratio_len_di,'\
        'ratio_len_tri,ratio_len_tetra,ratio_len_penta,ratio_len_hexa,ratio_len_total,ex_num,ex_den,insert_pos\n')
        with open(ssr_all_file+ALL_mer+inname,"r") as fp:
            fp.readline()  
            is_end=False
            while True:
                line=fp.readline()  
                if not line:  
                    is_end=True
                    print(inname,"   "+"Invalid data")
                    break    
                datas=line.split(',')
                ratio_len_total=int(datas[16])
                #last_pos = float(datas[1])
                front_pos = int(datas[1]) 
                if ratio_len_total<90:
                    continue  
                else:         
                    lastline=line
                    front=int(datas[1])  
                    last_ratio_len_total=ratio_len_total
                    last_type=datas[0]
                    if ratio_len_total in HMDP_zoom:
                        flag_HML = 'HMDP'
                    elif ratio_len_total in MMDP_zoom:
                        flag_HML = 'MMDP'
                    elif ratio_len_total in LMDP_zoom:
                        flag_HML = 'LMDP'
                    else:
                        print("ERROR---1")
                    break
              
            while not is_end:  
            
                lines=[lastline]        
                #print(lines)
                max_length=last_ratio_len_total
                min_length = last_ratio_len_total
                
                while True: 
                    line=fp.readline()
                    if not line: 
                        is_end=True
                        break
                        
                    datas=line.split(',')
                    #print(datas[16])
                    ratio_len_total=int(datas[16])  
                    if ratio_len_total<90:
                        front_pos = int(datas[1]) 
                        continue             
                        
                    now=int(datas[1])  
                    if now - front_pos < 1000 and now - front_pos >0:
                        npeak_e= str(front_pos)
                        npeak_e = npeak_e[-3:]
                        if npeak_e != '000':
                            #print(now)
                            #print(front_pos)
                            #print("-----------------------------------------"+"\n")
                            insert_pos = str(front_pos+1) 
                            flag_pos = front_pos
                            front_pos = now
                            flag_ins = 1 
                    if ratio_len_total in HMDP_zoom:
                        flag_HMLs = 'HMDP'
                    elif ratio_len_total in MMDP_zoom:
                        flag_HMLs = 'MMDP'
                    elif ratio_len_total in LMDP_zoom:
                        flag_HMLs = 'LMDP'
                    else:
                        print("Error---2")
                    if now-front<=1000 and datas[0]==last_type and flag_HML == flag_HMLs:
                        lines.append(line)   
                        #print(lines)
                        max_length=max(max_length,ratio_len_total)  
                        min_length =min(min_length,ratio_len_total) 
                        front=now
                        last_ratio_len_total=ratio_len_total
                        last_type=datas[0]
                    else:
                        lastline=line
                        front=now
                        last_ratio_len_total=ratio_len_total
                        last_type=datas[0]
                        if ratio_len_total in HMDP_zoom:
                            flag_HML = 'HMDP'
                        elif ratio_len_total in MMDP_zoom:
                            flag_HML = 'MMDP'
                        elif ratio_len_total in LMDP_zoom:
                            flag_HML = 'LMDP'
                        else:
                            print("Error---3")
                        break
                
                
                if max_length<150:
                    LP+=1
                    label='LMDP'
                    if min_length <90:
                        print("Peak classification error---LMDP")
                        print(str(LP)+"------"+outname)
                        #print(outname)
                        print("---------------------------------------")
                elif max_length<300:
                    MP+=1
                    label='MMDP'
                    if min_length <150:
                        print("Peak classification error---MMDP")
                        print(str(MP)+"------"+outname)
                        #print(outname)
                        print("---------------------------------------")
                else:
                    sHP+=1
                    label='HMDP'
                    if min_length <300:
                        print("Peak classification error---HMDP")
                        print(str(sHP)+"------"+outname)
                        print("---------------------------------------")
                """
                elif max_length<425:
                    HP+=1
                    label='HP'
                else:
                    sHP+=1
                    label='sHP'
                """
                    
                for line in lines:  
                    datas=line.split(',')
                    line_pos = int(datas[1])
                    #print(line_pos)
                    #print(datas[0])
                    name=(datas[0].split(': '))[1] 
                    if flag_name  == 1:
                        name = peak_name_new
                    #print(name)
                    if flag_ins == 0:
                        if  label=='LMDP':
                            line=f'{name}-{label}{LP},{name}-{label},'+line
                        elif label=='MMDP':
                            line=f'{name}-{label}{MP},{name}-{label},'+line       
                        else:
                            line=f'{name}-{label}{sHP},{name}-{label},'+line
                    elif flag_ins == 1:
                        if flag_insl == 0:                            
                            if label == 'LMDP':
                                line = f'{name}-{label}{LP},{name}-{label},' + line  
                                peak_count = f'{name}-{label}{LP}'
                            elif label == 'MMDP':
                                line = f'{name}-{label}{MP},{name}-{label},' + line  
                                peak_count = f'{name}-{label}{MP}'
                            else:
                                line = f'{name}-{label}{sHP},{name}-{label},' + line
                                peak_count = f'{name}-{label}{sHP}'
                            flag_insl = 1
                        elif flag_insl == 1:
                            if label=='LMDP':
                                if line_pos -flag_pos < 1000 and line_pos -flag_pos > 0:
                                    line=f'{name}-{label}{LP},{name}-{label},'+line.split('\n')[0]+","+insert_pos+"\n"  
                                    flag_ins = 0
                                    flag_insl = 0
                                else:
                                    line = f'{name}-{label}{LP},{name}-{label},' + line  
                            elif label=='MMDP':
                                if line_pos -flag_pos < 1000 and line_pos -flag_pos > 0:
                                    line=f'{name}-{label}{MP},{name}-{label},'+line.split('\n')[0]+","+insert_pos+"\n"       
                                    flag_ins = 0
                                    flag_insl = 0
                                else:
                                    line = f'{name}-{label}{MP},{name}-{label},' + line  
                            else:
                                if line_pos -flag_pos < 1000 and line_pos -flag_pos > 0:
                                    #print("test1")
                                    line=f'{name}-{label}{sHP},{name}-{label},'+line.split('\n')[0]+","+insert_pos+"\n"
                                    flag_ins = 0
                                    flag_insl = 0
                                else:
                                    #print("test2")
                                    line = f'{name}-{label}{sHP},{name}-{label},' + line
                                    #print(line)
                        #print(line)
                        else:
                            print("---Error inserting position data---" + "\n")
                    else:
                        print("---Error inserting position data---"+"\n")
                    pf.write(line)
    print("\n"+outname+"---"+"Generated")

def txt_merge(txt_name,txt_cm):
    with open(outfile+txt_name+".txt","r",encoding="utf-8") as fp:
        #fp.readline()
        data=fp.readlines()
    datas=data[2:]
    with open('in_outfile\\'+txt_name.split('-S')[0]+"_all.txt","a+",encoding="utf-8") as pf:
        if txt_cm == 0:
            pf.write(data[0])
            pf.write(data[1])
        for temp_data in datas:
            pf.write(temp_data)
        #pf.write(datas)
    #print('%d: '%txt_cm,txt_name+".txt","Successfully processed")
 

def multi_proc(name_ct):
    peak_name = name_ct+"_MP+"+".csv"
    ssr_name = name_ct+"_all-classify"+".txt"
    peak_in = open(infile+peak_name, "r") 
    ssr_in = open(infile+ssr_name, "r") 
    peak_all = peak_in.readlines()
    ssr_all = ssr_in.readlines()
    #ssr_row=len(ssr_all)
    #print(ssr_row)
    peak_in.close()
    ssr_in.close()
    txt_data = pd.read_csv(infile+name_ct+"_all-classify"+".txt",sep='\t')
    check_data = txt_data
    #print(txt_dict) 
    ssr_type = []
    ssr_count=[]
    peak_counts = []
    peak_cpl = []
    peak_suml = []
    type_all = {}
    check_type = {}
    dict_pse = {} 
    type_name = peak_name.split(".")[0] + "_type.csv"
    type_out = open(outfile+type_name, "w")
    type_out.write('Peak_name,Peak_type,Type,position,label,ratio_len_total,motif_type,ratio(%),size_sum,motif_abundance,abundance_ratio(%),count_sum,count_sum_all\n')
    inx = 1  
    flag_start = 0
    flag_check = 0
    for i in peak_all[1:]:
        i = i.strip().split(',')  
        #print(i)
        peak_e = int(i[3])  
        #print(peak_e)
        #peak_s = peak_e-1000  
        npeak_e=i[3]
        npeak_e = npeak_e[-3:]
        if npeak_e != '000':
           lnpeak_e=i[3]
           e_len = len(lnpeak_e)
           lnpeak_e=lnpeak_e[:-3]
           peak_s=int(lnpeak_e.ljust(e_len,'0'))+1 
           inx -= 1
        else:
            peak_s = peak_e - 999  
        #print(peak_s)
        #print(peak_e)
        if len(i) == 22:
            #print(i)
            if i[21] != '':
                peak_s = int(i[21])
        t_l = float(i[18])     
        dict_pse[peak_s] = peak_e
        ssr_dict = {}
        dssr_len={}
        peak_count = {}
        peak_cpd = {} 
        #ssr_countp ={}
        #flag_check = 0
        for j in ssr_all[inx:]: 
            ssr = j.strip().split() 
            ssr_n = ssr[0]  
            ssr_s = int(ssr[3]) 
            ssr_e = int(ssr[4])
            ssr_l = int(ssr[2])
            if ssr_e < peak_s: 
                inx += 1
                flag_count = 0
            elif ssr_s > peak_s and ssr_e <= peak_e: 
                if ssr_n not in ssr_dict.keys(): 
                    #print(ssr_s)
                    ssr_dict[ssr_n] = ssr_l/t_l  
                    dssr_len[ssr_n] = ssr_l
                    if ssr_n not in peak_count.keys():
                        peak_count[ssr_n] = 1
                    else:
                        peak_count[ssr_n] += 1
                else:
                    ssr_dict[ssr_n] += ssr_l/t_l 
                    dssr_len[ssr_n] += ssr_l
                    if ssr_n not in peak_count.keys():
                        peak_count[ssr_n] = 1
                    else:
                        peak_count[ssr_n] += 1
                inx += 1 
                flag_count = 0

            elif ssr_s <= peak_s and ssr_e <= peak_e: #
                ov_len =ssr_e -peak_s + 1
                com_value = peak_s - ssr_s 
                if ssr_n not in ssr_dict.keys(): 
                    #print(ssr_s)
                    ssr_dict[ssr_n] = ov_len/t_l  
                    dssr_len[ssr_n] = ov_len
                    if ssr_n not in peak_count.keys():
                        if ov_len > com_value:
                            peak_count[ssr_n] = 1
                    else:
                        if ov_len > com_value:
                            peak_count[ssr_n] += 1
                else:
                    ssr_dict[ssr_n] += ov_len/t_l 
                    dssr_len[ssr_n] += ov_len
                    if ssr_n not in peak_count.keys():
                        if ov_len > com_value:
                            peak_count[ssr_n] = 1
                    else:
                        if ov_len > com_value:
                            peak_count[ssr_n] += 1
                inx += 1 
                flag_count = 0
            elif  ssr_e > peak_e:  # ssr_e k_end 
                if ssr_s <= peak_e and peak_s <= ssr_s:
                    ov_len = ssr_e - peak_e
                    k_len = ssr_l-ov_len  
                    com_value = peak_e - ssr_s +1
                    #print(k_len)
                    if ssr_n not in ssr_dict.keys(): 
                        #print(ssr_s)
                        ssr_dict[ssr_n] = k_len/t_l  
                        dssr_len[ssr_n] = k_len
                        if ssr_n not in peak_count.keys():
                            if com_value >= ov_len:
                                peak_count[ssr_n] = 1
                        else:
                            if com_value >= ov_len:
                                peak_count[ssr_n] += 1
                    else:
                        #print(ssr_s)
                        ssr_dict[ssr_n] += k_len/t_l 
                        dssr_len[ssr_n] += k_len
                        if ssr_n not in peak_count.keys():
                            if com_value >= ov_len:
                                peak_count[ssr_n] = 1
                        else:
                            if com_value >= ov_len:
                                peak_count[ssr_n] += 1
                elif ssr_s <= peak_e and ssr_s < peak_s:
                    ov_len = ssr_e - peak_e
                    len_temp = peak_s - ssr_s
                    ov_len = ov_len + len_temp
                    k_len = ssr_l-ov_len
                    com_value = ssr_e - ssr_s  
                    if ssr_n not in ssr_dict.keys(): 
                        #print(ssr_s)
                        ssr_dict[ssr_n] = k_len/t_l  
                        dssr_len[ssr_n] = k_len
                        if ssr_n not in peak_count.keys():
                            if flag_count == 0:
                                peak_count[ssr_n] = 1
                                flag_count = 1
                        else:
                            if flag_count == 0:
                                peak_count[ssr_n] += 1
                                flag_count = 1
                    else:
                        #print(ssr_s)
                        ssr_dict[ssr_n] += k_len/t_l 
                        dssr_len[ssr_n] += k_len
                        if ssr_n not in peak_count.keys():
                            if flag_count == 0:
                                peak_count[ssr_n] = 1
                                flag_count = 1
                        else:
                            if flag_count == 0:
                                peak_count[ssr_n] += 1
                                flag_count = 1
                break
        ssr_type.append(ssr_dict) 
        ssr_count.append(dssr_len)
        peak_counts.append(peak_count)
        peak_sum= sum(peak_count.values())
        peak_suml.append(peak_sum)
        for t_p in peak_count:
            peak_cpd[t_p] = peak_count[t_p]/peak_sum
        peak_cpl.append(peak_cpd)
        if len(peak_count) != 0: 
            for t_p in peak_count:
                if t_p not in type_all:
                    type_all[t_p] = peak_count[t_p]
                else:
                    type_all[t_p] += peak_count[t_p]
        sum_c=0
        for temp_s in dssr_len.keys():
            sum_c+=dssr_len[temp_s]
        if sum_c != t_l:
            print("Incorrect data statistics："+"\n")
            print("sum_c:  "+str(sum_c)+"--------"+"ratio_len_total: "+str(t_l)+"------"+"Peak_name:  "+str(i[0])+"\n")
            print("peak_e :"+str(peak_e)+"\n")
            print("ssr_s:  "+str(ssr_s)+"\n")
            print("inx:  "+str(inx)+"\n")
            print("----------------------------------------"+"\n")
            err_1 = "sum_c:  "+str(sum_c)+"--------"+"ratio_len_total: "+str(t_l)+"------"+"Peak_name:  "+str(i[0])+"\n"
            err_2 = "peak_e :"+str(peak_e)+"\n"
            err_3 = "ssr_s:  "+str(ssr_s)+"\n"
            err_4 = "inx:  "+str(inx)+"\n"
            err_5 = "----------------------------------------"+"\n"
            with open(error_log+name_ct+"-len_total_err1.txt", "a+") as fp_err:
                fp_err.write("Data error"+"\n")
                fp_err.write(err_1)
                fp_err.write(err_2)
                fp_err.write(err_3)
                fp_err.write(err_4)
                fp_err.write(err_5)
        if t_l>1000:
            print("Ratio_ Len_ Total greater than 1000，Data error"+"----"+str(ssr_s)+"  "+str(i[0])+"  "+str(i[3])+"\n")
            err_1 = "Ratio_ Len_ Total greater than 1000，Data error"+"----"+str(ssr_s)+"  "+str(i[0])+"  "+str(i[3])+"\n"
            with open(error_log+name_ct+"-len_total1001_err2.txt","a+") as pf_err:
                pf_err.write(err_1)

    for t in range(len(ssr_type)):
        peak_line = peak_all[1+t].strip().split(',') 
        sort_type = sorted(ssr_type[t].items(), key=lambda d: d[1], reverse=True)  
        sort_count = sorted(ssr_count[t].items(), key=lambda d: d[1], reverse=True) 
        sort_peak = sorted(peak_counts[t].items(), key=lambda d: d[1], reverse=True) 
        sort_peakp = sorted(peak_cpl[t].items(), key=lambda d: d[1], reverse=True)  
        new_ratio=[]
        new_num=[]
        last_type = []  
        new_count=[]
        new_peakt = [] 
        new_peakc = [] 
        new_peakp = [] 
        for temp_r in sort_type:
            temp_num=round(temp_r[1],4) 
            temp_num="%.1f%%"%(temp_num * 100)
            temp_num=str(temp_num).split('%')
            new_ratio.append(temp_r[0])
            new_num.append(temp_num[0])
        for temp_c in  sort_count:
             new_count.append(str(temp_c[1]))
             last_type.append(str(temp_c[0]))
        for temp_cc in sort_peak:
            new_peakt.append(str(temp_cc[0]))
            new_peakc.append(str(temp_cc[1]))
        for temp_p in sort_peakp:
            temp_num=round(temp_p[1],4) 
            temp_num="%.1f%%"%(temp_num * 100)
            temp_num=str(temp_num).split('%')
            new_peakp.append(str(temp_num[0]))
        if peak_suml[t] == 0:
            new_peakc.append('0')
            new_peakp.append('0')
        type_out.write(",".join(peak_line[:5]) + "," + peak_line[18] + "," + ";".join(last_type) +"," +";".join(new_num)+","+";".join(new_count)+","+";".join(new_peakt)+","+";".join(new_peakp)+","+";".join(new_peakc)+","+str(peak_suml[t])+"\n")
    type_out.close()
    print("%s"%peak_name+"与%s Processing completed"%ssr_name)
    print("\n")
    #print("Task finished!")



def classify_pro(peak_namef):
    type_name = peak_namef.split(".")[0] + "_last.csv"
    with open(outfile+type_name, "w") as pf:
        pf.write('Peak_name,Peak_type,Type,position,label,pD1RD sum,motif_type,ratio(%),Motif size sum,motif_abundance,abundance_ratio(%),count_sum,count_sum_all,pD1RD\n')
        with open(outfile+peak_namef, "r",encoding=get_encoding(outfile+peak_namef)) as fp:
            fp.readline()  
            is_end=False
            start_line = fp.readline()
            start_datas = start_line.split(',')
            peak_name1 = start_datas[0]
            while not is_end:
                while True:
                    line=fp.readline()  
                    if not line:  
                        is_end=True
                        pf.write(start_line)
                        print("Scan ended---1")
                        break     
                    datas=line.split(',')
                    peak_name=datas[0]
                    if peak_name1 != peak_name: 
                        pf.write(start_line)
                        start_line=line
                        peak_name1=datas[0]
                    elif peak_name1 == peak_name:  
                        last_peak_name=peak_name
                        lines=[start_line]
                        lines.append(line)
                        #print(lines)
                        break 
                    else:
                        print(peak_name+"-----"+"Classification processing error")
                while not is_end:  
                    while True: 
                        new_line=fp.readline()
                        if not new_line: 
                            is_end=True
                            print("Scan End ---2: Same peak at the end_ Name, processed")
                            break
                        new_datas=new_line.split(',')
                        new_peak_name=new_datas[0]  
                        if last_peak_name != new_peak_name: 
                            peak_name1=new_peak_name
                            start_line=new_line
                            break                               
                        elif last_peak_name == new_peak_name: 
                            lines.append(new_line)
                            last_peak_name=new_peak_name
                            #print(lines)
                        else:
                            print("Processing error")
                            break
                    #last_len=len(lines)
                    #print(last_len)
                    #print(lines)
                    dict_type={}
                    dict_typea={}
                    dict_type_ratio={}
                    dict_typea_ratio ={}
                    dict_type_count={}
                    type_num=0
                    type_numa = 0
                    all_ratio_len_total = 0
                    all_mer_count_sum = 0
                    for line_typet in lines: 
                        temp_linet = line_typet.split(',')
                        all_ratio_len_total=int(temp_linet[5])+all_ratio_len_total
                        all_mer_count_sum = int(temp_linet[12])+all_mer_count_sum
                    #print(all_ratio_len_total)
                    count_m=len(lines)
                    for line_type in lines: 
                        temp_line=line_type.split(',')
                        ratio_len_total=int(temp_line[5])
                        motif_type=temp_line[6].split(';')
                        count_s=temp_line[8].strip().split(';')
                        motif_abundance = temp_line[9].split(';')
                        count_suma = temp_line[11].strip().split(';')
                        #print(ratio)
                        n=0
                        for temp_motif_type in motif_type:
                            temp_motif_type=temp_motif_type.lstrip() 
                            if type_num==0:  
                                dict_type[temp_motif_type]=int(count_s[n])
                            elif type_num!=0 and temp_motif_type not in dict_type: 
                                dict_type[temp_motif_type] = int(count_s[n])
                            elif type_num!=0 and temp_motif_type in dict_type: 
                                dict_type[temp_motif_type] = int(count_s[n]) + dict_type[temp_motif_type]
                            n+=1
                        type_num +=1
                        k=0
                        for temp_at in motif_abundance:
                            temp_at=temp_at.lstrip() 
                            if type_numa==0:  
                                dict_typea[temp_at]=int(count_suma[k])
                            elif type_numa!=0 and temp_at not in dict_typea: 
                                dict_typea[temp_at] = int(count_suma[k])
                            elif type_numa!=0 and temp_at in dict_typea: 
                                dict_typea[temp_at] = int(count_suma[k]) + dict_typea[temp_at]
                            k+=1
                        type_numa +=1
                    #print(dict_type)
                    
                    for temp_dict in dict_type.keys(): 
                        dict_type_ratio[temp_dict]=dict_type[temp_dict]/all_ratio_len_total
                    
                    for temp_p in dict_typea.keys():
                        dict_typea_ratio[temp_p]=dict_typea[temp_p]/all_mer_count_sum
                    
                    sort_type_ratio = sorted(dict_type_ratio.items(), key=lambda d: d[1], reverse=True) 
                    sort_type_count = sorted(dict_type.items(), key=lambda d: d[1], reverse=True) 
                    sort_at_ratio =  sorted(dict_typea_ratio.items(), key=lambda d: d[1], reverse=True) 
                    sort_at_count =  sorted(dict_typea.items(), key=lambda d: d[1], reverse=True) 
                    #print(sort_type_ratio)
                    #p=0.00
                    type_new=[]
                    type_a_new =[]
                    ratio_new=[]
                    ratio_a_new=[]
                    count_new=[]
                    count_a_new =[]
                    dict_ratio={}
                    dict_a_ratio ={}
                    for r in sort_type_ratio:
                        #p+=r[1]
                        type_new.append(str(r[0]))  
                        dict_ratio[r[0]] = r[1]     
                        #temp_r=r[1]
                        #temp_r="%.1f%%" % (temp_r * 100)
                        temp_rs   = round(r[1],4)
                        temp_rs="%.1f%%"%(temp_rs * 100)
                        temp_rs=str(temp_rs).split('%')
                        ratio_new.append(temp_rs[0])
                    for ra in sort_at_ratio:
                        #p+=r[1]
                        type_a_new.append(str(ra[0]))  
                        dict_a_ratio[ra[0]] = ra[1]     
                        #temp_r=r[1]
                        #temp_r="%.1f%%" % (temp_r * 100)
                        temp_rs   = round(ra[1],4)
                        temp_rs="%.1f%%"%(temp_rs * 100)
                        temp_rs=str(temp_rs).split('%')
                        ratio_a_new.append(temp_rs[0])
                    for c in sort_type_count:
                        count_new.append(str(c[1]))                    
                    for ca in sort_at_count:
                        count_a_new.append(str(ca[1]))                    
                    #print(lines[0])
                    line_p = lines[0].split(',')
                    pos_m = line_p[3]
                    list_len=[]
                    for line_t in lines:
                        line_ta =line_t.split(',')
                        list_len.append(int(line_ta[5])) 
                        #print(int(line_ta[5]))
                    len_max = max(list_len)
                    #print(list_len)
                    pf.write(temp_line[0]+"m"+str(count_m)+","+",".join(temp_line[1:5])+","+str(all_ratio_len_total)+","+";".join(type_new)+","+";".join(ratio_new)+","+";".join(count_new)+","+";".join(type_a_new)+","+";".join(ratio_a_new)+","+";".join(count_a_new)+","+str(all_mer_count_sum)+","+str(len_max)+"\n")
                    break
        print(type_name+"---"+"Processing completed")




def count_sum_pro(peak_namel):
    if len(seg_list) != 0:
        global data_seg
        gap_flag = '1'
    else:
        gap_flag = '0'

    type_name = peak_namel.split("MP")[0] + "LP+_count.csv"
    HP_name = peak_namel.split("MP")[0] + "HP_count.csv"
    MP_name = peak_namel.split("MP")[0] + "MP_count.csv"
    LP_name = peak_namel.split("MP")[0] + "LP_count.csv"
    with open(lastfile+all_file+type_name, "w") as pf:
        pf.write('Peak_name,Peak_type,S_position,E_position,label,Motif type,Motif size sum,Percentage(%),Integrated motif type,pD1RD,pD1RD sum,motif_abundance,count_sum,abundance_ratio(%),count_sum_all,ucsc_name\n')
        with open(outfile+peak_namel, "r",encoding=get_encoding(outfile+peak_namel)) as fp:
            lines = fp.readlines()[1:]
            for line in lines:
                data = line.split(',')
                chr_type=data[1].split('-')
                peak_type=chr_type[2]  
                n_peak = data[0]
                if peak_type == 'HMDP' or peak_type == 'MMDP' or peak_type == 'LMDP':
                #if peak_type == 'LMDP':
                    new_data=line.strip().split(',') 
                    name_peak = new_data[0] # peak_name
                    type_peak = new_data[1] # peak_type
                    e_position = new_data[3] #E_position
                    peak_lable = new_data[4] # lable
                    type_motif = new_data[6].split(';')  
                    a_type_motif = new_data[9].split(';') 
                    sum_count =  new_data[8].split(';')  
                    a_sum_count = new_data[11].split(';')  
                    type_ratio = new_data[7].split(';')  
                    a_type_ratio = new_data[10].split(';')  
                    ratio_len_total=int(new_data[5])      
                    a_ratio_len_total = int(new_data[12]) 
                    peak_n = name_peak.split('-')
                    peak_chr = peak_n[0]+'-'+peak_n[1]
                    if gap_flag == '1':
                        data_seg_temp = data_seg[data_seg['filename'].str.contains(peak_chr)]
                    #m_peak = name_peak.split('P')[1]
                    if 'm' in str(name_peak):
                        #print(peak_chr)
                        npeak_e = str(e_position)[-3:]
                        if npeak_e != '000':
                            e_len = len(str(e_position))
                            lnpeak_e = str(e_position)[:-3]
                            new_e_position = int(lnpeak_e.ljust(e_len, '0'))  
                            s_position = int(new_e_position) - (int(peak_n[2].split('m')[1])-1)*1000
                            s_position = s_position + 1
                            #print(s_position)
                        else:
                            s_position = int(e_position) - int(peak_n[2].split('m')[1]) * 1000
                            s_position = s_position + 1
                        if gap_flag == '1':
                            #print(data_seg)
                            for rows in data_seg_temp.index:
                                if peak_chr in data_seg_temp.loc[rows]['filename']:
                                    #print(type(data_seg.loc[rows]['from']))
                                    if 0 < data_seg_temp.loc[rows]['from'] - s_position <1000:
                                        if  int(e_position) > data_seg_temp.loc[rows]['from']:
                                            s_position = data_seg_temp.loc[rows]['from']
                                        #print(s_position)
                                else:
                                    break
                    else:
                        npeak_e = str(e_position)[-3:]
                        if npeak_e != '000':
                            e_len = len(str(e_position))
                            lnpeak_e = str(e_position)[:-3]
                            s_position = int(lnpeak_e.ljust(e_len, '0'))  
                            s_position = s_position + 1
                        else:
                            s_position = int(e_position) - 999
                        if gap_flag == '1':
                            #print(data_seg)
                            for rows in data_seg_temp.index:
                                #print(peak_chr)
                                #print(data_seg.loc[row]['filename'])
                                if peak_chr in data_seg_temp.loc[rows]['filename']:
                                    #print(data_seg)
                                    if 0< data_seg_temp.loc[rows]['from'] - s_position <1000:
                                        if  int(e_position) > data_seg_temp.loc[rows]['from']:
                                            s_position = data_seg_temp.loc[rows]['from']
                                        #print(s_position)
                                else:
                                    break
                    if len(new_data)== 13:
                        len_max = new_data[5]
                        len_sum = new_data[5]   
                    elif len(new_data) == 14:
                        #print(new_data)
                        #print(new_data[9])
                        if new_data[13] == '':
                            len_max = new_data[5]
                            len_sum = new_data[5]  
                            #print(len_max)
                        else:
                            len_max = new_data[13] 
                            len_sum = new_data[5] 
                    else:
                        print(len(new_data))
                        print(new_data)
                        print("Data error - please check the data")
                    dict_motif_type = {} 
                    a_dict_motif_type = {} 
                    #list_type=[]
                    new_type=[]
                    a_new_type =[]                    
                    new_ratio=[]
                    a_new_ratio=[]
                    type_classfy=[]
                    new_tr=[]
                    new_st = []
                    n =0
                    for temp_m in type_motif:
                        temp_m=temp_m.lstrip()
                        dict_motif_type[temp_m] = int(sum_count[n])/ratio_len_total
                        #list_type.append(dict_motif_type[temp_m])
                        n+=1
                    k =0
                    for temp_ma in a_type_motif:
                        temp_ma=temp_ma.lstrip()
                        a_dict_motif_type[temp_ma] = int(a_sum_count[k])/a_ratio_len_total
                        #list_type.append(dict_motif_type[temp_m])
                        k+=1    
                    list_type = sorted(dict_motif_type.items(), key=lambda d: d[1], reverse=True)
                    a_list_type = sorted(a_dict_motif_type.items(), key=lambda d: d[1], reverse=True)
                    #p=0.00
                    flag_a = 0
                    AG_p = 0.00
                    TC_p = 0.00
                    AG_ap = 0.00
                    TC_ap = 0.00
                    flag_size_AG = 0
                    flag_size_TC = 0
                    flag_count_AG = 0
                    flag_count_TC = 0
                    #print(list_type)
                    for temp_c in list_type:
                        type_e = temp_c[0]
                        type_p = temp_c[1]
                        if 'T' not in type_e and 'C' not in type_e:
                            AG_p += type_p
                        if 'A' not in type_e and 'G' not in type_e:
                            TC_p += type_p
                    if AG_p>=0.5:
                        # print(list_type)
                        # print('AG_P')
                        flag_size_AG = 1
                    if TC_p>=0.5:
                        # print(list_type)
                        # print('TC_P')
                        flag_size_TC = 1
                     
                    for temp_ca in a_list_type:
                        atype_e = temp_ca[0]
                        atype_p = temp_ca[1]
                        if 'T' not in atype_e and 'C' not in atype_e:
                            AG_ap += atype_p
                        if 'A' not in atype_e and 'G' not in atype_e:
                            TC_ap += atype_p
                    if AG_ap>=0.5:
                        flag_count_AG = 1
                        # print(a_list_type)
                        # print('AG_AP')
                    if TC_ap>=0.5:
                        flag_count_TC = 1
                        # print(a_list_type)
                        # print('TC_AP')
                        
                    for temp_a in a_list_type:
                        c_type = temp_a[0]
                        c_type_p = temp_a[1]
                        if c_type_p >= 0.667:
                            a_temp_type = f'[{c_type}]h'
                            flag_a = 1
                            break
                        elif c_type_p >= 0.5 and c_type_p < 0.667:
                            a_temp_type = f'[{c_type}]m'
                            flag_a = 1
                            break
                        elif c_type_p >= 0.334 and c_type_p < 0.5:
                            a_temp_type = f'[{c_type}]l'
                            flag_a = 1
                            break
                        elif flag_count_AG == 1:
                            a_temp_type = '[AGnt]'
                            flag_a = 1
                            break
                        elif flag_count_TC == 1:
                            a_temp_type = '[TCnt]'
                            flag_a = 1
                            break
                    for temp_r in list_type:
                         p = temp_r[1]
                         type_name_new = temp_r[0]
                         if p>=0.667:
                             #print("p>=0.667----------------------------------------")
                             #print(type_name_new)
                             temp_type = f'[{type_name_new}]h'
                             type_classfy.append(temp_type)
                             break
                             #print("------------------------------------------------")
                         elif p>=0.5 and p<0.667:
                             #print("p>=0.5----------------------------------------")
                             #print(type_name_new)
                             temp_type = f'[{type_name_new}]m'
                             type_classfy.append(temp_type)
                             break
                             #print("------------------------------------------------")
                         elif p >=0.334 and p<0.5:
                             temp_type = f'[{type_name_new}]l'
                             type_classfy.append(temp_type)
                             break
                         elif flag_size_AG == 1:
                             type_classfy.append('[AGnt]')
                             break
                         elif flag_size_TC == 1:
                             type_classfy.append('[TCnt]')
                             break
                         elif flag_a == 1:
                             type_classfy.append(a_temp_type)
                             break
                         else:
                             type_classfy.append('[MOTIFmix]')
                             break
                    m = 0
                    motif_new=[]
                    motif_sum=[]
                    #print(name_peak)
                    #print(type_ratio)
                    type_ratio = list(map(float,type_ratio)) 
                    #print(number)
                    type_ratio = sorted(type_ratio,key=float, reverse=True) 
                    for temp_ty in type_motif:
                        temp_ty=temp_ty.lstrip()
                        motif_new.append(temp_ty)
                        motif_sum.append(str(int(sum_count[m])))
                        ty_some = f'{temp_ty}({type_ratio[m]})' 
                        #print(ty_some)
                        new_ratio.append(str(type_ratio[m]))
                        new_tr.append(ty_some)
                        m+=1
                    ne =0
                    for new_te in  new_type:
                        tr_some =f'{new_te}({type_ratio[ne]})'
                        ne+=1
                        new_st.append(tr_some)
                    pf.write(name_peak+","+type_peak+","+str(s_position)+","+str(e_position)+","+peak_lable+","+";".join(motif_new)+","+";".join(motif_sum)+","+";".join(new_ratio)+","+";".join(type_classfy)+","+str(len_max)+","+str(len_sum)+","+";".join(a_type_motif)+","+";".join(a_sum_count)+","+";".join(a_type_ratio)+","+str(a_ratio_len_total)+","+name_peak+"("+";".join(type_classfy)+")"+"\n")
                   #Peak_name,E_position,Motif type,Motif size sum,Percentage(%),Integrated motif type,pD1RD,pD1RD sum

    data_all = pd.read_csv(lastfile+all_file+type_name)
    #print(data_all)
    """
    #Save HMDP
    """
    data_HP = data_all[~data_all['Peak_type'].str.contains('-MMDP|-LMDP')]
    data_HP = data_HP.drop(columns=['Peak_type'])
    #print(data_HP)
    #df1 = df.drop(labels='d', axis=1)
    data_HP.to_csv(lastfile+hpfile+HP_name,sep=',',index = False)
    print(HP_name + "---" + "File generated")
    """
    #Save MMDP
    """
    data_MP = data_all[~data_all['Peak_type'].str.contains('-HMDP|-LMDP')]
    data_MP = data_MP.drop(columns=['Peak_type'])
    #df1 = df.drop(labels='d', axis=1)
    data_MP.to_csv(lastfile+mpfile+MP_name,sep=',',index = False)
    print(MP_name+"---"+"File generated")
    """
    #Save LMDP
    """
    data_LP = data_all[~data_all['Peak_type'].str.contains('-MMDP|-HMDP')]
    data_LP = data_LP.drop(columns=['Peak_type'])
    #df1 = df.drop(labels='d', axis=1)
    data_LP.to_csv(lastfile+lpfile+LP_name,sep=',',index = False)
    print(LP_name+"---"+"File generated")
    
    
    

def xDMA_bedf(filen):
    pd_df = pd.read_csv(xDMAfile+filen)
    #print(pd_df)
    pd_df = pd_df.drop(['label', 'Motif type','Motif size sum','Percentage(%)','pD1RD sum','pD1RD','motif_abundance','count_sum','abundance_ratio(%)','count_sum_all','ucsc_name'],axis=1) #删除列
    
    dict_peakc = {} #chr name
    dict_peakt = {} # type colour
    for row in pd_df.index:
        peak_name = pd_df.loc[row]['Peak_name']
        motif_type = pd_df.loc[row]['Integrated motif type']
        peak_n = peak_name.split('-')

        if '0' in peak_n[1] and peak_n[1] != '10' and peak_n[1] != '20':
            chr_name = 'chr'+peak_n[1].split('0')[1]
        else:
            chr_name = 'chr' + peak_n[1]
        motif_n = motif_type.split(']')[0].split('[')[1]

        if 'nt' in motif_n or 'mix' in motif_n:
            if motif_n == 'AGnt': 
                type_colour = '255,0,0'
            elif motif_n == 'TCnt': 
                type_colour = '128,0,128'
            elif motif_n == 'MOTIFmix':
                type_colour = '169,169,169'
            else:
                print("Data error---1")
        else:
            if len(motif_n) == 2: 
                type_colour = '0,0,0'
            elif len(motif_n) == 3: 
                type_colour = '255,128,0'
            elif len(motif_n) == 1:
                type_colour = '0,255,255'
            elif len(motif_n) == 4: 
                type_colour = '255,255,0'
            elif len(motif_n) == 5: 
                type_colour = '0,255,0'
            elif len(motif_n) == 6:
                type_colour = '0,0,255'
            else:
                print("Data error---2")
        dict_peakc.setdefault('chr_name', []).append(chr_name)
        dict_peakt.setdefault('type_colour', []).append(type_colour)
    df_c = pd.DataFrame(dict_peakc)
    df_t = pd.DataFrame(dict_peakt)
    pd_df = pd.concat([pd_df,df_c,df_t],axis=1)
    pd_df['score'] = 0
    pd_df['strand'] = '.'
    pd_df['S_position'] = pd_df['S_position'] - 1
    pd_df['rs_position'] = pd_df['S_position']
    pd_df['rE_position'] = pd_df['E_position']
    pd_df['Integrated motif type'] = pd_df['Integrated motif type'].apply(lambda x: f"({x})")
    pd_df['peak_XDMA'] = pd_df['Peak_name']+pd_df['Integrated motif type']
    order = ['chr_name','S_position', 'E_position','peak_XDMA','score','strand','rs_position','rE_position','type_colour'] 
    pd_df = pd_df[order]
    pd_df.to_excel(bed_outfile+ filen.split('.csv')[0]+"_ucsc_bed.xlsx",index=False,engine="openpyxl")
    pd_df.to_csv(bed_outfile+ filen.split('.csv')[0]+".bed", sep='\t', index=False,header=False)
    print(filen.split('.csv')[0]+"_ucsc_bed.xlsx")
    print("------Ucsc bed file data has been generated------"+"\n")






def ucsc_regionf(file_name):
    xDMA_data = pd.read_csv(xDMAfile+file_name)
    xDMA_data['chr_name'] = xDMA_data['Peak_name']
    xDMA_data['s_pos'] = xDMA_data['S_position'] - 50000
    xDMA_data['e_pos'] = xDMA_data['S_position'] + 49999
    xDMA_data['check_value']  = xDMA_data["s_pos"] < 0
    list_reg = xDMA_data[(xDMA_data['check_value'] == True)].index.tolist() 
    if len(list_reg) != 0:
        for i_reg in list_reg:
            xDMA_data.loc[i_reg,('s_pos')] = 1
            xDMA_data.loc[i_reg,('e_pos')] = 100000
    xDMA_data['chr']= xDMA_data['chr_name'].str.split('-', expand=True)[0]
    xDMA_data['num'] = xDMA_data['chr_name'].str.split('-', expand=True)[1]
    xDMA_data['peak'] = xDMA_data['chr_name'].str.split('-', expand=True)[2]
    #xDMA_data = pd.concat([xDMA_data, xDMA_data['chr_name'].str.split('-', expand=True)], axis=1).drop('chr_name',axis=1)
    #print(xDMA_data.columns.values)
    xDMA_data.drop(['peak','check_value'],axis=1,inplace=True)
    xDMA_data['chr'] = 'chr'
    #print(xDMA_data)
    xDMA_data = xDMA_data.astype({'num': 'str'})
    xDMA_data['num'] = xDMA_data['num'].str.replace(r'^0','', regex=True) 
    #print(xDMA_data)
    xDMA_data['chr_name'] = xDMA_data['chr']+xDMA_data['num']
    xDMA_data = xDMA_data.drop(columns=['chr','num'])
    order = ['chr_name','s_pos', 'e_pos','Peak_name','S_position','E_position','label','Motif type','Motif size sum','Percentage(%)','Integrated motif type','pD1RD','pD1RD sum','motif_abundance','count_sum','abundance_ratio(%)','count_sum_all']
    xDMA_data = xDMA_data[order]
    xDMA_data.to_csv(ucsc_regfile+file_name.split('.csv')[0]+"_ucsc_region.csv", index=False)
    print(file_name.split('.csv')[0]+"_ucsc_region.csv"+"---generated")
        
    
    
peak_nc = config_control[0]["peak_name"] #Generate peak
type_dc = config_control[1]["motif_type"]   #Statistic type
type_mc = config_control[1]["type_mer"]  # Merge Peak
type_ic = config_control[1]["Integrated_type"]    #Type Naming
lp_nc = config_control[2]["lp_num"]  # Used to control the xMDP counting program
xdma_bc = config_control[2]["xDMA_bed"] #Used to control whether to generate bed files
ucsc_rc = config_control[2]["ucsc_region"] #Used to control whether to generate interval plots of xMDP peaks on ucsc
file_mc = config_control[2]["file_mer"] #Used for file merging, output to merged file


num_cores = multiprocessing.cpu_count()-4




if peak_nc == '1':
    Parallel(n_jobs=num_cores)(delayed(csv_peak)(temp_i["name"]) for temp_i in config_peak)
    print("\n"+"Peak initial name successful-1")

print("-------------------------------------------------------"+"\n")





if type_dc == "1":  #Statistic type
    Parallel(n_jobs=num_cores)(delayed(multi_proc)(peak_nt["name"]) for peak_nt in config_peak)
    print("\n"+"Type statistics completed---2")
    print("-----------------------------------------------------------------------------------------------------")
    print("\n")
        



if type_mc == "1": # Merge Peak
    Parallel(n_jobs=num_cores)(delayed(classify_pro)(peak_n["name"]+"_MP+"+"_type.csv") for peak_n in config_peak)
    print("\n"+"Merge peak completed---3")
    print("-----------------------------------------------------------------------------------------------------")
    print("\n")
 

if type_ic == "1": #Type Naming
    Parallel(n_jobs=num_cores)(delayed(count_sum_pro)(peak_np["name"]+"_MP+"+"_type_last.csv") for peak_np in config_peak)
    print("\n"+"Type naming completed---4")
    print("-----------------------------------------------------------------------------------------------------")
    print("\n")
    count = 0
    for peak_n in config_peak: 
        peak_tn=peak_n["name"]+"_MP+"+"_type_last.csv"
        name_hp = peak_tn.split("MP")[0] + "HP_count.csv"
        name_mp = peak_tn.split("MP")[0] + "MP_count.csv"
        name_lp = peak_tn.split("MP")[0] + "LP_count.csv"
        name_all = peak_tn.split("MP")[0] + "LP+_count.csv"
        species_name = peak_n["name"].split('-')[0]
        #print(name_all)
        all_namehp = species_name+"-"+"peak-HMDP_all.csv"
        all_namemp = species_name+"-"+"peak-MMDP_all.csv"
        all_namelp = species_name+"-"+"peak-LMDP_all.csv"
        all_name = species_name+"-"+"xDMA-_all.csv"
        data_mer_hp = pd.read_csv(lastfile+hpfile+name_hp,encoding=get_encoding(lastfile+hpfile+name_hp))
        data_mer_mp = pd.read_csv(lastfile+mpfile+name_mp,encoding=get_encoding(lastfile+mpfile+name_mp))
        data_mer_lp = pd.read_csv(lastfile+lpfile+name_lp,encoding=get_encoding(lastfile+lpfile+name_lp))
        data_mer_all = pd.read_csv(lastfile+all_file+name_all,encoding=get_encoding(lastfile+all_file+name_all))
        if flag_replace == '1':
            data_mer_hp['Peak_name'] = data_mer_hp['Peak_name'].str.replace(old_name, new_name)
            #print(data_mer_hp)
            data_mer_mp['Peak_name'] =data_mer_mp['Peak_name'].str.replace(old_name, new_name)
            data_mer_lp['Peak_name'] =data_mer_lp['Peak_name'].str.replace(old_name, new_name)
            data_mer_all['Peak_name'] =data_mer_all['Peak_name'].str.replace(old_name, new_name)
        if count == 0:
            if (os.path.exists(lastfile+xDMA+all_namehp)):
                os.remove(lastfile+xDMA+all_namehp)
            if (os.path.exists(lastfile+xDMA+all_namemp)):
                os.remove(lastfile+xDMA+all_namemp)
            if (os.path.exists(lastfile+xDMA+all_namelp)):
                os.remove(lastfile+xDMA+all_namelp)
            if (os.path.exists(lastfile+all_file+all_name)):
                os.remove(lastfile+all_file+all_name)
            data_mer_hp.to_csv(lastfile+xDMA+all_namehp,sep=',',index = False,mode ='a+')
            data_mer_mp.to_csv(lastfile+xDMA+all_namemp,sep=',',index = False,mode ='a+')
            data_mer_lp.to_csv(lastfile+xDMA+all_namelp,sep=',',index = False,mode ='a+')
            data_mer_all.to_csv(lastfile+all_file+all_name,sep=',',index = False,mode ='a+')
        else:
            data_mer_hp.to_csv(lastfile+xDMA+all_namehp,sep=',',index = False,mode ='a+',header = False)
            data_mer_mp.to_csv(lastfile+xDMA+all_namemp,sep=',',index = False,mode ='a+',header = False)  
            data_mer_lp.to_csv(lastfile+xDMA+all_namelp,sep=',',index = False,mode ='a+',header = False)
            data_mer_all.to_csv(lastfile+all_file+all_name,sep=',',index = False,mode ='a+',header = False)
        count=count+1
    print("\n"+"Merging xMDP files completed---5")
    print("-----------------------------------------------------------------------------------------------------")
    print("\n")
    
    

    
shp_nd =pd.DataFrame() 
hp_nd =pd.DataFrame() 
mp_nd =pd.DataFrame()   
lp_nd =pd.DataFrame() 

dict_shp = {}    
dict_hp ={} 
dict_mp ={} 
dict_lp = {}    


if lp_nc == '1': #HMDP count statistics
    lp_file = os.listdir(outfile)
    for i_file in lp_file:
        if '_MP+_type_last.csv' in i_file:
            chr_n = i_file.split('_MP+_type_last.csv')[0]
            data = pd.read_csv(outfile+i_file,encoding=get_encoding(outfile+i_file))
            shp_cn = data[data['Peak_type'].str.contains('-HMDP')].shape[0]
            #hp_cn = data[data['Peak_type'].str.contains('-HP')].shape[0]
            mp_cn = data[data['Peak_type'].str.contains('-MMDP')].shape[0]
            lp_cn = data[data['Peak_type'].str.contains('-LMDP')].shape[0]
            #print(lp_cn)
            dict_shp[chr_n] = shp_cn
            #dict_hp[chr_n] = hp_cn
            dict_mp[chr_n] = mp_cn
            dict_lp[chr_n] = lp_cn
    shp_nd = pd.DataFrame(list(dict_shp.items()),columns=['chr_name', 'HMDP_num'])
    #hp_nd = pd.DataFrame(list(dict_hp.items()),columns=['chr_name', 'HP_num'])
    mp_nd = pd.DataFrame(list(dict_mp.items()),columns=['chr_name', 'MMDP_num'])
    lp_nd = pd.DataFrame(list(dict_lp.items()),columns=['chr_name', 'LMDP_num'])
    shp_nd = pd.concat([shp_nd,mp_nd,lp_nd],axis=1)
    shp_nd.to_excel(lastfile +all_file+ chr_n.split('-')[0]+"_xDMA_COUNT-all.xlsx", index=False,engine="openpyxl")
    print("\n"+"xMDP count statistics completed---6")
    print("-----------------------------------------------------------------------------------------------------")
    print("\n")
    #print(lp_nd)
    
    
if xdma_bc == '1':
    xDMA_list = os.listdir(xDMAfile)  
    Parallel(n_jobs=num_cores)(delayed(xDMA_bedf)(XDMA_f) for XDMA_f in xDMA_list)
    print("\n"+"xMDP_ Ucsc bed file completed---7")
    print("-----------------------------------------------------------------------------------------------------")
    print("\n")




if ucsc_rc == '1':  #Used to generate ucsc image intervals
    xDMA_lf = os.listdir(xDMAfile)  #Obtain the file name under xDMAfile and return it as a list
    Parallel(n_jobs=num_cores)(delayed(ucsc_regionf)(temp_f) for temp_f in xDMA_lf)
    print("\n"+"xDMA_ucsc——Region file completed---8")
    print("-----------------------------------------------------------------------------------------------------")
    print("\n")
    
    
if file_mc == '1': #Used for merging files and summarizing views
    peak_lf = os.listdir(infile)
    motif_typef = os.listdir(outfile)
    head_namet = motif_typef[0].split('-')[0]
    for temp_fc in peak_lf:
        if '.csv' in temp_fc:
            head_namep = temp_fc.split('-')[0]
            break
    peak_file = head_namep + '-Microsatellites_density_peaks_ALL.csv'
    motif_tf = head_namet +'-peaks_motif_type_ALL.csv'
    if (os.path.exists(file_mers+peak_file)):
        os.remove(file_mers+peak_file)
    if (os.path.exists(file_mers+motif_tf)):
        os.remove(file_mers+motif_tf)
    count_pc = 0
    for temp_fc in peak_lf:
        if '.csv' in temp_fc:
            data_peaks = pd.read_csv(infile+temp_fc)
            if count_pc == 0:
               data_peaks.to_csv(file_mers+peak_file,sep=',',index = False,mode ='a+')
               count_pc +=1
            else:
               data_peaks.to_csv(file_mers+peak_file,sep=',',index = False,mode ='a+',header = False)
    print(peak_file+"---"+"merge completed"+"\n")
    count_mc = 0
    for temp_ft in motif_typef:
        if'_type.csv' in temp_ft:
            data_type = pd.read_csv(outfile+temp_ft,encoding=get_encoding(outfile+temp_ft))
            if count_mc == 0:
               data_type.to_csv(file_mers+motif_tf,sep=',',index = False,mode ='a+')
               count_mc +=1
            else:
               data_type.to_csv(file_mers+motif_tf,sep=',',index = False,mode ='a+',header = False)
    print(motif_tf + "---" + "File merge completed" + "\n")
    print("\n"+"File merge completed---9")
    print("-----------------------------------------------------------------------------------------------------")
    print("\n")

print("\n"+"Processing completed ----- Software version：V1.0")
input()