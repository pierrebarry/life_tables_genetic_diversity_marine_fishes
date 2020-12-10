# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 15:15:46 2019

@author: pierr
"""

import numpy as np
import json
import sys

## ALL Samples
SPECIES=['Spilc','Sspra','Eencr','Hgutt','Mmerl','Scabr','Dlab','Msurm',
         'Lmorm','Dpunt','Peryt','Cjuli','Ssard','Pminu','Cgale','Aboye',
         'Afall','Aminu','Gnige','Scine','Lbude','Styph','Scant','Tdela']
LOCA=['Li','Mu','Fa','Ga']
NUMBER=['1','2','3','4','5','6']

SAMPLES=[]
for i in range(24):
    for j in range(4):
        for k in range(6):
            SAMPLES+=[SPECIES[i]+LOCA[j]+NUMBER[k]]

with open('/share/tycho_poolz1/pagagnaire/COGEDIV/pbarry/'+str(sys.argv[1])+'/fastp/fastp_report_'+str(sys.argv[2])+'.json') as json_fastp:
    a = json.load(json_fastp)

#with open(str(sys.argv[3]) as json_fastp:
#     a = json.load(json_fastp

result=[]
# get sample
#c=a['command']
#for i in range(24*4*6):
#    if c.find(SAMPLES[i])>0:
#        sample=c[73:81]
#        break

#if sample[0:1]=="Dl":
#    c=a['command']
#    for i in range(24*4*6):
#        if c.find(SAMPLES[i])>0:
#            sample=c[73:80]
 #           break

result+=[str(sys.argv[2])]
## total reads_before_filtering
result+=[a['summary']['before_filtering']['total_reads']]
result+=[a['summary']['before_filtering']['total_bases']]
result+=[a['summary']['before_filtering']['gc_content']]
result+=[a['summary']['before_filtering']['q20_rate']]
result+=[a['summary']['before_filtering']['q30_rate']]

## filtering result

reads_inital=a['summary']['before_filtering']['total_reads']
base_inital=a['summary']['before_filtering']['total_bases']
result+=[a['filtering_result']['passed_filter_reads']/reads_inital*100]
result+=[a['filtering_result']['corrected_reads']/reads_inital*100]
result+=[a['filtering_result']['corrected_bases']/base_inital*100]
result+=[a['filtering_result']['low_quality_reads']/reads_inital*100]
result+=[a['filtering_result']['too_many_N_reads']/reads_inital*100]
result+=[a['filtering_result']['too_short_reads']/reads_inital*100]
result+=[a['filtering_result']['low_complexity_reads']/reads_inital*100]

## total after_filtering

result+=[a['summary']['after_filtering']['total_reads']]
result+=[a['summary']['after_filtering']['total_bases']]
result+=[a['summary']['after_filtering']['gc_content']]
result+=[a['summary']['after_filtering']['q20_rate']]
result+=[a['summary']['after_filtering']['q30_rate']]
result+=[a['duplication']['rate']]

with open('fastp/Summary_fastp.txt',"a") as f:
	 np.savetxt(f, [result],newline='\n', fmt='%s', delimiter=";")
     
         
    
    
    
