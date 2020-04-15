#coding:utf-8
# e-nco-ding: shi-ft_jis
#c-od-i-n-:-ut-f---8
#!/usr/bin/env python
#時系列データのグラフ描画プログラム、各評価値のstep推移をグラフ化する。（evaluate2ALL2.0用）
#パーティクル平均の値を取得する
#Akira Taniguchi (2017/02/23-2018/02/22-2018/09/13-2019/05/19)
#保存ファイル名（プログラム最後の方）と読み込みファイル名の指定に注意
#AUROのリバイズのためのコード（比較手法数を減らしたバージョン）

import sys
import string
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from __init__ import *
#import math

#sns.set(style="darkgrid")
sns.set_style("whitegrid", {'grid.linestyle': '--'})
current_palette = sns.color_palette()
#sns.set_palette("muted")
#sns.color_palette("muted")
sns.color_palette(current_palette)

step = data_step_num
#DATA = ['(A) SpCoSLAM','(B) SpCoSLAM with AW+WS','(C) 2.0 (10 FLR)','(D) 2.0 (10 FLR + RS)','(E1) 2.0 (1 FLR + SBU)','(E2) 2.0 (10 FLR + SBU)','(E3) 2.0 (20 FLR + SBU)']
DATA = ['(A)','(B)','(C)','(D)','(E1)','(E2)','(E3)']

data1 ='SpCoA'
HYOUKA = ['NMIc','NMIi',"The number of spatial concepts","The number of position distributions","The number of segmented words","PAR (Word)","PAR (Sentence)","ARIc","ARIi","Calculation time [sec.]","PRR","EAR $L$","EAR $K$"] #"Estimation error of $L$","Estimation error of $K$"]#
datasetNUM = 0
datasetname = "" #datasets[int(datasetNUM)]
learningdata = ["A3LDK","B3LDK","C3LDK","D3LDK","E13LDK","E23LDK","E33LDK"]
#["p30a20g10sfix","p30a20g10nfnlsfix","p1a20g10s","batch100a20g10s"]


#trialname = raw_input("data_name?(**_m???_NUM) > ")
hs = raw_input("Evaluation? [NMIc(0),NMIi(1),L(2),K(3),SEG(4),PARwMI(5),PARs(6),ARIc(7),ARIi(8),Time(9),PRR(10),EAR_L(11),EAR_K(12)]> ")
h = int(hs)
data_num1 = raw_input("data_start_num?(DATA_**) > ")
data_num2 = raw_input("data_end_num?(DATA_**) > ")
N = int(data_num2) - int(data_num1) +1
#filename = raw_input("Read_Ct_filename?(.csv) >")
S = int(data_num1)

ARIc_M = [0.0 for c in xrange(len(DATA)*step*N)]
ARIi_M = [0.0 for c in xrange(len(DATA)*step*N)]
NMIc_M = [0.0 for c in xrange(len(DATA)*step*N)]
NMIi_M = [0.0 for c in xrange(len(DATA)*step*N)]
PARs_M = [0.0 for c in xrange(len(DATA)*step*N)]
PARw_M = [0.0 for c in xrange(len(DATA)*step*N)]
EL_M   = [0.0 for c in xrange(len(DATA)*step*N)]
TL_M   = [0.0 for c in xrange(step)]
EK_M   = [0.0 for c in xrange(len(DATA)*step*N)]
TK_M   = [0.0 for c in xrange(step)]
ESEG_M = [0.0 for c in xrange(len(DATA)*step*N)]
TSEG_M = [0.0 for c in xrange(step)]
PSEG_M = [0.0 for c in xrange(step)]
PRR_M = [0.0 for c in xrange(len(DATA)*step*N)]
TIME_M = [0.0 for c in xrange(len(DATA)*step*N)]
EAR_L_M   = [0.0 for c in xrange(len(DATA)*step*N)]
EAR_K_M   = [0.0 for c in xrange(len(DATA)*step*N)]
#PAR = [0.0]*len(DATA)*step*N  #np.random.normal(0.5, 0.5, 500)
i = 0
#d = 0
for d in xrange(len(DATA)):
  trialname = learningdata[d]
  for s in range(N):
    i = 0
    #if (d == 0 or 2):
    #  data_name = 'TAMD2_sig'
    for line in open(datafolder + trialname + str(s+1).zfill(3) + '/' + trialname + str(s+1).zfill(3) +'_meanEvaluation2ALL.csv', 'r'):
      itemList = line[:-1].split(',')
      if (i != 0) and (itemList[0] != '') and (i <= step):
        #print i,itemList
        ARIc_M[d*N + (i-1)*len(DATA)*N + s] = float(itemList[0])
        ARIi_M[d*N + (i-1)*len(DATA)*N + s] = float(itemList[1])
        NMIc_M[d*N + (i-1)*len(DATA)*N + s] = float(itemList[2])
        NMIi_M[d*N + (i-1)*len(DATA)*N + s] = float(itemList[3])
        PARs_M[d*N + (i-1)*len(DATA)*N + s] = float(itemList[10])
        PARw_M[d*N + (i-1)*len(DATA)*N + s] = float(itemList[4])
        #PARw_M[s] = PARw_M[s] + [float(itemList[5])]
        EL_M[d*N + (i-1)*len(DATA)*N + s]   = float(itemList[6])
        EK_M[d*N + (i-1)*len(DATA)*N + s]   = float(itemList[7])
        TSEG_M[i-1] = float(itemList[8])
        ESEG_M[d*N + (i-1)*len(DATA)*N + s] = float(itemList[9])
        #PAR[d*N + (i-1)*len(DATA)*N + s] =  float(itemList[3])
      i = i + 1
    """
    i = 0
    for line in open(datafolder + trialname + str(s+1).zfill(3) + '/' + trialname + str(s+1).zfill(3) +'_EvaluationPRR2.csv', 'r'): #_meanEvaluationPRR2.csv
      itemList = line[:-1].split(',')
      if (i != 0) and (itemList[0] != '') and (i <= step):
        PRR_M[d*N + (i-1)*len(DATA)*N + s] = float(itemList[0])
      i = i + 1
    """
    if(h == 9):
      i = 1
      for line in open(datafolder + trialname + str(s+1).zfill(3) + '/time_step.txt', 'r'): #_meanEvaluationPRR2.csv
        itemList = line[:-1].split(',')
        if (itemList[0] != '') and (i <= step):
          TIME_M[d*N + (i-1)*len(DATA)*N + s] = float(itemList[1])
        i = i + 1
    """
    i = 0
    for line in open(datafolder + '/Evaluation/PSEG.csv', 'r'):
      itemList = line[:-1].split(',')
      if (i != 0) and (itemList != '') and (i <= step):
        PSEG_M[i-1] = float(itemList[0])
      i = i + 1
    """
    


H_M = [NMIc_M,NMIi_M,EL_M,EK_M,ESEG_M,PARw_M,PARs_M,ARIc_M,ARIi_M,TIME_M,PRR_M,EAR_L_M,EAR_K_M]
#h = 0
#for hyouka in HYOUKA:
if (1):
  hyouka = HYOUKA[int(h)]
  iteration = []
  for i in range(step):
    iteration = iteration + [i+1 for j in range(N*len(DATA))] #+ [2 for i in range(10*len(DATA))] + [3 for i in range(10*5)] + [4 for i in range(10*5)] + [5 for i in range(10*5)] + [6 for i in range(10*5)] + [7 for i in range(10*5)] + [8 for i in range(10*5)] + [9 for i in range(10*5)] + [10 for i in range(10*5)]
  method = []#[DATA[i] for k in range(N) for i in range(len(DATA)) for j in range(step)]
  for p in range(step):
    for d in range(len(DATA)):
      method = method + [DATA[d] for i in range(N)] #+ [DATA[1] for j in range(N)] + [DATA[2] for k in range(N)] #+ [DATA[3] for l in range(10)] #+ ['F' for l in range(10)]
  
  subject = [i for i in range(N)]*len(DATA)*step
  
  
  """
  ####################
  ##SpCoAの結果はstep=50のみ(未完成、読み込み未実装)
  iteration = iteration + [int(len(DATA)+1)]
  method = method + [data1]
  subject = subject + [i for i in range(N)]
  ####################
  """
  
  
  if  (h == 2): #L
    #TL_Mを読み込む
    i = 0
    for line in open( datasetfolder + datasetname + 'Lnum.csv', 'r'):
      itemList = line[:].split(',')
      #print itemList
      TL_M = [int(itemList[j]) for j in range(step)]
      i = i + 1
    
    H_M[h] = H_M[h] + TL_M
    iteration = iteration + [i+1 for i in range(step)]
    method = method + ["True" for i in range(step)]
    subject = subject + [0 for i in range(step)]
    
  elif(h == 3): #K
    #TK_Mを読み込む
    i = 0
    for line in open( datasetfolder + datasetname + 'Knum.csv', 'r'):
      itemList = line[:].split(',')
      TK_M = [int(itemList[j]) for j in range(step)]
      i = i + 1
    
    H_M[h] = H_M[h] + TK_M
    iteration = iteration + [i+1 for i in range(step)]
    method = method + ["True" for i in range(step)]
    subject = subject + [0 for i in range(step)]
    
  elif(h == 4): #SEG
    #TSEG_M,PSEG_Mは読み込んである
    H_M[h] = H_M[h] + [TSEG_M[j] for j in range(len(TSEG_M))] + [PSEG_M[j] for j in range(len(PSEG_M))]
    
    iteration = iteration + [i+1 for i in range(step)] + [i+1 for i in range(step)]
    method = method + ["Morpheme" for i in range(step)] + ["Phrase" for i in range(step)]
    subject = subject + [0 for i in range(step)] + [0 for i in range(step)]
    
  elif(h == 11):
    #TL_Mを読み込む
    i = 0
    for line in open( datasetfolder + datasetname + 'Lnum.csv', 'r'):
      itemList = line[:].split(',')

      #print itemList
      TL_M = [int(itemList[j]) for j in range(step)]
      i = i + 1
    print len(TL_M),TL_M
    #EARの計算
    for i in range(step):
      for s in range(N):
        for d in range(len(DATA)):
          EAR_L_M[d*N + (i)*len(DATA)*N + s] =max( 1.0 - abs( EL_M[d*N + (i)*len(DATA)*N + s] - TL_M[i] )/TL_M[i], 0)
    #EL_M[d*N + (i-1)*len(DATA)*N + s] abs( EL_M[d*N + (i)*len(DATA)*N + s] - TL_M[i] ) #
    H_M[h] = EAR_L_M
  
  elif(h == 12):
    #TK_Mを読み込む
    i = 0
    for line in open( datasetfolder + datasetname + 'Knum.csv', 'r'):
      itemList = line[:].split(',')
      TK_M = [int(itemList[j]) for j in range(step)]
      i = i + 1
    #EARの計算
    for i in range(step):
      for s in range(N):
        for d in range(len(DATA)):
          EAR_K_M[d*N + (i)*len(DATA)*N + s] =max( 1.0 - abs( EK_M[d*N + (i)*len(DATA)*N + s] - TK_M[i] )/TK_M[i], 0)
    #EL_M[d*N + (i-1)*len(DATA)*N + s] abs( EK_M[d*N + (i)*len(DATA)*N + s] - TK_M[i] ) #
    H_M[h] = EAR_K_M
  
  data = {'step':iteration, 'method':method, 'subject':subject, hyouka:H_M[h]}
  #data2 = {'timepoint':iteration, 'ROI':method, 'subject':subject,'BOLD signal':PAR}
  
  df2 = pd.DataFrame(data)
  
  # Plot the response with standard error
  #sns.tsplot(data=gammas, time="timepoint", unit="subject",condition="ROI", value="BOLD signal"), estimator=np.median
  #sns.tsplot(data=df22, time="timepoint", unit="subject",condition="ROI", value="BOLD signal"), err_style="ci_bars", err_kws={"alpha": .3}
  markers = ['^','v','o','D','s','H','>','<','d','X']
  #print len(iteration),len(subject),len(method),len(H_M[h])
  #print iteration
  #print subject
  #print method
  #print H_M[h]
  AAA = sns.tsplot(data=df2, time="step", unit="subject",condition="method", value=hyouka)  
  plt.subplots_adjust(bottom=0.15, top=0.9, wspace=None, hspace=None)
  #left=0.60, , right=0.85, top=0.9800
  
  #A = np.array([[0,1,0,1,0,1,0,1], [1,0,1,0,1,0,1,0], [.5,.5,.5,.5,.5,.5,.5,.5]]), ci=10
  for i in range(len(DATA)):
    AAA.lines[i].set_marker(markers[i])
  location = 'upper left'
  if (h == 11 or h == 12):
    location = 'lower right'
  AAA.legend(loc=location,ncol=1,fontsize=16,labelspacing = 0.4) #prop={'size':10})#title='method','upper left'
  #sns.tsplot(A, ci=50)
  #print np.std([0,1,0.5])
  if (h == 2 or h == 3):
    plt.ylim([0.0,14])#max(H_M[h])])
  elif (h == 4 or h == 9):
    plt.ylim([0.0,max(H_M[h])])
  elif (h == 5):
    plt.ylim([0.0,0.5])
  elif (h == 6):
    plt.ylim([0.0,1.0])
  elif (h == 11 or h == 12):
    plt.ylim([0.0,1.0]) #5
  else:
    plt.ylim([0.0,1.0])
  plt.xticks(fontsize=16)
  plt.yticks(fontsize=16)
  plt.ylabel(hyouka,fontsize=18)
  plt.xlabel("step",fontsize=18)
  
  print df2
  #df2.to_csv("./text"+ hyouka+".csv")
  
  ######type 1 font#####
  plt.rcParams['ps.useafm'] = True
  plt.rcParams['pdf.use14corefonts'] = True
  #plt.rcParams['text.usetex'] = True
  
  #fig = AAA.get_figure()
  #plt.savefig(datafolder + '/Evaluation/' + trialname + '_' + data_num1 + '_' + data_num2 +  '_' + hyouka + '2_ALL2.0.eps', dpi=300)#, transparent=True #eps wa hanntoumei ga muri.
  plt.savefig(datafolder_evaluation + '/Evaluation/' + trialname + '_' + data_num1 + '_' + data_num2 +  '_' + hyouka + '2_ALL2.0sAURO.pdf', dpi=300)#, transparent=True
  plt.savefig(datafolder_evaluation + '/Evaluation/' + trialname + '_' + data_num1 + '_' + data_num2 +  '_' + hyouka + '2_ALL2.0sAURO.png', dpi=300)#, transparent=True
  #h = h+1
  
#plt.show()
print "close."

#fp.close()
