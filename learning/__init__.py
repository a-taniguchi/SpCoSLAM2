#coding:utf-8
#パラメータ設定ファイル
#Akira Taniguchi 2017/01/18-2018/02/11-
import numpy as np

####################パラメータ####################
R = 30 #100        #パーティクル数(run_gmapping.shの値と合わせる必要あり)

UseFT = 1#1       #画像特徴を使う場合（１）、使わない場合（０）
UseLM = 1#1       #言語モデルを更新する場合（１）、しない場合（０）

CNNmode = 5     #CNN最終層1000次元(1)、CNN中間層4096次元(2)、PlaceCNN最終層205次元(3)、SIFT(0)

if CNNmode == 0:
  Descriptor = "SIFT_BoF"
  DimImg = 100  #画像特徴の次元数
elif CNNmode == 1:
  Descriptor = "CNN_softmax"
  DimImg = 1000 #画像特徴の次元数
elif CNNmode == 2:
  Descriptor = "CNN_fc6"
  DimImg = 4096 #画像特徴の次元数（高次元）
elif CNNmode == 3:
  Descriptor = "CNN_Place205"
  DimImg = 205  #画像特徴の次元数
elif CNNmode == 4:
  Descriptor = "hybridCNN"
  DimImg = 1183  #画像特徴の次元数
elif CNNmode == 5:
  Descriptor = "CNN_Place365"
  DimImg = 365  #画像特徴の次元数

dimx = 2                  #xtの次元数（x,y）

##初期(ハイパー)パラメータ
alpha0 = 20.00              #場所概念のindexのCRPハイパーパラメータ
gamma0 = 0.100              #位置分布のindexのCRPハイパーパラメータ
beta0 = 0.1 #0.2               #場所の名前Wのハイパーパラメータ
chi0  = 0.1 #0.2               #画像特徴のハイパーパラメータ
k0 = 1e-3                  #μのパラメータ
m0 = np.zeros(dimx)     #μのパラメータ
V0 = np.eye(dimx)*2 #*1000              #Σのパラメータ
n0 = 3.0 #2.0 #3.0                    #Σのパラメータ（＞次元数）
k0m0m0 = k0*np.dot(np.array([m0]).T,np.array([m0]))

#L = 100                  #場所概念の数50#
#K = 100                  #位置分布の数50#

#latticelmパラメータ
knownn = 3#2 #[2,3,4] #          #言語モデルのn-gram長 (3)
unkn   = 3#2 #[3,4] #            #綴りモデルのn-gram長 (3)
annealsteps  = 10#3 #[3, 5,10]    #焼き鈍し法のステップ数 (3)
anneallength = 15#5 #[5,10,15]    #各焼き鈍しステップのイタレーション数 (5)
burnin   = 100 #100     #burn-inのイタレーション数 (20)
samps    = 100 #100     #サンプルの回数 (100)
samprate = 100 #100     #サンプルの間隔 (1, つまり全てのイタレーション)
ramdoman = 0 #5 #0#1 #焼きなましパラメータをランダムにする（0:しない、0以外：最大値：各パラメータ値＋randoman）

#SpCoSLAM 2.0 追加要素
wic = 1  #1:wic重みつき、0:wic重みなし
LAG = 10 #10 #固定ラグ活性化のラグ値(it,ct) (LAG>=1; SpCoSLAM1.0:1)
LMLAG = 10 #LAG #固定ラグ活性化のラグ値(St) (LMLAG>=1), 固定ラグ活性化しない（LMLAG==0）
tyokuzen = 1 #直前のステップの言語モデルで音声認識（１）、ラグ値前の言語モデルで音声認識（0）
tyokuzenNPYLM = 0 ##NPYLMも固定ラグにする（バグがあるので使わない）
LMtype = "lattice_learn_NPYLM"  #"lattice" #latticelm:"lattice", NPYLM:"NPYLM"(未実装), lattice→learn→NPYLM:lattice_learn_NPYLM
LMweight = "WS"  #"weight" #wf*ws="weight", P(S{1:t}|c{1:t-1},α,β)/p(S{1:t}|β) = "WS"

#multiCPU for latticelm (a number of CPU)
multiCPU = 1 #2  #max(int(multiprocessing.cpu_count()/2)-1,2) )

##rosbag data playing speed (normal = 1.0)
rosbagSpeed = 1.0 #0.5#2

#latticelm or NPYLM (SpCoSLAM_NPYLM.shを実行)
#latticeLM = 1     #latticelm(1), NPYLM(0)
NbestNum = 10 #n-bestのｎをどこまでとるか（n<=10）


#パーティクルのクラス（構造体）
class Particle:
  def __init__(self,id,x,y,theta,weight,pid):
    self.id = id
    self.x = x
    self.y = y
    self.theta = theta
    self.weight = weight
    self.pid = pid
    #self.Ct = -1
    #self.it = -1


#Juliusパラメータ
#Juliusフォルダのsyllable.jconf参照
JuliusVer = "v4.4" #"v.4.3.1" #
HMMtype = "DNN"  #"GMM"
lattice_weight = "AMavg"  #"exp" #音響尤度(対数尤度："AMavg"、尤度："exp")
wight_scale = -1.0

if (JuliusVer ==  "v4.4"):
  Juliusfolder = "/home/akira/Dropbox/Julius/dictation-kit-v4.4/"
else:
  Juliusfolder = "/home/akira/Dropbox/Julius/dictation-kit-v4.3.1-linux/"

if (HMMtype == "DNN"):
  lang_init = 'syllableDNN.htkdic' 
else:
  lang_init = 'web.000.htkdic' # 'trueword_syllable.htkdic' #'phonemes.htkdic' # 初期の単語辞書（./lang_mフォルダ内）
lang_init_DNN = 'syllableDNN.htkdic' #なごり


####################ファイル####################
#パスはUbuntu使用時とWin使用時で変更する必要がある。特にUbuntuで動かすときは絶対パスになっているか要確認。
#win:相対パス、ubuntu:絶対パス
datafolder   = "/home/akira/Dropbox/SpCoSLAM/data/" #"./../datadump/" #

speech_folder = "/home/akira/Dropbox/Julius/directory/SpCoSLAM/*.wav" #*.wav" #音声の教示データフォルダ(Ubuntuフルパス)
speech_folder_go = "/home/akira/Dropbox/Julius/directory/SpCoSLAMgo/*.wav" #*.wav" #評価用の音声データフォルダ(Ubuntuフルパス)
lmfolder = "/home/akira/Dropbox/SpCoSLAM/learning/lang_m/"

datasetfolder = "/home/akira/Dropbox/SpCoSLAM/rosbag/"
dataset1 = "albert-b-laser-vision/albert-B-laser-vision-dataset/"
bag1 = "albertBimg.bag"
#dataset2 = "MIT_Stata_Center_Data_Set/"   ##用意できてない
#datasets = {"albert":dataset1,"MIT":dataset2}
datasets = [dataset1] #[dataset1,dataset2]
bags = [bag1] #run_rosbag.pyにて使用
scantopic = ["scan", "base_scan _odom_frame:=odom_combined"]
#datasetname = ""        #データセットフォルダ名（データセットを変更する場合、未使用）

#評価用正解データファイル
correct_Ct = 'Ct_correct.csv'  #データごとの正解のCt番号
correct_It = 'It_correct.csv'  #データごとの正解のIt番号
correct_data = 'SpCoSLAM_human.csv'  #データごとの正解の文章（単語列、区切り文字つき）(./data/)
correct_data_SEG = 'SpCoSLAM_SEG.csv'  #データごとの正解の文章（単語列、区切り文字つき）(./data/)
correct_name = 'name_correct.csv'  #データごとの正解の場所の名前（音素列）

N_best_number = 10  #PRR評価用のN-bestのN
margin = 10*0.05 # 地図のグリッドと位置の値の関係が不明のため(0.05m/grid)*margin(grid)=0.05*margin(m)
