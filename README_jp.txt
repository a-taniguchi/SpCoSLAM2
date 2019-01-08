# SpCoSLAM2.0

SpCoSLAM2.0の実装 (オンライン場所概念獲得、地図生成、語彙獲得の高精度化・高速化)  
SpCoSLAM2.0の論文の実験において使用したプログラム  

We propose an online learning algorithm based on a Rao-Blackwellized particle filter for spatial concept acquisition and mapping. We have proposed a nonparametric Bayesian spatial concept acquisition model (SpCoA). We propose a novel method (SpCoSLAM) integrating SpCoA and FastSLAM in the theoretical framework of the Bayesian generative model. The proposed method can simultaneously learn place categories and lexicons while incrementally generating an environmental map. 

【実行環境】  
Ubuntu　14.04  
Pythonバージョン2.7.6  
ROS indigo  
画像特徴抽出: Caffe (リファレンスモデル: Places-365 resnet152 推奨)  
音声認識: Julius dictation-kit-v4.4 GMM-HMM/DNN-HMM (日本語音節辞書使用、ラティス出力)  
未知語の語彙獲得ありの場合：latticelm 0.4, OpenFST  

IROS2017では、albert-B-laser-vision-datasetのrosbagファイルを使用  

【実行準備】  
・学習データセットのパス指定、トピック名を合わせるなど（__init__.py、run_gmapping.sh）  
・学習データセットの時刻情報より、教示時刻を保存したファイルを作成  
・音声データファイルを用意しておく。__init__.pyでファイルパス指定。  
・学習プログラム実行の前に、CNN_place.pyを実行。画像特徴のファイルフォルダを作成しておく。  
・パーティクル数の指定は、__init__.pyとrun_gmapping.shの両方変更する必要がある。  

【実行方法】  
cd ~/SpCoSLAM/learning  
./SpCoSLAM2.0.sh  
->trialname?(output_folder) >*output_folder_name* を入力

【注意事項】  
・run_rosbag.pyでたまにgflag関係のエラーがでることがあるが、ファイル読み込み失敗が原因。ほっておけば再読み込みしてくれて、動くので問題ない。  
・低スペックPCでは、gmappingの処理が追いつかずに地図がうまくできないことがある。  


このリポジトリにはgmappingが含まれます。  
./catkin_ws/src/フォルダ以下のファイルはgmappingのオリジナルバージョンのライセンス（License: CreativeCommons-by-nc-sa-2.0）に従います。  

---
このプログラムを使用したものを公開される場合は、必ず引用情報を明記してください。

Reference:　　
Akira Taniguchi, Yoshinobu Hagiwara, Tadahiro Taniguchi, and Tetsunari Inamura, "An Improved and Scalable Online Learning of Spatial Concepts and Language Models with Mapping", arXiv:1803.03481, 2018. (Preprint submitted)  

Original paper:  
https://arxiv.org/abs/1803.03481  


Sample video:  
https://youtu.be/TN081g15G84  
https://youtu.be/_6S-mNtjn44  


2019/01/08  Akira Taniguchi
