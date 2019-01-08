# SpCoSLAM2.0

SpCoSLAM2.0の実装 (オンライン場所概念獲得、地図生成、語彙獲得の高精度化・高速化)  
SpCoSLAM2.0の論文の実験において使用したプログラム  

【Abstract】
We propose a novel online learning algorithm, called SpCoSLAM 2.0, for spatial concepts and lexical acquisition with high accuracy and scalability. 
Previously, we proposed SpCoSLAM as an online learning algorithm based on unsupervised Bayesian probabilistic model that integrates multimodal place categorization, lexical acquisition, and SLAM. 
However, our previous algorithm had limited estimation accuracy owing to the influence of the early stages of learning, and increased computational complexity with added training data. 
Therefore, we introduce techniques such as fixed-lag rejuvenation to reduce the calculation time while maintaining an accuracy higher than that of the previous algorithm. 
The results show that, in terms of estimation accuracy, the proposed algorithm exceeds the previous algorithm and is comparable to batch learning. 
In addition, the calculation time of the proposed algorithm does not depend on the amount of training data and becomes constant for each step of the scalable algorithm. 
Our approach will contribute to the realization of long-term spatial language interactions between humans and robots.

【実行環境】  
Ubuntu　14.04  
Pythonバージョン2.7.6  
ROS indigo  
画像特徴抽出: Caffe (リファレンスモデル: Places-365 resnet152 推奨)  
音声認識: Julius dictation-kit-v4.4 GMM-HMM/DNN-HMM (日本語音節辞書使用、ラティス出力)  
未知語の語彙獲得ありの場合：latticelm 0.4, OpenFST  

IROS2017と同様に、albert-B-laser-vision-datasetのrosbagファイルを使用  

【実行準備】  
・ 学習データセットのパス指定、トピック名を合わせるなど（__init__.py、run_gmapping.sh）  
・ 学習データセットの時刻情報より、教示時刻を保存したファイルを作成  
・ 音声データファイルを用意しておく。__init__.pyでファイルパス指定。  
・ 学習プログラム実行の前に、CNN_place.pyを実行。画像特徴のファイルフォルダを作成しておく。  
・ パーティクル数の指定は、__init__.pyとrun_gmapping.shの両方変更する必要がある。  

【実行方法】  
cd ~/SpCoSLAM/learning  
./SpCoSLAM2.0.sh  
->trialname?(output_folder) > *output_folder_name* を入力

【注意事項】  
・run_rosbag.pyでたまにgflag関係のエラーがでることがあるが、ファイル読み込み失敗が原因。ほっておけば再読み込みしてくれて、動くので問題ない。  
・低スペックPCでは、gmappingの処理が追いつかずに地図がうまくできないことがある。  


このリポジトリにはgmappingが含まれます。  
./catkin_ws/src/フォルダ以下のファイルはgmappingのオリジナルバージョンのライセンス（License: CreativeCommons-by-nc-sa-2.0）に従います。  

---
このプログラムを使用したものを公開される場合は、必ず引用情報を明記してください。

Reference:　　
Akira Taniguchi, Yoshinobu Hagiwara, Tadahiro Taniguchi, and Tetsunari Inamura, 
"An Improved and Scalable Online Learning of Spatial Concepts and Language Models with Mapping", 
arXiv:1803.03481. (Preprint submitted)  

Original paper:  
https://arxiv.org/abs/1803.03481  


Sample video:  
https://youtu.be/TN081g15G84  
https://youtu.be/_6S-mNtjn44  


2019/01/08  Akira Taniguchi
