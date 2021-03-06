

【フォルダ】
/lamg_m/：Julus用の単語辞書ファイル（日本語音節辞書など）が入っているフォルダ
/torioki/：とりおきファイル


【ファイル】
CNN_place.py：CNNの特徴量をそのまま保存する.(PlaceCNN版)
Julius1best_gmm.py：主にPRR評価用に使用。Juliusで音声認識し、n-best音声認識結果を得る
JuliusLattice_gmm.py：学習プログラムと連携して、Juliusで音声認識しWFST音声認識結果を得る(GMM decoding)
JuliusLattice_dnn.py：学習プログラムと連携して、Juliusで音声認識しWFST音声認識結果を得る(DNN decoding)
JuliusLattice_gmm_LM.py：(パーティクルごとにLMを持つバージョン)
README.txt：このファイル
SpCoSLAM.sh：オンライン学習実行用基本ファイル
SpCoSLAM_LM.sh：(パーティクルごとにLMを持つバージョン)
SpCoSLAM_visualization.sh：自動で学習結果を逐次描画するためのシェルスクリプト（実際はこちらはあまり使ってない）・（下記参照[rviz上に位置分布を描画する方法]）
SpCoSLAMp1．sh：SpCoA側のパーティクル数１のとき用。(gmapping用のパーティクル数はまた別)
__init__.py：ハイパーパラメータなどの初期値やファイルパス指定用のファイル
artistanimation.py：CNNの結果のヒストグラムを可視化したものを動画化するためのプログラム（お試し用）
autovisualization.py：自動で学習結果を逐次描画するためのプログラム(保存は別でしなければならない)
collectmapclocktime.py：動画作成のために、地図の時刻をひとつのファイルにまとめて生成するためのもの。

evaluate.py:学習結果のデータを評価するプログラム（ARIc,ARIi,NMIc,NMIi,PARs,PARw,L,K,TSEG,ESEG）
evaluate2_MI.py:学習結果のデータを評価するプログラム（基本的にこれを使う）
（ARI(0埋め),NMI（0埋め）,PARw（MI重みづけあり）,PARw(MIなし),L,K）,分割単語数なし
evaluate2_ALL.py:学習結果のデータを評価するプログラム（基本的にこれを使う）全ての評価値を算出
（ARI(0埋め),NMI（0埋め）,PARw（MI重みづけあり）,PARw(MIなし),L,K,分割単語数,PARs）
evaluateMI.py:evaluate2_MI.pyの改良前（？）のプログラム。本質的には同じだが一部異なる（PARs使用、PARｗ重み付けなし不使用）
evaluatePARw.py:PARwのMIありなし比較用（evaluate2_MI.pyの下位互換）
evaluateSEG.py:(SEG：単語分割数比較用：場所の名前とそれ以外)
evaluateWORD.py：おそらく未使用
evaluate_PRR.py：PRR評価用プログラム（範囲はデータセットから決まる版）
evaluate_PRR2.py：PRR評価用プログラム（範囲指定版、IROS2017ではこちらを使用）
evaluate_PRR_MI.py：PRR評価用プログラム（範囲はデータセットから決まる版）+URによるMI効果のある重みづけ計算
evaluate_p(stlit).py：P(St|it)を推定するプログラム（語彙の定性的評価用）

gmapping.sh：FastSLAM実行用のシェルスクリプト
learnFastSLAM.py：FastSLAM用の必要最低限のデータファイル書き出しプログラム（？）
learnSpCoSLAM3.2.py：SpCoSLAM online learning program for IROS2017 (無駄なコメントアウトコードを省いたバージョン+bugfix)、かつ、SpCoAのオンライン版（SpCoSLAMから言語モデル更新と画像特徴を除けるバージョン）
learnSpCoSLAM3.py：(for one particle)
learnSpCoSLAM_LM.py：パーティクルごとにLM
learnSpCoSLAM2.0.py：SpCoSLAM online learning program for IROS2018
map_saver.py：地図を逐次的に保存していくプログラム（rospy使用）

meanEvaluate.py:評価値のファイルをまとめて平均と標準偏差を計算して保存するプログラム(0MI)
meanEvaluate2_MI.py：評価値のファイルをまとめて平均と標準偏差を計算して保存するプログラム(2MI)
meanEvaluate2_ALL2.0.py：評価値のファイルをまとめて平均と標準偏差を計算して保存するプログラム2.0(2ALL)
meanEvaluate_PRR.py：PRRの値をまとめて保存するプログラム
meanEvaluate_tsplot.py:時系列データのグラフ描画プログラム、各評価値のstep推移をグラフ化する。（evaluate2_MI用）パーティクル平均の値を取得する
(meanEvaluation0MI.csv->PARs,TSEG,ESEG, meanEvaluation2MI.csv->ARIc,ARIi,NMIc,NMIi,PARw,L,K meanEvaluationPRR_MI.csv->PRR, meanEvaluationSEG.csv->PSEG)
meanEvaluate_tsplot2.0.py:時系列データのグラフ描画プログラム、各評価値のstep推移をグラフ化する。（evaluate2_MI用）パーティクル平均の値を取得する 2.0
meanEvaluate_tsplot2.py:時系列データのグラフ描画プログラム、各評価値のstep推移をグラフ化する。(なにかのために仮使用したやつ)
meanEvaluate_tsplot3.py:時系列データのグラフ描画プログラム、各評価値のstep推移をグラフ化する。最大尤度のパーティクルの評価値を取得する
new_place_draw.py：学習した場所領域のサンプルをrviz上に可視化するプログラム(石伏（サンプルプロット）→磯辺（ガウス概形）→彰)
new_place_draw_online.py：オンライン可視化用
new_position_draw_online.py：ロボットの自己位置描画用
pcolor_small_W.py：学習結果の単語分布のカラーマップ描画
pcolor_small_WMI.py：学習結果の単語分布のカラーマップ描画（MI重み付け）
pcolor_small_phi2.py：位置分布のindexの多項分布のカラーマップ描画

run_FastSLAM.py：FastSLAMを実行するための子プログラム
run_SpCoSLAM.py：SpCoSLAMを実行するための子プログラム
run_SpCoSLAM_LM.py：SpCoSLAMを実行するための子プログラム(パーティクルごとにLMを持つバージョン)
run_SpCoSLAMp1.py：SpCoSLAMを実行するための子プログラム（パーティクル数１のとき用）
run_gmapping.sh：gmappingを実行するための子プログラム
run_mapviewer.py：map_serverコマンドを実行する小プログラム（未使用？）
run_mapviewer.sh：run_mapviewer.pyを実行するためのシェルスクリプト（未使用？）
run_rosbag.py：rosbagを実行するための子プログラム
run_roscore.sh：roscoreを実行するための子プログラム
saveSpCoMAP.rviz：rvizファイル
saveSpCoMAP_online.rviz：rvizファイル（オンライン描画用）



【実行準備】
・学習データセットのパス指定、トピック名を合わせるなど（__init__.py、run_gmapping.sh）
・学習データセットの時刻情報より、教示時刻を保存したファイルを作成（）
・音声データファイルを用意しておく。__init__.pyでファイルパス指定。
・学習プログラム実行の前に、CNN_place.pyを実行。画像特徴のファイルフォルダを作成しておく。
・パーティクル数の指定は、__init__.pyとrun_gmapping.shの両方変更する必要がある。

【実行方法】
cd ~/Dropbox/SpCoSLAM/learning
./SpCoSLAM.sh

【注意事項】
・run_rosbag.pyでtimespeedたまにgflag関係のエラーがでることがあるが、ファイル読み込み失敗が原因。callback関数なのでほっておけば再読み込みしてくれて、動くので問題ない。
・低スペックPCでは、gmappingの処理が追いつかずに地図がうまくできないことがある。
・learnSpCoSLAM**.pyの方のweightはws*wfのみ（単語辞書更新も）（wzはパーティクル近似）。gmappingの方でこれを読み込み、wzにかけている。

-----
[rviz上に位置分布を描画する方法]
Googleドライブのファイル参照。
roscore
rviz -d ~/Dropbox/SpCoSLAM/learning/saveSpCoMAP_online.rviz 
python ./autovisualization.py p30a20g10sfix008

個別指定の場合
cd ./data/
rosrun map_server map_server ./p30a20g10sfix008/map/map361.yaml
python ./new_place_draw.py p30a20g10sfix008 50 23 

-------------------------------------------------
更新日時
2017/02/12 Akira Taniguchi
2017/03/12 Akira Taniguchi
2018/02/07 Akira Taniguchi
