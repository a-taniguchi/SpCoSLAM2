#! /bin/sh
#Akira Taniguchi 2017/02/03-2018/11/25
#The path setting of the output data folder is necessary.

DATAPATH='/home/akira/Dropbox/SpCoSLAM/data/'

echo -n "trialname?(output_folder) >"
read trialname

#echo -n "dataset_number?(0:albert,1:MIT) >"
#read datasetNUM

#trialname=test3
datasetNUM=0

echo $trialname > $DATAPATH'trialname.txt'

echo $DATAPATH$trialname
mkdir $DATAPATH$trialname
mkdir $DATAPATH$trialname'/particle'
mkdir $DATAPATH$trialname'/weight'
mkdir $DATAPATH$trialname'/map'
mkdir $DATAPATH$trialname'/img'

SCAN=scan
gnome-terminal --command './run_roscore.sh'
sleep 5
gnome-terminal --command './run_gmapping.sh '$SCAN
gnome-terminal --command 'python ./map_saver.py '$trialname
gnome-terminal --command 'python ./run_rosbag.py '$trialname' '$datasetNUM
#cd ./learning
gnome-terminal --command 'python ./run_SpCoSLAM2.0.py '$trialname' '$datasetNUM

#python ./learning/learnSpCoSLAM2.py $trialname $datasetNUM


