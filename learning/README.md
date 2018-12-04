# Learning programs

## 【Folder】  
/lamg_m/: A folder including word dictionary files (Japanese syllable dictionary)  


## 【Files】  
- `CNN_place.py`: Save image features extracting by CNN (ver. Places-CNN)  
- `CNN_place_multiCPU.py`: Save image features extracting by CNN (ver. Places-CNN)  (using multiple CPU)   
- `Julius1best_gmm.py`: Using for evaluation of place recognition rate (PRR). The speech files are recognized by the speech recognition Julius. It can get n-best speech recognition results.  (version HMM-GMM)  
- `Julius1best_dnn.py`: Using for evaluation of place recognition rate (PRR). The speech files are recognized by the speech recognition Julius. It can get n-best speech recognition results.  (version HMM-DNN)  
- `JuliusLattice_gmm.py`: The speech files are recognized by the speech recognition Julius. It can get WFST speech recognition results for learning programs.  (version HMM-GMM)  
- `JuliusLattice_dnn.py`: The speech files are recognized by the speech recognition Julius. It can get WFST speech recognition results for learning programs.  (version HMM-DNN)  
- `JuliusLattice_dec.py`: The speech files are recognized by the speech recognition Julius. It can get WFST speech recognition results for learning programs. (version changeable to HMM-GMM or HMM-DNN by setting of __init__.py)  
- `README.txt`: This file  
- `SpCoSLAM2.0.sh`: A main execution script for online learning of SpCoSLAM 2.0  
- `__init__.py`: A file for setting file paths and initial hyper-parameters  
- `autovisualization.py`: A program for automatically drawing learning results sequentially
(Save can be done with screenshots etc.)
- `collectmapclocktime.py`: The times on the generated map files in one file collectively. For creating a movie.
- `gmapping.sh`: Shell script for FastSLAM
- `learnSpCoSLAM2.0.py`: SpCoSLAM 2.0 online learning program (By setting of `__init__.py`, using language model update, using image features, imporved algorithm, and scalable algorithm can be chainged in SpCoSLAM 2.0)
- `map_saver.py`: Saving a environmental map (using rospy)
- `new_place_draw.py`: Visualization of position distributions (Gaussian distributions) on rviz 
- `new_place_draw_online.py`: For online visualization of learning result
- `new_position_draw_online.py`: For visualization of the robot position
- `run_SpCoSLAM2.0.py`: Sub-program for performing SpCoSLAM 2.0  
- `run_gmapping.sh`: Sub-program for performing gmapping
- `run_mapviewer.py`: Sub-program for performing map_server command (not used?)
- `run_mapviewer.sh`: Shell script for performing run_mapviewer.py (not used?)
- `run_rosbag.py`: Sub-program for performing rosbag
- `run_roscore.sh`: Sub-program for performing roscore
- `saveSpCoMAP.rviz`: Rviz setting file
- `saveSpCoMAP_online.rviz`: Rviz setting file for online visualization


-----
[How to visualize the position distributions on rviz]  
`roscore`  
`rviz -d ./*/SpCoSLAM/learning/saveSpCoMAP_online.rviz `  
`python ./autovisualization.py p30a20g10sfix008`  

In case of individual specification  
`rosrun map_server map_server ./p30a20g10sfix008/map/map361.yaml`  
`python ./new_place_draw.py p30a20g10sfix008 50 23 `  
 
-------------------------------------------------  
Updated date  
2017/02/12, 03/12, 
2018/01/12, 04/24 Akira Taniguchi  (SpCoSLAM)  
2018/05/19 Akira Taniguchi  (SpCoSLAM 2.0)  
2018/12/04 Akira Taniguchi  (update)  
