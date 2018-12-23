# SpCoSLAM 2.0

Implementation of SpCoSLAM 2.0 (An Improved and Scalable Online Learning of Spatial Concepts and Language Models with Mapping)  
This repository includes the source codes used for the experiments in our paper.  

【Other repositories】  
 [SpCoSLAM_Lets](https://github.com/EmergentSystemLabStudent/SpCoSLAM_Lets): Wrapper of SpCoSLAM for mobile robots (Recommended)  
 [SpCoSLAM](https://github.com/a-taniguchi/SpCoSLAM): Implementation of SpCoSLAM (Online Spatial Concept and Lexical Acquisition with Simultaneous Localization and Mapping)   
 [SpCoSLAM_evaluation](https://github.com/a-taniguchi/SpCoSLAM_evaluation): The codes for the evaluation or the visualization in our paper  

## Abstract of SpCoSLAM 2.0
We propose a novel online learning algorithm, SpCoSLAM 2.0 for spatial concepts and lexical acquisition with higher accuracy and scalability.
In previous work, we proposed SpCoSLAM as an online learning algorithm based on the Rao-Blackwellized particle filter.
The proposed method can simultaneously learn place categories and lexicons while incrementally generating an environmental map. 
However, this conventional algorithm had problems such as the decrease of the estimation accuracy due to the influence of the early stages of learning as well as the increase of the computational complexity with the increase of the training data.
Therefore, we first developed an improved algorithm by introducing new techniques such as rejuvenation.
Next, we developed a scalable algorithm to reduce the calculation time while maintaining a higher accuracy than the conventional algorithm.

Figure: The graphical model of SpCoSLAM   
<img src="https://github.com/a-taniguchi/SpCoSLAM/blob/master/img/graphicalmodel02.jpg" width="520px">


## 【Execution environment】  
- Ubuntu 14.04  
- Python 2.7.6  
- ROS indigo  
- CNN feature extracter: Caffe (Reference model: [Places-365 resnet152](http://places.csail.mit.edu/))  
- Speech recognition system: Julius dictation-kit-v4.4 GMM-HMM/DNN-HMM (Using Japanese syllabary dictionary, lattice output)  
- If you perform the lexical acquisition (unsupervised word segmentaiton): [latticelm 0.4](http://www.phontron.com/latticelm/) and OpenFST  

In our paper of IROS2018, we used a rosbag file of open-dataset [albert-B-laser-vision-dataset](https://dspace.mit.edu/handle/1721.1/62291).

## 【Preparation for execution】  
- Path specification of training dataset, matching ros topic name etc (`__init__.py` and `run_gmapping.sh`)
- Create a file that stores the teaching time from the time information of the training dataset
- Prepare speech data files. Specify the file path in `__init__.py`  
- Start `CNN_place.py` before running the learning program  
  Create a folder for files of image features  
- To specify the number of particles, you need to change both `__ init__.py` and `run_gmapping.sh`  
- Change the path of the folder name in `/catkin_ws/src/openslam_gmapping/gridfastslam/gridslamprocessor.cpp`    
  We changed this file only.  
  [Note] If the original `gmapping` has already been installed on your PC, you need to change the uninstallation or path setting of `gmapping`.

## 【Execution procedure】
`cd ~/SpCoSLAM2-master/learning `  
`./SpCoSLAM2.0.sh `  
`->trialname?(output_folder) >*output_folder_name* `  

## 【Notes】
- Sometimes `gflag`-related errors sometimes appear in `run_rosbag.py`. 
  It is due to file reading failure. 
  It will reload and it will work so it will not be a problem.
- On low spec PCs, processing of gmapping can not catch up and maps can not be done well.

- This repository contains `gmapping`.
  The following files of `./catkin_ws/src/` folder follow the license of the original version of gmapping (License: CreativeCommons-by-nc-sa-2.0).

---
If you use this program to publish something, please describe the following citation information.

Akira Taniguchi, Yoshinobu Hagiwara, Tadahiro Taniguchi, and Tetsunari Inamura, "SpCoSLAM 2.0: An Improved and Scalable Online Learning of Spatial Concepts and Language Models with Mapping", arXiv:1803.03481, 2018. (Preprint submitted)  

Original paper:  
https://arxiv.org/abs/1803.03481  


Sample video:  
https://youtu.be/TN081g15G84  
https://youtu.be/_6S-mNtjn44  

2018/05/19  Akira Taniguchi  
2018/11/22  Akira Taniguchi (update)  
2018/12/23  Akira Taniguchi (update)  

