# SpCoSLAM 2.0

SpCoSLAM: Online Spatial Concept and Lexical Acquisition with Simultaneous Localization and Mapping  
Implementation of SpCoSLAM 2.0 (An Improved and Scalable Online Learning of Spatial Concepts and Language Models with Mapping)  
This repository includes the source codes used for the experiments in our paper [1].  
The source codes contains the original algorithm of SpCoSLAM [2].  

## Other repositories  
 [SpCoSLAM_Lets](https://github.com/EmergentSystemLabStudent/SpCoSLAM_Lets): ROS Wrapper of SpCoSLAM for real mobile robots  
 [SpCoSLAM](https://github.com/a-taniguchi/SpCoSLAM): Implementation of SpCoSLAM (Old version)   
 [SpCoSLAM_evaluation](https://github.com/a-taniguchi/SpCoSLAM_evaluation): The codes for the evaluation or the visualization in our paper  

## Abstract of SpCoSLAM
We propose an online learning algorithm based on a Rao-Blackwellized particle filter for spatial concept acquisition and mapping. We have proposed a nonparametric Bayesian spatial concept acquisition model (SpCoA). We propose a novel method (SpCoSLAM) integrating SpCoA and FastSLAM in the theoretical framework of the Bayesian generative model. The proposed method can simultaneously learn place categories and lexicons while incrementally generating an environmental map.   

## Abstract of SpCoSLAM 2.0
We propose a novel online learning algorithm, called SpCoSLAM 2.0, for spatial concepts and lexical acquisition with high accuracy and scalability. Previously, we proposed SpCoSLAM as an online learning algorithm based on unsupervised Bayesian probabilistic model that integrates multimodal place categorization, lexical acquisition, and SLAM. However, our previous algorithm had limited estimation accuracy owing to the influence of the early stages of learning, and increased computational complexity with added training data. Therefore, we introduce techniques such as fixed-lag rejuvenation to reduce the calculation time while maintaining an accuracy higher than that of the previous algorithm. The results show that, in terms of estimation accuracy, the proposed algorithm exceeds the previous algorithm and is comparable to batch learning. In addition, the calculation time of the proposed algorithm does not depend on the amount of training data and becomes constant for each step of the scalable algorithm. Our approach will contribute to the realization of long-term spatial language interactions between humans and robots.  

Figure: The graphical model of SpCoSLAM   
<img src="https://github.com/a-taniguchi/SpCoSLAM/blob/master/img/graphicalmodel02.jpg" width="520px">


## Execution environment
A new version will probably be available.  
- Ubuntu 14.04  
- Python 2.7.6  
- ROS indigo  
- CNN feature extracter: Caffe (Reference model: [Places-365 resnet152](http://places.csail.mit.edu/))  
- Speech recognition system: Julius dictation-kit-v4.4 GMM-HMM/DNN-HMM (Using Japanese syllabary dictionary, lattice output)  
- If you perform the lexical acquisition (unsupervised word segmentaiton): [latticelm 0.4](http://www.phontron.com/latticelm/) and OpenFST  

In our paper, we used a rosbag file of open-dataset [albert-B-laser-vision-dataset](https://dspace.mit.edu/handle/1721.1/62291).

## Preparation for execution  
- Path specification of training dataset, matching ros topic name etc (`__init__.py` and `run_gmapping.sh`)
- Create a file that stores the teaching time from the time information of the training dataset
- Prepare speech data files. Specify the file path in `__init__.py`  
- Start `CNN_place.py` before running the learning program  
  Create a folder for files of image features  
- To specify the number of particles, you need to change both `__ init__.py` and `run_gmapping.sh`  
- Change the path of the folder name in `/catkin_ws/src/openslam_gmapping/gridfastslam/gridslamprocessor.cpp`    
  We changed this file only.  
  [Note] If the original `gmapping` has already been installed on your PC, you need to change the uninstallation or path setting of `gmapping`.

## Execution procedure
`cd ~/SpCoSLAM2-master/learning `  
`./SpCoSLAM2.0.sh `  
`->trialname?(output_folder) >*output_folder_name* `  

## Notes
- Sometimes `gflag`-related errors sometimes appear in `run_rosbag.py`. 
  It is due to file reading failure. 
  It will reload and it will work so it will not be a problem.
- On low spec PCs, processing of gmapping can not catch up and maps can not be done well.

- This repository contains `gmapping`.
  The following files of `./catkin_ws/src/` folder follow the license of the original version of gmapping (License: CreativeCommons-by-nc-sa-2.0).

---
## Reference
If you use this program to publish something, please describe the following citation information.

[1] Akira Taniguchi, Yoshinobu Hagiwara, Tadahiro Taniguchi, and Tetsunari Inamura, "An Improved and Scalable Online Learning of Spatial Concepts and Language Models with Mapping", Autonomous Robots, Feb. 2020. (Open Access) [[LINK]](https://link.springer.com/article/10.1007/s10514-020-09905-0) [[arXiv]](https://arxiv.org/abs/1803.03481)  
[2] Akira Taniguchi, Yoshinobu Hagiwara, Tadahiro Taniguchi, and Tetsunari Inamura, "Online Spatial Concept and Lexical Acquisition with Simultaneous Localization and Mapping", IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS), 2017. [[LINK]](http://ieeexplore.ieee.org/document/8202243/) [[arXiv]](https://arxiv.org/abs/1704.04664)  

Sample videos:  
[Ver. SpCoSLAM(IROS2017)]  
 - https://youtu.be/z73iqwKL-Qk  
 
[Ver. SpCoSLAM 2.0]   
 - https://youtu.be/H5yztfmxGbc  
 - https://youtu.be/_6S-mNtjn44  


2018/05/19  Akira Taniguchi  
2020/04/10  Akira Taniguchi (update)  

=======
# SpCoSLAM 2.0

SpCoSLAM: Online Spatial Concept and Lexical Acquisition with Simultaneous Localization and Mapping  
Implementation of SpCoSLAM 2.0 (An Improved and Scalable Online Learning of Spatial Concepts and Language Models with Mapping)  
This repository includes the source codes used for the experiments in our paper [1].  
The source codes contains the original algorithm of SpCoSLAM [2].  

## Other repositories  
 [SpCoSLAM_Lets](https://github.com/EmergentSystemLabStudent/SpCoSLAM_Lets): ROS Wrapper of SpCoSLAM for real mobile robots  
 [SpCoSLAM](https://github.com/a-taniguchi/SpCoSLAM): Implementation of SpCoSLAM (Old version)   
 [SpCoSLAM_evaluation](https://github.com/a-taniguchi/SpCoSLAM_evaluation): The codes for the evaluation or the visualization in our paper  

## Abstract of SpCoSLAM
We propose an online learning algorithm based on a Rao-Blackwellized particle filter for spatial concept acquisition and mapping. We have proposed a nonparametric Bayesian spatial concept acquisition model (SpCoA). We propose a novel method (SpCoSLAM) integrating SpCoA and FastSLAM in the theoretical framework of the Bayesian generative model. The proposed method can simultaneously learn place categories and lexicons while incrementally generating an environmental map.   

## Abstract of SpCoSLAM 2.0
We propose a novel online learning algorithm, called SpCoSLAM 2.0, for spatial concepts and lexical acquisition with high accuracy and scalability. Previously, we proposed SpCoSLAM as an online learning algorithm based on unsupervised Bayesian probabilistic model that integrates multimodal place categorization, lexical acquisition, and SLAM. However, our previous algorithm had limited estimation accuracy owing to the influence of the early stages of learning, and increased computational complexity with added training data. Therefore, we introduce techniques such as fixed-lag rejuvenation to reduce the calculation time while maintaining an accuracy higher than that of the previous algorithm. The results show that, in terms of estimation accuracy, the proposed algorithm exceeds the previous algorithm and is comparable to batch learning. In addition, the calculation time of the proposed algorithm does not depend on the amount of training data and becomes constant for each step of the scalable algorithm. Our approach will contribute to the realization of long-term spatial language interactions between humans and robots.  

Figure: The graphical model of SpCoSLAM   
<img src="https://github.com/a-taniguchi/SpCoSLAM/blob/master/img/graphicalmodel02.jpg" width="520px">


## Execution environment
A new version will probably be available.  
- Ubuntu 14.04  
- Python 2.7.6  
- ROS indigo  
- CNN feature extracter: Caffe (Reference model: [Places-365 resnet152](http://places.csail.mit.edu/))  
- Speech recognition system: Julius dictation-kit-v4.4 GMM-HMM/DNN-HMM (Using Japanese syllabary dictionary, lattice output)  
- If you perform the lexical acquisition (unsupervised word segmentaiton): [latticelm 0.4](http://www.phontron.com/latticelm/) and OpenFST  

In our paper, we used a rosbag file of open-dataset [albert-B-laser-vision-dataset](https://dspace.mit.edu/handle/1721.1/62291).

## Preparation for execution  
- Path specification of training dataset, matching ros topic name etc (`__init__.py` and `run_gmapping.sh`)
- Create a file that stores the teaching time from the time information of the training dataset
- Prepare speech data files. Specify the file path in `__init__.py`  
- Start `CNN_place.py` before running the learning program  
  Create a folder for files of image features  
- To specify the number of particles, you need to change both `__ init__.py` and `run_gmapping.sh`  
- Change the path of the folder name in `/catkin_ws/src/openslam_gmapping/gridfastslam/gridslamprocessor.cpp`    
  We changed this file only.  
  [Note] If the original `gmapping` has already been installed on your PC, you need to change the uninstallation or path setting of `gmapping`.

## Execution procedure
`cd ~/SpCoSLAM2-master/learning `  
`./SpCoSLAM2.0.sh `  
`->trialname?(output_folder) >*output_folder_name* `  

## Notes
- Sometimes `gflag`-related errors sometimes appear in `run_rosbag.py`. 
  It is due to file reading failure. 
  It will reload and it will work so it will not be a problem.
- On low spec PCs, processing of gmapping can not catch up and maps can not be done well.

- This repository contains `gmapping`.
  The following files of `./catkin_ws/src/` folder follow the license of the original version of gmapping (License: CreativeCommons-by-nc-sa-2.0).

---
## Reference
If you use this program to publish something, please describe the following citation information.

[1] Akira Taniguchi, Yoshinobu Hagiwara, Tadahiro Taniguchi, and Tetsunari Inamura, "An Improved and Scalable Online Learning of Spatial Concepts and Language Models with Mapping", Autonomous Robots, Feb. 2020. (Open Access) [[LINK]](https://link.springer.com/article/10.1007/s10514-020-09905-0) [[arXiv]](https://arxiv.org/abs/1803.03481)  
[2] Akira Taniguchi, Yoshinobu Hagiwara, Tadahiro Taniguchi, and Tetsunari Inamura, "Online Spatial Concept and Lexical Acquisition with Simultaneous Localization and Mapping", IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS), 2017. [[LINK]](http://ieeexplore.ieee.org/document/8202243/) [[arXiv]](https://arxiv.org/abs/1704.04664)  

Sample videos:  
[Ver. SpCoSLAM(IROS2017)]  
 - https://youtu.be/z73iqwKL-Qk  
 
[Ver. SpCoSLAM 2.0]   
 - https://youtu.be/H5yztfmxGbc  
 - https://youtu.be/_6S-mNtjn44  


2018/05/19  Akira Taniguchi  
2020/04/10  Akira Taniguchi (update)  

