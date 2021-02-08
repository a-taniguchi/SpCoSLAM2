## Online SpCoSLAM for SIGVerse

These codes are online learning of SpCoSLAM for SIGVerse simulator.

[Note] The gmapping (SLAM part) is separated. It needs to prepare the position data estimated by self-localization or SLAM in advance.


## Execution procedure

` python ./learnSpCoSLAM2.0_SIGVerse.py trialname datasetNUM `

`__init__.py`: It is necessary to set the file path of the relevant part in advance. 

Original dataset is here: https://github.com/a-taniguchi/SpCoDataset/tree/master/SIGVerseV3/similar  
[Note] Speech signal data is not included in the dataset.


