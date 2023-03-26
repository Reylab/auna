**AUNA** (**AU**dio and **N**erual **A**ctivity) is a graphic tool for marking audio events and analizing the neural data simultaneously recorded. It is currently an **alpha** version.

 
### Requirements  
 
- MATLAB >= 9.6 (R2019a)  
- Signal Processing Toolbox  
- Image Processing Toolbox  
- Statistics and Machine Learning Toolbox  
 
### Installation  
Add folder and subfolders to the Matlab path.  
_MAtNWB_ has to be installed to open NWB files. This can be done by running:  
 
```  
!git clone [https://github.com/NeurodataWithoutBorders/matnwb.git](https://github.com/NeurodataWithoutBorders/matnwb.git)  
cd matnwb  
addpath(genpath(pwd));  
generateCore();  
```

 
### Instructions  
Configure the parameters in the par_auna.m file, go to the folder with the data and then run `auna` in the command window.
The events and zones are saved insiede the file  **auna_data.mat**.

### Files
#### Audio file
AUNA loads the audio from a binary file **NSX[channel_number].NC5** with values in int16. The metadata will be loaded from from a file named **NSX_TimeStamps.mat** with the variables _sr_ (sampling rate in Hz) and _lts_ (number of samples) .
If a file NSx.mat is found in the path, AUNA will try to load the channel using the metadata table in that file. 

#### Units files
If par.neuroformat = **WC**: The units activity will be loaded from the files: **times_[channel_number].mat**
If par.neuroformat = **NWB**: The units activity will be loaded from the files: **'chan_[channel_number].nwb**

#### Session configuration and events
By default the annotations are saved inside the file **auna_data.mat**. The zones structure has a field for each zone with its limits and pre/post times to create rasters.
The concepts struct has a field for each concept with times in ms for each annotation.

### Shortcuts  
 
- numerical minus('-'): maximum timescale (zoom out icon)  
- space bar: pause/resume (||)  
- right_arrow: next data block (>>)  
- left_arrow: prev data block (<<)  
- return: Add mark  
- 1-9: choose the event number
  
