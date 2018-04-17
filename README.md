# Features-extraction-based-on-Beta-Elliptic-Model-and-Fuzzy-Elementary-Perceptual-Codes
The presented code consists of the preprocessing and the segmentation of online handwriting into a sequence of Beta strokes. After that, from each Beta stroke, we extract a set of static and dynamic features using four features extraction techniques based on the Beta-Elliptic model and the Fuzzy Elementary Perceptual Codes. Then, all the segments which are composed of N consecutive Beta strokes are categorized into groups and subgroups according to their position and their geometric characteristics.      
- The first features extraction technique using  Beta-Elliptic model called the Advanced Overlapped Beta Strategy (AOBS).  
- The second features extraction technique using  Beta-Elliptic model called the Simplified Beta Strategy (SBS).  
- The third features extraction technique using  combination between Advanced Overlapped Beta Strategy and Fuzzy Elementary Perceptual Codes (AOBSFEPC).  
- The fourth features extraction technique using  combination between Simplified Beta Strategy and Fuzzy Elementary Perceptual Codes (SBSFEPC).    

All details found in our paper https://arxiv.org/abs/1804.05661

Getting started

Run the script with Matlab: Run_BEMFEPC.m

What should the person wishing to run this locally do?  

The person wishing to run this locally should do these steps:
1) Create one folder for each writer in the folder data (for example the folder "1" containing the files .inkml of the writer 1, the folder "2" containing the files .inkml of the writer 2, ...). You can use the folder data for the first Run.
2) Modify line 10 : path_folder_Samples = ['../data/',t,'/'];
Line 10 must contain the path of the folder named "data". 
3) Modify line 11 : path_folder_Results = ['../results/'];
The results will be saved also in the created folder named "results"

Dataset 

We're using public IBM_UB_1 dataset 
You can download the dataset by using this link: https://cubs.buffalo.edu/research/hwdata
