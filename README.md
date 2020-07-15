# Online Robust PCA
An implementation of online robust PCA (ORPCA) described in the paper: 
* [Online Robust PCA via Stochastic Optimization](https://papers.nips.cc/paper/5131-online-robust-pca-via-stochastic-optimization) by Feng et al. [NIPS 2013]     

However, it is modified slightly since missing data points cannot be reconstructed in this paper. (The reconstruction of missing data were not performed here, but it can be performed by chaning the missing percentage parameter 'perc_miss' in the m-file 'main_orpca.m'.)

## Dataset
[VIRAT Video Dataset](https://viratdata.org) was used here (not in the paper), which is designed to be realistic, natural and challenging for video surveillance domains in terms of its resolution, background clutter, diversity in scenes, and human activity/event categories than existing action recognition datasets. It has become a benchmark dataset for the computer vision community.

## Results
### Foreground and Background
