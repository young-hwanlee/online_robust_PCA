# Online Robust PCA
An implementation of online robust PCA (ORPCA) described in the paper: 
* [Online Robust PCA via Stochastic Optimization](https://papers.nips.cc/paper/5131-online-robust-pca-via-stochastic-optimization) by Feng et al. [NIPS 2013]     

However, it is modified slightly since missing data points cannot be reconstructed in this paper. (The reconstruction of missing data were not performed here, but it can be performed by chaning the missing percentage parameter 'perc_miss' in the m-file 'main_orpca.m'.)

## Dataset
[VIRAT Video Dataset](https://viratdata.org) was used here (not in the paper), which is designed to be realistic, natural and challenging for video surveillance domains in terms of its resolution, background clutter, diversity in scenes, and human activity/event categories than existing action recognition datasets. It has become a benchmark dataset for the computer vision community.

## Results
### Original frame size
Since the frame size of the original video files is large, they were cropped when running the simulation. The original and gray frame are shown below:  
<table align='center'>
<tr align='center'>
<td> Original frame </td>
<td> Gray frame </td>
</tr>
<td><img width="576" alt="original_figure" src="https://user-images.githubusercontent.com/67979833/87625647-a6bfd400-c6f8-11ea-9b7b-7c90083cdf33.png">
<td><img width="576" alt="gray_figure" src="https://user-images.githubusercontent.com/67979833/87625645-a6273d80-c6f8-11ea-98c6-66b3906742ce.png">
</tr>
</table>

### Foreground and Background
After running the algorithm, the results (original, background, and foreground frames) are displayed below 
where
* Foreground: moving object is captured by the sparse matrix E
* Background: non-moving object is captred by the low rank matrix X     
(If fine-tuning the parameters, running more frames (or increasing the epoch), the sparse matrix E may capture only the moving object (not the background.))

<table align='center'>
<tr align='center'>
<td> Original 1 </td>
<td> Original 2 </td>
<td> Original 3 </td>
<td> Original 4 </td>
</tr>
<tr>
<td><img width="576" alt="data_Z_1" src="https://user-images.githubusercontent.com/67979833/87624863-c229df80-c6f6-11ea-85c3-b68e448fc8ca.png">
<td><img width="576" alt="data_Z_2" src="https://user-images.githubusercontent.com/67979833/87624864-c229df80-c6f6-11ea-9f7a-53fab3f3a1d5.png">
<td><img width="576" alt="data_Z_3" src="https://user-images.githubusercontent.com/67979833/87624867-c229df80-c6f6-11ea-83ec-ab7e5610b6c7.png">
<td><img width="576" alt="data_Z_4" src="https://user-images.githubusercontent.com/67979833/87624868-c2c27600-c6f6-11ea-95d1-1f870a879d1f.png">
</tr>
<tr align='center'>
<td> Background 1 </td>
<td> Background 2 </td>
<td> Background 3 </td>
<td> Background 4 </td>
</tr>
<tr>
<td><img width="576" alt="low_rank_X_1" src="https://user-images.githubusercontent.com/67979833/87624870-c2c27600-c6f6-11ea-84d9-51189cadcb69.png">
<td><img width="576" alt="low_rank_X_2" src="https://user-images.githubusercontent.com/67979833/87624871-c2c27600-c6f6-11ea-9c8c-220d9683ff48.png">
<td><img width="576" alt="low_rank_X_3" src="https://user-images.githubusercontent.com/67979833/87624872-c2c27600-c6f6-11ea-8760-030bcfa8a8f2.png">
<td><img width="576" alt="low_rank_X_4" src="https://user-images.githubusercontent.com/67979833/87624873-c2c27600-c6f6-11ea-91c7-dc86f7074961.png">
</tr>
<tr align='center'>
<td> Foreground 1 </td>
<td> Foreground 2 </td>
<td> Foreground 3 </td>
<td> Foreground 4 </td>
</tr>
<tr>
<td><img width="576" alt="sparse_E_1" src="https://user-images.githubusercontent.com/67979833/87624874-c35b0c80-c6f6-11ea-9083-52d6e0e709b3.png">
<td><img width="576" alt="sparse_E_2" src="https://user-images.githubusercontent.com/67979833/87624875-c35b0c80-c6f6-11ea-9774-6b8da772387f.png">
<td><img width="576" alt="sparse_E_3" src="https://user-images.githubusercontent.com/67979833/87624876-c35b0c80-c6f6-11ea-8578-1884429359d4.png">
<td><img width="576" alt="sparse_E_4" src="https://user-images.githubusercontent.com/67979833/87624877-c35b0c80-c6f6-11ea-8b68-dc6ad9bb09f0.png">
</tr>
</table>
