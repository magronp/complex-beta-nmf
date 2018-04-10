# Complex NMF with the beta-divergence GitHub repository

Code for complex NMF with the beta-divergence

Here, you will find the code related to complex NMF with the beta-divergence. It is applied to the task of harmonic/percussive source separation.

If you use any of the things existing in this repository, please cite the [corresponding paper](https://hal.archives-ouvertes.fr/hal-01741278). 


### Dataset set-up

To reproduce the experiments conducted in our paper, you will need to download the [Dexmixing Secret Database (DSD100)](http://www.sisec17.audiolabs-erlangen.de) and to place its content in the `dataset`.

If you use this dataset, you will end up with the proper directory structure and file names, as used in the `functions/get_data_DSD_cbetanmf` function.

If you want to use a different dataset, then you have two options: 
- either you format your file names and directory structure to match the one from DSD100;
- or you modify the file reading function `get_data_DSD_cbetanmf` to suit your needs.


### Experiments

The two main scripts to reproduce the experiments are in the `scripts` folder:

`learning_beta_cbnmf` will perform a grid search for beta and plot the separation quality (SDR, SIR and SAR) as a function of this parameter. This apply for Complex NMF and NMF.

`separation_cbetanmf` will perform the benchmark conducted in the paper, that is, comparing NMF and complex NMF with Euclidean, IS and beta divergences. The audio files will be record in the `audio_files` folder and the script `comopute_score_cbetanmf` computes the SDR, SIR and SAR for this benchmark.
