% Set the settings used in the experiments

Fs = 44100;

% STFT parameters
Nfft = 4096;
Nw = 4096;
hop = Nw/4;
wtype = 'hann';

% Paths
dataset_path = 'datasets/DSD100/';
audio_path = 'complex-beta-nmf/audio_files/';
metrics_path = 'complex-beta-nmf/metrics/';

% Data
Nsongs = 50;
J = 2;          %number of sources

% Algorithms
algos ={'NMF';'cEuNMF';'cISNMF';'cbetaNMF'};
Nalgo = length(algos);

% NMF
K = 50;
Ktot = K*J;
iter_nmf = 100;

% Parameters
kappa = 1;
beta_cnmf = 0.5;
beta_nmf = 0.75;
