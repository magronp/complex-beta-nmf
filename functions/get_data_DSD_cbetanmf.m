function [sm,x,Sm,X,Vdico] = get_data_DSD_cbetanmf(dataset_path,testordev,num_piece,Fs,Nfft,Nw,hop,wtype)

if nargin<8
    wtype = 'hann';
end

% Time limits for the dictionary learning and separation segments
time1 = 70; time2=80; time3=90;

%%% Get files path
dir_path = strcat(dataset_path,'Sources/',testordev);
aux = dir(dir_path);
song_name = aux(num_piece+2).name;
list_instr = {'bass','drums','other','vocals'};
L = strcat(dir_path,'/',song_name,'/',list_instr,'.wav');


%%% Get dicos signals STFT

dico_percu = audioread(L{2},[Fs*time1 Fs*time2-1]);
dico_percu=mean(dico_percu,2)';

dico_bass = audioread(L{1},[Fs*time1 Fs*time2-1]);
dico_other = audioread(L{3},[Fs*time1 Fs*time2-1]);
dico_vocals = audioread(L{4},[Fs*time1 Fs*time2-1]);
dico_harmo = mean(dico_bass + dico_other + dico_vocals,2)';

s_dico = [dico_percu;dico_harmo];
Vdico = abs(STFT(s_dico,Nfft,hop,Nw,wtype)).^2;


%%% Get the test signals
test_percu = audioread(L{2},[Fs*time2 Fs*time3-1]);
test_percu=mean(test_percu,2)';

test_bass = audioread(L{1},[Fs*time2 Fs*time3-1]);
test_other = audioread(L{3},[Fs*time2 Fs*time3-1]);
test_vocals = audioread(L{4},[Fs*time2 Fs*time3-1]);
test_harmo = mean(test_bass + test_other + test_vocals,2)';

sm = [test_percu;test_harmo];

% STFT
Sm = STFT(sm,Nfft,hop,Nw,wtype);
sm = iSTFT(Sm,Nfft,hop,Nw,wtype);

% Mixture
x = sum(sm,1); X = sum(Sm,3);

end