clear all; close all; clc;
test_or_dev = 'Test';
set_settings_cbetanmf;

score = zeros(J,3,Nalgo,Nsongs);
for ind=1:Nsongs
    
    clc; fprintf('BSS - Song %d / %d \n',ind,Nsongs);
     
    % Load original sources
    s1 = audioread(strcat(audio_path,int2str(ind),'_percu_orig.wav'));
    s2 = audioread(strcat(audio_path,int2str(ind),'_harmo_orig.wav'));
    sm = [s1 s2]';
    
    % Loop over the algorithms
   for al=1:Nalgo
       
       % Load estimated sources
       s_estim1 =  audioread(strcat(audio_path,int2str(ind),'_percu_estim_',algos{al},'.wav'));
       s_estim2 =  audioread(strcat(audio_path,int2str(ind),'_harmo_estim_',algos{al},'.wav'));
       se = [s_estim1 s_estim2]';
       
       % BSS eval
       [sdr,sir,sar] = GetSDR(se,sm);
       score(:,:,al,ind) = [sdr sir sar];
   end
   
end

% Record score
li = squeeze(isnan(score(1,1,1,:)));
score(:,:,:,li) = [];
save(strcat(metrics_path,'bss_complex-beta-nmf.mat'),'score');

wo = mean(score,4);
squeeze(wo(1,:,:))'
