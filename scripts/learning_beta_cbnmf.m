clc; clear all; close all;
test_or_dev = 'Dev';
set_settings_cbetanmf;

% Beta-divergence range
B = 0:0.1:2;
Nb = length(B);
score =zeros(J,3,Nb,Nsongs,2);

for ind=1:Nsongs

    % Source generation
    [sm,x,Sm,X,Vdico] = get_data_DSD_cbetanmf(dataset_path,test_or_dev,ind,Fs,Nfft,Nw,hop);
    [F,T] = size(X);
    V =abs(X).^2;
    
    % Dictionnary learning on isolated sources
    clc; fprintf('Song %d / %d \n Dico learning \n',ind,Nsongs);
    W = cell(1,J); W_matrix = zeros(F,Ktot);
    for j=1:J
        Vj = Vdico(:,:,j);
        
        % normalize (l2-norm) the spectrogram in each time frame
        Vjnorm = Vj ./ repmat(sqrt(sum(Vj.^2))+eps,[F 1]);
        
        % k-means clustering
        [~,waux] = kmeans(Vjnorm',K);
        
        % Normalize the learned dictionary
        waux = waux' ./ repmat(sqrt(sum(waux'.^2))+eps,[F 1]);
        
        W{j}=waux+eps;
        W_matrix(:,(j-1)*K+1:j*K) = W{j};
    end
    
    % Initial mu and same H for all algos / beta
    muini = repmat(angle(X),[1 1 J]);
    Hini = cell(1,J); Hini_matrix = zeros(Ktot,T);
    for j=1:J
        Hini{j} = rand(K,T); Hini_matrix((j-1)*K+1:j*K,:) = Hini{j};
    end
    
    % Loop over beta
    for b=1:Nb
        beta = B(b);
        clc; fprintf('Song %d / %d \n beta %d / %d  \n',ind,Nsongs,b,Nb);
        
        % NMF
        Se = zeros(F,T,J);
        [~,H] = NMF(abs(V),W_matrix,Hini_matrix,iter_nmf,beta,0,0,ones(F,T),0);
        aux =  wiener(X,W_matrix,H);
        Se(:,:,1) = sum(aux(:,:,1:K),3); Se(:,:,2) = sum(aux(:,:,K+1:end),3);
        se = real(iSTFT(Se,Nfft,hop,Nw,wtype));
        [sdr,sir,sar] = GetSDR(se,sm);
        score(:,:,b,ind,1) = [sdr sir sar];
        
        % Complex beta-NMF
        mpost = complex_betaNMF(X,W,Hini,muini,kappa,iter_nmf,beta,0);
        se = real(iSTFT(mpost,Nfft,hop,Nw,wtype));
        [sdr,sir,sar] = GetSDR(se,sm);
        score(:,:,b,ind,2) = [sdr sir sar];
    end

end

% Record scores
save(strcat(metrics_path,'learning_beta_cbnmf.mat'),'score','B');

% Plot results
li = isnan(score(1,1,1,:,1));
score(:,:,:,li,:) = [];
SDR = mean(score(:,1,:,:,:),4); SDR = squeeze(mean(SDR,1));
SIR = mean(score(:,2,:,:,:),4); SIR = squeeze(mean(SIR,1));
SAR = mean(score(:,3,:,:,:),4); SAR = squeeze(mean(SAR,1));

% NMF
figure;
subplot(1,3,1); plot(B,SDR(:,1)); xlabel('\beta','fontsize',16); title('SDR (dB)');
subplot(1,3,2); plot(B,SIR(:,1)); xlabel('\beta','fontsize',16); title('SIR (dB)');
subplot(1,3,3); plot(B,SAR(:,1)); xlabel('\beta','fontsize',16); title('SAR (dB)');

% Complex beta-NMF
figure;
subplot(1,3,1); plot(B,SDR(:,2)); xlabel('\beta','fontsize',16); title('SDR (dB)');
subplot(1,3,2); plot(B,SIR(:,2)); xlabel('\beta','fontsize',16); title('SIR (dB)');
subplot(1,3,3); plot(B,SAR(:,2)); xlabel('\beta','fontsize',16); title('SAR (dB)');
