clc; clear all; close all;
test_or_dev = 'Test';
set_settings_cbetanmf;

for ind=1:Nsongs

    % Source generation
    [sm,x,Sm,X,Vdico] = get_data_DSD_cbetanmf(dataset_path,test_or_dev,ind,Fs,Nfft,Nw,hop);
    [F,T] = size(X);
    V = abs(X).^2;
    
    % Dictionnary learning on isolated sources
    clc; fprintf('Song %d / %d \n',ind,Nsongs);
    fprintf('Dico learning \n');
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
    
    % Initial mu and same initial H for all algos
    muini = repmat(angle(X),[1 1 J]);
    Hini = cell(1,J); Hini_matrix = zeros(Ktot,T);
    for j=1:J
        Hini{j} = rand(K,T);
        Hini_matrix((j-1)*K+1:j*K,:) = Hini{j};
    end
    
    %%% Separation
    Se = zeros(F,T,J,Nalgo);
     
    % NMF
    fprintf('NMF \n');
    [~,H] = NMF(V,W_matrix,Hini_matrix,iter_nmf,beta_nmf,0,0,ones(F,T),0); aux =  wiener_nmf(X,W_matrix,H);
    Se(:,:,1,1) = sum(aux(:,:,1:K),3); Se(:,:,2,1) = sum(aux(:,:,K+1:end),3);
    
    % Complex EuNMF
    fprintf('Complex EuNMF \n');
    aux=cnmf_wisa(X,W,Hini,iter_nmf);  Se(:,:,:,2) = aux;
    
    % Complex ISNMF
    fprintf('Complex ISNMF \n');
    aux = complex_ISNMF_unconstrained(X,W,Hini,muini,kappa,iter_nmf,0); Se(:,:,:,3) = aux;
    
    % Complex beta-NMF
    fprintf('Complex betaNMF \n');
    aux = complex_betaNMF(X,W,Hini,muini,kappa,iter_nmf,beta_cnmf,0); Se(:,:,:,4) = aux;

    
   % Synthesis
    s_estim = zeros(J,length(sm),Nalgo);
    for al=1:Nalgo
        s_estim(:,:,al) = real(iSTFT(Se(:,:,:,al),Nfft,hop,Nw,'hann'));
    end
    
    % Record
    audiowrite(strcat(audio_path,int2str(ind),'_percu_orig.wav'),sm(1,:),Fs);
    audiowrite(strcat(audio_path,int2str(ind),'_harmo_orig.wav'),sm(2,:),Fs);
    for al = 1:Nalgo
        audiowrite(strcat(audio_path,int2str(ind),'_percu_estim_',algos{al},'.wav'),squeeze(s_estim(1,:,al)),Fs);
        audiowrite(strcat(audio_path,int2str(ind),'_harmo_estim_',algos{al},'.wav'),squeeze(s_estim(2,:,al)),Fs);
    end

end
