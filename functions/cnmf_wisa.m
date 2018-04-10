% Complex NMF with phase unwrapping constraints and intra-source additivity
%
% Please refer to: Complex NMF under phase constraints based on signal modeling: application to audio source separation
% P. Magron, R. Badeau and B. David, in Proc. of IEEE ICASSP 2016
%
% Inputs:
%     X : F*T mixture's STFT
%     Wini : a cell with J initial dictionaries (F*K)
%     Hini : a cell with J initial activations matrix (K*T)
%     Niter : number of iterations
%     upW : =1 to update the dictionary of 0 if not
%
% Outputs:
%     Xhat : F*T*J estimated sources
%     Xmix : F*T estimated mixture
%     W,H : dictionnary and activations matrices
%     Phi : exp(i*phi) phases matrices (F*T*J)
%     err : cost function over iterations

function [Xhat,Xmix,W,H,Phi,err]=cnmf_wisa(X,Wini,Hini,Niter,upW)

[F,T] = size(X);
J = length(Hini);

if nargin<5
    upW = 0;  %don't update the dictionary by default
end

% Initialization
W = Wini;
H = Hini;
Phi = repmat(X./(abs(X)+eps),[1,1,J]);

% init auxiliary variables
[Xhat,Xmix,B,Vhat,G] = update_aux(X,W,H,Phi);

% error
err = zeros(2,Niter+1);
err(:,1) = norm(X-Xmix).^2;

for it=1:Niter

    % Update Phi
    Phi = exp(1i*angle(B));
        
    % up aux param
    beta = abs(B);
    
    for j=1:J
	% Update H
        H{j} = H{j} .* (  W{j}' * (beta(:,:,j)./(G(:,:,j)+eps)) ) ./ ( W{j}'*( (Vhat(:,:,j))./(G(:,:,j)+eps) )  + eps);
        H{j} = max(H{j},eps);

        % Update W
        if upW
		W{j} = W{j} .* ( (beta(:,:,j)./(G(:,:,j)+eps))*H{j}' ) ./ ( ( (Vhat(:,:,j))./(G(:,:,j)+eps)*H{j}' )  + eps);
	   	W{j} = max(W{j},eps);

            	sumWj = sqrt(sum(W{j}).^2);
            	W{j} = W{j} ./repmat(sumWj,[F 1]) ;
            	H{j} = repmat(sumWj',[1 T]) .* H{j};           
        end


    end
    
    % update auxiliary variables
    [Xhat,Xmix,B,Vhat,G] = update_aux(X,W,H,Phi);
    
    %error
    err(:,it+1) = norm(X-Xmix).^2;
end


end


% Update auxiliary variables
function [Xhat,Xmix,B,Vhat,G] = update_aux(X,W,H,Phi)
    J = size(Phi,3);
    Vhat = computeVk(W,H);
    Xhat = Vhat .* Phi;
    Xmix = sum(Xhat,3);
    G = Vhat ./ repmat(sum(Vhat,3)+eps,[1 1 J]);
    B = repmat(X-Xmix,[1,1,J]).*G + Xhat;
end

% Compute Vk=WkHk (tensor of size F*T*J)
function V = computeVk(W,H)

    J = length(W); F = size(W{1},1); T = size(H{1},2);
    V = zeros(F,T,J);
    for j=1:J
        V(:,:,j) = W{j}*H{j};
    end

end