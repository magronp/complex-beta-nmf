%  Multiplicative algorithms for the weighted NMF - with beta divergence
% 
% Inputs :
%     V : nonnegative matrix F*T
%     Wini : initial dictionary F*K
%     Hini : initial activation matrix K*T
%     Niter : number of iterations
%     beta : value of beta for the divergence
%     sigma_s : sparsity weight
%     p : sparsity parameter (p-norm)
%     mask : F*T weighting matrix
%     upW : update W (1) or not (0)
% 
% Outputs :
%     W : dictionary matrix
%     H : activation matrix
%     err : cost function
%     time : computation time

function [W,H,err,time] = NMF(V,Wini,Hini,Niter,beta,sigma_s,p,mask,upW)

verbose =0;

tinit=tic;
[F,T] = size(V); K=size(Wini,2);
V = V+eps;

if nargin<9
    upW=1;
end

if nargin<8
    mask = ones(F,T);
end

if nargin<7
    p=1;
end

if nargin<6
    sigma_s = (norm(V)^2)/(K^(1-p/2))*10^(-5);
end

if nargin<5
	beta = 1; % Defaut: KL divergence
end

% initialization
W=Wini; H=Hini;
WH=W*H+eps;

err = zeros(2,Niter+1);
time = zeros(1,Niter+1);
err(1,1) = beta_div(V,WH,beta,mask);
err(2,1) = sigma_s * sum(H(:).^p);

% algorithm
for k=1:Niter
    
    if verbose
        fprintf('iter %d ',k);
    end
    
    
    % Update H
    H = H.*   (W'*(V.*mask.*WH.^(beta-2)))  ./  (W'*(mask.*WH.^(beta-1)) + sigma_s*p*(H+eps).^(p-1) +eps);
    WH=W*H+eps;
    
    
    % Update W
    if upW
        W = W.*   ((V.*mask.*WH.^(beta-2))*H')  ./  ((mask.*WH.^(beta-1))*H'+eps);
        sumW = sqrt(sum(W).^2);
        W = W ./(repmat(sumW,[F 1])+eps) ;
        H = (repmat(sumW',[1 T])+eps) .* H;
        WH=W*H+eps;
    end
    

    % L2 normalization

    % error
    err(1,k+1) = beta_div(V,WH,beta,mask);
    err(2,k+1) = sigma_s * sum(H(:).^p);
    
    time(k+1) = toc(tinit);
end

end


% Beta-divergence
function D = beta_div(X,Y,beta,mask)

if nargin<4
    mask = ones(size(X));
end

switch beta
    case 0
        E = X./Y - log(X./Y)-1;
    case 1
        E = X .* log(X./Y)+Y-X;
    otherwise
        E =1/(beta*(beta-1)) * (X.^beta + (beta-1)*Y.^beta - beta*X.*Y.^(beta-1));
end

Emask = E .* mask;
D = sum(Emask(:));

end