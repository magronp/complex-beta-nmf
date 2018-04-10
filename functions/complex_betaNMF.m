function [m_post,W,H,mu,cost] = complex_betaNMF(X,Wini,Hini,muini,kappa,Niter,beta,upW)

if nargin<8
    upW=1;
end

[F,T,J] = size(muini);

% Anisotropy parameters
lambda = besseli(1,kappa) ./ besseli(0,kappa) *sqrt(pi)/2;
rho=besseli(2,kappa)./besseli(0,kappa) - lambda.^2;

% Initial parameters
W = Wini; H = Hini;
v = zeros(size(muini));
for j=1:J
   v(:,:,j) = W{j}*H{j}; 
end
mu = muini;

% Cost
cost = zeros(1,Niter+1);

% EM algo
for iter=1:Niter
    
    % --------------------------------- E step

    % ----- sources moments
    gamma = (1-lambda.^2).*v;
    c = rho.*v .* exp(2i*mu);

    % ----- mixture moments
    gamma_X = sum(gamma,3);      
    c_X = sum(c,3);
    detGamma = gamma_X.^2 - abs(c_X).^2;

    % ----- posterior moments
    m_post = ( (gamma.*gamma_X - c.*conj(c_X)).*(X) + (c.*gamma_X - gamma.*c_X).* conj(X)  ) ./ (detGamma+eps) ;
    gamma_post = gamma - ( gamma_X .* (gamma.^2+abs(c).^2) - 2 * gamma .* real(c.*conj(c_X))  )  ./ (detGamma+eps);
    c_post = c - (2*gamma.*gamma_X.*c - gamma.^2 .* c_X - c.^2 .* conj(c_X)  )  ./ (detGamma+eps);

    % ----- phase-corrected posterior power
    p = ((1-lambda.^2)*(gamma_post+abs(m_post).^2)-rho.*real(exp(-2*1i*mu).*(c_post+m_post.^2))  )/((1-lambda.^2)^2-rho.^2)  ;
    
    
    % Cost
    cost(iter) = bdiv(p,v,beta) ;
    
    
    % --------------------------------- M step

    % ----- NMF parameters
    for j=1:J
        
        % Get current source parameters
        vj = v(:,:,j)+eps;
        pj = p(:,:,j);

        % Update H
        H{j} =  H{j} .* ( W{j}' *(pj.*vj.^(beta-2)))  ./ (  W{j}' * ( vj.^(beta-1)  ) ) ;
        vj= W{j}*H{j}+eps;
        
        % Update W
        if upW
            W{j} =  W{j} .*  ( (pj.*vj.^(beta-2))* H{j}' )  ./ (  (vj.^(beta-1))* H{j}'  )  ;
            sumWj = sqrt(sum(W{j}).^2);
            W{j} = W{j} ./repmat(sumWj,[F 1]) ;
            H{j} = repmat(sumWj',[1 T]) .* H{j};           
        end
        
        v(:,:,j)= W{j}*H{j};
    end
    
        
    % ----- phase parameters
    %mu = angle(c_post+m_post.^2)/2;
    mu = (angle(c_post+m_post.^2)-pi)/2;
    
end

% Finally, one E step for having the most up to date posterior mean

% ----- update moments
gamma = (1-lambda.^2).*v;
c = rho.*v .* exp(2i*mu);
gamma_X = sum(gamma,3);      
c_X = sum(c,3);
detGamma = gamma_X.^2 - abs(c_X).^2;
m_post = ( (gamma.*gamma_X - c.*conj(c_X)).*(X) + (c.*gamma_X - gamma.*c_X).* conj(X)  ) ./ (detGamma+eps);
    
% Cost
p = ((1-lambda.^2)*(gamma_post+abs(m_post).^2)-rho.*real(exp(-2*1i*mu).*(c_post+m_post.^2))  )/((1-lambda.^2)^2-rho.^2)  ;
cost(iter+1) = bdiv(p,v,beta) ;
    
end

% Beta-divergence
function [co] = bdiv(p,v,beta) 

J = size(p,3);
co = 0;
for j=1:J
    co = co + beta_div(p,v,beta);
end
   
end

function D = beta_div(X,Y,beta,mask)

X=X+eps; Y=Y+eps;

if nargin<4
    mask = ones(size(X));
end

switch beta
    case 0
        E = X./Y - log(X./Y) -1;
    case 1
        E = X .* log(X./Y)+Y-X;
    otherwise
        E =1/(beta*(beta-1)) * (X.^beta + (beta-1)*Y.^beta - beta*X.*Y.^(beta-1));
end

Emask = E .* mask;
D = sum(Emask(:));

end