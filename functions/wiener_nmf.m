% Wiener filtering from NMF decomposition
%
% Inputs :
%     X : complex mixture's STFT
%     W, H : spectral templates and temporal activation
% 
% Outputs :
%     C : complex rank-1 separated components

function [C] = wiener_nmf(X,W,H)

F = size(X,1);
[K,T] = size(H);

WkHk = zeros(F,T,K);
for k=1:K
    WkHk(:,:,k) = W(:,k)*H(k,:);
end

C = WkHk .* X ./ (W*H+eps);

end
