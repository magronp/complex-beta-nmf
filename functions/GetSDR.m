function [SDR,SIR,SAR,perm] = GetSDR(se,s)
% se are the estimated signals
% s are the reference signals


%%% Errors %%%
if nargin<2, error('Not enough input arguments.'); end

% Make sure we have column vectors
if(size(se,1)<size(se,2))
    se=se';
end
if(size(s,1)<size(s,2))
    s=s';
end

se = bsxfun(@minus,se,mean(se));
s  = bsxfun(@minus,s,mean(s));

[nsampl,nsrc]=size(se);
[nsampl2,nsrc2]=size(s);
if nsrc2~=nsrc, error('The number of estimated sources and reference sources must be equal.'); end
if nsampl2~=nsampl, error('The estimated sources and reference sources must have the same duration.'); end

%%% Performance criteria %%%
% Computation of the criteria for all possible pair matches
SDR=zeros(nsrc,nsrc);
SIR=zeros(nsrc,nsrc);
SAR=zeros(nsrc,nsrc);
for jest=1:nsrc,
    for jtrue=1:nsrc,
        [SDR(jest,jtrue),SIR(jest,jtrue),SAR(jest,jtrue)] = compute_measures(se(:,jest),s,jtrue);
    end
end
% Selection of the best ordering
perm=perms(1:nsrc);
nperm=size(perm,1);
meanSIR=zeros(nperm,1);
for p=1:nperm,
    meanSIR(p)=mean(SIR((0:nsrc-1)*nsrc+perm(p,:)));
end
[~,popt]=max(meanSIR);
perm=perm(popt,:).';
SDR=SDR((0:nsrc-1).'*nsrc+perm);
SIR=SIR((0:nsrc-1).'*nsrc+perm);
SAR=SAR((0:nsrc-1).'*nsrc+perm);
end

function [SDR,SIR,SAR] = compute_measures(se,s,j)

Rss = s' * s; % nsrc x nrsc
this_s = s(:,j);
% Get the scaling factor for the clean sources
a=this_s' * se / Rss(j,j);
e_true = a*this_s;
e_res = se-a*this_s;
Sss=sum((e_true).^2);
Snn=sum((e_res).^2);

SDR=10*log10(Sss/Snn);

% Get the SIR
Rsr = s' * e_res; % nsrc x 1
b = Rss \ Rsr;

e_interf = s * b;
e_artif = e_res - e_interf;

SIR=10*log10(Sss/sum((e_interf).^2));
SAR=10*log10(Sss/sum((e_artif).^2));

end


