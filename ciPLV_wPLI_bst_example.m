clear;clc
load('test.mat');
%HA=HA(1,:);
[nA,ntime]=size(HA);
[nB,ntime]=size(HB);
% commented below, since we found out that repmat-ing HA directly is faster (see
% below)
%iA = repmat(1:nA, 1, nB)';
%iB = reshape(repmat(1:nB, nA, 1), [], 1);
R_wPLI=zeros(nA,nB);


% PLV
tic
phaseA = HA ./ abs(HA);
phaseB = HB ./ abs(HB);
csd=phaseA*phaseB';
R_PLV=abs(csd/ntime);
t=toc;
disp(['PLV, ' num2str(t) ' seconds']);

% ciPLV
tic
phaseA = HA ./ abs(HA);
phaseB = HB ./ abs(HB);
csd=phaseA*phaseB';
R_ciPLV=abs((imag((csd))/ntime)./sqrt(1-(real((csd))/ntime).^2));
t=toc;
disp(['ciPLV, ' num2str(t) ' seconds']);

% wPLI
tic
R_wPLI=zeros(nA,nB);
for itime=1:ntime
    phaseA_t = HA(:,itime) ./ abs(HA(:,itime));
    phaseB_t = HB(:,itime) ./ abs(HB(:,itime));
    csd=phaseA_t*phaseB_t';
    cdi = imag(csd);
    R_wPLI=R_wPLI+(abs(cdi).*sign(cdi))'./abs(cdi)';
end
wPLI=abs(R_wPLI/ntime);
t=toc;
disp(['wPLI, ' num2str(t) ' seconds']);

% wPLI debiased
tic
HA=HA';HB=HB';
sin_pd=sin(angle(repmat(HA,[1 nA])./repelem(HB,1, nA)));
a=abs(mean(sin_pd)./mean(abs(sin_pd)));
R_wPLI(1:length(a))=a;
t=toc;
disp(['wPLI phase difference, ' num2str(t) ' seconds']);

% faster implementation of wPLI (debiased), identical to the previous one
tic
num = abs(imag(phaseA*phaseB'));
den = zeros(nA,nB);
for t = 1:ntime
    den = den + abs(imag(phaseA(:,t) * phaseB(:,t)'));
end
wPLI_db = num./den;
t=toc;
disp(['wPLI db summing terms, ' num2str(t) ' seconds']);

if(any([nA,nB]==1))
    c=corr(wPLI_db(:),R_wPLI(:), 'rows','complete');
else
    c=compareconn(wPLI_db,R_wPLI);
end
disp(['comparison between the two wPLIs = ' num2str(c)])

if(any([nA,nB]==1))
    c=corr(wPLI_db(:),wPLI(:), 'rows','complete');
else
    c=compareconn(wPLI_db,wPLI);
end
disp(['comparison between wPLI db and wPLI = ' num2str(c)])

if(any([nA,nB]==1))
    c=corr(wPLI_db(:),R_ciPLV(:), 'rows','complete');
else
    c=compareconn(wPLI_db,R_ciPLV);
end
disp(['comparison between wPLI db and ciPLV = ' num2str(c)])
