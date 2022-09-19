clear;clc
% uncomment/comment accordingly below.
% test_channels is 274 variables, 361 time points
% test_vertices is 4002 variables, 52 time points
% load('test_channels.mat');
load('test_vertices.mat')
HA=HA(1,:); % this is for 1xN scenario, comment for NxN
[nA,~]=size(HA);
[nB,nt]=size(HB);


% PLV
tic
phaseA = HA ./ abs(HA);
phaseB = HB ./ abs(HB);
csd=phaseA*phaseB';
R_PLV=abs(csd/nt);
t=toc;
disp(['PLV, ' num2str(t) ' seconds']);

% ciPLV
tic
phaseA = HA ./ abs(HA);
phaseB = HB ./ abs(HB);
csd=phaseA*phaseB';
R_ciPLV=abs((imag((csd))/nt)./sqrt(1-(real((csd))/nt).^2));
t=toc;
disp(['ciPLV, ' num2str(t) ' seconds']);

% wPLI
tic
R_wPLI=zeros(nB,nA);
for itime=1:nt
    phaseA_t = HA(:,itime) ./ abs(HA(:,itime));
    phaseB_t = HB(:,itime) ./ abs(HB(:,itime));
    csd=phaseA_t*phaseB_t';
    cdi = imag(csd);
    R_wPLI=R_wPLI+(abs(cdi).*sign(cdi))'./abs(cdi)';
end
wPLI=abs(R_wPLI/nt);
t=toc;
disp(['wPLI, ' num2str(t) ' seconds']);

% wPLI debiased based on sine differences
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
for t = 1:nt
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
disp(['comparison between the two debiased wPLIs = ' num2str(c)])

if(any([nA,nB]==1))
    c=corr(wPLI_db(:),wPLI(:), 'rows','complete');
else
    c=compareconn(wPLI_db,wPLI);
end
disp(['comparison between wPLI debiased and wPLI = ' num2str(c)])

if(any([nA,nB]==1))
    c=corr(wPLI_db(:),R_ciPLV(:), 'rows','complete');
else
    c=compareconn(wPLI_db,R_ciPLV);
end
disp(['comparison between wPLI debiased and ciPLV = ' num2str(c)])

if(any([nA,nB]==1))
    c=corr(R_ciPLV(:),wPLI(:), 'rows','complete');
else
    c=compareconn(R_ciPLV,wPLI);
end
disp(['comparison between ciPLV and wPLI = ' num2str(c)])

[nr,nc]=size(R_ciPLV);
if nr==nc
    N=max(nr,nc);
    Isubdiag = find(tril(ones(N),-1));
    R_ciPLV=R_ciPLV(Isubdiag);
    wPLI_db=wPLI_db(Isubdiag);
end
scatter(R_ciPLV,wPLI_db);xlim([-.05 1.05]);ylim([-.05 1.05])
