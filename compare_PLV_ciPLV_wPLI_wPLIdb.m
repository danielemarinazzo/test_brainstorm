clear;clc
% uncomment/comment accordingly below.
% test_channels is 274 variables, 361 time points
% test_vertices is 4002 variables, 52 time points
load('test_channels.mat');
%load('test_vertices.mat')
%HA=HA(1,:); % this is for 1xN scenario, comment for NxN
[nA,~]=size(HA);
[nB,nt]=size(HB);


% PLV
tic
phaseA = HA ./ abs(HA);
phaseB = HB ./ abs(HB);
csd=phaseA*phaseB';
PLV=abs(csd/nt);
t=toc;
disp(['PLV, ' num2str(t) ' seconds']);

% ciPLV
tic
phaseA = HA ./ abs(HA);
phaseB = HB ./ abs(HB);
csd=phaseA*phaseB';
ciPLV=abs((imag((csd))/nt)./sqrt(1-(real((csd))/nt).^2));
t=toc;
disp(['ciPLV, ' num2str(t) ' seconds']);

% wPLI sign cdi
tic
num=zeros(nB,nA);
den=zeros(nB,nA);
phaseA = HA ./ abs(HA);
phaseB = HB ./ abs(HB);
for t=1:nt
    cdi=imag(phaseA(:,t) * phaseB(:,t)');
    num=num+(abs(cdi).*sign(cdi))';
    den=den+abs(cdi)';
end
wPLI_sc=abs(num/nt)./(den/nt);
t=toc;
disp(['wPLI sign cdi, ' num2str(t) ' seconds']);

% wPLI based on sin differences
tic
data=real(HB);
[nc,ns]=size(HB);
phs=angle(HB);
tpli_num = complex(zeros(nc));
tpli_den = complex(zeros(nc));
for s=1: ns
    dphs=bsxfun(@minus, phs(:, s), phs(:, s)');
    tpli_den=tpli_den+abs(sin(dphs));
    tpli_num=tpli_num+sin(dphs);
end
wPLI_sindiff=abs(tpli_num/ns)./(abs(tpli_den)/ ns);
t=toc;
disp(['wPLI sin differences, ' num2str(t) ' seconds']);

% wPLI ratio imag csd
tic
phaseA = HA ./ abs(HA);
phaseB = HB ./ abs(HB);
num = imag(phaseA*phaseB');
den = zeros(nA,nB);
for t = 1:nt
    den = den + abs(imag(phaseA(:,t) * phaseB(:,t)'));
end
wPLI_csdrat = abs(num./den);
t=toc;
disp(['wPLI ratio imag csd, ' num2str(t) ' seconds']);

% wPLI based on sine of a ratio without loop
tic
HA=HA';HB=HB';
sin_pd=sin(angle(repmat(HA,[1 nA])./repelem(HB,1, nA)));
a=abs(mean(sin_pd)./mean(abs(sin_pd)));
N = numel(a) ;
wPLI_sinrat = NaN(ceil(sqrt(N))) ;
wPLI_sinrat(1:N) = a ;
t=toc;
disp(['wPLI sine ratio, ' num2str(t) ' seconds']);


% wPLI debiased ratio imag csd
tic
num = imag(phaseA*phaseB');
den = zeros(nA,nB);sqd = zeros(nA,nB);
for t = 1:nt
    den = den + abs(imag(phaseA(:,t) * phaseB(:,t)'));
    sqd = sqd + imag(phaseA(:,t)*phaseB(:,t)').^2;
end
wPLI_db_csdrat = (num.^2-sqd)./(den.^2-sqd);
t=toc;
disp(['wPLI debiased ratio imag csd, ' num2str(t) ' seconds']);

% fieldtrip wPLI
tic
num = zeros(nA,nB);
den = zeros(nA,nB);
for t = 1:nt
    num = num + imag(phaseA(:,t)*phaseB(:,t)');
    den = den + abs(imag(phaseA(:,t) * phaseB(:,t)'));
end
wPLI_ft = abs(num./den);
t=toc;
disp(['wPLI fieldtrip, ' num2str(t) ' seconds']);

% fieldtrip db wPLI
tic
num = zeros(nA,nB);
den = zeros(nA,nB);
sqd = zeros(nA,nB);
for t = 1:nt
    num = num + imag(phaseA(:,t)*phaseB(:,t)');
    den = den + abs(imag(phaseA(:,t) * phaseB(:,t)'));
    sqd = sqd + imag(phaseA(:,t)*phaseB(:,t)').^2;
end
wPLI_db_ft = (num.^2-sqd)./(den.^2-sqd);
t=toc;
disp(['wPLI db fieldtrip, ' num2str(t) ' seconds']);


measures={'PLV','ciPLV','wPLI_sc','wPLI_sindiff','wPLI_sinrat',...
    'wPLI_csdrat','wPLI_ft','wPLI_db_csdrat','wPLI_db_ft'};
comp_meas=ones(length(measures));
for imeas=1:length(measures)
    for jmeas=1:imeas
        if imeas~=jmeas
            eval(['A=',measures{imeas},';']);
            eval(['B=',measures{jmeas},';']);
            [nr,nc]=size(A);
            if(any([nr,nc]==1))
                comp_meas(imeas,jmeas)=corr(A(:),B(:), 'rows','complete');
                comp_meas(jmeas,imeas)=comp_meas(imeas,jmeas);
            else
                comp_meas(imeas,jmeas)=compareconn(A,B);
                comp_meas(jmeas,imeas)=comp_meas(imeas,jmeas);
            end
        end
    end
end
figure
set(groot,'defaultAxesTickLabelInterpreter','none');
imagesc(comp_meas);%colormap gray
set(gca,'XTick',1:length(measures),'XTickLabel',measures)
set(gca,'YTick',1:length(measures),'YTickLabel',measures)
xtickangle(45);colorbar;

figure
[nr,nc]=size(ciPLV);
if nr==nc
    N=max(nr,nc);
    Isubdiag = find(tril(ones(N),-1));
    ciPLV_vec=ciPLV(Isubdiag);
    wPLI_db_ft_vec=wPLI_ft(Isubdiag);
end
scatter(ciPLV_vec,wPLI_db_ft_vec);xlim([-.05 1.05]);ylim([-.05 1.05])
xlabel('ciPLV');ylabel('debiased wPLI')