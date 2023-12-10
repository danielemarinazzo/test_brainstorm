% this script compares different connectivity measures in the static case
% (variables and time points, single trial), starting from raw data (68
% regions, 5801 time points), sampling frequency 600 Hz

clear;clc

load('raw_data.mat')
data=data';
[b, a] = butter(3, [13, 30] / (fs/ 2), 'bandpass'); % design a filter in beta band
%%
% either we first filter then hilbert
filtered_beta=filtfilt(b, a, data);
HA=hilbert(filtered_beta);HA=HA';
HB=HA;

% or filter and hilbert in fourier space
%HA = my_filtfilt ( b, a, data, 'true');HA=HA';
%HB=HA;

%%
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
for t=1:nt
    cdi=imag(HA(:,t) * HB(:,t)');
    num=num+(abs(cdi).*sign(cdi))';
    den=den+abs(cdi)';
end
wPLI_sc=abs(num/nt)./(den/nt);
t=toc;
disp(['wPLI sign cdi, ' num2str(t) ' seconds']);

% wPLI ratio imag csd
tic
num = imag(HA*HB');
den = zeros(nA,nB);
for t = 1:nt
    den = den + abs(imag(HA(:,t) * HB(:,t)'));
end
wPLI_csdrat = abs(num./den);
t=toc;
disp(['wPLI ratio imag csd, ' num2str(t) ' seconds']);

% wPLI debiased ratio imag csd
tic
num = imag(HA*HB');
den = zeros(nA,nB);sqd = zeros(nA,nB);
for t = 1:nt
    den = den + abs(imag(HA(:,t) * HB(:,t)'));
    sqd = sqd + imag(HA(:,t)*HB(:,t)').^2;
end
wPLI_db_csdrat = (num.^2-sqd)./(den.^2-sqd);
t=toc;
disp(['wPLI debiased ratio imag csd, ' num2str(t) ' seconds']);

% "fieldtrip-like" wPLI (there they do it for many trials though)
% https://github.com/fieldtrip/fieldtrip/blob/master/connectivity/ft_connectivity_wpli.m
tic
num = zeros(nA,nB);
den = zeros(nA,nB);
for t = 1:nt
    num = num + imag(HA(:,t)*HB(:,t)');
    den = den + abs(imag(HA(:,t) * HB(:,t)'));
end
wPLI_ft = abs(num./den);
t=toc;
disp(['wPLI fieldtrip, ' num2str(t) ' seconds']);

% "fieldtrip-like" debiased wPLI (there they do it for many trials though)
% https://github.com/fieldtrip/fieldtrip/blob/master/connectivity/ft_connectivity_wpli.m
tic
num = zeros(nA,nB);
den = zeros(nA,nB);
sqd = zeros(nA,nB);
for t = 1:nt
    num = num + imag(HA(:,t)*HB(:,t)');
    den = den + abs(imag(HA(:,t) * HB(:,t)'));
    sqd = sqd + imag(HA(:,t)*HB(:,t)').^2;
end
wPLI_db_ft = (num.^2-sqd)./(den.^2-sqd);
t=toc;
disp(['wPLI db fieldtrip, ' num2str(t) ' seconds']);


measures={'PLV','ciPLV','wPLI_sc',...
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