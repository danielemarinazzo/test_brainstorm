% this script compares ciPLV and wPLI starting from raw data (68
% regions, 5801 time points), sampling frequency 600 Hz, for different
% window lengths.

clear;clc
figure;
load('raw_data.mat')
datatot=data;
window_lengths=[300, 600, 1200, 4800];
for i_length=1:length(window_lengths)
    data=datatot(:,1:window_lengths(i_length))';
    Ws=[4, 8]; % define band pass
    [b, a] = butter(3, Ws / (fs/ 2), 'bandpass'); % design a filter
    hilbert_in_freq=1; % whether you want to Hilbert in frequency domain, together with filtering

    if hilbert_in_freq
        HA = my_filtfilt ( b, a, data, 'true');HA=HA';
        HB=HA;
    else
        filtered=filtfilt(b, a, data);
        HA=hilbert(filtered);HA=HA';
        HB=HA;
    end


    %%
    [nA,~]=size(HA);
    [nB,nt]=size(HB);
    phaseA = HA ./ abs(HA);
    phaseB = HB ./ abs(HB);



    % ciPLV
    tic
    csd=phaseA*phaseB';
    ciPLV=abs((imag((csd))/nt)./sqrt(1-(real((csd))/nt).^2));
    t=toc;
    disp(['ciPLV, ' num2str(t) ' seconds']);


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



    subplot(2,2,i_length)
    [nr,nc]=size(ciPLV);
    if nr==nc
        N=max(nr,nc);
        Isubdiag = find(tril(ones(N),-1));
        ciPLV_vec=ciPLV(Isubdiag);
        wPLI_db_csdrat_vec=wPLI_db_csdrat(Isubdiag);
    end
    scatter(ciPLV_vec,wPLI_db_csdrat_vec);xlim([-.05 1.05]);ylim([-.05 1.05])
    xlabel('ciPLV');ylabel('debiased wPLI')
    title([num2str(window_lengths(i_length)) ' points,' num2str(window_lengths(i_length)/fs) ' seconds'])
end