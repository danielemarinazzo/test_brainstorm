load('raw_data.mat')
data=data(1,2001:3200);
time=0:1/fs:length(data)/fs-1/fs;

Ws=[8, 12]; %bandpass
[b, a] = butter(3, Ws / (fs/ 2), 'bandpass'); % design a filter

% filter than Hilbert
filtered=filtfilt(b, a, data);
HF=abs(hilbert(filtered));

% Hilbert in frequency
FH=abs(my_filtfilt ( b, a, data', 'true')); FH=FH';

plot(time,FH,'r','DisplayName','FH');hold on;
plot(time,HF,'b','DisplayName','HF');
plot(time,data,'Color',[.7 .7 .7],'DisplayName','data');hold off;
legend