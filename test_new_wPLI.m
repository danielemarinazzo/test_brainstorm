clear
load('HA_HB_iA_iB.mat')
[nA,~]=size(HA);
[nB,ntime]=size(HB);
R = reshape(mean(sin(angle(HA(iA,:)')-angle(HB(iB,:)')))' ./ mean(abs(sin(angle(HA(iA,:)')-angle(HB(iB,:)'))))',[],nB);
R=abs(R);
phaseA = HA ./ abs(HA);
phaseB = HB ./ abs(HB);
num = abs(imag(phaseA*phaseB'));
den = zeros(nA,nB);
for t = 1:ntime
    den = den + abs(imag(phaseA(:,t) * phaseB(:,t)'));
end
wPLI_db = num./den;
c=corr(wPLI_db(:),abs(R(:)), 'rows','complete')