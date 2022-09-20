% this script tests several connectivity measures in a time-resolved way,
% i.e. computing phase consistency across trials per each time point.
% The simulated data represent coupling starting at 1/3 of the trial,
% either with instantaneous coupling, or with lagged one. 

clear
load test_PLI_timevarying

% choose here "delay" or "nodelay". datah is already filtered and hilbert transformed
 datah=datah_simulated_delay; %lagged influence, also corrected measures should change
% datah=datah_simulated_nodelay; %zero lag influence, only noncorrected measures should change
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[nc, ns, nt]=size(datah);
nA=2;nB=2;
iA = repmat(1:nA, 1, nB)';
iB = reshape(repmat(1:nB, nA, 1), [], 1);

ndat=datah ./ abs(datah);

Rs=complex(zeros(nt,nA*nB),0);ciRs=Rs;wRs=Rs;wRsnum=Rs;wRsden=Rs;wRsdiff=Rs;wRsdb=Rs;
% Brainstorm-style implementation
tic
for itrials=1:ns
    HA=squeeze(datah(:,itrials,:));HB=squeeze(datah(:,itrials,:));
    phaseA = HA(iA,:) ./ abs(HA(iA,:));
    phaseB = HB(iB,:) ./ abs(HB(iB,:));
    Rs=Rs+(phaseA .* conj(phaseB))';
    ciRs=ciRs+(imag((phaseA .* conj(phaseB))'))./(real((phaseA .* conj(phaseB))')/nt).*conj((real(phaseA .* conj(phaseB))')/nt);
    cdd = phaseA .* conj(phaseB);
    cdi = imag(cdd);
    wRs=wRs+(abs(cdi).*sign(cdi))'./abs(cdi)';
    wRsnum=wRsnum+(abs(cdi).*sign(cdi))';
    wRsdiff=wRsdiff-(((abs(cdi).*sign(cdi))').^2)/nt;
    wRsden=wRsden+abs(cdi)';
end
toc

PLV=abs(Rs(:,2)/ns);
ciPLV=abs(ciRs(:,2)/ns);
wPLI=abs(wRs/ns);
wPLI_db=sqrt((wRsnum.^2-wRsdiff)./((wRsden.^2-wRsdiff)));
figure
plot(PLV)
hold on
plot(ciPLV,'r')
plot(wPLI(:,2),'g')
plot(wPLI_db(:,2),'m')
ylim([0 1]);
legend('PLV', 'ciPLV','wPLI','wPLI\_db');

