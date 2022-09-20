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

R_PLV=complex(zeros(nt,nA*nB),0);R_ciPLV=R_PLV;R_WPLI=R_PLV;
R_dbWPLI_num=R_PLV;R_dbWPLI_den=R_PLV;R_dbWPLI_sqd=R_PLV;
% Brainstorm-style implementation
tic
for itrials=1:ns
    HA=squeeze(datah(:,itrials,:));HB=squeeze(datah(:,itrials,:));
    phaseA = HA(iA,:) ./ abs(HA(iA,:));
    phaseB = HB(iB,:) ./ abs(HB(iB,:));
    R_PLV=R_PLV+(phaseA .* conj(phaseB))';
    R_ciPLV=R_ciPLV+(imag((phaseA .* conj(phaseB))'))./(real((phaseA .* conj(phaseB))')/nt).*conj((real(phaseA .* conj(phaseB))')/nt);
    cdd = phaseA .* conj(phaseB);
    cdi = imag(cdd);
    R_WPLI=R_WPLI+(abs(cdi).*sign(cdi))'./abs(cdi)';
    R_dbWPLI_num=R_dbWPLI_num+(abs(cdi).*sign(cdi))';
    R_dbWPLI_sqd=R_dbWPLI_sqd-(((abs(cdi).*sign(cdi))').^2)/nt;
    R_dbWPLI_den=R_dbWPLI_den+abs(cdi)';
end
toc

% here just abs and divide by the number of trials, except for the debiased
% wPLI. We should find a solution for this latter part for implementation in
% Brainstorm 

PLV=abs(R_PLV(:,2)/ns);
ciPLV=abs(R_ciPLV(:,2)/ns);
wPLI=abs(R_WPLI/ns);
wPLI_db=sqrt((R_dbWPLI_num.^2-R_dbWPLI_sqd)./((R_dbWPLI_den.^2-R_dbWPLI_sqd)));
figure
plot(PLV)
hold on
plot(ciPLV,'r')
plot(wPLI(:,2),'g')
plot(wPLI_db(:,2),'m')
ylim([0 1]);
legend('PLV', 'ciPLV','wPLI','wPLI debiased');

