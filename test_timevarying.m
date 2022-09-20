% this script tests several connectivity measures in a time-resolved way,
% i.e. computing phase consistency across trials per each time point.
% The simulated data represent coupling starting at 1/3 of the trial,
% either with instantaneous coupling, or with lagged one.

clear
load test_PLI_timevarying

% choose here "delay" or "nodelay". datah is already filtered and hilbert transformed
datah=datah_simulated_delay; %lagged influence, also corrected measures should change
%datah=datah_simulated_nodelay; %zero lag influence, only noncorrected measures should change
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[nc, ns, nt]=size(datah); %channels per trials per timepoints
nA=2;nB=2;
iA = repmat(1:nA, 1, nB)';
iB = reshape(repmat(1:nB, nA, 1), [], 1);

% Brainstorm-style implementation - initialize the "R" structures
R_PLV=complex(zeros(nt,nA*nB),0);R_ciPLV=R_PLV;
R_WPLI=zeros(nc,nc,nt);R_WPLI_db=zeros(nc,nc,nt);


tic
for itrials=1:ns
    HA=squeeze(datah(:,itrials,:));HB=squeeze(datah(:,itrials,:));
    phaseA = HA(iA,:) ./ abs(HA(iA,:));
    phaseB = HB(iB,:) ./ abs(HB(iB,:));
    R_PLV=R_PLV+(phaseA .* conj(phaseB))';
end
PLV=abs(R_PLV(:,2)/ns);
t=toc;disp(['PLV, ' num2str(t) ' seconds']);
tic
for itrials=1:ns
    HA=squeeze(datah(:,itrials,:));HB=squeeze(datah(:,itrials,:));
    phaseA = HA(iA,:) ./ abs(HA(iA,:));
    phaseB = HB(iB,:) ./ abs(HB(iB,:));
    R_ciPLV=R_ciPLV+(imag((phaseA .* conj(phaseB))'))./(real((phaseA .* conj(phaseB))')/nt).*conj((real(phaseA .* conj(phaseB))')/nt);
end
ciPLV=abs(R_ciPLV(:,2)/ns);
t=toc;disp(['ciPLV, ' num2str(t) ' seconds']);
tic
for i=1:nc
    for j=1:i
        if i~=j
            HA=squeeze(datah(i,:,:));HB=squeeze(datah(j,:,:));
            cdi = imag(HA .* conj(HB));
            R_WPLI(i,j,:)=abs( mean( abs(cdi).*sign(cdi) ,1) )./mean(abs(cdi),1);
            R_WPLI(j,i,:)=R_WPLI(i,j,:);
        end
    end
end
wPLI=squeeze(R_WPLI(1,2,:));
t=toc;disp(['wPLI, ' num2str(t) ' seconds']);


tic
for i=1:nc
    for j=1:i
        if i~=j
            HA=squeeze(datah(i,:,:));HB=squeeze(datah(j,:,:));
            cdi = imag(HA .* conj(HB));
            imagsum      = sum(cdi,1);
            imagsumW     = sum(abs(cdi),1);
            debiasfactor = sum(cdi.^2,1);
            R_WPLI_db(i,j,:)= (imagsum.^2 - debiasfactor)./(imagsumW.^2 - debiasfactor);
            R_WPLI_db(j,i,:)=R_WPLI_db(i,j,:);
        end
    end
end
wPLI_db=squeeze(R_WPLI_db(1,2,:));
t=toc;disp(['wPLI debiased, ' num2str(t) ' seconds']);


figure
plot(PLV,'b')
hold on
plot(ciPLV,'r')
plot(wPLI,'g')
plot(wPLI_db,'m')
ylim([0 1]);
legend('PLV', 'ciPLV','wPLI','wPLI debiased');

