% Build matrix for SVD of subspace detection

function [events,A,ind3,stack,stack_d,stack_d2] = buildmtx(events,C0,L0,ind2,dt,cutl)

% Choose a master event
[value,index] = max(mean(C0(ind2,ind2))); % Choose event which is, in average, the most correlated to all the others
ts=max(abs(L0(ind2(index),ind2))); % Get lags corresponding to master event, then get the maximum one

% Do a quick quality control check - if abs(lag) > 1 sec. to master, throughout
% to avoid side lobe problems / early S wave alignment issues
for i = 1:length(ind2) % If more than >1 then remove them from further analysis (svd)
    if abs(L0(ind2(index),ind2(i))) >= 1
        events{ind2(i)}.hasSeis = 2;
    end
end

ind3=[]; % Get events that will be selected for SVD
for i=1:length(events)
    if events{i}.hasSeis == 3
        ind3 = [ind3 i];
    end
end

master = ind2(index); % Event ID# of master event

% Do a reality check of the events, decide if some should be removed:
figure; hold on
for i=1:length(ind3)
    plot((1:length(events{ind3(i)}.signal.d))+(L0(master,ind3(i))/dt),events{ind3(i)}.signal.d./max(abs(events{ind3(i)}.signal.d))+(1.0*i-1),'k','linewidth',1)
end
title('Events included'); xlabel('Time (samples relative to master event)'); ylabel('Events')

% find the maximum lag and add to "signal" length of A
aft = max(L0(master,ind3)); % find the max length of lag to master event
bef = abs(min(L0(master,ind3))); % find the min lag to master event
if bef==0; bef=dt; end % In case master event is the first one min(L0)=0 (generates error later)

aft_dt = aft/dt; aft_dt = round(aft_dt);
bef_dt = bef/dt; bef_dt = round(bef_dt);
k=bef_dt+1;
signal_length = cutl+(bef_dt+aft_dt);
A=zeros(cutl+aft_dt+bef_dt,length(ind3));
%another value to fill with zeros the column of A with missing values
%for solving the problem with the K04 station
lr =cutl+aft_dt+bef_dt;
for i = 1:length(ind3)
    l = L0(master,ind3(i))/dt; l = round(l);
    m = l+k;
    if m == 0
        ll = length(events{ind3(i)}.signal.d);
        A(1:ll,i) = events{ind3(i)}.signal.d(:)'; clear ll
    else
        l1 = length(zeros(1,m-1));
        l2 = length(events{ind3(i)}.signal.d);
        if isempty(l1)==1; l1=0; end
        A(1:l1+l2,i) = [zeros(1,m-1) events{ind3(i)}.signal.d(:)']; clear l1 l2
    end
    clear m l
end
clear k
A = A(bef_dt:end-aft_dt-1,:);

% Compute the stack of these events for later use
stack = sum(A')./length(ind3);

% % Incorporate the phase-weighted stack (Thurber et al., 2014)
% dhilb = hilbert(A); iph = angle(dhilb);
% spha = sum(iph,2)/size(iph,2); %spha = spha.^2;
% statmp = stack.*spha';
% statmp = statmp/max(abs(statmp));
% stack_p = statmp; clear dhilb iph spha statmp

% Compute the time derivatives of the stack
stack_d = diff(stack);
stack_d2 = diff(stack,2);
