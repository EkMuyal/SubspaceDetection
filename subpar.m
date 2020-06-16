% Design of the empirical subspace paraneters

function [s_mat,s,U1,U2,W] = subpar(U,S,stack,stack_d,A,nv,f,dt)

[~,na] = size(A);
if nv > na 
    nv = na;
    disp(['The value of nv was to high, the new value of nv = ' num2str(na)])
end
clear na

% s = STACK
% U2 = nv SINGULAR VECTORS
% s_mat = EMPIRICAL SUBSPACE
% Take each section of the continuous record and scan the subspace detector
% through the sections - normalize all the amplitudes!

U1 = U(:,1);%./max(abs(U(:,1)));
U2 = U(:,1:nv);%max(max(abs(U(:,1:2))));

% Plot of stack, U1 and other singular vectors
figure; hold on;
plot((1:length(stack))*dt,stack/max(abs(stack)),'r','linewidth',1.5)
plot((1:length(U1))*dt,U1/max(abs(U1)) +1,'k','linewidth',1)

for jj=2:size(U2,2)
    plot((1:length(U2))*dt,U2(:,jj)/max(abs(U2(:,jj))) +jj,'b','linewidth',1)
end
title('Stack and sing. vectors of SVD decomp'); xlabel('Time (s)'); 
set(gca,'YTick',0:2); set(gca,'YTickLabel',{'stack' 'U1' 'U2,2'});

s = stack./norm(stack); % Normalized stack
s_d = stack_d./norm(stack_d); % Normalized 1st deriv of stack
%  s_p = stack_p./norm(stack_p); % Normalized phase-weighted stack
s_mat = [s; [s_d 0]];
s = s';
s_mat = s_mat';
signal_length = length(U1);
% Prepare the Taper
window = hanning(2*floor(signal_length*f),'periodic');
lwin = length(window);
half_win = window(1:lwin/2);
half_win2 = fliplr(half_win');
flat_win = ones(signal_length-length(window),1);
W = [half_win; flat_win; half_win2'];

% Weight the vectors:
s1=S(1,1)/(S(1,1)+S(2,2));
s2=S(2,2)/(S(1,1)+S(2,2));
sm=0; for ii=1:nv; sm = sm+S(ii,ii); end
for ii=1:nv
    st = S(ii,ii)/sm; 
    U2(:,ii) = U2(:,ii)*st; % Singular vectors
    U2(:,ii) = U2(:,ii)./norm(U2(:,ii)); % JB: replaced U by U2
    clear st
end

display(['Representativeness of the selected ' num2str(nv)...
    ' SVD vectors is ' num2str(sm/sum(diag(S(1:size(S,2),:)))*100)])

s_mat(:,1) = s_mat(:,1)*s1; % stack
s_mat(:,2) = s_mat(:,2)*s2; % 1st deriv of stack
% s_mat(:,3) = s_mat(:,3)*s1; % phase weighted stack

s_mat(:,1) = s_mat(:,1)./norm(s_mat(:,1));
s_mat(:,2) = s_mat(:,2)./norm(s_mat(:,2));
% s_mat(:,3) = s_mat(:,3)./norm(s_mat(:,3));
clear sm s1 s2
