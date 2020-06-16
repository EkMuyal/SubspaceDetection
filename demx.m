% Routine to de-multiplex a signal multiplexed
% Input
% sig: multiplexed signal
% nc: number of components in the multiplexed signal
% Output
% data: signal demultiplexed

function data = demx(sig,nc)

% Check if it can be de-multiplexed
test = length(sig)/nc;
test = test - round(test);
if test ~= 0; 
    disp('Signal length not a factor of nc. Cannot be de-multiplexed');
    data = [];
    return
end

% De-multiplexing
kk = 1;
data = zeros(nc,length(sig)/nc);
for ii=1:nc:length(sig)
    for jj = 1:nc
        data(jj,kk) = sig(ii+jj-1);
    end
    kk = kk + 1;
end

