% Routine to multiplex a signal 3-C component
% Input
% sig: signal with 3 components
% Output
% data: signal multiplexed

function data = mx(sig)

% Check size of data
[n,m] = size(sig);
if n ~= 3 && m ~= 3; disp('Signal isn t a 3-C signal. Stop'); return; end
if n ~= 3; sig = sig'; end
[n,m] = size(sig); % n: 3 comp; m: length of signal

% Multiplexing
kk = 1;
data = zeros(3*m,1);
for ii=1:m
    for jj = 1:n
        data(kk,1) = sig(jj,ii);
        kk = kk + 1;
    end
end

