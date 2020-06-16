% Do the cross-correlation between parent events to work out potential
% clusters
% cc: cut-off for linkage between events, the smaller it is, the closer you 
% are (so you are in the same cluster)
% opt: 'l' for complete signal (signal variable), 's' for short signal
% (clip variable)
% hs: hasSeis label for further processing
% v: choose set of events (1 is the main cluster with most events, and then
% the other clusters with less events increasing this number)

function [events,C,L,ind,cutl] = cclink(events,dt,cc,opt,hs,v)

% Cross-correlation part
if strcmp(opt,'l') == 1; cutl = length(events{1}.signal.d); end
if strcmp(opt,'s') == 1; cutl = length(events{1}.clip.d); end

if mod(cutl,2)~=0; n = cutl+1; else n=cutl; end

if strcmp(opt,'l') == 1; halfpt = (n/2); end
if strcmp(opt,'s') == 1; halfpt = (n/2)-1; end

tmax = (halfpt+1)*dt;
tplus = dt*(0:halfpt);
tmins = tplus(2:halfpt) - tmax;

if strcmp(opt,'l') == 1; lags = [tplus tmins]; end
if strcmp(opt,'s') == 1; lags = [tplus tmins -dt 0]; end % JB: added 0

C=zeros(length(events),length(events)); % xcorr matrix
L=zeros(length(events),length(events)); % lags matrix


for i = 1:length(events)
     
    if events{i}.hasSeis == 1
        if strcmp(opt,'l') == 1; 
            X = fft(events{i}.signal.d,cutl);
            x = events{i}.signal.d;
        end
        if strcmp(opt,'s') == 1; 
            X = fft(events{i}.clip.d,cutl);
            x = events{i}.clip.d;
        end
        
        for j = 1:length(events)            
            if events{j}.hasSeis == 1
                if strcmp(opt,'l') == 1;
                    Y = fft(events{j}.signal.d,cutl);
                    y = events{j}.signal.d;
                end
                if strcmp(opt,'s') == 1;
                    Y = fft(events{j}.clip.d,cutl);
                    y = events{j}.clip.d;
                end
                
                c = real(ifft(X .* conj(Y)))/sqrt((sum(x.*x)*sum(y.*y)));                
                [cmax,maxlag] = max(c);                
                C(i,j) = cmax;                
                L(i,j) = lags(maxlag);                
                clear cmax maxlag Y y c
            end
        end
        clear X x
    end
end

% Work out the linkage
disp('Computing the Linkage')
D = 1 - C;
for i = 1:length(events)
    D(i,i) = 0;
end
Y = squareform(D); % Convert NxN matrix to vector of lower half matrix values
Z = linkage(Y); % Do a matrix of hierarchical clusters of the rows of Y

T = cluster(Z,'cutoff',cc,'criterion','distance'); % Clustering: use distance as a criterion
counts=histc(T,1:length(T));
[~,val]=sort(counts,'descend');
% v=1:5; % pull out the five largest sets
[ind]=find(T==val(v(1))); % Find events in first set of events (largest set)

disp(['Number of events in this group = ' num2str(length(ind))])

% Relabel these hasSeis
for i = 1:length(ind)
    events{ind(i)}.hasSeis = hs;
end

