% Modified by Camilo De La Hoz Lozano - JB Tary, 2019
% Subspace Detection of microseismicity from Sea of Marmara
% First step of subspace detection for the Marmara 2011 dataset
% Modified from Greg Beroza files
% This file is for building templates from detected/catalog events and
% doing the subspace detection
% One can use station specific parameters (filtering and clustering), as
% long as the singular vectors are calculated with the same ones, and that
% continuous data also use the same parameters
% Optimum params G01: f [5 15] c 0.2 -> 27 in group <- all # might be wrong
% Optimum params G03: f [5 15] c 0.2 -> 16 in group
% Optimum params G05: f [5 20] c 0.4 -> 11 in group (to recut A)
% Optimum params G06: f [5 15] c 0.4 -> 6 in group
% Optimum params K04: f [5 15] c 0.25 -> 11 in group
% Optimum params G07: f [5 20] c 0.5 -> 6 in group

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                 PROCESS EVENT WAVEFORMS                             %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP: 1

% Events should be identified prior to this script.
clear
path = 'ParentEvents';
parents = dir([path '/*.mat']);
stan = 'G01'; % For a start, pick a station

display('Reading in parent event Waveforms')
for i=1:length(parents)
    %%% THIS PATH NEED TO LEAD TO THE WAVEFORM DIRECTORIES
    pathtmp = [path '/' parents(i).name];
    events{i}.seis = load(pathtmp);
    events{i}.hasSeis = 1;
    events{i}.dt = 1/events{i}.seis.fsmaster;
    clear pathtmp
end

% Save pristine events variable for further processing
evsv = events;

% Process and Filter waveforms
opt = 'Z';
nPoles = 4;
fcl = 5; fch = 15;
winP = [11 586]; % Shorter window at beginning of parent events containing the P wave
[events,dt,b,a] = prepro(events,stan,nPoles,fcl,fch,winP(1),winP(2),opt);

% Do the long signal cross correlation and linkage to find a candidate group
c = 0.2; % Loose cut-off: 0.5
opt2 = 'l'; hs = 3; % Or hs=2; but only those with hs=3 will be selected for SVD
cluster = 1; % Subset (cluster) of events to select (1 is the main one)
[events,C,L,ind,cutl] = cclink(events,dt,c,opt2,hs,cluster);

% % Redo the cross correlation, alignment and linkage to improve the group
% % Used for the SVD
% c = 0.2; % Stricter cut-off (0.2)
% opt2 = 's'; hs = 3;
% [events,C0,L0,ind2] = cclink(events,dt,c,opt2,hs); clear opt2

figure; 
subplot(2,2,1); imagesc(1:length(C),1:length(C),C); title('C')
subplot(2,2,2); imagesc(1:length(L),1:length(L),L); title('L')
% subplot(2,2,3); imagesc(1:length(C0),1:length(C0),C0); title('C0')
% subplot(2,2,4); imagesc(1:length(L0),1:length(L),L0); title('L0')

% Choose a master event and build matrix for SVD

[events,A,ind3,stack,stack_d,stack_d2] = buildmtx(events,C,L,ind,dt,cutl);

figure; hold on
for i=1:length(ind3)
    plot((1:length(A))*dt,A(:,i)./max(abs(A(:,i)))+(1.0*i-1),'k','linewidth',1)
end
i=i+1; plot((1:length(stack))*dt,stack./max(abs(stack))+(1.0*i-1),'r','linewidth',1.5)
i=i+1; plot((1:length(stack_d))*dt,stack_d./max(abs(stack_d))+(1.0*i-1),'b','linewidth',1.5)
i=i+1; plot((1:length(stack_d2))*dt,stack_d2./max(abs(stack_d2))+(1.0*i-1),'b','linewidth',1)
% i=i+1; plot((1:length(stack_p))*dt,stack_p./max(abs(stack_p))+(1.0*i-1),'g','linewidth',1.5)
title('Events included for stacking / stack / 1st deriv / 2nd deriv'); 
xlabel('Time (sec.)'); ylabel('Events and stacks')

%% Selection of events for all stations in network
% events variable will stay the "master" variable, we will now build the
% same variables for other stations
% STEP: 5

stans = ['G03';'G05';'G06';'K04';'G07'];

for ii = 1:size(stans,1)
    % Pre-process the parent events
    [evtmp,dttmp,~,~] = prepro(evsv,stans(ii,:),nPoles,fcl,fch,winP(1),winP(2),opt);
    [evtmp,Ctmp,Ltmp,~,cutltmp] = cclink(evtmp,dttmp,1,opt2,1,cluster);
    for jj = 1:length(ind)
        if ~isfield(evtmp{ind(jj)},'signal'); continue; end
        evtmp{ind(jj)}.hasSeis = hs;
    end
    
    % Get the variables for SVD
    [evtmp,Atmp,~,stacktmp,stack_dtmp,~] = buildmtx(evtmp,Ctmp,Ltmp,ind,dttmp,cutltmp);
    
    eval(['events' num2str(ii) ' = evtmp;']);
    eval(['dt' num2str(ii) ' = dttmp;']);
    eval(['A' num2str(ii) ' = Atmp;']);
    eval(['stack' num2str(ii) ' = stacktmp;']);
    eval(['stack_d' num2str(ii) ' = stack_dtmp;']);
    
    figure; hold on
    plot((1:length(stacktmp))*dt,stacktmp./max(abs(stacktmp))+1,'r','linewidth',1.5)
    plot((1:length(stack_dtmp))*dt,stack_dtmp./max(abs(stack_dtmp))+2,'b','linewidth',1.5)
    title('Stack (red) / 1st deriv (blue)');
    xlabel('Time (sec.)'); ylabel('Stacks')
    
    clear *tmp
end

%% Do the SVD
% STEP: 2

% MATRIX A IS ALL THE WAVEFORMS
% U are the left singular vectors ** these are the important ones
% S are the singular values
% V are the right singular vectors
[U,S,V] = svd(A);

% Design of the empirical subspace detector 
nv = 2; % Number of singular vectors to select
f = 0.05; % input portion of signal to be tapered
[s_mat,s,U1,U2,W] = subpar(U,S,stack,stack_d,A,nv,f,dt); 

save(['svd_marmara_' stan '.mat'],'-v7.3');

%% Do the SVD for the other stations of the network
% STEP: 6

for ii = 1:size(stans,1)
    
    eval(['Atmp = A' num2str(ii) ';']);
    eval(['stacktmp = stack' num2str(ii) ';']);
    eval(['stack_d = stack_d' num2str(ii) ';']);
    eval(['dttmp = dt' num2str(ii) ';']);
    
    % Design of the empirical subspace detector
    nv = 2; % Number of singular vectors to select
    f = 0.05; % input portion of signal to be tapered
    [Utmp,Stmp,Vtmp] = svd(Atmp);
    [s_mattmp,stmp,U1tmp,U2tmp,Wtmp] = subpar(Utmp,Stmp,stacktmp,stack_d,Atmp,nv,f,dttmp);
    
    eval(['U_' num2str(ii) ' = Utmp;']);
    eval(['S' num2str(ii) ' = Stmp;']);
    eval(['V' num2str(ii) ' = Vtmp;']);
    eval(['s_mat' num2str(ii) ' = s_mattmp;']);
    eval(['s' num2str(ii) ' = stmp;']);
    eval(['U1' num2str(ii) ' = U1tmp;']);
    eval(['U2' num2str(ii) ' = U2tmp;']);
    eval(['W' num2str(ii) ' = Wtmp;']);
    
    clear *tmp
end

clear events evsv
save('svd_marmara_allsta.mat','-v7.3');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%            APPLY THE DETECTOR TO CONTINUOUS RECORDS                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP: 3

% Read and filter the continous data using the same parameters as we used 
% to form the design set
% Continuous data inside directories with station names (ex: 'G01')
pathcont = 'Data';
clear events evsv parents % This variable just take too much space

stan = 'G01'; stans = ['G03';'G05';'G06';'K04';'G07'];
stanc = [stan;stans];
for ll = 1:size(stanc,1)
for ii = 25:31 % Loop on the days
    mkdir([num2str(ii) '_07_2011'])
    cd([num2str(ii) '_07_2011'])
    
    datemin = datenum([2011,07,ii,0,0,0]); % Start date of data for analysis
    datemax = datenum([2011,07,ii,23,0,0]); % End date of data for analysis
    
    hours = preprocont(pathcont,stanc(ll,:),datemin,datemax,dt,b,a,opt);
    hours_x = preprocont(pathcont,stanc(ll,:),datemin,datemax,dt,b,a,'X');
    hours_y = preprocont(pathcont,stanc(ll,:),datemin,datemax,dt,b,a,'Y');
    
    % This is how many samples you skip between each correlation
    dstep = 10; % Increase this number to speed things up
    
    % Save all variables except 'events', but including hours inside a new variable
    save(['marmara_svd_knw_' stanc(ll,:) '.mat'],'-v7.3');
    
    clear datemin datemax hours*
    cd ..
end
end

%% Do the detection on the other stations of the network
% STEP: 7

clear events* % This variable just take too much space
% Continuous data inside directories with station names (ex: 'G01')
pathcont = 'Data';
dstep = 10; % Increase this number to speed things up

for ii = 25:31 % Loop on the days
    datemin = datenum([2011,07,ii,0,0,0]); % Start date of data for analysis
    datemax = datenum([2011,07,ii,23,0,0]); % End date of data for analysis
    
    if ~exist([num2str(ii) '_07_2011'],'dir'); mkdir([num2str(ii) '_07_2011']); end
    cd([num2str(ii) '_07_2011'])
    
    for jj = 1:size(stans,1)
        eval(['dttmp = dt' num2str(jj) ';']);
        hours = preprocont(pathcont,stans(jj,:),datemin,datemax,dttmp,b,a,opt);
        hours_x = preprocont(pathcont,stans(jj,:),datemin,datemax,dttmp,b,a,'X');
        hours_y = preprocont(pathcont,stans(jj,:),datemin,datemax,dttmp,b,a,'Y');
        
        save(['marmara_svd_knw_' stans(jj,:) '.mat'],...
            'hours','hours_x','hours_y','a','b','dstep',...
            ['U1' num2str(jj)],['W' num2str(jj)],...
            ['U2' num2str(jj)],['s' num2str(jj)],...
            ['s_mat' num2str(jj)],'-v7.3');
        clear *tmp hours*
    end
    
    cd ..
end

%% Detection statistics part using stack, subspace and empirical subspace
% STEP: 8

clear
stan = 'G01'; % Reference station
for ii = 26 % Days used, example use day 26
    cd([num2str(ii) '_07_2011'])
    
    load(['marmara_svd_knw_' stan '.mat'],'hours',...
        'a','b','dstep','U1','W','U2','s','s_mat')
    
    disp('Scanning ...');
    for jj=1:length(hours)
        
        display(hours{jj}.name);
        if hours{jj}.hasSeis == 1
            N = hours{jj}.signal.d;
            for q=1:floor(((length(N)-length(U1))/dstep)+1)
                y = N(((q-1)*dstep)+1:((q-1)*dstep+length(U1))); % the cut
                y1 = detrend(y,'constant'); % remove the mean
                y2 = filtfilt(b,a,y1); % the filter
                y2 = y2(:);
                y3 = y2.*W; % taper
                y4 = y3./norm(y3); % fix the norm
                z(q) = (y4'*s*s'*y4)/(y4'*y4); % run the stack
                z1(q) = (y4'*U2*U2'*y4)/(y4'*y4); % run the U1/U2
                z2(q) = (y4'*s_mat*s_mat'*y4)/(y4'*y4); % run the stack and 1st deriv
                clear y y1 y2 y3 y4
            end
            disp('Saving ...')
            save(strcat(hours{jj}.name,'_knwn.mat'),'z','z1','z2','dstep');
            clear z z1 z2 N
        end
    end
    
    clear hours a b dstep U1 U2 W s s_mat
    cd ..
end
disp('Done')

%% Detection statistics for the other stations of the network
% STEP: 9

clear
stans = ['G03';'G05';'G06';'K04';'G07'];
for ii = 26 % Days used, example use day 26
    cd([num2str(ii) '_07_2011'])
    
    for ll = 1:size(stans,1)
        load(['marmara_svd_knw_' stans(ll,:) '.mat'],'hours',...
            'a','b','dstep',['U1' num2str(ll)],['W' num2str(ll)],...
            ['U2' num2str(ll)],['s' num2str(ll)],...
            ['s_mat' num2str(ll)])
        
        eval(['U1 = U1' num2str(ll) ';']);
        eval(['W = W' num2str(ll) ';']);
        eval(['U2 = U2' num2str(ll) ';']);
        eval(['s = s' num2str(ll) ';']);
        eval(['s_mat = s_mat' num2str(ll) ';']);
        
        disp('Scanning ...');
        for jj=1:length(hours)
            
            display(hours{jj}.name);
            if hours{jj}.hasSeis == 1
                N = hours{jj}.signal.d;
                for q=1:floor(((length(N)-length(U1))/dstep)+1)
                    y = N(((q-1)*dstep)+1:((q-1)*dstep+length(U1))); % the cut
                    y1 = detrend(y,'constant'); % remove the mean
                    y2 = filtfilt(b,a,y1); % the filter
                    y2 = y2(:);
                    y3 = y2.*W; % taper
                    y4 = y3./norm(y3); % fix the norm
                    z(q) = (y4'*s*s'*y4)/(y4'*y4); % run the stack
                    z1(q) = (y4'*U2*U2'*y4)/(y4'*y4); % run the U1/U2
                    z2(q) = (y4'*s_mat*s_mat'*y4)/(y4'*y4); % run the stack and 1st deriv
                    clear y y1 y2 y3 y4
                end
                disp('Saving ...')
                save(strcat(stans(ll,:),'_sta',num2str(ll),'_',hours{jj}.name,'_knwn.mat'),'z','z1','z2','dstep');
                clear z z1 z2 N
            end
        end
        
        clear hours a b dstep U1* U2* W* s_mat*
        clearvars s* -except stans
    end
    
    cd ..
end

disp('Done')

%% Estimate the average moveout for the other stations
% STEP: 10

clear
path = 'ParentEvents';
parents = dir([path '/*.mat']);
stan = 'G01'; % Reference station

display('Reading in parent event Waveforms')
for i=1:length(parents)
    pathtmp = [path '/' parents(i).name];
    events{i}.seis = load(pathtmp);
    events{i}.hasSeis = 1;
    events{i}.dt = 1/events{i}.seis.fsmaster;
    clear pathtmp
end

% Process and Filter waveforms
opt = 'Z';
nPoles = 4;
fcl = 5; fch = 15;
winP = [11 586]; % Shorter window at beginning of parent events containing the P wave
[events,dt,~,~] = prepro(events,stan,nPoles,fcl,fch,winP(1),winP(2),opt);

% Do the long signal cross correlation and linkage to find a candidate group
c = 0.2; % Loose cut-off: 0.5
opt2 = 'l'; hs = 3; % Or hs=2; but only those with hs=3 will be selected for SVD
cluster = 1; % Subset (cluster) of events to select (1 is the main one)
[events,~,~,ind,~] = cclink(events,dt,c,opt2,hs,cluster);
stans = ['G03';'G05';'G06';'K04';'G07'];
% Option to use only events with high CC on other stations to calculate timings
% 'withcc' or 'withoutcc'
option = 'withcc';

switch option
    case 'withcc'
        disp('Using CC with threshold 0.7 on other stations to get moveout')
        for ii = 1:length(ind)
            events2{ii} = events{ind(ii)};
            events2{ii}.hasSeis = 1;
        end
        
        for ii = 1:size(stans,1)
            [events3,dt,~,~] = prepro(events2,stans(ii,:),nPoles,fcl,fch,winP(1),winP(2),opt);
            [events3,~,~,ind2,~] = cclink(events3,dt,0.7,opt2,hs,cluster);
            
            kk = 1;
            for jj = 1:length(ind2)
                evtmp = events3{ind2(jj)};
                for ll = 1:size(evtmp.seis.stap); statmp{ll} = evtmp.seis.stap(ll,:); end
                [~,~,ib] = intersect(stans(ii,:),statmp);
                if isempty(ib); clear *tmp ib ll; continue; end
                if length(evtmp.seis.ppicks) ~= length(statmp); clear *tmp ib ll; continue; end
                
                picks(kk) = evtmp.seis.ppicks(ib); kk = kk + 1;
                clear *tmp ib ll
            end
            
            picks_final(ii,1) = median(picks);
            picks_final(ii,2) = length(picks);
            clear events3 ind2
        end
        
    case 'withoutcc'
        disp('Not using CC on other stations to get moveout, all events included')
        for ii = 1:size(stans,1)
            kk = 1;
            for jj = 1:length(ind)
                evtmp = events{ind(jj)};
                for ll = 1:size(evtmp.seis.stap); statmp{ll} = evtmp.seis.stap(ll,:); end
                [~,~,ib] = intersect(stans(ii,:),statmp);
                if isempty(ib); clear *tmp ib ll; continue; end
                if length(evtmp.seis.ppicks) ~= length(statmp); clear *tmp ib ll; continue; end
                
                picks(kk) = evtmp.seis.ppicks(ib); kk = kk + 1;
                clear *tmp ib ll
            end
            
            picks_final(ii,1) = median(picks);
            picks_final(ii,2) = length(picks);
            clear picks
        end
end

kk = 1;
for jj = 1:length(ind)
    evtmp = events{ind(jj)};
    for ll = 1:size(evtmp.seis.stap); statmp{ll} = evtmp.seis.stap(ll,:); end
    [~,~,ib] = intersect(stan,statmp);
    if isempty(ib); clear *tmp ib ll; continue; end
    if length(evtmp.seis.ppicks) ~= length(statmp); clear *tmp ib ll; continue; end
    
    picksref(kk) = evtmp.seis.ppicks(ib); kk = kk + 1;
    clear *tmp ib ll
end
sv = length(picksref); 
picksref = median(picksref); picksref(1,2) = sv; clear sv
mo = picks_final(:,1) - (ones(size(picks_final(:,1)))*picksref(1));
dstep = 10;

clear events* parents
save('moveout')

%% Sum the detection statistics for all stations
% STEP: 11

clear
drec = dir('*07_2011'); % All directories, one per day
load('moveout','mo','stan','stans','dt','dstep'); % Need moveout relative to reference station and sampling interval
tshift = round(mo*(1/dt)/dstep); % Time shifts relative to reference station in samples within z* variables

for ii = 1:length(drec)
    files = dir([drec(ii).name '/obs01*']); % Get files per hour
    disp(['Working on ' drec(ii).name])
    
    for jj = 1:length(files)
        % Get detection statistics (stack: z, subspace: z1, empirical subspace: z2)
        load([drec(ii).name '/' files(jj).name]); % From reference station (G01)
        
        ztmp(1,:) = z;
        z1tmp(1,:) = z1;
        z2tmp(1,:) = z2; clear z z1 z2
                
        for ll = 1:size(stans,1) % For the other stations
            filesta = dir([drec(ii).name '/*sta' num2str(ll) '*']); % For the other stations
            load([drec(ii).name '/' filesta(jj).name]);
            % Cut all time series to be from beginning to largest length of
            % U1 for concatenation (see code on detection above)
            if length(z) ~= length(ztmp)
                ztmp = ztmp(:,1:min([length(z) length(ztmp)]));
                z1tmp = z1tmp(:,1:min([length(z1) length(z1tmp)]));
                z2tmp = z2tmp(:,1:min([length(z2) length(z2tmp)]));
                z = z(1,1:min([length(z) length(ztmp)]));
                z1 = z(1,1:min([length(z1) length(z1tmp)]));
                z2 = z(1,1:min([length(z2) length(z2tmp)]));
            end
            ztmp(ll+1,:) = z;
            z1tmp(ll+1,:) = z1;
            z2tmp(ll+1,:) = z2; clear z z1 z2 filesta
        end        
        
        % Shift time series before summing (to correct for moveout of event)
        znew = zeros(size(ztmp,1),size(ztmp,2)+max(tshift));
        z1new = zeros(size(z1tmp,1),size(z1tmp,2)+max(tshift));
        z2new = zeros(size(z2tmp,1),size(z2tmp,2)+max(tshift));
        
        znew(1,1:length(ztmp)) = ztmp(1,:); % Reference station
        z1new(1,1:length(z1tmp)) = z1tmp(1,:); 
        z2new(1,1:length(z2tmp)) = z2tmp(1,:); 
        for yy = 2:length(tshift)+1
            znew(yy,tshift(yy-1)+1:tshift(yy-1)+size(ztmp,2)) = ztmp(yy,:);
            z1new(yy,tshift(yy-1)+1:tshift(yy-1)+size(z1tmp,2)) = z1tmp(yy,:);
            z2new(yy,tshift(yy-1)+1:tshift(yy-1)+size(z2tmp,2)) = z2tmp(yy,:);
        end
        
        zsum = sum(znew,1); % Sum all stations
        zmad = mad(zsum,1); % Median absolute deviation of xsum
        zsum = zsum(1,1:size(ztmp,2)); % Take only the part consistent with reference station (assuming tshift always > 0)
        
        z1sum = sum(z1new,1);
        z1mad = mad(z1sum,1);
        z1sum = z1sum(:,1:size(z1tmp,2));
        
        z2sum = sum(z2new,1);
        z2mad = mad(z2sum,1);
        z2sum = zsum(:,1:size(z2tmp,2));
                
        save([drec(ii).name '/net_' files(jj).name(7:end)],'zsum','zmad',...
            'z1sum','z1mad','z2sum','z2mad','znew','z1new','z2new','stan','stans','dstep');
        
        clear z* dstep yy ll
    end
    
    clear jj files
end

%% Sum the detection statistics for all stations - 2nd option
% Use the max in a given time window instead of moveout
% STEP: 11

clear
drec = dir('*07_2011'); % All directories, one per day
wins = 1; % Time window in sec. /2
dt = 0.008;
stan = 'G01'; % Reference station
stans = ['G03';'G05';'G06';'K04';'G07'];

for ii = 1:length(drec)
    files = dir([drec(ii).name '/obs01*']); % Get files per hour
    disp(['Working on ' drec(ii).name])
    
    for jj = 1:length(files)
        disp(['Working on ' files(jj).name])
        % Get detection statistics (stack: z, subspace: z1, empirical subspace: z2)
        load([drec(ii).name '/' files(jj).name]); % From reference station (G01)
        
        znew(1,:) = z/max(z); % Normalize and give double weights
        z1new(1,:) = z1/max(z);
        z2new(1,:) = z2/max(z); z2new(1,z2new(1,:)<0.25) = 0;
        z2new(1,:) = 2*(z2/max(z));
        clear z z1 z2
                
        for ll = 1:size(stans,1) % For the other stations
            % Order stations: 'OBS-3','OBS-5','OBS-6','KOERI-4','OBS-7'
            filesta = dir([drec(ii).name '/*sta' num2str(ll) '*']); % For the other stations
            load([drec(ii).name '/' filesta(jj).name]);
            % Cut all time series to be from beginning to largest length of
            % U1 for concatenation (see code on detection above)
            if length(z) ~= length(znew)
                znew = znew(:,1:min([length(z) length(znew)]));
                z1new = z1new(:,1:min([length(z1) length(z1new)]));
                z2new = z2new(:,1:min([length(z2) length(z2new)]));
                z = z(1,1:min([length(z) length(znew)]));
                z1 = z(1,1:min([length(z1) length(z1new)]));
                z2 = z(1,1:min([length(z2) length(z2new)]));
            end
            znew(ll+1,:) = z/max(z); % Normalized
            z1new(ll+1,:) = z1/max(z1);
            z2new(ll+1,:) = z2/max(z2); z2new(ll+1,z2new(ll+1,:)<0.25) = 0; % z2new(ll+1,:) = z2new(ll+1,:).^2;
            clear z z1 z2 filesta
        end
        
        % Do the sliding sum on the detection statistics
        for ll = 1:size(znew,2)
            is = [1 2 4 5]; % Station indexes to sum
            witmp = round(wins/dt/dstep);
            bs = ll-witmp; be = ll+witmp;
            if bs<1; bs=1; end % Beginning
            if be>size(znew,2); be=size(znew,2); end % End
            
            zsum(ll) = sum(max(znew(is,bs:be),[],2)); % Sum stations
            z1sum(ll) = sum(max(z1new(is,bs:be),[],2));
            z2sum(ll) = sum(max(z2new(is,bs:be),[],2));
            
            clear witmp bs be is
        end
        
        for ll = 1:size(znew,2)
            witmp = round(wins/dt/dstep);
            bs = ll-witmp; be = ll+witmp;
            if bs<1; bs=1; end % Beginning
            if be>size(znew,2); be=size(znew,2); end % End
            zsc(ll) = (median(z2sum(1,bs:be)) - mean(z2sum))/std(z2sum); % Type z-score
            % prs = prctile(z2new(:,ll),[25 75]); zsc(ll) = prs(2) - prs(1); % IQR metric
            
            clear witmp prs bs be
        end
        % zmad = mad(zsum,1); % Median absolute deviation of zsum
        % z1mad = mad(z1sum,1); % Not in use for now
        % z2mad = mad(z2sum,1);
                
        save([drec(ii).name '/net2_' files(jj).name(7:end)],'zsum','zsc',...
            'z1sum','z2sum','znew','z1new','z2new','stan','stans','dstep','dt');
        
        clear z* dstep yy ll
    end
    
    clear jj files
end
disp('Done !')
