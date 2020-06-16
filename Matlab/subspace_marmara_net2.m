% Modified by Camilo De La Hoz Lozano - JB Tary 2019
% Second step of subspace detection for the Marmara 2011 dataset
% Modified from Greg Beroza files
% This file is for applying a threshold to time series after subspace
% detection to find detections
clear; close all

% Global parameters
stan = 'G01'; % For a start, pick a station, in this case it is going to be the main station
w_l = 15; % number of seconds to search for catalog match
% Length of the signals to save in samples before and after detection 
% (check length of singular vectors of all stations and potential moveout to determine them)
l1 = 250; % 2 sec. before
l2 = 1875; % 15 sec. after
% Define a threshold for each method
st_threshold = 0.5; % Stack Threshold
sb_threshold = 0.5; % Subspace Threshold
em_threshold = 0.5; % Empirical Subspace Threshold
% Check that the template P isn't matching with detector S (check fail saves to adjust it)
stdS = 0.00001;

% Load in the existing catalog data
disp('Getting catalog information')
path = 'ParentEvents';
parents = dir([path '/*.mat']);

DATE=NaN(length(parents),1); 
LON=DATE; 
LAT=DATE; 
DEP=DATE; 
MAG=DATE; 
TT=DATE;
for i=1:length(parents)
    load([path '/' parents(i).name],'newhdr','comp','ppicks')
    [~,iv] = intersect(comp,stan);
    if isempty(iv)==1
        DATE(i,1)=-999; LON(i,1)=-999; LAT(i,1)=-999; DEP(i,1)=-999;
        MAG(i,1)=-999; TT(i,1)=-999;
        continue;
    end
    DATE(i,1)=datenum(newhdr(9),newhdr(10),newhdr(11),...
            newhdr(12),newhdr(13),newhdr(14));
    LON(i,1)=newhdr(3);
    LAT(i,1)=newhdr(2);
    DEP(i,1)=newhdr(4);
    MAG(i,1)=newhdr(15);
    TT(i,1)=ppicks(iv);
    clear newhdr comp ppicks
end

% Change to median of parent events used for the stack
TT = median(TT); % For now take median of travel time (search over several values interesting?)

for hh = 26 % Days used, example use day 26
cd([num2str(hh) '_07_2011'])

count_stack=1;
count_empsub=1;
count_subspace=1;

% Load the hours of the stations, for latter use in the ccrepick code
load(['marmara_svd_knw_' stan '.mat'],...
    'W','hours','U1','U2','dstep','stack','dt'); % Missing variables (JB)
signal_length = length(stack); % Added, JB

% Create time (in samples) vector of correlation indexes (start indexes)
N = hours{1}.signal.d; q=1:floor(((length(N)-length(U1))/dstep)+1);
xb = ((q-1)*dstep)+1; clear N q dstep

% This part of the code is for loading the data that will be stored in B 
% for latter use in subspace_ccrepick code
svdfiles = dir('*svd*knw*'); % Get all files with svds
for ii = 1:length(svdfiles)
    T = load([svdfiles(ii).name],'hours','hours_x','hours_y','stan');
%     T = load([svdfiles(ii).name],'hours','stan');
    % Load the data of all of the stations
    hours_list{1,ii} = T.hours; % The 3 components
    hours_list{2,ii} = T.hours_x;
    hours_list{3,ii} = T.hours_y;
    stan_list{ii} = svdfiles(ii).name(end-6:end-4);
    clear T
end
comp_order = 'ZXY';
clear ii svdfiles

%%
disp('Building catalog')
for p = 1:length(hours) % Run through each hour of data
    tic
    disp(['Working on ' hours{p}.name])
    % Fix the event time for matches with templates
    % Set the current hour

    currenthr = datenum([hours{p}.signal.nz(1),1,1,hours{p}.signal.nz(3),...
        hours{p}.signal.nz(4),hours{p}.signal.nz(5)+1]); % Add 1 sec.: totally artificial linked to file name
    currenthr = currenthr+hours{p}.signal.nz(2)-1;
 
    % Match the list to current hour (find events within this hour)
    date_active = DATE-currenthr;
    index=find((date_active > 0) & (date_active < 1/24));
    datelist = date_active(index);
    if isempty(index)==1; datelist=NaN; end; clear index
    
    datelist = datelist*(1/dt)*86400; % Convert to seconds/samples
    datelist = datelist + (TT*(1/dt)); % Add the theoretical travel time (s)
    
    % Load the stored information from the scan for this hour for main
    % station
    load(strcat(hours{p}.name,'_knwn.mat'));
    % Normalize the z function, avoid problems latter
    z = z/max(z);
    z1 = z1/max(z1);
    z2 = z2/max(z2);
    
    % Make a new directory for everything this hour to go into
    mkdir(strcat(hours{p}.name))
    cd(strcat(hours{p}.name))
    
    set(gcf,'visible','off')
    subplot(3,1,1); plot(xb,z,'k'); xlabel('Time (samples)'); ylabel('z')
    hold on; plot([min(xb) max(xb)],[em_threshold em_threshold],'--r','linewidth',1.5)
    subplot(3,1,2); plot(xb,z1,'k'); xlabel('Time (samples)'); ylabel('z1')
    subplot(3,1,3); plot(xb,z2,'k'); xlabel('Time (samples)'); ylabel('z2')
    saveas(gcf,[hours{p}.name(1:end-4),'Zs','.png'],'png');
    close all
    
    fid=fopen('catalog.txt','a'); % open the stack catalog file for write
    
    %% STACK
    count_stack = subdetec('stack',stack,U1,U2,stan_list,stan,...
        hours_list,z,xb,st_threshold,stdS,DATE,MAG,TT,count_stack,datelist,...
        currenthr,dt,w_l,p,signal_length,fid,l1,l2,comp_order);
    
    %% Empirical Subspace
    count_empsub = subdetec('empirsubspace',stack,U1,U2,stan_list,stan,...
        hours_list,z2,xb,em_threshold,stdS,DATE,MAG,TT,count_empsub,datelist,...
        currenthr,dt,w_l,p,signal_length,fid,l1,l2,comp_order);
    
    %% Subspace
    count_subspace = subdetec('subspace',stack,U1,U2,stan_list,stan,...
        hours_list,z1,xb,sb_threshold,stdS,DATE,MAG,TT,count_subspace,datelist,...
        currenthr,dt,w_l,p,signal_length,fid,l1,l2,comp_order);
    
    %% Remaining part
    clear date_active datelist currenthr Q z z1 z2
    
    [~,~,ib] = intersect(stan,stan_list);
    
    % Figures to check the detections
    cd stack
    files=dir('*.mat');
    for ii=1:length(files)
        load(files(ii).name,'Bamp')
        temp = zeros(1,l1+l2); % To avoid problems with signals of different sizes
        temp(1:length(Bamp{ib})) = Bamp{ib}(1,:)/max(abs(Bamp{ib}(1,:)));
        arr1(ii,:) = temp;
        clear Bamp temp
    end
    cd ..
    
    cd empirsubspace
    files=dir('*.mat');
    for ii=1:length(files)
        load(files(ii).name,'Bamp')
        temp = zeros(1,l1+l2); 
        temp(1:length(Bamp{ib})) = Bamp{ib}(1,:)/max(abs(Bamp{ib}(1,:)));
        arr2(ii,:) = temp;
        clear Bamp temp
    end
    cd ..
    
    cd subspace
    files=dir('*.mat');
    for ii=1:length(files)
        load(files(ii).name,'Bamp')
        temp = zeros(1,l1+l2); 
        temp(1:length(Bamp{ib})) = Bamp{ib}(1,:)/max(abs(Bamp{ib}(1,:)));
        arr3(ii,:) = temp;
        clear Bamp temp
    end
    cd ..
    
    set(gcf, 'visible', 'off')
    subplot(1,3,1)
    if exist('arr1','var')==1
        imagesc((0:(l1+l2)-1)*dt,1:size(arr1,1),arr1); colormap(gray)
        xlabel('Time (s)'); ylabel('Events')
        title('Events using the stack')
    end
    
    subplot(1,3,2)
    if exist('arr2','var')==1
        imagesc((0:(l1+l2)-1)*dt,1:size(arr2,1),arr2); colormap(gray)
        xlabel('Time (s)'); ylabel('Events')
        title('Events using the empirical subspace')
    end
    
    subplot(1,3,3)
    if exist('arr3','var')==1
        imagesc((0:(l1+l2)-1)*dt,1:size(arr3,1),arr3); colormap(gray)
        xlabel('Time (s)'); ylabel('Events')
        title('Events using the subspace')
    end
    
    if exist('arr1','var')==1 || exist('arr2','var')==1 || exist('arr3','var')==1
        saveas(gcf,[hours{p}.name(1:end-4),'_Detections','.png'],'png');
    end
    clear arr* i ib
    close all
    cd ..
    
    toc
end

fclose all;
clear W hours U1 U2 dstep stack dt signal_length xb svdfiles hours_list stan_list
clear comp_order currenthr
cd ..
end

disp('Done')
