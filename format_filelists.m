% Create complete list of sac files/dates for a given station
% 
% INPUT:
% dirgen - path to directory of the station's component directories
% sta: station name (ex: 'G01')
% comp - component letters/numbers corresponding to directory's name (ex:
% 'Z' for 'G01Z')
% datemin - start date of dataset for complete list
% datemax - end date of dataset for complete list
% 
% OUTPUT:
% infoc - complete list of sac files
% datec - corresponding dates
% 
% Example:
% dirgen = '/media/jbtary/OSDisk/Users/jb.tary/Data/Marmara/2011/';
% datemin = datenum([2011,04,15,12,0,0]); % 15th of April, 12pm
% datemax = datenum([2011,07,31,23,0,0]); % 31st of July, 23am
% sta = 'G01'; comp = 'Z';
% infoc = format_filelists(dirgen,datemin,datemax,sta,comp);

function [infoc,datec] = format_filelists(dirgen,datemin,datemax,sta,comp) 
    
% Complete list of sac files for each station
info = dir([dirgen sta '/' sta comp '/*.sac']);

day = arrayfun(@(x) str2num(x.name(7:8)),info,'UniformOutput',false);
month = arrayfun(@(x) str2num(x.name(10:11)),info,'UniformOutput',false);
yr = arrayfun(@(x) str2num(x.name(13:16)),info,'UniformOutput',false);
hr = arrayfun(@(x) str2num(x.name(18:19)),info,'UniformOutput',false);
mm = arrayfun(@(x) str2num(x.name(21:22)),info,'UniformOutput',false);

day = cell2mat(day); month = cell2mat(month); yr = cell2mat(yr);
hr = cell2mat(hr); mm = cell2mat(mm);
hrm = round(hr+(mm/60));
datetmp = datenum([yr,month,day,hrm,zeros(size(hrm,1),1),zeros(size(hrm,1),1)]);
[~,ia] = unique(datetmp);

info = info(ia,:);
datetmp2 = datetmp(ia,:);
% To check which files were removed: setdiff(datetmp,datetmp2)
clear day month yr hr mm hrm ia

% Exhaustive list of dates to compare with actual lists from OBSs

date(1,1) = datemin*24; datetmp = datemin*24; % To avoid roundoff errors
kk = 2;
while datetmp < datemax*24
    date(kk,1) = date(kk-1,1)+1;
    datetmp = date(kk,1);
    kk = kk+1;
end
date = date/24;

% Fit list of each OBS in complete file        
[~,ia,~] = intersect(datetmp2,date);
if isempty(ia) == 1
    disp('No files found in the date range provided')
    infoc = [];
    datec = [];
else
    name1 = extractfield(info,'name')';
    name1 = name1(ia,1);
    infoc = name1;
    datec = datetmp2(ia,1);
end

clear datetmp* info ia ib name1

