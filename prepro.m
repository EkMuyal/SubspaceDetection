% Pre-process parent event data, multiplexed or not
% Then cut the data for further processing
% opt: 4 choices 'X' for x comp, 'Y' for y comp, 'Z' for z comp, 'M' for
% 3-C multiplexed signal

function [events,dt,b,a] = prepro(events,stan,nPoles,fcl,fch,winPa,winPb,opt)

% Process and Filter waveforms
display('Processing and Filtering the Waveforms')
dt = events{1}.dt;
nyq = 1/(2*dt);
Wn =[fcl fch]/nyq; % Wn must be between 0 and 1, where 1 = nyquist frequency
[b,a] = butter(nPoles,Wn,'bandpass');
winP = [winPa winPb]; % Short window at beginning of parent events containing the P wave

for i = 1:length(events)
    
    % Find station in file (or not)
    stalist = events{i}.seis.comp;
    [~,iv] = intersect(stalist,stan);
    
    if isempty(iv) ~= 1
        
        if strcmp(opt,'X') == 1 || strcmp(opt,'M') == 1
            x = events{i}.seis.dataxs{iv}; % Take horizontal x component
            x = x - mean(x); % Remove the mean
            x = filtfilt(b,a,x); % Filter
            xlong.d = x'; xshort.d = x(1,winP(1):winP(2))'; % Make cuts
            xlong = staper(xlong,0.05); % Taper_x
            xshort = staper(xshort,0.05);
            % RSSQ - FIX to Unit Amplitude
            rssq1x = max(abs(xlong.d));
            xlong.d = xlong.d / rssq1x;
            rssq2x = max(abs(xshort.d));
            xshort.d = xshort.d / rssq2x;
        end
        
        if strcmp(opt,'Y') == 1 || strcmp(opt,'M') == 1
            y = events{i}.seis.datays{iv}; % Take horizontal y component
            y = y - mean(y);
            y = filtfilt(b,a,y);
            ylong.d = y'; yshort.d = y(1,winP(1):winP(2))';
            ylong = staper(ylong,0.05); % Taper_y
            yshort = staper(yshort,0.05);
            
            rssq1y = max(abs(ylong.d));
            ylong.d = ylong.d / rssq1y;
            rssq2y = max(abs(yshort.d));
            yshort.d = yshort.d / rssq2y;
        end
        
        if strcmp(opt,'Z') == 1 || strcmp(opt,'M') == 1
            z = events{i}.seis.datazs{iv}; % Take vertical z component
            z = z - mean(z);
            z = filtfilt(b,a,z);
            zlong.d = z'; zshort.d = z(1,winP(1):winP(2))';
            zlong = staper(zlong,0.05); % Taper_z
            zshort = staper(zshort,0.05);
            
            rssq1z = max(abs(zlong.d));
            zlong.d = zlong.d / rssq1z;
            rssq2z = max(abs(zshort.d));
            zshort.d = zshort.d / rssq2z;
        end
        
        % Save processed copies and multiplexed ones if needed
        % This step is for one station multiplexing
        if strcmp(opt,'M') == 1
            ml(1,:) = xlong.d;
            ml(2,:) = ylong.d;
            ml(3,:) = zlong.d;
            [xyzl]=mx(ml);
            xyzlong.d = xyzl';
            
            ms(1,:) = xshort.d;
            ms(2,:) = yshort.d;
            ms(3,:) = zshort.d;
            [xyzs]=mx(ms);
            xyzshort.d = xyzs';
            
            % Test on moving-avg for multiplexed: end-up with signal very
            % similar to largest amplitude component (Z here)
            % xyzlong.d = movsum(xyzlong.d,3)/3;
            
            events{i}.signal = xyzlong;
            events{i}.clip = xyzshort;
        else
            if strcmp(opt,'X') == 1; siglong = xlong; sigshort = xshort; end
            if strcmp(opt,'Y') == 1; siglong = ylong; sigshort = yshort; end
            if strcmp(opt,'Z') == 1; siglong = zlong; sigshort = zshort; end
            
            events{i}.signal = siglong;
            events{i}.clip = sigshort;
        end
        
        events{i}.signal = srmean(events{i}.signal);
        events{i}.clip = srmean(events{i}.clip);
        
        clear x xlong xshort rssq1x rssq2x
        clear y ylong yshort rssq1y rssq2y
        clear z zlong zshort rssq1z rssq2z
        clear xyzlong xyzshort xyzs xyzl ml ms siglong sigshort
        
    else
        events{i}.hasSeis = 0;
    end
    
    clear stalist iv
end

% Re-cut long signals for cross-correlation
% JB KO4
if strcmp(stan(1),'K') == 1
    kk = 1;
    for i=1:length(events)
        if events{i}.hasSeis == 0; continue;
        else
            siz(kk) = length(events{i}.signal.d); kk = kk + 1;
        end
    end
    sizm = min(siz);
    
    for i=1:length(events)
        if events{i}.hasSeis == 0; continue;
        else
            events{i}.signal.d = events{i}.signal.d(winP(1):sizm);
        end
    end
    clear siz sizm i kk
else
    
    % for G stations
    kk = 1;
    for i=1:length(events)
        if events{i}.hasSeis == 0; continue; end
        siz(kk) = length(events{i}.signal.d); kk = kk+1;
    end
    sizm = min(siz); clear siz kk
    
    jj = 1;
    for i=1:length(events)
        if events{i}.hasSeis == 0 ; continue; end
        events{jj}.signal.d = events{i}.signal.d(winP(1):sizm); jj = jj+1;
    end
    clear siz sizm i kk jj
    
end
