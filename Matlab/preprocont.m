% Pre-process continuous data, multiplexed or not
% opt: 4 choices 'X' for x comp, 'Y' for y comp, 'Z' for z comp, 'M' for
% 3-C multiplexed signal

function hours = preprocont(pathcont,stan,datemin,datemax,dt,b,a,opt)

if strcmp(stan(1),'G') == 1 % Ifremer stations
    switch opt
        case 'M'
            [infoc_x,~] = format_filelists([pathcont '/' ],datemin,datemax,...
                stan,'X');
            filist_x = infoc_x;
            [infoc_y,~] = format_filelists([pathcont '/' ],datemin,datemax,...
                stan,'Y');
            filist_y = infoc_y;
            [infoc_z,~] = format_filelists([pathcont '/' ],datemin,datemax,...
                stan,'Z');
            filist_z = infoc_z;
            infoc = [filist_x filist_y filist_z];
            clear filist*
        case 'X'
            [infoc,~] = format_filelists([pathcont '/' ],datemin,datemax,...
                stan,'X');
        case 'Y'
            [infoc,~] = format_filelists([pathcont '/' ],datemin,datemax,...
                stan,'Y');
        case 'Z'
            [infoc,~] = format_filelists([pathcont '/' ],datemin,datemax,...
                stan,'Z');
        otherwise
            disp('Unknown station or component')
    end
    
else % Koeri stations
    switch opt
        case 'M'
            [infoc_x,~] = format_filelists([pathcont '/' ],datemin,datemax,...
                stan,'E');
            filist_x = infoc_x;
            [infoc_y,~] = format_filelists([pathcont '/' ],datemin,datemax,...
                stan,'N');
            filist_y = infoc_y;
            [infoc_z,~] = format_filelists([pathcont '/' ],datemin,datemax,...
                stan,'Z');
            filist_z = infoc_z;
            infoc = [filist_x filist_y filist_z];
            clear filist*
        case 'X'
            [infoc,~] = format_filelists([pathcont '/' ],datemin,datemax,...
                stan,'E');
        case 'Y'
            [infoc,~] = format_filelists([pathcont '/' ],datemin,datemax,...
                stan,'N');
        case 'Z'
            [infoc,~] = format_filelists([pathcont '/' ],datemin,datemax,...
                stan,'Z');
        otherwise
            disp('Unknown station or component')
    end
    
end

disp('Reading in continuous data')

for ii = 1:size(infoc,1)
    if strcmp(stan(1),'G') == 1 % Ifremer stations
        switch opt
            case 'M'
                
                pathx = [pathcont '/' stan '/' stan 'X/' infoc{ii,1}];
                pathy = [pathcont '/' stan '/' stan 'Y/' infoc{ii,2}];
                pathz = [pathcont '/' stan '/' stan 'Z/' infoc{ii,3}];
                
                hours{ii}.seis_x = rsac(pathx);
                hours{ii}.seis_y = rsac(pathy);
                hours{ii}.seis_z = rsac(pathz);
                hours{ii}.hasSeis = 1;
                hours{ii}.name = infoc{ii,3};
                
                clear pathx pathy pathz
            otherwise
                path = [pathcont '/' stan '/' stan opt '/' infoc{ii,1}];
                hours{ii}.seis = rsac(path);
                hours{ii}.hasSeis = 1;
                hours{ii}.name = infoc{ii,1};
                
                clear path
        end
        
    else % Koeri stations
        switch opt
            case 'M'
                
                pathx = [pathcont '/' stan '/' stan 'E/' infoc{ii,1}];
                pathy = [pathcont '/' stan '/' stan 'N/' infoc{ii,2}];
                pathz = [pathcont '/' stan '/' stan 'Z/' infoc{ii,3}];
                
                hours{ii}.seis_x = rsac(pathx);
                hours{ii}.seis_y = rsac(pathy);
                hours{ii}.seis_z = rsac(pathz);
                hours{ii}.hasSeis = 1;
                hours{ii}.name = infoc{ii,3};
                
                clear pathx pathy pathz
            otherwise
                if strcmp(opt,'X'); path = [pathcont '/' stan '/' stan 'E/' infoc{ii,1}]; end
                if strcmp(opt,'Y'); path = [pathcont '/' stan '/' stan 'N/' infoc{ii,1}]; end
                if strcmp(opt,'Z'); path = [pathcont '/' stan '/' stan 'Z/' infoc{ii,1}]; end
                hours{ii}.seis = rsac(path);
                hours{ii}.hasSeis = 1;
                hours{ii}.name = infoc{ii,1};
                
                clear path
        end
    end
end

disp('Processing and Filtering the Waveforms')

for ii=1:length(hours)
    
    switch opt
        case 'M'
            x = hours{ii}.seis_x;
            y = hours{ii}.seis_y;
            z = hours{ii}.seis_z;
            
            if hours{ii}.seis_z.dt ~= dt % Resampling for stations with different sampling intervals
                [p,q] = rat(dt/hours{ii}.seis_z.dt);
                x.d= resample(x.d,q,p);
                y.d= resample(y.d,q,p);
                z.d= resample(z.d,q,p);
                clear p q
            end
            
            rssq = max(abs(x.d)); x.d = x.d./rssq; clear rssq % Normalize
            rssq = max(abs(y.d)); y.d = y.d./rssq; clear rssq
            rssq = max(abs(z.d)); z.d = z.d./rssq; clear rssq
            
            x.d = filtfilt(b,a,x.d); % Filter
            y.d = filtfilt(b,a,y.d);
            z.d = filtfilt(b,a,z.d);
            
            x = srmean(x); % Remove the mean
            y = srmean(y);
            z = srmean(z);
            
            % Do the multiplexing
            h(1,:) = x.d; h(2,:) = y.d; h(3,:) = z.d;
            xyzh=mx(h); clear h
            h = z;
            h.d = xyzh';
            h.t = 0:dt:((length(h.d)-1)*dt);
            h.t = h.t';
            
            hours{ii}.signal = h;
            clear x y z h
            
        otherwise
            x = hours{ii}.seis;
            
            if hours{ii}.seis.dt ~= dt % Resampling for stations with different sampling intervals
                [p,q] = rat(dt/hours{ii}.seis.dt);
                x.d = resample(x.d,q,p);
                x.dt = dt;
                x.t = 0:dt:((length(x.d)-1)*dt);
                x.t = x.t';
                clear p q
            end
            rssq = max(abs(x.d)); x.d = x.d./rssq; clear rssq % Normalize
            x.d = filtfilt(b,a,x.d); % Filter
            x = srmean(x); % Remove the mean
            hours{ii}.signal = x;
            clear x
    end

end

