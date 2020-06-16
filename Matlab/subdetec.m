% Function to check for detections for subspace detector (stack, subspace,
% empirical subspace)
% Some inputs:
%   signame: type of detection (ex: 'stack', 'empirsubspace' or 'subspace')
%   stack, U1, U2: stack and svd time series
%   stan_list: list of stations used for detection
%   stan: main station used for detection (ex: 'G01')
%   hours_list: hours of continuous data
%   dd: time series used for detection
%   xb: time vector corresponding to 'dd' time series
%   thrs: threshold for detection (ex: 0.2)
%   count: count for current detection time series
%   p: index of current hour to analyze in hours_list
%   l1: length of signal to save before detection
%   l2: length of signal to save after detection
%   comp_order: order of the components for further processing

function count = subdetec(signame,stack,U1,U2,stan_list,stan,...
    hours_list,dd,xb,thrs,stdS,DATE,MAG,TT,count,datelist,...
    currenthr,dt,w_l,p,signal_length,fid,l1,l2,comp_order)

% Derivative
Z=diff(dd);

% Open the write file and CD to the 'signame' directory
mkdir(signame) % Make a new directory for 'signame' info
cd(signame)

% Do the stack run
disp(['Building the ' signame ' catalog'])
for q=1+10:length(dd)-10 % Added +10 and -10
    
    % Check for a peak in the CC function
    % Check the value is above a threshold
    % Check the first derivative for a sign change
    if dd(q) > thrs && sign(Z(q-1)) == 1 && sign(Z(q)) == -1 && sign(Z(q+1)) == -1
        
        % Check this is the maximum in the neighborhood of 0.25 seconds
        % Eliminate other matches within this 10 samples range
        vect = [dd(q-10:q-1) dd(q+1:q+10)];
        if dd(q) > vect
            
            % Process the template match to save the waveform
            smp=xb(q); % Sample number in time series
            
            % Put the name of the stations in a cell array
            for qq =1:length(stan_list)
                B{1,qq} = stan_list{qq};
            end
            clear qq
            
            % Get the index of the main station
            [~,~,ib] = intersect(stan,B);
            
            for qq =1:length(stan_list)
                h1{1,qq} = hours_list{1,qq};
                h2{1,qq} = hours_list{2,qq};
                h3{1,qq} = hours_list{3,qq};
            end
            
            % Load the data and cut the d vector
            for qq =1:length(stan_list)
                if smp+l2-1 > length(h1{qq}{p}.signal.d)
                    sig(1,:) = h1{qq}{p}.signal.d(smp-l1:end); % Cut
                    sig(2,:) = h2{qq}{p}.signal.d(smp-l1:end); % Cut
                    sig(3,:) = h3{qq}{p}.signal.d(smp-l1:end); % Cut
                    B{2,qq} = sig; clear sig
                else if smp-l1 < 1
                        sig(1,:) = h1{qq}{p}.signal.d(1:smp+l2-1); % Cut
                        sig(2,:) = h2{qq}{p}.signal.d(1:smp+l2-1); % Cut
                        sig(3,:) = h3{qq}{p}.signal.d(1:smp+l2-1); % Cut
                        B{2,qq} = sig; clear sig
                    else
                        sig(1,:) = h1{qq}{p}.signal.d(smp-l1:smp+l2-1); % Cut
                        sig(2,:) = h2{qq}{p}.signal.d(smp-l1:smp+l2-1); % Cut
                        sig(3,:) = h3{qq}{p}.signal.d(smp-l1:smp+l2-1); % Cut
                        B{2,qq} = sig; clear sig
                    end
                end
            end
            clear qq h1 h2 h3
            
            starttimeseries = smp-l1;
            
            % Processing of the data
            for qq = 1:length(stan_list)
                Bamp{qq} = B{2,qq}; % save a waveform with true amplitude
                sig = Bamp{qq}; maxamp = max(max(abs(sig)));
                sig = sig/maxamp;
                B{2,qq} = sig; % change to unit amplitude
                clear sig
            end
            clear qq
            
            mainb = B(2,ib); mainbamp = Bamp{ib}; clear ib
            
            B1 = mainb{1}./norm(mainb{1}); % norm - save as B1
            
            %                 B.G01 = detrend(B.G01,'constant'); % detrend
            %                 B.G01 = filter(b,a,B.G01); % filter
            %                 Bamp.G01 = B.G01; % save a waveform with true amplitude
            %                 B1 = B.G01./norm(B.G01); % norm - save as B1
            %                 B.G01 = B.G01./max(abs(B.G01)); % change to unit amplitude
            
            % Graphics
            set(gcf, 'visible', 'off')
            stemp = mainb{1};
            if l1+signal_length <= length(stemp)
                stemp = stemp(:,l1:l1+signal_length-1);
                plot((0:signal_length-1)*dt,stemp); hold on;
                plot((0:signal_length-1)*dt,stack/max(abs(stack))+1,'r')
                plot((0:signal_length-1)*dt,U1/max(abs(U1))+2,'b')
                plot((0:signal_length-1)*dt,U2(:,2)/max(abs(U2(:,2)))+3,'g')
                title([q dd(q) xb(q)])
                set(gca,'YTick',0:3); ylim([-1 4])
                set(gca,'YTickLabel',{'Event' 'stack' 'U1' 'U2'});
                xlabel('seconds'); clear stemp
            end
            
            % Check that the template P isn't matching with detector S
            % Note: the threshold here will likely need to be
            % changed for specific examples (check in fail saves)
            
            if std(mainb{1}(1:300)) > stdS
                count = count + 1;
                S_amp=max(abs(mainbamp(300:end)));
                
                iv = find(isnan(datelist)==1);
                if isempty(iv)==1
                    tmp = abs(datelist-q);
                    [junk,idx] = min(tmp);
                    closest = datelist(idx);
                    closest_dt = abs(closest-q);
                else % If NaN then pass to the next loop
                    closest_dt = (1/dt)*w_l + 10;
                end
                clear iv
                
                if closest_dt < (1/dt)*w_l
                    % Then there is a match to a catalog (so not good)
                    tmp=abs(DATE-(currenthr+((closest/(1/dt)-TT)/86400)));
                    [junk,indx] = min(tmp);
                    M = MAG(indx);
                    % Print the event information to a text file
                    fprintf(fid,'%s,%7d,%s,%s,%4.2f,%f,%s\n',num2str(count-1),q,datestr(currenthr+(q/(1/dt)-TT)/86400,31),...
                        datestr(currenthr+((smp/(1/dt))/86400),'yyyy-mm-dd HH:MM:SS.FFF'),M,S_amp,'1,1');
                else
                    % Then it is a new event (so good)
                    fprintf(fid,'%s,%7d,%s,%s,%s,%f,%s\n',num2str(count-1),q,datestr(currenthr+(q/(1/dt)-TT)/86400,31),...
                        datestr(currenthr+((smp/(1/dt))/86400),'yyyy-mm-dd HH:MM:SS.FFF'),'NaN',S_amp,'2,1');
                end
                clear tmp junk idx
                % Save the image use png for quick views - eps for figure making
                saveas(gcf,[num2str(count-1),'_',num2str(q),'.png'],'png');
                save([num2str(count-1),'_',num2str(q)],'Bamp','B','currenthr','smp','dt','starttimeseries','comp_order')
                %                     print('-depsc2','-r300',[num2str(count-1),'_',num2str(q),'.eps']);
                
                % Save the P/S mismatches to check
            else
                saveas(gcf,['fail_',num2str(q),'.png'],'png');
            end
        end
        clear B S_amp Bamp B1 h closest closest_dt smp mainb mainbamp
        clear starttimeseries
        % Clear gcf info - very important for memory
        clf
        
        close all
        close all hidden
        close force all
        
    end
end

cd ..

