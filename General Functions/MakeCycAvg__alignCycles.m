function [type,starts,ends,stims] = MakeCycAvg__alignCycles(info,Fs,ts,stim)
    if contains(info.dataType,'Translation')||contains(info.goggle_ver,'Moogles')
        type = 1;
        starts = 1;
        ends = [];
        toggle = stim{2};
%         %%
%         for i = 2:length(toggle)
%             if toggle(i) ~= toggle(i-1) % new cycle
%                 if i < length(toggle)
%                     ends = [ends i-1];
%                     starts = [starts i];
%                 else
%                     ends = [ends i];
%                 end 
%             end
%         end
%         snip_len = median(ends-starts);
%         ends = starts + snip_len - 1;
%         stims = stim{1};
        
        %%
        trig = abs(diff(toggle));
        starts = find(trig==1);
        snip_len = round(median(diff(starts)),0);
        ends = starts + snip_len - 1;
        starts(ends>length(ts)) = [];
        ends(ends>length(ts)) = [];
        stims = stim{1};
        stims = stims(starts(1):ends(1));
    elseif contains(info.dataType,{'RotaryChair','aHIT'})||contains(info.goggle_ver,'Moogles') %Align based on real/virtual motion traces
        if contains(info.dataType,'Sine')
            type = 1;
            fparts = split(info.dataType,'-');
            freqs = fparts(contains(fparts,'Hz'));
            freq = zeros(1,length(freqs));
            for i = 1:length(freqs)
                tempfreq = strrep(freqs(i),'Hz','');
                tempfreq = strrep(tempfreq,'p','.');
                freq(i) = str2double(tempfreq);
            end
            amp = str2double(strrep(fparts{contains(fparts,'dps')},'dps','')); % but what about translations?
%             if contains(info.dataType,'Rotation')
%                 amp = str2double(strrep(fparts{contains(fparts,'dps')},'dps',''));
%             elseif contains(info.dataType,'Translation')
%                 tempamp = strrep(fparts{contains(fparts,'mpsq')},'mpsw','');
%                 tempamp = strrep(tempamp,'p','.');
%                 amp = str2double(tempamp);
%             end
            snip_len = floor(Fs/min(freq));
            template = zeros(length(freq),snip_len);  
            for i = 1:length(freq)
                template(i,:) = sin(2*pi*freq(i)*ts(1:snip_len));
            end
            template = sum(template,1);
            template = amp/max(template)*template;
            if size(stim,2)==1
                template = template';
            end
            %Find the mismatch between signal and template
            errors = NaN(1,1000);
            sub_i = floor(linspace(1,length(stim)-snip_len,1000));
            for i = 1:1000  
                errors(i) = sum((template-stim(sub_i(i):sub_i(i)+snip_len-1)).^2);
            end  
            errors = (errors - min(errors))/(max(errors)-min(errors));
            errors2 = errors;
            errors2(errors>0.1)=NaN;
            poss_val = find(~isnan(errors2));
            region_s = [poss_val(1),poss_val([false,diff(poss_val)>1])];
            region_e = [poss_val([diff(poss_val)>1,false]),poss_val(end)];
            starts = region_s;
            for i = 1:length(starts)
                temp_err = NaN*errors;
                temp_err(region_s(i):region_e(i)) = errors(region_s(i):region_e(i));
                [~,m_ind] = min(temp_err);
                starts(i) = sub_i(m_ind);
            end        
            %Now find zeros crossings on the stimulus trace itself
            neg_stim = find(stim < 0);
            pos_stim = find(stim > 0);
            for i = 1:length(starts)
                if stim(starts(i)) < 0 %look for next positive value
                    l_i = pos_stim(find(pos_stim>starts(i),1,'first'));
                    if isempty(l_i)
                        l_i = length(starts);
                    end 
                    poss_i = [l_i-1 l_i];
                else %look for previous negative value
                    l_i = neg_stim(find(neg_stim<starts(i),1,'last'));
                    if isempty(l_i)
                        l_i = 1;
                    end 
                    poss_i = [l_i l_i+1];
                end
                [~,ind] = min(abs(stim(poss_i)));
                starts(i) = poss_i(ind); 
            end
            starts = unique(starts); %remove duplicates
            if length(starts) > 1
                snip_len1 = min(diff(starts));
                snip_len2 = length(stim) - starts(end);
                if abs(snip_len2-snip_len)/snip_len < 0.01 %Less than 1% off of the expected cycle length
                    snip_len = min([snip_len1 snip_len2]);
                else
                    snip_len = snip_len1; 
                end
            end 
            ends = starts + snip_len - 1;
            %Delete incomplete cycles
            starts(ends>length(stim)) = [];
            ends(ends>length(stim)) = [];
            all_stim = zeros(snip_len,length(starts));
            for i = 1:length(starts)
                all_stim(:,i) = stim(starts(i):ends(i));
            end
            %Remove any obvious erroneous motion traces
            tol = 0.2; %Amplitude can be 20% wrong and still be tolerates
            rm_tr = abs(max(all_stim)-amp)/amp > tol | abs(min(all_stim)+amp)/amp > tol;
            starts(rm_tr) = [];
            ends(rm_tr) = [];
            all_stim(:,rm_tr) = [];
            stims = mean(all_stim,2); 
        elseif contains(info.dataType,'Step')
            type = 2;
            stims = stim;
            starts = 1;
            ends = length(ts);
        elseif contains(info.dType,'Impulse')
            %ADD CODE HERE
            type = 1;
        else
            error('Unknown Data Type (RotaryChair/aHIT)')
        end 
        
    elseif contains(info.dataType,'eeVOR') %align using the trigger signal
        
        if contains(info.dataType,{'65Vector','MultiVector'})
            type = 1;
            %The trigger is actually showing when the trapezoids start and end. There
            %are only 20 cycles of the stimulus applied and there are 40 trigger
            %toggles. To make the stimulus trace, I assumed the trigger was high when
            %the virtual velocity was linear and that the max virtual velocity is 50
            %dps (estimated from PJB's 2019 manuscript).
            %The above was confirmed by the the VOG code.
            %I also want both excitation and inhibition phases of the stimulus during
            %alignemnt.
            %Find time window for alignment
            trig = diff(stim);
            starts = find(trig==1)-1;
            starts = starts(1:2:end);
            snip_len = min(diff(starts));
            ends = starts + snip_len;
            %Create model stimulus trace
            stims = stim(starts(1):ends(1));
            ind = find(stims==1);
            end1 = ind(diff(ind)>1);
            start2 = ind(find(diff(ind)>1)+1);
            stims(ind(1):end1) = linspace(0,50,length(ind(1):end1));
            stims(end1+1:start2-1) = 50*ones(length(end1+1:start2-1),1);
            stims(start2:ind(end)) = linspace(50,0,length(start2:ind(end)));
        elseif contains(info.dataType,'Sine') %sine (toggle = new cycle)
            type=1;
            trig = abs(diff(stim));
            starts = find(trig==1);
            snip_len = round(median(diff(starts)),0);
            ends = starts + snip_len - 1;
            fparts = split(info.dataType,'-');
            tempfreq = strrep(fparts{contains(fparts,'Hz')},'Hz','');
            tempfreq = strrep(tempfreq,'p','.');
            freq = str2double(tempfreq);
            amp = str2double(strrep(fparts{contains(fparts,'dps')},'dps',''));
            stims = amp*sin(2*pi*freq*(ts(1:snip_len)));
        elseif contains(info.dataType,{'PulseTrain','Autoscan'}) % (high = on, low = off)
            type=1;
            trig = (diff(stim));
            starts = find(trig==1);
            snip_len = min(diff(starts));
            ends = starts + snip_len - 1;
            stims = 50*stim(starts(1):ends(1));
        elseif contains(info.dataType,'Activation') % (low = dark, high = light but don't cycle average here)
            type=2;
            stims = stim;
            starts = 1;
            ends = length(ts);
        else
            error('Unknown Data Type (eeVOR)')
        end
        if ends(end) > length(ts)
           starts(end) = [];
           ends(end) = [];
        end
    else
        error('Unknown Data Type')
    end 
end