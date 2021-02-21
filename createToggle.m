function [toggle] = createToggle(trigger,Fs,freq,method)
% % creates toggle from sinusoid trigger


switch method
    case 1
        % Method 1, which just chooses crossing points of 0 (close to it) but produces cycles of different lengths
        toggle = 0;
        tog = 0;
        for i = 2: length(trigger)
            if (trigger(i)*trigger(i-1) < 0) && (trigger(i) > 0)
                if tog == 0
                    tog = 1;
                elseif tog == 1
                    tog = 0;
                end
            else
                toggle(i) = tog;
            end
        end
        
        
    case 2
        % Method 2, which uses sine frequency and sampling frequency to
        % create equally sized cycles (except for maybe an incomplete cycle
        % at the end)
        T = 1/freq;
        i = 1;
        timepts = 1;
        while length(trigger) >= i*Fs*T + 1 % while we are still not at the end of time vector
            timepts = [timepts, i*Fs*T+1];
            i = i+1;
        end
        tog = 0;
        toggle = zeros(1,length(time));
        for j = 1:length(timepts)-1
            toggle(timepts(j):timepts(j+1)) = tog;
            if tog == 0
                tog = 1;
            elseif tog == 1
                tog = 0;
            end
        end
        
        % toggle(timepts(end):length(toggle)) = tog;
        toggle(timepts(end):length(toggle)-1) = tog;
        if tog == 0
            tog = 1;
        elseif tog == 1
            tog = 0;
        end
        toggle(length(toggle)) = tog;
        
    case 3
        % Method 3, a combination of the two methods above, finds first
        % value going from negative to positive and then starts the toggle
        % from there based on T and Fs

        i = 2;
        while (trigger(i) < 0) || (trigger(i)*trigger(i-1) > 0)
            i = i + 1;
            if i + 1 > length(trigger)
                return
            end
        end
        % when exiting while loop, i will be the first positive coefficient
        % of the trigger pulse after crossing 0, so here we start toggling
        T = 1/freq;
        n = 1;
        timepts = i;
        while length(trigger) >= n*Fs*T + i % while we are still not at the end of time vector
            timepts = [timepts, n*Fs*T+i];
            n = n+1;
        end
        tog = 1;
        toggle = zeros(1,length(trigger));
        for j = 1:length(timepts)-1
            toggle(timepts(j):timepts(j+1)) = tog;
            if tog == 0
                tog = 1;
            elseif tog == 1
                tog = 0;
            end
        end
        
        toggle(timepts(end):length(toggle)) = tog;

end