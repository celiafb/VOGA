%% MakeCycAvg.m
% This function was created to automate the data analysis pipeline for the
% LDVOG data sets as much as possible.
% Cycle average before filtering

function done = MakeCycAvg(path,Seg_Path,Cyc_Path,Experimenter,version,exp_types)
clc;       
load('VNELcolors.mat','colors')
% Fill colors for cycle selection
%Color of cycle on long graph
colors.cyc_keep = [0.75 1 0.75]; 
colors.cyc_bold_k = [0 0.75 0];
colors.cyc_rm = [1 1 1];
colors.cyc_bold_r = [0.85 0.85 0.85];
%% Load in data
progress_tab = assessProgress(path);
progress_i = [find(~progress_tab{:,2}&~progress_tab{:,3});find(progress_tab{:,2}|progress_tab{:,3})]; %put unanalyzed files at the top
if ~isempty(exp_types)
    progress_i(ismember(progress_i,find(~contains(progress_tab{:,1},exp_types)))) = [];
end
set(0,'units','pixels')  
Pix_SS = get(0,'MonitorPositions');
pix_wid = Pix_SS(1,3);
pix_height = Pix_SS(1,4);
f = figure;
set(f,'Units','normalized','Position',[1 - 500/pix_wid, 0, 500/pix_wid, 1])
uit = uitable(f,'Data',table2cell(progress_tab(progress_i,:)),'ColumnName',progress_tab.Properties.VariableNames);
uit.Units = 'normalized';
uit.Position = [20/500, 20/pix_height, 1 - 40/500, 1 - 40/pix_height];
uit.ColumnWidth = [{320},{55},{50}];
[indx,tf] = nmlistdlg('PromptString','Select a file to analyze:',...
                       'SelectionMode','single',...
                       'ListSize',[350 300],...
                       'ListString',table2cell(progress_tab(progress_i,1)));  
close(f)
%If the user selects cancel
if ~tf
    disp('Operation Ended.')
    done = true;
    return; 
end
done = false;
a = dir([Seg_Path,filesep,'*.mat']);
a = {a.name}';
a = a(progress_i);
In_FileName = a{indx};
load([Seg_Path,filesep,In_FileName],'Data');
%% Initialize Figure
fig = figure(1);
delete(findall(gcf,'type','annotation')) %in case there are leftover anotations
fig.Units = 'inches';
fig.Position = [0 0 11 10];
%Title
annotation('textbox',[0 .9 1 .1],'String',strrep(strrep(strrep(In_FileName,'_',' '),'-',' '),'.mat',''),'FontSize',14,...
    'HorizontalAlignment','center','EdgeColor','none');
%% Extract raw position data
info = Data.info;
info.Analyzer = Experimenter;
info.ver = version;
info.colors = colors;
Fs = Data.Fs;
if contains(info.dataType,'Activation')
    %preserve the time because this will be used to rejoin them
    te = Data.Time_Eye;
    ts = Data.Time_Stim; 
else
    te = Data.Time_Eye - Data.Time_Eye(1);
    ts = Data.Time_Stim - Data.Time_Stim(1);
end
%Assign Trigger
if contains(info.goggle_ver,'Moogles') %MOOG room coil system
    stim = Data.Trigger;
else %NKI or LDVOG Trigger = internal gyro or comes from the PCU for eeVOR
    if contains(info.dataType,'RotaryChair')
        if isfield(Data,'HeadVel_Z')
            stim = Data.HeadVel_Z;
        else
            stim = Data.HeadMPUVel_Z; 
        end   
    elseif contains(info.dataType,'aHIT')
        if contains(info.dataType,'LHRH')
            stim = Data.HeadVel_Z;
        elseif contains(info.dataType,'LARP')
            stim = (Data.HeadVel_X - Data.HeadVel_Y)/sqrt(2);
        elseif contains(info.dataType,'RALP')
            stim = (Data.HeadVel_X + Data.HeadVel_Y)/sqrt(2);
        else
            stim = Data.Trigger;
        end    
    else
        stim = Data.Trigger; 
    end       
end
%Fix huge number of NaN values in torsion traces of NKI traces
if contains(info.goggle_ver,'NKI')
    if sum(isnan(Data.LE_Position_X)) > 0.9*length(te) %less than 10% data integrity
        Data.LE_Position_X = zeros(length(te),1); %set to 0 so no torsion
    else
        Data.LE_Position_X = spline(te(~isnan(Data.LE_Position_X)),Data.LE_Position_X(~isnan(Data.LE_Position_X)),te);
    end    
    if sum(isnan(Data.RE_Position_X)) > 0.9*length(te)
        Data.RE_Position_X = zeros(length(te),1);
    else
        Data.RE_Position_X = spline(te(~isnan(Data.RE_Position_X)),Data.RE_Position_X(~isnan(Data.RE_Position_X)),te);
    end    
end
%% Cycle Align
%Defines a variable type that determines what type of filtering and
%analysis needs to be done. type = 1: cycle average, 2: full trace (vel
%step/activation)
[type,starts,ends,stims] = MakeCycAvg__alignCycles(info,Fs,ts,stim);
%Remove any unnecessary trace the start and end
chk = 1;
%keep_inds = zeros(ends(1)-starts(1)+1,length(starts));
%keep_inds = zeros(end)
for i = 1:length(starts)
    keep_inds(:,i) = starts(i):ends(i);
end
keep_inds = keep_inds - starts(1)+1;
if ends(end)-length(te) > 0         % CFB edits since after making all the cycles the same length, ends was longer than te
    ts = [ts zeros(1,ends(end)-length(te))];
    Data.LE_Position_X = [Data.LE_Position_X ;zeros(1,ends(end)-length(te))'];
    Data.LE_Position_Y = [Data.LE_Position_Y ;zeros(1,ends(end)-length(te))'];
    Data.LE_Position_Z = [Data.LE_Position_Z ;zeros(1,ends(end)-length(te))'];
    Data.RE_Position_X = [Data.RE_Position_X ;zeros(1,ends(end)-length(te))'];
    Data.RE_Position_Y = [Data.RE_Position_Y ;zeros(1,ends(end)-length(te))'];
    Data.RE_Position_Z = [Data.RE_Position_Z ;zeros(1,ends(end)-length(te))'];
    stim{1} = [stim{1} zeros(1,ends(end)-length(te))];
    stim{2} = [stim{2} zeros(1,ends(end)-length(te))];
    te = [te zeros(1,ends(end)-length(te))];
end
if contains(info.dataType,'Activation')
    te = te(starts(1):ends(end));
    ts = ts(starts(1):ends(end));
    t_snip = reshape(ts(1:length(stims)),1,[]);
else
    te = te(starts(1):ends(end)) - te(starts(1));
    ts = ts(starts(1):ends(end)) - ts(starts(1));
    t_snip = reshape(ts(1:length(stims))-ts(1),1,[]);
end
stim = stim(starts(1):ends(end));
Data.LE_Position_X = Data.LE_Position_X(starts(1):ends(end));
Data.LE_Position_Y = Data.LE_Position_Y(starts(1):ends(end));
Data.LE_Position_Z = Data.LE_Position_Z(starts(1):ends(end));
Data.RE_Position_X = Data.RE_Position_X(starts(1):ends(end));
Data.RE_Position_Y = Data.RE_Position_Y(starts(1):ends(end));
Data.RE_Position_Z = Data.RE_Position_Z(starts(1):ends(end));
%% Set some Defaults
%Filter Traces with intial guesses
if type == 1
    if contains(info.goggle_ver,'NKI')
        pos_med = [21;21;3;3;3;3];
        pos_spline = [1;1;0.999995;0.999995;0.9999995;0.9999995];
    else
        pos_med = [11;11;3;3;3;3];
        pos_spline = [0.99995;0.99995;0.999995;0.999995;0.9999995;0.9999995];
    end
    pos_sgolay = [2*ones(6,1);5*ones(6,1)];
    if contains(info.dataType,'Sine')
        vel_smooth = round(t_snip(end)*40);
    else
        vel_smooth = 30;
    end
    vel_spline = NaN;
    vel_acc = NaN; 
    vel_med = NaN;
    keep_tr = true(1,length(starts));
else %velstep and activation
    pos_med = [11;11;5;5;3;3];
    pos_spline = NaN(6,1);
    pos_sgolay = [2*ones(6,1);5*ones(6,1)];
    vel_smooth = NaN;
    vel_spline = NaN;
    vel_acc = 10;
    vel_med = Fs+1-mod(Fs,2);
    keep_tr = true(1,length(starts));
end
line_wid.norm = 0.5;
line_wid.bold = 2;
%% Once analyzeable, here is the while loop they stay in until saving or exiting
%You can change the order/existance of these options without ruining
%anything because the comparrisons are all string based
opts = {'Start Over','Filter Position','Filter Velocity','Select Cycles','Set Y-axis limits','Save'};
ind = 1; %Run the start procedure first
while ~strcmp(opts{ind},'Save') %Run until it's ready to save or just hopeless
    if strcmp(opts{ind},'Start Over')
        filt_params_p = [pos_med;pos_spline;pos_sgolay];
        filt_params_v = [vel_smooth;vel_spline;vel_acc;vel_med];        
        YLim_Pos = [-30 30];
        YLim_Vel = [-100 100];
        [filt,Data_calc,LE_V,RE_V,Data_cal,Data_In] = MakeCycAvg__filterTraces(type,filt_params_p,filt_params_v,te,ts,Data,keep_inds); 
        CycAvg = MakeCycAvg__makeStruct(LE_V,RE_V,keep_tr,Data,Fs,t_snip,stims,info,filt,In_FileName);
        ha = MakeCycAvg__plotFullCycAvg([],type,colors,line_wid,YLim_Pos,YLim_Vel,te,ts,t_snip,stim,stims,Data,Data_In,Data_cal,Data_calc,LE_V,RE_V,CycAvg,keep_inds,keep_tr);
        % Determine if data are analyzeable
        analyze = nmquestdlg('Are these data analyzeable?','','Yes','No','Yes');
        if strcmp(analyze,'No') 
            save([Cyc_Path,filesep,'NotAnalyzeable_',In_FileName],'Data')
            %If the file was previously designated as analyzeable, remove that file
            d1 = dir([Cyc_Path,filesep,'*.mat']);
            if ~isempty(d1)
                d1 = {d1.name};
                if ismember(['CycAvg_',In_FileName],d1)
                    delete([Cyc_Path,filesep,'CycAvg_',In_FileName])
                end
            end
            return;
        end
    elseif strcmp(opts{ind},'Filter Position')
        %Get new parameter values
        prompt = {['Median',newline,'L X:'],'R X:','L Y:','R Y:','L Z:','R Z:',...
            ['Spline',newline,'L X:'],'R X:','L Y:','R Y:','L Z:','R Z:',...
            ['Sav-Gol 1',newline,'L X:'],'R X:','L Y:','R Y:','L Z:','R Z:',...
             ['Sav-Gol 2',newline,'L X:'],'R X:','L Y:','R Y:','L Z:','R Z:'};
        dlgtitle = 'Filter position';
        definput = strrep(cellfun(@(x) num2str(x,10),num2cell(filt_params_p),'UniformOutput',false),'NaN','');
        temp_filt_params_p = cellfun(@str2double,inputdlgcol(prompt,dlgtitle,[1 10],definput,'on',4,[11 6.5 3.5 4]));
        if ~isempty(temp_filt_params_p) %Didn't hit cancel
            filt_params_p = temp_filt_params_p;
            [filt,Data_calc,LE_V,RE_V,Data_cal,Data_In] = MakeCycAvg__filterTraces(type,filt_params_p,filt_params_v,te,ts,Data,keep_inds); 
            CycAvg = MakeCycAvg__makeStruct(LE_V,RE_V,keep_tr,Data,Fs,t_snip,stims,info,filt,In_FileName);
            ha = MakeCycAvg__plotFullCycAvg(ha,type,colors,line_wid,YLim_Pos,YLim_Vel,te,ts,t_snip,stim,stims,Data,Data_In,Data_cal,Data_calc,LE_V,RE_V,CycAvg,keep_inds,keep_tr);
        end
    elseif strcmp(opts{ind},'Filter Velocity')
        %Get new parameter values
        prompt = {'Irrlsmooth','Spline','Accel Thresh:','Median Filter:'};
        dlgtitle = 'Filter velocity';
        definput = strrep(cellfun(@(x) num2str(x,10),num2cell(filt_params_v),'UniformOutput',false),'NaN','');
        temp_filt_params_v = cellfun(@str2double,inputdlgcol(prompt,dlgtitle,[1 20],definput,'on',1,[11 7.25 2.25 3.25]));
        if ~isempty(temp_filt_params_v) %Didn't hit cancel
            filt_params_v = temp_filt_params_v;
            [filt,Data_calc,LE_V,RE_V,Data_cal,Data_In] = MakeCycAvg__filterTraces(type,filt_params_p,filt_params_v,te,ts,Data,keep_inds); 
            CycAvg = MakeCycAvg__makeStruct(LE_V,RE_V,keep_tr,Data,Fs,t_snip,stims,info,filt,In_FileName);
            ha = MakeCycAvg__plotFullCycAvg(ha,type,colors,line_wid,YLim_Pos,YLim_Vel,te,ts,t_snip,stim,stims,Data,Data_In,Data_cal,Data_calc,LE_V,RE_V,CycAvg,keep_inds,keep_tr);
        end
    elseif strcmp(opts{ind},'Select Cycles')
        keep_tr = MakeCycAvg__selectCycles(type,keep_tr);
        CycAvg = MakeCycAvg__makeStruct(LE_V,RE_V,keep_tr,Data,Fs,t_snip,stims,info,filt,In_FileName);
        ha = MakeCycAvg__plotFullCycAvg(ha,type,colors,line_wid,YLim_Pos,YLim_Vel,te,ts,t_snip,stim,stims,Data,Data_In,Data_cal,Data_calc,LE_V,RE_V,CycAvg,keep_inds,keep_tr);
        %[keep_tr,ha] = MakeCycAvg__selectCycles(ha,type,colors,line_wid,te,t_snip,stims,Fs,Data,info,filt,LE_V,RE_V,keep_inds,keep_tr,In_FileName);
    elseif strcmp(opts{ind},'Set Y-axis limits')
        %Get new parameter values
        prompt = {['Set Y-axis limits',newline,newline,'Position:',newline,newline,'Lower Limit:'],...
            'Upper Limit:',['Velocity:',newline,newline,'Lower Limit:'],'Upper Limit:'};
        dlgtitle = 'Y-axis Limits';
        definput = cellfun(@(x) num2str(x,10),num2cell([YLim_Pos,YLim_Vel]),'UniformOutput',false);
        out_nums = cellfun(@str2double,inputdlgcol(prompt,dlgtitle,[1 18],definput,'on',2,[11 8.25 3 2.25]));
        if ~isempty(out_nums)
            %Check to make sure they aren't reversed
            if out_nums(1)>out_nums(2) 
                YLim_Pos = [out_nums(2),out_nums(1)];
            else
                YLim_Pos = [out_nums(1),out_nums(2)];
            end
            %Check to make sure they aren't reversed
            if out_nums(3)>out_nums(4) 
                YLim_Vel = [out_nums(4),out_nums(3)];
            else
                YLim_Vel = [out_nums(3),out_nums(4)];
            end
        end
        if type==1
            set(ha(1),'YLim',YLim_Pos)
            set(ha(2),'YLim',YLim_Vel)
            set(ha(3),'YLim',YLim_Vel)
            set(ha(4),'YLim',YLim_Vel)
            set(ha(5),'YLim',YLim_Vel)           
        elseif type==2
            set(ha(1),'YLim',YLim_Pos)
            set(ha(2),'YLim',YLim_Vel)
        end
    else
        disp('This case not yet been defined.')
    end
    [ind,tf2] = nmlistdlg('PromptString','Select an action:',...
                       'SelectionMode','single',...
                       'ListSize',[100 100],...
                       'ListString',opts,...
                       'Position',[11,7.75,2,2.75]);  
    if tf2 == 0 %Treat this like an exit
        return;
    end
end
%% Save
CycAvg = MakeCycAvg__makeStruct(LE_V,RE_V,keep_tr,Data,Fs,t_snip,stims,info,filt,In_FileName);
savefig([Cyc_Path,filesep,In_FileName(1:end-4),'.fig'])
close;
%Save if you reached this point
MakeCycAvg__saveCycAvg(Cyc_Path,In_FileName,CycAvg);
end