%% VOGA.m
% This script should be used to segment, filter data, select cycles, and
% parameterize batches of files. 
%
% To initialize after MATLAB has restarted, navigate to the folder with this
% function and type either "VOGA" or "allfoldpath" into the command line 
% to add all functions to the path. If you typed "VOGA", hit cancel now.
%
% To analyze data, navivate to the path with the Raw File, Segmented Files,
% and Cycle Average folders and then run "VOGA."
%
% Click cancel to end the loop
%
% Written by Andrianna Ayiotis
% Updated on 2021-01-31

opts = {'Initialize','Make Folders','Segment','Cycle Average',...
    'Summary Table','Generate Figures','Set Version'};
ind = 1; %Run the start procedure first
tf1 = 1;
while tf1
    if strcmp(opts{ind},'Initialize')        
        %code_Path = [userpath,filesep,'VOGA'];
        code_Path = 'R:\Chow\MATLAB\VOGA';
        addpath(genpath(code_Path));
        path_folders = extractfield(dir,'name',find(extractfield(dir,'isdir')));
        expected = [any(contains(path_folders,'Raw Files')),any(contains(path_folders,'Segmented Files')),any(contains(path_folders,'Cycle Averages'))];
        if ~all(expected)
            uiwait(msgbox('This directory is missing one or more expected folders. Run "Make Folders" to make the missing folders.'))
        end
    elseif strcmp(opts{ind},'Make Folders') 
        %% Make Folders
        %See if the folders already exist or need to be renamed/created
        path_folders = extractfield(dir,'name',find(extractfield(dir,'isdir')));
        if any(contains(path_folders,'Raw LD VOG Files'))
            movefile('Raw LD VOG Files','Raw Files')
        elseif ~any(contains(path_folders,'Raw Files'))
            mkdir('Raw Files')
        end
        if any(contains(path_folders,'CycAvg')) 
            movefile('CycAvg','Cycle Averages')
        elseif any(contains(path_folders,'Cyc_Avg')) 
            movefile('Cyc_Avg','Cycle Averages')
        elseif ~any(contains(path_folders,'Cycle Averages'))
            mkdir('Cycle Averages')
        end
        if ~any(contains(path_folders,'Segmented Files'))
            mkdir('Segmented Files')
        end
    elseif strcmp(opts{ind},'Segment') 
        %% Segment
        Raw_Path = [cd,filesep,'Raw Files'];
        Seg_Path = [cd,filesep,'Segmented Files'];
        %Transfer NKI Raw Files from their subfolders if they exist
        moveNKIfiles(Raw_Path)
        %Detect and process log files and austoscan files
        logtoNotes(Raw_Path)
        % Select files to segment
        all_files = extractfield(dir(Raw_Path),'name',find(~extractfield(dir(Raw_Path),'isdir')));
        VOG_files = all_files(((contains(all_files,'SESSION')&contains(all_files,'.txt'))|contains(all_files,'.dat'))...
            &~contains(all_files,'-Notes.txt'));
        if isempty(VOG_files) 
            uiwait(msgbox('No VOG files found in the Raw Files folder.'))
        else
            [indx,tf] = nmlistdlg('PromptString','Select files to segment:','ListSize',[300 300],'ListString',VOG_files);
            if tf == 1
                sel_files = VOG_files(indx);
                for i = 1:length(sel_files)
                    In_Path = [Raw_Path,filesep,sel_files{i}];
                    Segment(In_Path,Seg_Path)
                end
            end
        end
    elseif strcmp(opts{ind},'Cycle Average')
        %% Filter and select cycles
        path = cd; % Assumes you are in the right directory already
        Seg_Path = [path,filesep,'Segmented Files'];
        Cyc_Path = [path,filesep,'Cycle Averages'];
        % Get version and experimenter info from the file
        if ~any(contains(extractfield(dir(code_Path),'name'),'VerInfo.txt'))
            writeInfoFile(code_Path);
        end
        data = readtable([code_Path,filesep,'VerInfo.txt'],'ReadVariableNames',false);
        version = data{1,2}{:};
        Experimenter = data{2,2}{:};
        done = false;
        % Experiment types
        progress_tab = assessProgress(path);
        all_exp_names = progress_tab{:,1};
        for i = 1:length(all_exp_names)
            dash = strfind(all_exp_names{i},'-');
            all_exp_names{i} = all_exp_names{i}(1:dash(3)-1);
        end
        exp_names = unique(all_exp_names);
        if length(exp_names)>1
            [indx,tf] = nmlistdlg('PromptString','Select experiment types to analyze:',...
                           'SelectionMode','multiple',...
                           'ListSize',[350 300],...
                           'ListString',exp_names); 
        else
            tf = 0;
        end
        if tf == 0
            exp_types = {};
        else
            exp_types = exp_names(indx);
        end
        while(~done) %run until the user hits cancel on analyzing a file
            done = MakeCycAvg(path,Seg_Path,Cyc_Path,Experimenter,version,exp_types);
        end
    elseif strcmp(opts{ind},'Summary Table')
        path = cd; % Assumes you are in the right directory already
        Cyc_Path = [path,filesep,'Cycle Averages'];
        MakeCycleSummaryTable(path,Cyc_Path)
    elseif strcmp(opts{ind},'Generate Figures')
        path = cd; % Assumes you are in the right directory already
        Seg_Path = [path,filesep,'Segmented Files'];
        Cyc_Path = [path,filesep,'Cycle Averages'];
        % Get version and experimenter info from the file
        if ~any(contains(extractfield(dir(code_Path),'name'),'VerInfo.txt'))
            writeInfoFile(code_Path);
        end
        data = readtable([code_Path,filesep,'VerInfo.txt'],'ReadVariableNames',false);
        version = data{1,2}{:};
        Experimenter = data{2,2}{:};
        return;
    elseif strcmp(opts{ind},'Set Version')
        writeInfoFile(code_Path);
    end
    %% Poll for new reponse
    [ind,tf1] = nmlistdlg('PromptString','Select an action:',...
                       'SelectionMode','single',...
                       'ListSize',[100 100],...
                       'ListString',opts,...
                       'Position',[0,7.75,2,2.75]);    
end
disp('VOGA instance ended.')