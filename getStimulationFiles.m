function [coilFileCells, stimParamCells, counterForFiles, stimE, refE, stim2E, ref2E]=getStimulationFiles(directory,endOrgan,typeOfStim)

%{
INPUT:
    directory - string containing the directory folder
    endOrgan - which end organ are you stimulation? Utricle, Saccule, LH, LA, LP
    typeOfStim - what type of file? PulseTrain, StaticTilt, VirtualTranslation, VirtualRotation, StaticBaseline
    
OUTPUT:
    coilFileCells - cell array of the coil file names
    stimParamCells - cell array of the txt file names that contain stim parameters
    counterForFiles - total number of files that are sent
    stimE - the stimulating electrode number
    refE - the reference electrode number
    stim2E - if the file used dual bipolar
    ref2E - if the file used dual bipolar

%}


%make a structure of the list of files in the directory
filesList=dir(directory);
counterForFiles=1;
%first find the .coil files
idxCoilFile=~cellfun('isempty',strfind({filesList.name},'.coil'));

%then find the files for endOrgan we want
idxEndOrganFile=~cellfun('isempty',strfind({filesList.name},endOrgan));

%then files for type of stim we want.
idxStimTypeFile=~cellfun('isempty',strfind({filesList.name},typeOfStim));

%multiply the arrays to find the subset of files we want
idxToList=idxCoilFile.*idxEndOrganFile.*idxStimTypeFile;



%Now use the idxToList to figure out which file names to compile and return
for i=1:length(idxToList)
    if (idxToList(i)==1)
        coilFileCells(counterForFiles)=cellstr(filesList(i).name);
        stimParamCells(counterForFiles)=regexprep(coilFileCells(counterForFiles),'.coil','_ExptDetail.txt');
        
        stimE(counterForFiles)=regexp(coilFileCells(counterForFiles),'((?<=Stim)\d*)','match');
        %stimE(counterForFiles)=regexprep(searchForFreq{counterForFiles},'p','.');
            
        %searchForRefE(counterForFiles)=regexp(coilFileCells(counterForFiles),'(?<=-)[^"]+(?=(dps-))','match');
        refE(counterForFiles)=regexp(coilFileCells(counterForFiles),'((?<=Ref)\d*)','match');
        
        counterForFiles=counterForFiles+1;
    end
end

counterForFiles=counterForFiles-1;




end