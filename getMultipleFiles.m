function [mpuFileCells, coilFileCells, freq, amplitude, counterForFiles, tiltAxis, tiltDeg]=getMultipleFiles(directory,movementDirection,theta,phi)

%{
INPUT:
    directory - string containing the directory folder
    direction - the direction of movement that you want to get files for
    theta - only used if doing oblique angle files - should be angle in
    horizontal plane, with 0deg out R ear, and 180 out L ear, 90 out nose
    phi - angle from vertical. Up = 0. Out nose=+90
OUTPUT:
    mpuFileCells - cell array of the mpu file names
    coilFileCells - cell array of the coil file names
    freq - the frequencies of each file in order
    amp - the amplitudes of each file in order
    counterForFiles - total number of files that are sent

NOTE: if you get error with this function that says freqCell not found, there
weren't any files found so check your directory name and the angle numbers
you are putting into the function

%}


%make a structure of the list of files in the directory
filesList=dir(directory);
counterForFiles=1;
%first find the .coil files
idxCoilFile=~cellfun('isempty',strfind({filesList.name},'.coil'));

%then find the files for direction we want
idxDirectionFile=~cellfun('isempty',strfind({filesList.name},movementDirection));

%take out any impulse files -MRC edit 20160927
idxImpulseFile = cellfun('isempty',strfind({filesList.name},'Impulse'));


%multiply the arrays to find the subset of files we want
idxToList=idxCoilFile.*idxDirectionFile.*idxImpulseFile;


if(~strcmp(movementDirection,'StaticTilt'))
    if (strcmp(movementDirection,'ObliqueAngle')||strcmp(movementDirection,'ObliqueAngleCombo')||strcmp(movementDirection,'ObliqueAngleHorizontal')||strcmp(movementDirection,'ObliqueAngleSaggital')||strcmp(movementDirection,'ObliqueAngleCoronal'))
        
        %because files are named slightly differently, have to change the 'string find' commands
        if(strcmp(movementDirection,'ObliqueAngleHorizontal')||strcmp(movementDirection,'ObliqueAngle'))
            idxAngle=~cellfun('isempty',strfind({filesList.name},strcat('-',num2str(theta),'deg')));
        else
            idxAngle=~cellfun('isempty',strfind({filesList.name},strcat('-Theta',num2str(theta),'deg')));
        end
        idxToList=idxToList.*idxAngle;
        
        if(phi~=90) %if phi=90, then its a Horizontal, and we didn't include phi in file name
            idxPhi=~cellfun('isempty',strfind({filesList.name},strcat(num2str(phi),'deg')));
            idxToList=idxToList.*idxPhi;
        end
    end
    
    
    %Now use the idxToList to figure out which file names to compile and return
    for i=1:length(idxToList)
        if (idxToList(i)==1)
            coilFileCells(counterForFiles)=cellstr(filesList(i).name);
            mpuFileCells(counterForFiles)=regexprep(coilFileCells(counterForFiles),'.coil','_MPUdata.txt');
            
            %if for translation (needs the mpsq)
            if (strcmp(movementDirection,'ObliqueAngle')||strcmp(movementDirection,'ObliqueAngleCombo')||strcmp(movementDirection,'Lateral')||strcmp(movementDirection,'ObliqueAngleHorizontal')||strcmp(movementDirection,'ObliqueAngleSaggital')||strcmp(movementDirection,'ObliqueAngleCoronal')||strcmp(movementDirection,'Surge')||strcmp(movementDirection,'Heave'))
                searchForFreq(counterForFiles)=regexp(coilFileCells(counterForFiles),'(?<=(q-))[^"]+(?=(hz))','match');
                freqCell(counterForFiles)=regexprep(searchForFreq{counterForFiles},'p','.');
                
                searchForAmp(counterForFiles)=regexp(coilFileCells(counterForFiles),'(?<=-)[^"]+(?=(mpsq-))','match');
                ampCell(counterForFiles)=regexprep(searchForAmp{counterForFiles},'p','.');
            else
                searchForFreq(counterForFiles)=regexp(coilFileCells(counterForFiles),'(?<=(s-))[^"]+(?=(hz))','match');
                freqCell(counterForFiles)=regexprep(searchForFreq{counterForFiles},'p','.');
                
                searchForAmp(counterForFiles)=regexp(coilFileCells(counterForFiles),'(?<=-)[^"]+(?=(dps-))','match');
                ampCell(counterForFiles)=regexprep(searchForAmp{counterForFiles},'p','.');
            end
            counterForFiles=counterForFiles+1;
        end
    end
    
    counterForFiles=counterForFiles-1;
    freq=str2double(freqCell);
    amplitude=str2double(ampCell);
    tiltAxis=[];
    tiltDeg=[];
    
elseif strcmp(movementDirection,'StaticTilt')
    for i=1:length(idxToList)
        if (idxToList(i)==1)
            coilFileCells(counterForFiles)=cellstr(filesList(i).name);
            mpuFileCells(counterForFiles)=regexprep(coilFileCells(counterForFiles),'.coil','_MPUdata.txt');
            
            %search for degree of tilt
            searchForTiltPhrase(counterForFiles)=regexp(coilFileCells(counterForFiles),'(?<=(-))[^"]+(?=(deg))','match');
            
            %search for Axis you Tilt around
            searchForTiltAxis(counterForFiles)=regexp(searchForTiltPhrase{counterForFiles},'((?<=Theta)\w*)','match'); %looks for phrase after Theta
            tiltAxisCell(counterForFiles)=regexprep(searchForTiltAxis{counterForFiles},'Neg','-'); %replace 'Neg' with a - sign
            
            
            searchForTiltDeg(counterForFiles)=regexp(searchForTiltPhrase{counterForFiles},'((?<=(-))\d*)','match');
            tiltDegCell(counterForFiles)=searchForTiltDeg{counterForFiles};
            
            %searchFor
            
            counterForFiles=counterForFiles+1;
        end
    end
    
    counterForFiles=counterForFiles-1;
    freq=[];
    amplitude=[];
    tiltDeg=str2double(tiltDegCell);
    tiltAxis=str2double(tiltAxisCell);
    
    
    
%    make3Dplot(false,true,true,53,sineR,sineL,movementDirection, freq(m), color, 4, 4);
    
    
end

end