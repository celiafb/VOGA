%% Load in data
clc; clear all;
directory='/Users/celia/Desktop/20210120-Ch156-normTVOR-REarGent/';
%directory = '/Volumes/labdata/Hageman/Moog/MoogChinchData-KH/20190523-Ch141-normTVOR/';
typeOfData='Movement'; %'Movement' data? Or 'ProsthesisOnly' Data?
 
afterProsthSyncDate=true; %after 10/19/2016?
gainsInFile=false; %are gains stored in file? If false - make sure you have a txt file with gain values for the script to read.
gainFile='20210120-Ch156-gains.txt';
haveOffsetFiles=true; %do you have offset files? if true, fill them in below. if false - fill in ZEROS below
coilOrientation1='20210120_095156_Ch156_offset_orientation1.coil';
coilOrientation2='20210120_095501_Ch156_offset_orientation2.coil';
 
 
 
Fs = 1000;
frame=0; %frame=0 is head frame of ref,  frame=1 is eye frame of ref
REF=2000; %point of reference for analysis - usually use the start of movemnet
%use 10000 for prosthesis static tilt. use 9000 for movement static tilt
%use 5000 for prosthesis sinusoids
%use 2000 for normal rotation sinusoids
 
positionFiles=true; %For analysis of positional data
velocityFiles=true; %For analysis of velocity data
 
%% Set Up Variables Unique to whether Motion or Prosthesis Files
switch(typeOfData)
    case 'Movement'
        movementDirection='Yaw'; %'Yaw' 'LARP' 'RALP' 'Pitch' 'Roll' 'ObliqueAngleHorizontal' 'ObliqueAngleSaggital' 'ObliqueAngleCoronal' 'ObliqueAngleCombo' 'Lateral' 'Surge' 'Heave' 'StaticTilt'
        %values only used for oblique angle
        theta=0;
        phi=90;
        
        [mpufileList,coilfileList,freqList,ampList,numOfFiles,tiltAxisList,tiltDegList]=getMultipleFiles(directory,movementDirection,theta,phi);
        
        fdr = strcat('Analyze_',gainFile(1:8),'_',gainFile(10:14),'_',typeOfData,'_',movementDirection);
        mkdir(directory, fdr); % here we will store the files to then analyze with VOGA
        
    case 'ProsthesisOnly'
        typeOfStim='VirtSineRotation';%'VirtSineTranslation'; %'PulseTrains' 'VirtSineRotation' 'StaticTilt' '
        endOrgan='LH';
        
        if(strcmp(typeOfStim,'PulseTrains'))
            freq=2;
        elseif (strcmp(typeOfStim,'VirtSineTranslation') || strcmp(typeOfStim,'VirtSineRotation'))
            freq=1;
        else
            freq=0; %for static tilts - this is just place holder
        end
        
        [coilfileList, stimParamCells, numOfFiles, stimE, refE]=getStimulationFiles(directory,endOrgan,typeOfStim);

        fdr = strcat('Analyze_',gainFile(1:8),'_',gainFile(10:14),'_',typeOfData,'_',typeOfStim,'_',endOrgan);
        mkdir(directory, fdr); % here we will store the files to then analyze with VOGA
        
end
%% Calculate/Load Gains
if (gainsInFile==true)
    tempGains=loadGains(directory,coilfileList{1}); %load gains from first file - assume all other files are same
        if tempGains == -1
            set(handles.messages,'String','Incorrect Coil File');
        else
            for i=1:3 %x, y, z
                ch0(i)=tempGains(i,1); %ch1
                ch1(i)=tempGains(i,2); %ch2
                ch2(i)=tempGains(i,3); %ch3
                ch3(i)=tempGains(i,4); %ch4
            end
        end
        GAINSR=[ch0(1:3),ch1(1:3)];
        GAINSL=[ch2(1:3),ch2(1:3)];
elseif(gainsInFile==false && strcmp(gainFile,'')) % fill these values in if you don't have a file with gains stores and gains arent saved in the data file itself
    GAINSR=[0.0926 -0.1702 -0.2310 -0.1162 0.2039 0.2897]; %For Ch133 on 5/11/17
    GAINSL=[-0.1453 0.2590 0.3367 0.1622 -0.2953 -0.3562]; %For Ch133 on 5/11/17
elseif (~strcmp(gainFile,'')) %if you have a txt file with gains saved
    [GAINSR, GAINSL]=getGains(directory,gainFile);
 
    GAINSR=[GAINSR(:,1);GAINSR(:,2)];
    GAINSL=[GAINSL(:,1);GAINSL(:,2)];
end
%% Calculate Offsets
if (haveOffsetFiles==true)
    if (afterProsthSyncDate)
        [ZEROS_R, ZEROS_L]=calcOffsets(strcat(directory,'/OffsetCheck/'),coilOrientation1,coilOrientation2,1);
    else
        [ZEROS_R, ZEROS_L]=calcOffsets(strcat(directory,'/OffsetCheck/'),coilOrientation1,coilOrientation2,0);
    end
else
    ZEROS_R=[7.58e7 -2.2e7 -1.15e8 2.85e7 1.65e7 -4.49e7];  
    ZEROS_L=[5.59e7 -2.5e6 -1.27e8 -1.02e8 4.49e7 2.08e8];
end
%% Movement Only - get MPU column
 
if (strcmp(typeOfData,'Movement'))
        % Column 3: Lateral (Left/Right - our +/- Y)
        % Column 4: Surge (Front/Back - our +/- X)
        % Column 5: Heave (Up/Down - our +/- Z)
        % Column 6: Roll (rotate around naso/occipital)
        % Column 7: Pitch (rotate around interaural)
        % Column 8: Yaw (rotate around earth vertical)
        if (strcmp(movementDirection,'LARP') || strcmp(movementDirection,'RALP'))
            sensorColumn=100;
        elseif (strcmp(movementDirection,'ObliqueAngle')||strcmp(movementDirection,'ObliqueAngleCombo')||strcmp(movementDirection,'ObliqueAngleHorizontal')||strcmp(movementDirection,'ObliqueAngleCoronal')||strcmp(movementDirection,'ObliqueAngleSaggital'))
            sensorColumn=200;
        elseif (strcmp(movementDirection,'StaticTilt'))
            sensorColumn=300;
        else
            sensorColumn=getSensorColumn(movementDirection);
        end
end
%% Load Data
cd(fullfile(directory,fdr));
for m = 1:numOfFiles
    
    clear mpuData coilData tempGains filteredRawData rotR rotL mpuAligned filteredRawAngVelR filteredRawAngVelL
    clear indexZeros newIndexZeros mpuShorten angVelRShorten angVelLShorten startCycleToCut endCycleToCut
    clear projectedCoilAngVelR projectedCoilAngVelL RPY projectedMPUAngVel filteredMPU
    clear coilDataWithSyncs coilData prosthSync TS_idx TS_time TS_interval TS_rate filteredRotrefL filteredRotrefR
    clear rotrefR rotrefL angVelR angVelL jointSineR jointSineL motionFit maxHeadMovementPointR filteredRawData
    clear maxHeadMovementPointL pshortenedFileR pshortenedFileL pshortenedMotion
    clear pmotionSync motionSync  startIdx vmotionSync
    
    coilFile=coilfileList{m};
    
    switch (typeOfData)
        case 'Movement'
            mpufile=mpufileList{m};
            
            if(~strcmp(movementDirection,'StaticTilt'))
                freq=freqList(m);
                
                %Use different NFilt Options based on the frequency of the sinusoidal
                %rotation - to avoid having a filter cutoff frequency too close to the
                %frequency of the file.
                if(isempty(freqList)||freq>0.5) %for frequencies above 0.5 Hz
                    NFilt=75;
                else
                    NFilt=100; %for frequencies less than or equal to 0.5 Hz
                end
            else
                tiltDeg=tiltDegList(m);
                tiltAxis=tiltAxisList(m);
                
                NFilt=75;
            end
            
            
            %Read in the mpu data and coil data
            mpuData=readmpu(directory,mpufile,0);
            if (afterProsthSyncDate==true)
                coilDataWithSyncs=readcoils(directory,coilFile,true);
                coilData(:,1:3)=coilDataWithSyncs(:,1:3);
                coilData(:,4:15)=coilDataWithSyncs(:,6:17);
            else
                coilData=readcoils(directory,coilFile,false);
            end
            
        case 'ProsthesisOnly'
            stimParamFile=stimParamCells{m};
            NFilt=50; %all sinusoidal data at this point is 1 Hz, so this filter value is fine for all stim files
            
            %Load coil data
            coilDataWithSyncs=readcoils(directory,coilFile,true);
            coilData(:,1:3)=coilDataWithSyncs(:,1:3);
            coilData(:,4:15)=coilDataWithSyncs(:,6:17);
            prosthSync(:,1:2)=coilDataWithSyncs(:,4:5);
            
            TS_idx=1+find(diff(prosthSync(:,1)));  %returns index values of sync signal for stimulation pulses
            TS_time= prosthSync(TS_idx, 1)-1 + prosthSync(TS_idx,2)/25000; %calculate actual time based on prosthesis sub-sample (column2 of sync). It is based on the 25000 samples avg'd for the 1kHz coil system sample
            TS_interval= diff(TS_time) / 1000; %for interval between pulses, divide by sampling frequency
            TS_rate= 1./TS_interval; %rate = 1/interval
            
    end
    %% Filter, Align, Rotational Kinematics
    
    %Filter the raw data before aligning it
    filteredRawData(:,1:3)=coilData(:,1:3);
    for i=4:15
        filteredRawData(:,i)=filtfilt(ones(1,NFilt)/NFilt,1,coilData(:,i));
    end
    
    switch(typeOfData)
        case 'Movement'
            %%Raw to Rot of the filtered raw data
            [filteredRotrefR,filteredRotrefL,mpuAligned]=align(mpuData,filteredRawData,REF,frame,[],GAINSR,GAINSL,ZEROS_R,ZEROS_L);
            
            filteredMPU(:,1:2)=mpuAligned(:,1:2);
            for i=3:8
                filteredMPU(:,i)=filtfilt(ones(1,NFilt)/NFilt,1,mpuAligned(:,i));
            end
            
            [unfilteredRotR,unfilteredRotL,unfilteredMpuAligned]=align(mpuData,coilData,REF,frame,[],GAINSR,GAINSL,ZEROS_R,ZEROS_L);
            
            %Rot2AngVel with the filtered data
            angVelRpreSpline=rot2angvel(filteredRotrefR)/pi*180 * 1000;
            angVelLpreSpline=rot2angvel(filteredRotrefL)/pi*180 * 1000;
            
            if (strcmp(movementDirection,'Yaw')||strcmp(movementDirection,'LARP') || strcmp(movementDirection,'RALP'))
                [angVelR, angVelL]=desacc_prosthCanal(angVelRpreSpline,angVelLpreSpline);
            else
                angVelR=angVelRpreSpline;
                angVelL=angVelLpreSpline;
            end
            
            
            mpuAligned(:,3:5)=mpuAligned(:,3:5)/65536*9.8*8; %For acceleration in m/s2
            filteredMPU(:,3:5)=filteredMPU(:,3:5)/65536*9.8*8; %For acceleration in m/s2
            mpuAligned(:,6:8)=mpuAligned(:,6:8)/65536*500; %For velocity
            filteredMPU(:,6:8)=filteredMPU(:,6:8)/65536*500; %For velocity
            
            %% For the quick analysis
            
            %Choose which mpu column or signal to use for analysis
            if (sensorColumn==200) %oblique angle
                %Need to extract the mag and direction of motion in direction of specific angle w/ unit vector (cos(ang),sin(ang)) - dot product
                
                mpuToUse = mpuAligned(:,3).*cosd(theta).*sind(phi) + mpuAligned(:,4).*sind(theta).*sind(phi)+mpuAligned(:,5).*cosd(phi);
                figure; plot(mpuAligned(:,3:5)); hold on; plot(mpuToUse); legend('Y','X','Z','4');
            elseif(sensorColumn==100) %LARP and RALP
                if(strcmp(movementDirection,'LARP'))
                    ang = (45-90)*pi/180;
                else
                    ang = (45)*pi/180;
                end
                
                RPY(:,1) = mpuAligned(:,7); %dps
                RPY(:,2) = mpuAligned(:,6);
                RPY(:,3) = mpuAligned(:,8);
                alpha = 0;
                beta = 0;
                gamma = 0;
                dt = 0.001;
                projectedAngVel = zeros(1,length(mpuAligned));
                for i = 1:length(RPY(:,1))
                    %projectedCoilAngVel(i,1) = rotR
                    alpha = alpha + dt*RPY(i,3)*pi/180;
                    beta = beta + dt*RPY(i,2)*pi/180;
                    gamma = gamma + dt*RPY(i,1)*pi/180;
                    R = rpyToMat(gamma,beta,alpha);
                    projectedAngVel(i) = (R*RPY(i,:)')'*[cos(ang);sin(ang);0];
                end
                
                z2=filtfilt(ones(1,50)/50,1,projectedAngVel);
                mpuToUse = projectedAngVel';
                
                
            elseif(strcmp(movementDirection,'StaticTilt'))
                %had to add 90 to tiltAxis value since defined differently than the horizontal plane oblique angles
                %tiltAxis=0=nasal ND. tiltAxis=90=LED. tiltAxis=-90=RED. tiltAxis=180=NU
                % IMPORTANT: MRC changed 6/25/2018 - Moog is now correct,
                % do not need to add 90 degrees
                mpuToUse = -filteredMPU(:,3).*cosd(tiltAxis) + filteredMPU(:,4).*sind(tiltAxis); %phi is 90 for all of these so I just removed those sind(90) since =1 and deleted cosd(90) term
                tiltTrace=abs(asind(mpuToUse./9.8));
                
            else %yaw pitch roll lateral surge heave
                mpuToUse=mpuAligned(:,sensorColumn);
            end
        case 'ProsthesisOnly'
            %Find rotation vectors of eye movements - don't need to align since sync signal is collected with the coil system thus already aligned
            [rotrefR]=gimbalproc(coilData(:,4:9), ZEROS_R, GAINSR, REF);
            [rotrefL]=gimbalproc(coilData(:,10:15), ZEROS_L, GAINSL ,REF);
            %^^for this, use rotrefR and rotrefL because that is in head reference frame
            [filteredRotrefR]=gimbalproc(filteredRawData(:,4:9), ZEROS_R, GAINSR, REF);
            [filteredRotrefL]=gimbalproc(filteredRawData(:,10:15), ZEROS_L, GAINSL ,REF);
            
            %Rot2AngVel with the filtered data
            angVelRpreSpline=rot2angvel(filteredRotrefR)/pi*180 * 1000;
            angVelLpreSpline=rot2angvel(filteredRotrefL)/pi*180 * 1000;
            
            if (strcmp(typeOfStim,'VirtSineRotation')||strcmp(typeOfStim,'VirtSineTranslation'))
                [angVelR, angVelL]=desacc_prosthCanal(angVelRpreSpline,angVelLpreSpline);
                %Meg edits 12/24/20
                if contains(coilFile,'hz')
                    ind = strfind(coilFile,'hz');
                    indMin = 1;
                    while (~isequal(coilFile(ind-indMin),'-'))
                        indMin = indMin + 1;
                    end
                    tempstr = coilFile(ind-indMin+1:ind-1);
                    if contains(tempstr,'p')
                        pind = strfind(tempstr,'p');
                        tempstr(pind) = '.';
                    end
                    freq = str2double(tempstr);
                end
            else
                angVelR=angVelRpreSpline;
                angVelL=angVelLpreSpline;
            end
            
            
    end
    %% Get rid of initial and final ramping time and any extra time
    if strcmp(typeOfData,'ProsthesisOnly') %don't need ramp cycles since sync pulse triggers start time
        motionSync=[TS_idx(2:end) 1./TS_interval./10];
        analysisType=typeOfStim;
        movementDirection='Prosthesis';
        
        startIdx=motionSync(1,1); %this triggers when the stimulation started
        stopIdx=motionSync(end,1); %this is the last sync signal meaning the pulses stopped
        
        pshortenedFileR=filteredRotrefR(startIdx:stopIdx,:); %
        pshortenedFileL=filteredRotrefL(startIdx:stopIdx,:);
        vshortenedFileR=angVelR(startIdx:stopIdx,:); %
        vshortenedFileL=angVelL(startIdx:stopIdx,:);
        
        pshortenedTT=1/1000:0.001:length(pshortenedFileR)/1000;
        pshortenedSyncs=zeros(length(pshortenedTT),1);
        pmotionSync(:,1)=motionSync(:,1)-startIdx+1;
        vshortenedTT=1/1000:0.001:length(vshortenedFileR)/1000;
        vshortenedSyncs=zeros(length(vshortenedTT),1);
        vmotionSync(:,1)=motionSync(:,1)-startIdx+1;
        
        for i=1:length(motionSync)
            if (motionSync(i,1)>0 && motionSync(i,1)<length(pshortenedTT))
                pshortenedSyncs(motionSync(i,1))=motionSync(i,2);
            end
        end
        pshortenedMotion=pshortenedSyncs;
        
        vshortenedMotion=pshortenedSyncs;
        
    else
        motionSync=mpuToUse;
        analysisType=movementDirection;
        numRamps=2;
        if(strcmp(movementDirection,'LARP')||strcmp(movementDirection,'RALP')||strcmp(movementDirection,'Yaw'))
            timePreAndPost=2; %seconds before motion/stimulation starts
        elseif(strcmp(direction,'Prosthesis'))
            timePreAndPost=10; %seconds before motion/stimulation starts
        else
            timePreAndPost=5; %seconds before motion/stimulation starts
        end
        
        extraTime=timePreAndPost*Fs+numRamps*Fs/freq;
        
        pshortenedFileR=filteredRotrefR(extraTime:(length(filteredRotrefR)-extraTime),:);
        pshortenedFileL=filteredRotrefL(extraTime:(length(filteredRotrefL)-extraTime),:);
        vshortenedFileR=angVelR(extraTime:(length(angVelR)-extraTime),:);
        vshortenedFileL=angVelL(extraTime:(length(angVelL)-extraTime),:);
        
        
        pshortenedMotion=motionSync(extraTime:(length(motionSync)-extraTime)); %don't need the -1 for position
        vshortenedMotion=motionSync(extraTime:(length(motionSync)-extraTime-1)); %Need the -1 for velocity
        
        headVel = filteredMPU(extraTime:(length(filteredMPU)-extraTime),6:8); % For VOGA compatibility
        headAccel = filteredMPU(extraTime:(length(filteredMPU)-extraTime),3:5); % For VOGA compatibility
        
        
        
        ptt=1/1000:0.001:length(filteredRotrefR)/1000;
        pshortenedTT=1/1000:0.001:length(pshortenedFileR)/1000;
        vtt=1/1000:0.001:length(angVelR)/1000;
        vshortenedTT=1/1000:0.001:length(vshortenedFileR)/1000;
        startIdx = extraTime;
        stopIdx = length(motionSync)-extraTime-1;
    end
    
    
    
    
    %% Finally, save in a Data structure and save as .mat file
    if strcmp(movementDirection,'Yaw')
        movementDirection = 'LHRH'; % VOGA convention
    end
    if strcmp(typeOfData,'Movement')
        typeOfData = 'MOOG'; % movement, VOGA convention
        if strcmp(movementDirection,'Yaw') || strcmp(movementDirection,'LARP') || strcmp(movementDirection,'RALP') || strcmp(movementDirection,'Pitch') || strcmp(movementDirection,'Roll')
            exptcond = 'Rotation';
            fname = sprintf(strcat(coilFile(17:21),'-','VisitNA','-',coilFile(1:8),'-',typeOfData,'-Sine-',exptcond,'-NoStimLight','-',movementDirection,'-',coilFile(dash_indices(1):dash_indices(end)-1),'.mat'));
            aaaa = 1;
        else
            exptcond = 'Translation';
        end
    else
        typeOfData = 'eeVOR'; % prosthesis, VOGA convention
    end
    
    pshortenedFileL = rot2fick(pshortenedFileL); % inputs XYZ and outputs ZYX
    pshortenedFileR = rot2fick(pshortenedFileR);
    
    var{m}.Fs = 1000;
    var{m}.Time_Eye = pshortenedTT;
    var{m}.Time_Stim = pshortenedTT;
    var{m}.Trigger = pshortenedMotion;
    var{m}.LE_Position_X = pshortenedFileL(:,3);
    var{m}.LE_Position_Y = pshortenedFileL(:,2);
    var{m}.LE_Position_Z = pshortenedFileL(:,1);
    var{m}.RE_Position_X = pshortenedFileR(:,3);
    var{m}.RE_Position_Y = pshortenedFileR(:,2);
    var{m}.RE_Position_Z = pshortenedFileR(:,1);
    if strcmp(typeOfData,'Movement')
        var{m}.HeadVel_X = headVel(:,1);
        var{m}.HeadVel_Y = headVel(:,2);
        var{m}.HeadVel_Z = headVel(:,3);
        var{m}.HeadAccel_X = headAccel(:,1);
        var{m}.HeadAccel_Y = headAccel(:,2);
        var{m}.HeadAccel_Z = headAccel(:,3);
    else
        var{m}.HeadVel_X = [];
        var{m}.HeadVel_Y = [];
        var{m}.HeadVel_Z = [];
        var{m}.HeadAccel_X = [];
        var{m}.HeadAccel_Y = [];
        var{m}.HeadAccel_Z = [];
    end
    var{m}.info.subject = '';
    var{m}.info.ear = '';
    var{m}.info.visit = '';
    var{m}.info.goggle_ver = '';
    var{m}.info.goggle_reorient_ang = '';
    var{m}.info.exp_date = gainFile(1:8);
    dir2 = dir(fullfile(directory,'*.docx'));
    var{m}.info.rawnotes = strcat(dir2.folder,'/',dir2.name);
    var{m}.info.rawfile = strcat(dir2.folder,'/',coilFile);
    var{m}.info.name = coilFile;
    under_indices = strfind(coilFile,'_');
    dash_indices = strfind(coilFile,'-');
    coil_indices = strfind(coilFile,'.coil');
    var{m}.info.dataType = coilFile(under_indices(end)+1:coil_indices(1)-1);
    
    Data  = var{m};
    save(fname,'Data');
    
end

fclose('all');
% So at this point we have, for both prosthesis or movement data, the
% actual data in pShortenedFileR/L, vShortenedFileR/L,
% p/vShortenedMotion (not sure which to use) with p/vShortenedTT