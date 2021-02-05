function [jointSineR, jointSineL, motionFit, maxHeadMovementPointR, maxHeadMovementPointL, shortenedFileR, shortenedFileL, shortenedMotion, startIdx, stopIdx]=getQuickFit(dataR,dataL,movementData,direction,freq,dataType,showPlot,customStatus,startCycleToCut,endCycleToCut,coilFile)
%dataType=0 for position, 1 for velocity

%INPUT VARIABLES:
%dataR - coil data (either position or velocity) of Right eye
%dataL - coil data (either position or velocity) of Left eye
%movementData - data from mpu or sync signal. If using sync signal, 1st column should be idxValues, 2nd column should be rate values.
%direction - used for actual motion to determine time pre and post since its different for rotation over translation. If stimulation: 'Prosthesis'
%freq - frequency for data fit
%dataType - 0 for position, 1 for velocity
%showPlot - true if you want to show the plot, false if you don't
%customStatus - true if you've cut cycles. false if you haven't cut cycles and want to plot the whole trace
%startCycleToCut - only need value if customStatus=true. how many start cycles do you want to cut?
%endCycleToCut - only need value if customStatus=false. how many end cycles do you want to cut?

Fs=1000;

if (strcmp(direction,'Prosthesis')) %don't need ramp cycles since sync pulse triggers start time
    
     if(customStatus==false) %if we aren't cutting any cycles
        startIdx=movementData(1,1); %this triggers when the stimulation started
        stopIdx=movementData(end,1); %this is the last sync signal meaning the pulses stopped         
     else %if we are cutting cycles
        startIdx=movementData(1,1)+(startCycleToCut*1000/freq); %this triggers when the stimulation started
        stopIdx=movementData(end,1)-(endCycleToCut*1000/freq); %this is the last sync signal meaning the pulses stopped     
     end  
         
        shortenedFileR=dataR(startIdx:stopIdx,:); %
        shortenedFileL=dataL(startIdx:stopIdx,:);

        shortenedTT=1/1000:0.001:length(shortenedFileR)/1000;
        shortenedSyncs=zeros(length(shortenedTT),1);
        movementData(:,1)=movementData(:,1)-startIdx+1;
        
        for i=1:length(movementData)
            if (movementData(i,1)>0 && movementData(i,1)<length(shortenedTT))
                shortenedSyncs(movementData(i,1))=movementData(i,2);
            end
        end
        
        shortenedMotion=shortenedSyncs;
        
else
    numRamps=2;
    if(strcmp(direction,'LARP')||strcmp(direction,'RALP')||strcmp(direction,'Yaw'))
        timePreAndPost=2; %seconds before motion/stimulation starts
    elseif(strcmp(direction,'Prosthesis'))
        timePreAndPost=10; %seconds before motion/stimulation starts
    else
        timePreAndPost=5; %seconds before motion/stimulation starts
    end
    
    if(customStatus==false)
        
        
        extraTime=timePreAndPost*Fs+numRamps*Fs/freq;
        
        shortenedFileR=dataR(extraTime:(length(dataR)-extraTime),:);
        shortenedFileL=dataL(extraTime:(length(dataL)-extraTime),:);
        
        if(dataType==0)%For position
            shortenedMotion=movementData(extraTime:(length(movementData)-extraTime)); %don't need the -1 for position
        elseif (dataType==1) %For velocity
            shortenedMotion=movementData(extraTime:(length(movementData)-extraTime-1)); %Need the -1 for velocity
        end
        
        tt=1/1000:0.001:length(dataR)/1000;
        shortenedTT=1/1000:0.001:length(shortenedFileR)/1000;
        startIdx = extraTime;
        stopIdx = length(movementData)-extraTime-1;

    else
        timePre=(timePreAndPost+startCycleToCut/freq)*Fs+numRamps*Fs/freq;
        timePost=(timePreAndPost+endCycleToCut/freq)*Fs+numRamps*Fs/freq;
        
        shortenedFileR=dataR(timePre:(length(dataR)-timePost),:);
        shortenedFileL=dataL(timePre:(length(dataL)-timePost),:);
        
        if(dataType==0)%For position
            shortenedMotion=movementData(timePre:(length(movementData)-timePost)); %don't need the -1 for position
        elseif(dataType==1)%For velocity
            shortenedMotion=movementData(timePre:(length(movementData)-timePost-1));
        end
        
        tt=1/1000:0.001:(length(dataR)/1000);
        shortenedTT=1/1000:0.001:length(shortenedFileR)/1000;
        startIdx = timePre;
    stopIdx = length(movementData)-timePost-1;
    end
    
    
end

%     %% desaccade routine, needs checking
%     derivVal = 0.05;
%     for i=1:6
%         if (i >3)
%             accel = diff(shortenedFileL(:,i-3));
%         else
%             accel = diff(shortenedFileR(:,i));
%         end
%         mask = abs(accel)>derivVal;
%         if (i > 3)
%             shortenedFileL(mask,i-3) = NaN;
%         else
%             shortenedFileR(mask,i) = NaN;
%         end
%     end

[xOffsetR,xMagR,xSinphaseR]=SingleFreqDFT(shortenedTT,shortenedFileR(:,1),freq);
[yOffsetR,yMagR,ySinphaseR]=SingleFreqDFT(shortenedTT,shortenedFileR(:,2),freq);
[zOffsetR,zMagR,zSinphaseR]=SingleFreqDFT(shortenedTT,shortenedFileR(:,3),freq);
[xOffsetL,xMagL,xSinphaseL]=SingleFreqDFT(shortenedTT,shortenedFileL(:,1),freq);
[yOffsetL,yMagL,ySinphaseL]=SingleFreqDFT(shortenedTT,shortenedFileL(:,2),freq);
[zOffsetL,zMagL,zSinphaseL]=SingleFreqDFT(shortenedTT,shortenedFileL(:,3),freq);

[motionOffset,motionMag,motionSinphase]=SingleFreqDFT(shortenedTT,shortenedMotion,freq);

%if you want to include offset - but we don't because only interested in the freq of the stimulus
% for i=1:10*Fs/freq
%     larpFitR(i)=larpOffsetR+larpMagR*sin(2*pi*freq*i/1000+larpSinphaseR);
%     ralpFitR(i)=ralpOffsetR+ralpMagR*sin(2*pi*freq*i/1000+ralpSinphaseR);
%     yawFitR(i)=yawOffsetR+yawMagR*sin(2*pi*freq*i/1000+yawSinphaseR);
%     larpFitL(i)=larpOffsetL+larpMagL*sin(2*pi*freq*i/1000+larpSinphaseL);
%     ralpFitL(i)=ralpOffsetL+ralpMagL*sin(2*pi*freq*i/1000+ralpSinphaseL);
%     yawFitL(i)=yawOffsetL+yawMagL*sin(2*pi*freq*i/1000+yawSinphaseL);
% end

for i=1:10*Fs/freq
    xFitR(i)=xMagR*sin(2*pi*freq*i/1000+xSinphaseR);
    yFitR(i)=yMagR*sin(2*pi*freq*i/1000+ySinphaseR);
    zFitR(i)=zMagR*sin(2*pi*freq*i/1000+zSinphaseR);
    xFitL(i)=xMagL*sin(2*pi*freq*i/1000+xSinphaseL);
    yFitL(i)=yMagL*sin(2*pi*freq*i/1000+ySinphaseL);
    zFitL(i)=zMagL*sin(2*pi*freq*i/1000+zSinphaseL);
    motionFit(i)=motionMag*sin(2*pi*freq*i/1000+motionSinphase);
end
%
jointSineR=[xFitR' yFitR' zFitR'];
jointSineL=[xFitL' yFitL' zFitL'];

[maxMotion,idxMaxMotion]=max(motionFit);

maxHeadMovementPointR=jointSineR(idxMaxMotion,:);
maxHeadMovementPointL=jointSineL(idxMaxMotion,:);

if(showPlot==1)
    
    figure;
    sgtitle(coilFile, 'Interpreter', 'none')
    set(gcf, 'Position', get(0, 'Screensize'));
    subplot(4,2,2)
    plot(shortenedFileR(:,1))
    hold on;
    plot(xFitR)
    title('Right Eye X')
    subplot(4,2,4)
    plot(shortenedFileR(:,2))
    hold on;
    plot(yFitR)
    title('Right Eye Y')
    subplot(4,2,6)
    plot(shortenedFileR(:,3))
    hold on;
    plot(zFitR)
    title('Right Eye Z')
    subplot(4,2,8)
    plot(shortenedMotion)
    hold on;
    plot(motionFit,'r')
    title('Motion/Sync Signal')
    
    subplot(4,2,1)
    plot(shortenedFileL(:,1))
    hold on;
    plot(xFitL)
    title('Left Eye X')
    subplot(4,2,3)
    plot(shortenedFileL(:,2))
    hold on;
    plot(yFitL)
    title('Left Eye Y')
    subplot(4,2,5)
    plot(shortenedFileL(:,3))
    hold on;
    plot(zFitL)
    title('Left Eye Z')
    subplot(4,2,7)
    plot(shortenedMotion)
    hold on;
    plot(motionFit,'r')
    title('Motion/Sync Signal')
    
end