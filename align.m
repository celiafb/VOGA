function [rotR,rotL,mpuAligned,coilsAligned,isCoilsAligned] = align(mpu, coils, REF, frame,handles,GAINSR,GAINSL,ZEROS_R,ZEROS_L)
%{
Align and analyze coil data
1. align the two files with the sync signal using dutycycle(X,Fs).
2. convert raw coil data to rotation vectors using analyzeCoilData() which calls gimbalproc().

INPUT:
  mpu - the mpu data containing sensor data and the saved sync signal
  coils - the eye coil data also containing the sync signal 
  REF - reference point used in gimbol proc.
  frame - 0 for head frame of reference, 1 for eye ball frame of reference
  handles - if using with the Data_GUI, the gains are contained in here
  GAINSR - if not using with the Data_GUI, send gains for Right eye (R eye coils connected to Ch1 and Ch2
     so gains should be: [Ch1X Ch1Y Ch1Z Ch2X Ch2Y Ch2Z]
  GAINSL - if not using with the Data_GUI, send gains for Left eye (L eye coils connected to Ch3 and Ch4
     so gains should be: [Ch3X Ch3Y Ch3Z Ch4X Ch4Y Ch4Z]
  ZEROS_R - if not using with the Data_GUI, send zeroes in for Right eye
  ZEROS_L - if not using with the Data_GUI, send zeroes in for Left eye

OUTPUT:
  rotR - the rotation vector for right eye movements in head frame of reference, after subtracting
    offsets, and using frame of reference, gains, and reference point as provided to gimbalproc();
  rotL - the rotation vector for left eye movements in head frame of reference, after subtracting
    offsets, and using frame of reference, gains, and reference point as provided to gimbalproc();
  mpuAligned - the mpu data with any access data points pre or post
    alignment deleted
  coilsAligned - the raw coil data with any access data points pre or post
    alignment deleted

returns [d, initcross,finalcross],
where d=data of dutycycle of each pulse,
init cross is vector elements correspond to mid-crossing of initial transition of each pulse
finalcross is vector elements that correspond to final-crossing
 takes parameter frame - right eye, left eye, or head frame to determine
 which rotation vector to return: 0=head, 1=eye
%}

if (~isempty(handles))
    GAINSR = [handles.ch0 handles.ch1];
    GAINSL = [handles.ch2 handles.ch3];
    ZEROS_R = handles.ZEROSR;
    ZEROS_L = handles.ZEROSL;
end

[mpuSYNC,initCrossMPU, finalCrossMPU]=dutycycle(mpu(30:length(mpu),2),1000);
[coilsSYNC,initCrossCoils, finalCrossCoils]=dutycycle(coils(30:length(coils),3),1000);

%pick the larger of the two starting values...then find that point where
%the larger value happened in both data sets and choose that as starting
%point.  You need to round data - chose the thousandth's place - cuz
%otherwise you won't find the matching number

isMpuAligned = true;
isCoilsAligned = true;
mpuSYNC=round(mpuSYNC*200)/200;
coilsSYNC=round(coilsSYNC*200)/200;

%if the difference in the first sync pulse  durations in the files is >0.2
if (abs(mpuSYNC(1)-coilsSYNC(1))>=0.2)
    %this means a wrap around occured before recording in one channel
    %so we want to start adding to the channel that wrapped around now,
    %and the other channel later once we see the wrap around
    if (mpuSYNC(1)<coilsSYNC(1))
        mpuSYNC=mpuSYNC+0.2;
        offsetAddMPU=0.2;
        offsetAddCoils=0;
    else
        coilsSYNC=coilsSYNC+0.2;
        offsetAddCoils=0.2;
        offsetAddMPU=0;
    end
else
    offsetAddMPU=0;
    offsetAddCoils=0;
end


for i=2:length(mpuSYNC)
    if (mpuSYNC(i)<(mpuSYNC(i-1)-offsetAddMPU))
        offsetAddMPU=offsetAddMPU+0.2;
    end
    mpuSYNC(i)=mpuSYNC(i)+offsetAddMPU;
end

for j=2:length(coilsSYNC)
    if (coilsSYNC(j)<(coilsSYNC(j-1)-offsetAddCoils))
        offsetAddCoils=offsetAddCoils+0.2;
    end
    coilsSYNC(j)=coilsSYNC(j)+offsetAddCoils;
end


if (mpuSYNC(1)<coilsSYNC(1))
    startingDutyValue=coilsSYNC(1);
    idxCOILS=round(initCrossCoils(1)*1000);
    idxMPU=round(initCrossMPU(find(mpuSYNC==startingDutyValue,1,'first'))*1000);
else
    startingDutyValue=mpuSYNC(1);
    idxCOILS=round(initCrossCoils(find(coilsSYNC==startingDutyValue,1,'first'))*1000);
    idxMPU=round(initCrossMPU(1)*1000);
end

if idxMPU == 0
    idxMPU = 1;
end
if idxCOILS == 0
    idxCOILS = 1;
end


%now that we have the starting index where the two files line up - lets pick an ending index
if ((length(mpu)-idxMPU)>(length(coils)-idxCOILS))
    totalLength=length(coils)-idxCOILS;
else
    totalLength=length(mpu)-idxMPU;
end


mpuAligned=mpu(idxMPU:(totalLength+idxMPU),:);
if isempty(mpuAligned)
    isMpuAligned = false;
end
if isMpuAligned
    %get rid of MPU offsets
    for i=3:size(mpuAligned,2)
        mpuAligned(:,i)=mpuAligned(:,i)-mpuAligned(REF,i);
    end

    
end

coilsAligned=coils(idxCOILS:(totalLength+idxCOILS),:);
if isempty(coilsAligned)
    isCoilsAligned = false;
end
if isCoilsAligned
    [rotRhead,rotLhead,rotRReye,rotLLeye,rotRref,rotLref] = analyzeCoilData(coilsAligned, REF, handles,GAINSR,GAINSL,ZEROS_R,ZEROS_L);
    switch frame
        case 0 %For head frame of reference
            rotR = rotRref; %CHANGE MADE BY KH 2/11/16 - should be rotRref and not rotRhead
            rotL = rotLref; %CHANGE MADE BY KH 2/11/16 - should be rotLref and not rotLhead
        case 1 %For eye socket frame of reference
            rotR = rotRReye;
            rotL = rotLLeye;
    end
else
    rotR =[];
    rotL =[];
end
