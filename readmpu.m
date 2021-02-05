%% Code to read in mpu file using directory and file name
function mpu = readmpu(directory, mpufile, monkeyVal)

%{
INPUT:
   directory - the directory of the mpu file you are opening
   mpufile - the name of the mpu .txt file you are opening
   monkeyVal - 0  if chinchilla on moog file. 1 if monkey on moog. this
      variable is becuase the MPUs are oriented differently on the Moog in the
      chinchilla vs monkey systems.

OUTPUT:
   mpu - data output with the proper +\- for the respective system (i.e.
   monkey or chinchilla)

The MPU File contains the following data in each column.
Column 1: FPGA Count
Column 2: Sync Signal
Column 3: Lateral (Left/Right - our +/- Y)    % For Monkey: Surge (our X)
Column 4: Surge (Front/Back - our +/- X)      % For Monkey: Lateral (our Y)
Column 5: Heave (Up/Down - our +/- Z)         % For Monkey: Heave (our Z)
Column 6: Pitch (rotate around interaural)    % For Monkey: Roll
Column 7: Roll (rotate around naso/occipital) % For Monkey: Pitch
Column 8: Yaw (rotate around earth vertical)  % For Monkey: Yaw
%}







f=fopen(strcat(directory,mpufile));
if f == -1
    mpu = [];
else
    mpu=fscanf(f,'%i',[8,inf]);
    mpu=mpu';
    fclose(f);
    
    if monkeyVal == 1 %If it is a monkey file
        dummyVarAccel = mpu(:,3);
        dummyVarRot = mpu(:,6);
        mpu(:,3) = mpu(:,4);
        mpu(:,4) = dummyVarAccel;
        mpu(:,6) = mpu(:,7);
        mpu(:,7) = dummyVarRot;
    else %if chinchilla file, just flip orientation of two columsn due to orientation of MPU
        mpu(:,3) = -1.*mpu(:,3); %due to orientation of the MPU and our definition of +/-
        mpu(:,6) = -1.*mpu(:,6); %due to orientation of the MPU and our definition of +/-
    end
end