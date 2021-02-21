
function [deSaccAngVelR, deSaccAngVelL]=desacc_prosthCanal(angVelR,angVelL)


accelSaccVal = 0.8;
samprate = 1000;

filteredRawAngVel(:,:,1)=angVelR;
filteredRawAngVel(:,:,2)=angVelL;

for k=1:2
 for j = 1:3
        clear accel findSac tempValStart tempValStop
        accel=diff(filteredRawAngVel(:,j,k));
        
        %if acceleration is greater than 0.3 degrees/ms^2, then mark that point as saccade
        for i=1:(length(filteredRawAngVel(:,j,k))-1)
            if abs(accel(i))>accelSaccVal %use 0.3 for 50dps
                findSac(i)=1;
            else
                findSac(i)=0;
            end
        end
        
        %plot to make sure everything looks ok
%         figure; hold on
%         plot(accel)
%         plot(findSac,'r')
        
        %find the index value for start and stop of each saccade
        if nnz(findSac)>0
            count=1;
            for i=1:(length(findSac(:)))
                %if first value is a saccade - set start point to index value 1
                if i==1
                    if findSac(i)==1
                        tempValStart(count)=1;
                    end
                    
                    %if increase - set start point to index value
                elseif (findSac(i)-findSac(i-1))==1
                    tempValStart(count)=i-1;
                    
                    %if decrease - set stop point to index value
                elseif (findSac(i)-findSac(i-1))==-1
                    tempValStop(count)=i+1;
                    count=count+1;
                end
                
                %if still in saccade when file ends - store stop val as last index value
                if i==length(findSac(:))
                    if findSac(i)==1
                        tempValStop(count)=i;
                    end
                end
            end
            

            
            %connect saccade index markers if it the gap between sacaddes is less than 10 ms
            count=1;
            if length(tempValStart)>1
                for i=2:length(tempValStart)
                    %if time between saccades is less than 15 ms, most likely all one
                    %saccade
                    if tempValStart(i)-tempValStop(count)<300
                        tempValStop(count)=tempValStop(i);
                        tempValStop(i)=NaN;
                        tempValStart(i)=NaN;
                        
                    else
                        count=i;
                    end
                    %if total detected saccade is less than 10ms, most likely not a saccade
                    %if tempValStop(i-1)-tempValStart(i-1)<10
                    %    tempValStop(i-1)=NaN;
                    %    tempValStart(i-1)=NaN;
                    %end
                end
            end
            
            %if starting point is in middle of saccade, set value to be
            count=1;
            for i=1:length(tempValStart)
                if tempValStart(i)>=0
                    range(count,1,j+3*(k-1))=tempValStart(i);
                    range(count,2,j+3*(k-1))=tempValStop(i);
                    count=count+1;
                end
            end
            
            %get rid of NaNs
            range(isnan(range(:,1,j+3*(k-1))),:) = [];
            
            %reset range values to be 20 ms before and after detected sacade index
            for i=1:nnz(range(:,1,j+3*(k-1)))
                if range(i,1,j+3*(k-1))>21
                    range(i,1,j+3*(k-1))=range(i,1,j+3*(k-1))-20;
                end
                
                if range(i,2,j+3*(k-1))<length(angVelL)-21
                    range(i,2,j+3*(k-1))=range(i,2,j+3*(k-1))+20;
                end
            end
        end
 end
    
end
 

if((sum(findSac))>6)
    %% Desaccade
    %create mask for desaccading
    desaccMask=ones(length(angVelL(:,1))-1,6);
    
    %range(a,b,c) --> a=# of saccade, b=2 (for start and end values of saccades)c=1 for Rz, 2 for Rl, 3 for Rr,
    for i=1:6
        for j=1:nnz(range(:,1,i))
            desaccMask(range(j,1,i):range(j,2,i),i)=0;
        end
    end
    
    maskedData=zeros(length(angVelL(:,1))-1,6);
    
    %plot the masked Data - showing flat lines where saccade was found
    try
        for i=0:5
            maskedData(:,i+1)=desaccMask(:,i+1).*filteredRawAngVel(1:end-1,mod(i,3)+1,round((i+1)/7)+1);
            %   figure(i+1)
            %   hold on
            %   plot(filteredRawAngVel(:,i+1),'b')
            %   plot(maskedData(:,i+1),'g')
        end
    catch
        return
    end
    
    splinefitData=filteredRawAngVel(:,:,:);
    
    %use spline function to smooth over saccades
    for i=1:6
        for j=1:nnz(range(:,1,i))
            k1=range(j,1,i);
            k2=range(j,2,i);
            %splinefitData(k1,i+1)=median(smoothdata((k1-15):(k1+5),i+1));
            %splinefitData(k2,i+1)=median(smoothdata((k2-5):(k2+15),i+1));
            %slope1=(median(smoothdata((k1-5):(k1+5),i+1))-median(smoothdata((k1-15):(k1-5),i+1)))/(10/samprate);
            %slope2=(median(smoothdata((k2+5):(k2+15),i+1))-median(smoothdata((k2-5):(k2+5),i+1)))/(10/samprate);
            
            slope1=maskedData(k1,i)/(10/samprate);
            slope2=maskedData(k2,i)/(10/samprate);
            
            splinefitData(k1:k2,mod((i-1),3)+1,(round((i)/7))+1)=spline([k1 k2],[slope1 splinefitData(k1,mod((i-1),3)+1,(round((i)/7))+1) splinefitData(k2,mod((i-1),3)+1,(round((i+1)/7))+1) slope2], k1:k2);
        end
    end
    
    %plot the splined data in black
    % for i=1:6
    %   figure(i)
    %   hold on
    %   plot(filteredRawAngVel(:,mod((i-1),3)+1,(round((i)/7))+1),'c','LineWidth',3)
    %   plot(splinefitData(:,mod((i-1),3)+1,(round((i)/7))+1),'k','LineWidth',3)
    %
    % end
    
    deSaccAngVelR=splinefitData(:,:,1);
    deSaccAngVelL=splinefitData(:,:,2);
else
    deSaccAngVelR=angVelR;
    deSaccAngVelL=angVelL;
end
%save(strcat(Pathname,filename,'_Desaccaded.mat'),'splinefitData');

