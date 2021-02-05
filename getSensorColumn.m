function [sensorColumn]=getSensorColumn(direction)

%{
INPUT:
direction - string variable indicating the direction of movement

OUTPUT:
sensorColumn - the column number for mpu data

%}



if (strcmp(direction,'Lateral'))
    sensorColumn=3;
elseif (strcmp(direction,'Surge'))
    sensorColumn=4;
elseif strcmp(direction,'Heave')
    sensorColumn=5;
elseif strcmp(direction,'Roll')
    sensorColumn=7;
elseif strcmp(direction,'Pitch')
    sensorColumn=6;
elseif strcmp(direction,'Yaw') 
    sensorColumn=8;
else %For Oblique Angle
    
end


end