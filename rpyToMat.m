%% returns the rotation matrix given rotations in roll, pitch, and yaw
function R = rpyToMat(r,p,y)
R = zeros(3,3);

sy = sin(y);
cy = cos(y);

sp = sin(p);
cp = cos(p);

sr = sin(r);
cr = cos(r);

Ryaw = [cy -sy 0;
      sy  cy 0;
      0    0 1];

Rpitch = [cp 0 sp;
           0 1  0;
          -sp 0 cp];
      
Rroll = [1 0 0;
         0 cr -sr;
         0 sr  cr];
     
R = Ryaw*Rpitch*Rroll;
end