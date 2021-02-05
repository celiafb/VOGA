function [temp]=loadGains(directory, coilfile)

%{
Function to load gains that are saved in a coil file

INPUT:
  directory - path to the directory where the coil file is located
  coilfile - name of the '.coil' file

OUTPUT:
  [temp] - gains for each X,Y,Z of each of the 4 coils


%}


f=fopen(strcat(directory,coilfile), 'rb');
if (f == -1)
    temp=-1;
else
    fseek(f,-96,'eof');
    temp = fread(f,[3,inf],'double');
end
