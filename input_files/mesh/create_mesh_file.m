%this code runs the triangle command to create  mesh file  for the model
% and then it converts the mesh files to new model_converted files for
% reading in the DG C code
%create mesh
command = './triangle slope_model.poly  -pqanA';
status = system(command)

%convert mesh files format
fid = fopen('slope_model.1.ele');
firstline = str2num(fgetl(fid));
num = firstline(1);
Triangle = zeros(num+1,5);
Triangle(1,1:3) = firstline;
Triangle(1,4:5) =0;

for ij = 2:size(Triangle,1)
    Triangle(ij,:) = str2num(fgetl(fid));
end
fclose(fid);

save slope_model_converted.1.ele Triangle -ascii 

fid = fopen('slope_model.1.node');
firstline = str2num(fgetl(fid));
num = firstline(1);
Point = zeros(num+1,4);
Point(1,:) = firstline;

for ij = 2:size(Point,1)
    Point(ij,:) = str2num(fgetl(fid));
end
fclose(fid);

save slope_model_converted.1.node Point -ascii


fid = fopen('slope_model.1.neigh');
firstline = str2num(fgetl(fid));
num = firstline(1);
Neighbor = zeros(num+1,4);
Neighbor(1,1:2) = firstline;
Neighbor(1,3:4) =0;

for ij = 2:size(Neighbor,1)
    Neighbor(ij,:) = str2num(fgetl(fid));
end
fclose(fid);

save model_converted.1.neigh Neighbor -ascii

