%this code does the mesh partition using metis

%read triangle mesh files
Triangle = load('slope_model_converted.1.ele');
triangle = [Triangle(2:end,2) Triangle(2:end,3) Triangle(2:end,4) ];

%step 1 write 'partition.mesh' file
fid=fopen('slope_partition.mesh','w+');
tic
fprintf(fid,'%d\n',size(triangle,1));
for i =1:size(triangle,1)
    fprintf(fid,'%d %d %d\n',triangle(i,:));
end
fclose(fid);
toc

Process_Num = 2;
% step 2: use metis to partition the mesh to the number of process desired 
command = "mpmetis slope_partition.mesh ";
command = strcat(command,num2str(Process_Num));
status = system(command);

