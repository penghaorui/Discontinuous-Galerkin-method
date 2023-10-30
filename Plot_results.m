%this code plots the ouput files of the DG code
clear all;
close all;

%subsampling interval for plotting wavefield snapshots
delta = 2;

%total model size in x and z direction
total_x = 2000;
total_z = 2000; 

All_node = load('output_files/All_node.txt');

%plot wavefield snapshots
Vz = load('output_files/Vz500'); Show_wavefield(All_node,total_x,total_z,delta,Vz);   title('Vz','fontsize',12); 
Vz = load('output_files/Vz1000'); Show_wavefield(All_node,total_x,total_z,delta,Vz);   title('Vz','fontsize',12); 
Vx = load('output_files/Vx1800'); Show_wavefield(All_node,total_x,total_z,delta,Vx);   title('Vx','fontsize',12); 
Vz = load('output_files/Vz1800'); Show_wavefield(All_node,total_x,total_z,delta,Vz);   title('Vz','fontsize',12); 
Vx = load('output_files/Vx2000'); Show_wavefield(All_node,total_x,total_z,delta,Vx);   title('Vx','fontsize',12); 
Vz = load('output_files/Vz2000'); Show_wavefield(All_node,total_x,total_z,delta,Vz);   title('Vx','fontsize',12); 

Vx = load('output_files/Vx2200'); Show_wavefield(All_node,total_x,total_z,delta,Vx);   title('Vx','fontsize',12); 
Vz = load('output_files/Vz2200'); Show_wavefield(All_node,total_x,total_z,delta,Vz);   title('Vx','fontsize',12); 



print(1,'output_files/slope_Vx_snapshot.jpeg','-djpeg','-r300');
print(2,'output_files/slope_Vz_snapshot.jpeg','-djpeg','-r300');

%plot point record
figure;
load('output_files/record_x')
load('output_files/record_z')
plot(1:size(record_x,1),record_x(:,1)),hold on
plot(1:size(record_z,1),record_z(:,1))


load('output_files/absorb_fun.txt'); 
Show_wavefield(All_node,total_w,total_h,delta,absorb_fun);


