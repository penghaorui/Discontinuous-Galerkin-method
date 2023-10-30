%this code precomputes the matrices used in DG theory, the functions used
%are from the given codes from the book: Nodal Discontinuous Galerkin Methods: Algorithms, Analysis, and Applications

clear all;

N = 2;
Np = (N+1)*(N+2)*0.5;
Fekete = Fekete_Degree2;

FeketeNew(:,1) = 2*Fekete(:,1) - 1;  %transform triangle
FeketeNew(:,2) = 2*Fekete(:,2) - 1;


V = Vandermonde2D(N,FeketeNew(:,1),FeketeNew(:,2));
inverV = inv(V);
[Vr,Vs] = GradVandermonde2D(N,FeketeNew(:,1),FeketeNew(:,2));

EdgePoint1 = [FeketeNew(6,:); FeketeNew(4,:); FeketeNew(1,:);];
EdgePoint2 = [FeketeNew(1,:); FeketeNew(2,:); FeketeNew(3,:)];
EdgePoint3 = [FeketeNew(6,:); FeketeNew(5,:); FeketeNew(3,:);];

Vb = Vandermonde1D(N,EdgePoint1(:,1));

inverVb = inv(Vb);

Vib1 = [V(6,:); V(4,:); V(1,:);];
Vib2 = [V(1,:); V(2,:); V(3,:);];
Vib3 = [V(6,:); V(5,:); V(3,:);];

S1 = V*(Vib1')*(inverVb')/Vb;
S2 = V*(Vib2')*(inverVb')/Vb;
S3 = V*(Vib3')*(inverVb')/Vb;

%save file for later computation
save input_files/V.txt V -ascii
save input_files/inverV.txt inverV -ascii
save input_files/Vb.txt Vb -ascii
save input_files/inverVb.txt inverVb -ascii
save input_files/Vr.txt Vr -ascii
save input_files/Vs.txt Vs -ascii
save input_files/FeketeNew.txt FeketeNew -ascii
save input_files/S1.txt S1 -ascii
save input_files/S2.txt S2 -ascii
save input_files/S3.txt S3 -ascii

