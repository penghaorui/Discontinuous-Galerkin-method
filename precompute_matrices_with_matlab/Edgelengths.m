function [L,Average] = Edgelengths(point,t)
%compute the lengths of the three edges 
L =zeros(size(t,1),3);
%three vertices
P1=zeros(1,2);
P2=zeros(1,2);
P3=zeros(1,2);

for i = 1:size(t,1)
    P1(1,1) = point(t(i,1),1); P1(1,2) = point(t(i,1),2);
    P2(1,1) = point(t(i,2),1); P2(1,2) = point(t(i,2),2);
    P3(1,1) = point(t(i,3),1); P3(1,2) = point(t(i,3),2);
    
    L(i,1)=sqrt((P1(1,1)-P2(1,1))^2 + (P1(1,2)-P2(1,2))^2);
    L(i,2)=sqrt((P2(1,1)-P3(1,1))^2 + (P2(1,2)-P3(1,2))^2);
    L(i,3)=sqrt((P1(1,1)-P3(1,1))^2 + (P1(1,2)-P3(1,2))^2);
end

Average = sum(L,1)/size(t,1);
return


