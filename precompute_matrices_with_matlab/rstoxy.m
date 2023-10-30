function [x,y] = rstoxy(tri,point,t,r,s)

% function [x,y] = rstoab(r,s)
% Purpose : Transfer from (r,s) -> (x,y) coordinates in triangle
Np = length(r);
P = zeros(3,2);
x = zeros(tri,Np);
y = zeros(tri,Np);

for num = 1: tri
    for i = 1:3
        P(i,1) = point(t(num,i),1);
        P(i,2) = point(t(num,i),2);
    end
    x(num,:) = -(r(:,1) + s(:,1))/2*P(1,1) + (r(:,1) + 1)/2*P(2,1) + ((s(:,1) + 1)/2*P(3,1));
    y(num,:) = -(r(:,1) + s(:,1))/2*P(1,2) + (r(:,1) + 1)/2*P(2,2) + ((s(:,1) + 1)/2*P(3,2));
end

return;
