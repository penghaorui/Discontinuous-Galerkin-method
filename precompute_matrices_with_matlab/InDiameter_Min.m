function [Lmin] = InDiameter_Min(Length)
 
 radius = zeros(size(Length,1),1);
 for i = 1:size(Length,1)
     a = Length(i,1);
     b = Length(i,2);
     c = Length(i,3);
     s = 0.5*(a+b+c);
     S = sqrt(s*(s-a)*(s-b)*(s-c));  %area of triangle
     radius(i,1) = 2*S/s;  %diameter of inscribed circle
 end
 Lmin = min(radius); %find the smallest

return
