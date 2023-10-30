function P = InCenter(x1,y1,x2,y2,x3,y3)
a = sqrt((x3 - x2)^2 + (y3 - y2)^2);
b = sqrt((x1 - x3)^2 + (y1 - y3)^2);
c = sqrt((x1 - x2)^2 + (y1 - y2)^2);
P(1,1) = (a*x1 + b*x2 + c*x3)/(a+b+c);
P(1,2) = (a*y1 + b*y2 + c*y3)/(a+b+c);
return