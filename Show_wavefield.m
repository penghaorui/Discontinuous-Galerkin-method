function Show_Square(All_node,total_w,total_h,dx,Data)

[x,y]=meshgrid(0:dx:total_w,0:dx:total_h);

z=griddata(All_node(:,1),All_node(:,2),Data(:,1),x,y);
figure;hold on
imagesc(z);figure(gcf);
colorbar

axis equal
 axis([0 total_w/dx 0 total_h/dx]);

