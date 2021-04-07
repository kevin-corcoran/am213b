%
clear
figure(1)
clf reset
axes('position',[0.15,0.15,0.75,0.75])
%
xa=[0:0.08:2.2]*pi;
y1=sin(xa)+xa/2;
y2=cos(xa);
y3=xa.*(2-xa/3)+1;
h1=plot(xa, y1, 'k-','linewidth',2.0);
hold on
h2=plot(xa, y2,'b-o','linewidth',1.0,'markerfacecolor','r','markersize',8);
h3=plot(xa, y3,'r-s','linewidth',1.0,'markerfacecolor','y','markersize',8);
%
set(gca,'fontsize',16)
xlabel('x')
ylabel('f(x)')
title('Comparison of several functions')
axis([0,7, -1.6, 6.4])
hL=legend([h1,h2,h3],'sin(x)+x/2','cos(x)','x*(2-x/3)+1',...
  'location','northwest');
set(hL,'fontsize',18)
%
%