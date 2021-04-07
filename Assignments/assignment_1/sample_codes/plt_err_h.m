%
clear
figure(1)
clf reset
axes('position',[0.15,0.15,0.75,0.75])
%
f =@(s) sin(s);
x=1;
h=2.^(-[0:0.5:10]);
fp1=(f(x+h)-f(x))./h;
fp2=(f(x+h)-f(x-h))./(2*h);
fpe=cos(x);
er1=abs(fp1-fpe);
er2=abs(fp2-fpe);
%
h1=loglog(h,er1,'b-o','linewidth',1.0,'markerfacecolor','r','markersize',8);
hold on
h2=loglog(h,er2,'r-s','linewidth',1.0,'markerfacecolor','y','markersize',8);
%
set(gca,'fontsize',16)
xlabel('h')
ylabel('e(h)')
title('Error of numerical differentiation')
axis([0.5e-3,1,0.5e-7,1e-0])
hL=legend([h1,h2],'First order method','Second order method',...
  'location','southeast');
set(hL,'fontsize',18)
%
%