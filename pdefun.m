function[c,f,s]=pdefun(x,t,u,ux)%建立偏微分方程函数
global yc;
global yc2
%t1=mod(t,365);
[y0,c0,~]= tyc(yc,x,u);
%[y0,c0]= tyc2(yc,x,u);
s=0;
y=y0*24*3600;%*(u<=0)+1.278*24*3600*(u>0);
f=ux;
%c=(1.879e6*(u<=0)+2.357e6*(u>0))/y;
c=c0/y;
end
