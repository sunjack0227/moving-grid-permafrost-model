function[pa,qa,pb,qb]=pdebc(xa,ua,xb,ub,t)   %建立偏微分方程的边界条件函数
global dwf;
global pd
global u
global tu %modis,8天10年
global sT
global ku
global pq %地热梯度
%pd=0; %设置运算模式，pd=0正常运算，%率定运算
%u=0;

if isempty(dwf)
%地表温度修正
T=sT+u;

%正弦拟合


ftype = fittype...
    ('a+b*sin(2*pi/365*x+c)', 'independent', 'x', 'dependent', 'y');
dwf = fit(tu',T,ftype,...
    'startpoint',[mean(T),max(abs(T))-abs(mean(T)),0]);
end


switch pd
    case 0
t1=t-365*300;
if t1<0
    t1=0;
end
%T=interp1(tu,dwf,mod(t,365*10),'PCHIP');
%pa=ua-(T+u);qa=0;


%pa=ua-dwf(t);qa=0;
%pa=ua-dwf(t)+u/36500*t1-u;qa=0;
%pa=ua-dwf(t)-ku/36500*t;qa=0;   %未来升温
pa=ua-dwf(t)+ku/2-ku/36500*t;qa=0; %历史升温
%pa=ua-dwf(t)+ku/2;qa=0; %50年前预热
 %pa=ua-(u+12*sin(2*pi/365*t-1.931));qa=0; 

pb=-pq;qb=1;
%pb=ub-u;qb=0;
%pb=0;qb=1;
%pb=-0.0196;qb=1; %唐古拉
%pb=-0.03;qb=1;
%pb=-0.04;qb=1;
%pb=-0.045;qb=1; %西大滩
%pb=-0.05;qb=1;
%pb=-0.065;qb=1;
%pb=-0.076;qb=1;
%pb=-0.035;qb=1; %五道梁
%pb=-0.02;qb=1; %风火山,北麓河
%pb=ub+0.8393;qb=0;
%pb=-0.01;qb=1; %昆仑山垭口
%pb=-0.0056;qb=1;%两道河

   
    case 1
    pa=ua-dwf(t)+ku/2;qa=0; %50年前
    pb=-pq;qb=1;
    
    case 2
    pa=ua-dwf(t)+ku/2-ku/36500*t;qa=0; %历史升温
    pb=-pq;qb=1;
       
    
    case 3
    pa=ua-dwf(t);qa=0; %当前
    pb=-pq;qb=1;
      
    case 4
    pa=ua-dwf(t)-ku/36500*t;qa=0;   %未来升温
     pb=-pq;qb=1;
 end
end