function sol =bxpdes(x,t,a,t0,ys,pd,yc,pq,ku,varargin) %地温计算模块
%a循环次数,t0初始年分

%x,t,a,t0,ys,pd,yc,pq,ku,dwf,kun,xs,sT,tu,u; 
persistent c;
if isempty(c)
 c=0;
end

try 
    dwf=varargin{1}; %地表温度函数
catch
    dwf=[];
end
try
    kun=varargin{2}; %假设历史模拟的年前预热（默认50年前）
catch
    kun=[];
end

try 
    xs=varargin{3}; %实测地温梯度中的深度
catch
    xs=[];
end

try
    sT=varargin{4}; %实测地表温度
catch
    sT=[];
end
try
    tu=varargin{5}; %实测地表温度时间
catch
    tu=[];
end
try
    u=varargin{6}; %地表温度修正参数
catch
    u=[];
end

%global xs
%global ys


%global yc2
%global yc
yc2=yc;



pdd0=1;
pdd1=0;
xs=x'; 
if isempty(ys)
ys=x;
 ys(:)=-1;
end
if (var(ys)-0)<1e-8
   % pdd0=0; %是否是随意赋初值
end
if pdd0==1
sol=ys;
end
m=0;
n=-4;
pdd=0;
i=1;
%ys
tic
while i<=a

%ts=t+t0*365+max(t)*(i-1);
ts=t+t0*365+ceil(max(t/365))*365*(i-1);

%ys=roundn(ys,-4);

if min(t)~=0  && i>1 % %t不从0开始,若直接影响最初值以及最初计算周期屏蔽&& i>1
    %ts=[t0*365+max(t)*(i-1),t+t0*365+max(t)*(i-1)];
    ts=[max(t+t0*365+ceil(max(t/365))*365*(i-2)),t+t0*365+ceil(max(t/365))*365*(i-1)];
    pdd1=1;
end

if pdd==0    
s=pdepe2(m,@bxpdefun,@bxpdeic,@bxpdebc,x,ts);
elseif pdd==1
  disp('降精度')
    s=pdepe3(m,@bxpdefun,@bxpdeic,@bxpdebc,x,ts);
end
% blm=['S',num2str(i)];
% eval([blm,'=','s',';']);
% lj=['C:\Users\Administrator\Desktop\m3\result\',blm,'.mat'];
% save(lj,blm);
[hang,l]=size(s);
l=length(ts);
if n>-1&pdd==1
    break;
elseif n>-1&pdd==0&pdd0~=0
     ys=sol(end,:);
     pdd=1;
     n=-4;
     continue;
  
end



if hang~=l%| max(max(abs(s)))>20;
    if pdd0==0
        if pdd==1
            return
        end
        pdd=1;
        continue;
    end
    %break;
      if n==-1
   ys(:)=ys+0.1*rand(1,length(ys));
   n=n+1
   continue;
      end
      ys=roundn(ys,n);
    
      n=n+1
      continue;
end

if i==1
   if pdd1==1 %t不从0开始,直接影响最初值 
     sol=s(2:end,:);
   else
     sol=s;
   end
  ys=s(end,:);
  pdd0=1;

else
    
   
    sol=[sol;s(2:end,:)];
    
     ys=s(end,:);
     
    
end
i=i+1;
pdd=0;
n=-4;
 c=c+1;
disp(c)
 toc
%totaltime=totaltime+toc;
%disp(['第' num2str(c) '次运行，总耗时' num2str(totaltime) '秒' ])
end

function[c,f,s]=bxpdefun(x,t,u,ux)%建立地温计算偏微分方程函数
%global yc;
%global yc2
%t1=mod(t,365);
if isempty(yc2)
    yc2=yc;
end
    
[y0,c0,~]= tyc(yc2,x,u);
%[y0,c0]= tyc2(yc,x,u);
s=0;
y=y0*24*3600;%*(u<=0)+1.278*24*3600*(u>0);
f=ux;
%c=(1.879e6*(u<=0)+2.357e6*(u>0))/y;
c=c0/y;
end

function[y]=bxpdeic(x)%建立偏微分方程的初始条件函数
   %global ys;
   %global xs;
   %yss=roundn(ys,-6);
 
   yx=xs(:,1);
  %y=interp1(yx,ys,x,'PCHIP');
  
   y=interp1(yx,ys,x);
end

function[pa,qa,pb,qb]=bxpdebc(xa,ua,xb,ub,t)   %建立偏微分方程的边界条件函数
% global dwf;
% global pd
% global u
% global tu %modis,8天10年
% global sT
% global ku
% global pq %地热梯度
% global kun %n年前初始条件
%pd=0; %设置运算模式，pd=0正常运算，%率定运算
%u=0;
if isempty(kun)
    kun=50;
end

if isempty(u)
    u=0;
end

if isempty(dwf)
%地表温度修正
T=sT+u;

%正弦拟合

z=find(isnan(T));
T(z)=[];
tu(z)=[];
ftype = fittype...
    ('a+b*sin(2*pi/365*x+c)', 'independent', 'x', 'dependent', 'y');
dwf = fit(tu',T,ftype,...
    'startpoint',[mean(T),max(abs(T))-abs(mean(T)),0]);
end


switch pd
    case 0  %实际地表温度驱动

T=interp1(tu,dwf,t,'PCHIP');
pa=ua-(T+u);qa=0;
pb=-pq;qb=1;


   
    case 1
    pa=ua-dwf(t)+ku/100*kun;qa=0; %kun年前
    pb=-pq;qb=1;
    
    case 2
    pa=ua-dwf(t)+ku/100*kun-ku/36500*t;qa=0; %历史升温
    pb=-pq;qb=1;
       
    
    case 3
    pa=ua-dwf(t);qa=0; %当前
    pb=-pq;qb=1;
      
    case 4
    pa=ua-dwf(t)-ku/36500*t;qa=0;   %未来升温
     pb=-pq;qb=1;
 end
end

end