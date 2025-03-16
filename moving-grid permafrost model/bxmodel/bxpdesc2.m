function [sol,Twsol,wsol,xx,yc2] =bxpdesc2(x,t,a,t0,yc2,ys,pd,yc,pq,ku,varargin) %DECP主程序温度模块 %a循环次数,t0初始年分
 %[sol1,Twsol1,wsol1,xx1] =bxpdesc2(xx(:,1)',t,1,0,[],sol(1,:),0,yc0,0.04,0,[],[],[],dwf,t) 
%varargin=dwf,kun,xs,sT,tu,u;
%global xs %当前含水量（2），含冰量（3）、饱和含水量（4）、过剩冰量（5）
%global ys
%global yc2
%y=[x',yc2(:,14),yc2(:,15)];
%yc2=[];

try 
    dwf=varargin{1};
catch
    dwf=[];
end
try
    kun=varargin{2};
catch
    kun=[];
end

try 
    xs=varargin{3};
catch
    xs=[];
end

try
    sT=varargin{4};
catch
    sT=[];
end
try
    tu=varargin{5};
catch
    tu=[];
end
try
    u=varargin{6};
catch
    u=[];
end


m=0;
ii=1;
xx=x';
soiltpye=[];
xs=x'; 
clear tyc2


if isempty(ys)
ys=x;
 ys(:)=-1;
end

if isempty(yc2)
tyc2(x,1);
else
tyc2(x,0); 
end

Tw(yc2);

h=waitbar(0,'计算中');
aH=[];
sol=[];
wsol=xs(:,2)'; %各层未冻水含量
fsol=xs(:,3)';%各层冰含量
Twsol=yc2(:,15);
   for i=2:length(x)
       ax(i)=x(i)-x(i-1);
   end
tic

while ii<=a
tk=t+t0*365+max(t)*(ii-1);


if min(t)~=0 && ii>1
    tk=[t0*365+max(t)*(ii-1),t+t0*365+max(t)*(ii-1)];
   
end

yt=tk(2)-tk(1);
tk=[tk tk(end)+yt];
for i=2:length(tk)-1
ts=[tk(i-1) tk(i) tk(i+1)];

s=pdepe2(m,@bxpdefun,@bxpdeic,@bxpdebc,x,ts);

ys0=s(1,:);
ys=s(2,:);

Tw(yc2);

%形变模块
afsol=xs(:,3)-fsol(end,:)';% 冰含量变化
if isempty(soiltpye)
try
    soiltpye=yc2(:,16);
catch
    soiltpye=x';
%土壤质地：0砂土；1细砾；2碎石土（粉粘粒含量<15%》）） ；3粉土；4黏土；默认砂土
    soiltpye(:)=0;
end
end
axh=Deformation(x,afsol,soiltpye);
%aH=[aH,aHH];

x=x+axh;
xs(:,1)=x'; 

%%%%%%%%%%%%
tyc2(x,2);
%%%%%%%%%%%%%
%Twsol0=(xs(:,5)-yc2(:,15)).*ax'; %排出融水量m³/m²
Twsol=[Twsol,yc2(:,15)];

xx=[xx x'];
wsol=[wsol;xs(:,2)'];
fsol=[fsol;xs(:,3)'];


if i==2
    ss=[ys0;ys];
      pd0=1;
    
else
       
    ss=[ss;ys];
      
end

%disp(tk(i))
      strh=['计算进度: ',num2str(tk(i)),'/',num2str(max(tk))]; 
     waitbar(i/(length(tk)-1),h,strh);
toc
end

if ii==1
    sol=ss;
    ys=ss(end,:);
  else
   
    sol=[sol;ss(2:end,:)];
      ys=ss(end,:);
   
    
end

ii=ii+1;
 ii-1
 
toc
end
clear tyc2

    function Tw(yc)
%w超限冰含量或参与相变总含水
%fw当前温度下含冰量

%global xs %当前含水量（2）、含冰量（3）、饱和（实际）含水量（4）、过剩冰含量（5）分布
%global ys
%  yc(:,14)=yc(:,9).*(yc(:,8)-yc(:,7))+yc(:,11).*(yc(:,10)-yc(:,8))+...
%        yc(:,13).*(yc(:,12)-yc(:,10)); 

[yci,~]=size(yc);
yx=xs(:,1);

for j=1:length(yx)
for k=1:yci

if yx(j)>=yc(k,1)&&yx(j)<yc(k,2)
 ysu=ys(j);
  xs(j,4)=yc(k,14);
  
 if yc(j,15)>xs(j,4) %存在过剩冰
   
    if ysu<yc(k,7)&&ysu<yc(k,8)
        xs(j,3)=xs(j,5);
         xs(j,2)=0;
    end
   if ysu>=yc(k,7)&&ysu<=yc(k,8)
          xs(j,2)=yc(k,9)*(ysu-yc(k,7));
          xs(j,3)=xs(j,5)-xs(j,2);
   end     
        
   if ysu>=yc(k,8)&&ysu<yc(k,10)
         xs(j,2)=yc(k,9)*(yc(k,8)-yc(k,7))+yc(k,11)*(ysu-yc(k,8));
         xs(j,3)=xs(j,5)-xs(j,2);
   end    
        
 if ysu>=yc(k,10)&&ysu<yc(k,12)
  xs(j,2)=yc(k,9)*(yc(k,8)-yc(k,7))+yc(k,11)*(yc(k,10)-yc(k,8))+...
       yc(k,13)*(ysu-yc(k,10)); 
   xs(j,3)=xs(j,5)-xs(j,2);
 end    
        
  if ysu>=yc(k,10)&&ysu>=yc(k,8)&&ysu>=yc(k,12)
       xs(j,2)=xs(j,5);
        xs(j,3)=0;
   end
     
if xs(j,3)>yc(j,15)
      xs(j,3)=yc(j,15);
      xs(j,2)=0;
 end

 else
     
    %  xs(j,4)=yc(k,9)*(yc(k,8)-yc(k,7))+yc(k,11)*(yc(k,10)-yc(k,8))+...
%        yc(k,13)*(yc(k,12)-yc(k,10));   %参与相变总含水含量
         
   if ysu<yc(k,7)&&ysu<yc(k,8)
       xs(j,3)=xs(j,4);
        xs(j,2)=0;
   
    
   elseif ysu>=yc(k,7)&&ysu<=yc(k,8)
          xs(j,2)=yc(k,9)*(ysu-yc(k,7));
          xs(j,3)=xs(j,4)-xs(j,2);
        
        
   elseif ysu>=yc(k,8)&&ysu<yc(k,10)
         xs(j,2)=yc(k,9)*(yc(k,8)-yc(k,7))+yc(k,11)*(ysu-yc(k,8));
         xs(j,3)=xs(j,4)-xs(j,2);
       
        
 elseif ysu>=yc(k,10)&&ysu<yc(k,12)
  xs(j,2)=yc(k,9)*(yc(k,8)-yc(k,7))+yc(k,11)*(yc(k,10)-yc(k,8))+...
       yc(k,13)*(ysu-yc(k,10)); 
   xs(j,3)=xs(j,4)-xs(j,2);
       
        
   elseif ysu>=yc(k,10)&&ysu>=yc(k,8)&&ysu>=yc(k,12)
       xs(j,2)=xs(j,4);
        xs(j,3)=0;
   end
 end
        break;
end
  
end
end
%fw=xs(j,3);
end

function  tyc2(x,pd) %制作具体各层参数
% global yc
% global yc2
% global xs
persistent yc20;

[yci2,li]=size(yc);

yc2(:,1)=x';
yc2(1:end-1,2)=x(2:end)';
yc2(1,1)=yc(1,1);
yc2(end,2)=yc(end,2);

if isempty(yc20)
for j=1:length(x)
for k=1:yci2
if x(j)>=yc(k,1)&&x(j)<yc(k,2)
 
    yc20(j,7:13)=yc(k,7:13);
    
   
    break;
end
end
end
end
    
if pd==0
xs(:,5)=yc2(:,15);
 fixw(yc20,x);  
end

if pd==1  %不存在yc2，初始化制作yc2
for j=1:length(x)
for k=1:yci2
if x(j)>=yc(k,1)&&x(j)<yc(k,2)
 
    yc2(j,3:li)=yc(k,3:li);
         
   
    break;
end
  
end
end
fixw(yc20,x);
xs(:,5)=yc2(:,15);
end

if pd==2   %已存在yc2，修改每层过剩冰值，水热参数不变
    yc2(:,15)=xs(:,3);

fixw(yc20,x);
       
end

end



function fixw(yc20,x) %计算地下冰量的变化引起热参数的变化
%global yc2
for j=1:length(x)   
if yc2(j,15)>yc2(j,14) 
if yc2(j,12)>yc2(j,10)
    yc2(j,13)=yc20(j,13)+(yc2(j,15)-yc2(j,14))/(yc2(j,12)-yc2(j,10));
elseif yc2(j,10)>yc2(j,8)
    yc2(j,11)=yc20(j,11)+(yc2(j,15)-yc2(j,14))/(yc2(j,10)-yc2(j,8));
elseif yc2(j,8)>yc2(j,7)
    yc2(j,9)=yc20(j,9)+(yc2(j,15)-yc2(j,14))/(yc2(j,8)-yc2(j,7));
end
else
    yc2(j,7:13)=yc20(j,7:13);
       
end
end 
end

function axh = Deformation(x,afsol,soiltype)        %DEFORMATION 地表形变计算
%   此处显示详细说明
%   ry:融沉系数 ry=0.1236;

%fry=@(w)0.6*(w-0.14); %0砂土0
% fry=@(w)0.5*(w-0.11); %1细砾1
 %fry=@(w)0.4*(w-0.11); %2碎石土（粉粘粒含量<15%》））
% fry=@(w)0.7*(w-0.18); %3粉土3
 %fry=@(w)0.6*(w-0.23); %4黏土4
   
%fy=0.027; %冻胀率
   %global xs;
   %global yc2
  Totalw=x';
%   j1=find(xs(:,4)>=yc2(:,15));%不存在过剩冰网格（活动层）
%   j2=find(xs(:,4)<yc2(:,15));%存在过剩冰网格（多年冻土）
  maxefsol=yc2(:,15)-xs(:,4); 
  maxefsol(maxefsol<0)=0; %最大过剩冰融化量
%   Totalw(j1)=xs(j1,4);
%   Totalw(j2)=xs(j2,5);
  Totalw=xs(:,4);
  afsol1=zeros(length(x),1);
  afsol2=zeros(length(x),1);
   ry=fry(Totalw,soiltype);
  
   
   for xi=1:length(x)
       if xi>1
       ax(xi)=x(xi)-x(xi-1);
       end
       
       if afsol(xi)<0          %融沉情况
       if abs(afsol(xi))>maxefsol(xi)
           afsol2(xi)=-maxefsol(xi);
        else
           afsol2(xi)=afsol(xi);
       end                    
       afsol1(xi)=-(abs(afsol(xi))-abs(afsol2(xi)));
       else                   %冻胀情况
          afsol1(xi)=afsol(xi);
          afsol2(xi)=0;
       end
   end
%     ah(j1)=0.5*afsol(j1)./Totalw(j1).*ry(j1).*ax(j1)';%活动层形变
%     ah(j2)=afsol(j2)./Totalw(j2).*ry(j2).*ax(j2)';%多年冻土融沉
    ah=afsol2.*ax'+afsol1./Totalw.*ry.*ax';
   for xi=2:length(x)
       axh(xi)=sum(ah(1:xi)); %各移动网格量
   end
   
   aH=sum(ah) ;%总形变量
end

function ry=fry(w,soiltype)
ry=w;
ry(:)=0;
i1=find(soiltype==0);
i2=find(soiltype==1);
i3=find(soiltype==2);
i4=find(soiltype==3);
i5=find(soiltype==4);
ry(i1)=0.6*(w(i1)-0.14);
ry(i2)=0.5*(w(i2)-0.11);
ry(i3)=0.4*(w(i3)-0.11);
ry(i4)=0.7*(w(i4)-0.18);
ry(i5)=0.6*(w(i5)-0.23);
ry=ry*0.5;
ry(ry<0)=0;
end


%偏微分方程条件与解
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
   y=interp1(yx,ys,x,'PCHIP');
 %u0=-1+0.1*rand;
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

T=interp1(tu,sT,t,'PCHIP');
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


function pwsol=pw(Twsol,x,yc20) %过剩冰融化排出的水量
   for i=2:length(x)
       ax(i)=x(i)-x(i-1);
        
   end
  [~,m]=size(Twsol);
   for i=1:m
    pwsol0(:,i)=Twsol(:,i)-yc20(:,14);
   end
   pwsol0(pwsol0<0)=0;
   
    pwsol0(:,:)=pwsol0(:,1)- pwsol0(:,:);
  
   pwsol=pwsol0.*ax'; %排出融水量m³/m²
   pwsol=sum(pwsol);
end

