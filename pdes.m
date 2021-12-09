function sol =pdes(x,t,a,t0) %a循环次数,t0初始年分
global xs
global ys
persistent c;
persistent totaltime;

if isempty(c)
 c=0;
 totaltime=0;
end

pdd0=1;
pdd1=0;
xs=x; 
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
s=pdepe2(m,@pdefun,@pdeic,@pdebc,x,ts);
elseif pdd==1
    '降精度'
    s=pdepe3(m,@pdefun,@pdeic,@pdebc,x,ts);
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
end