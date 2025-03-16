function [Pu,Hu,Hd,Hdu] = pua( sol,t,x,csu,csu1,uc,tc) %温度波幅，年温度变化深度,活动层深度,t0末几年
%Pu:年平均地温；Hu温度年变化深度；Hd活动层厚度；Hdu冻土上限年均温度
%csu1==0,不画图
%csu,1：计算末年值，2：计算初年值
[xn,~]=size(x);
if ~exist('csu')
    csu=1;
end
if ~exist('csu1')
    csu1=1;
end

t0=1;

    if size(t,1)==2
   % t00=tc(2,1);
    t=t(1,:);
    end
if nargin==7
    if size(tc,1)==2
   % t00=tc(2,1);
    tc=tc(1,:);
    end
    txt2='实测';
    pd2=0;
    global hc;
    h=hc;
    n2=ceil(max(tc)/365);
if csu<2
    b2=[n2-1 n2];
else
    b2=[0,1];
end
j2=find(tc>=b2(1)*365 & tc<b2(2)*365);
s2=uc(j2,:);
tm2=tc(j2)-b2(1)*365;
tm2=tm2';
pu2=nanmean(s2);
k2=length(h);
for i2=1:k2
  %s1=abs(s(:,i)-pu(i));
    %c=trapz(tm,s1);
    %a(i)=c/(365*(b(2)-b(1)))*pi/2;
    a2(i2)=(max(s2(:,i2))-min(s2(:,i2)))/2;
    if pd2==0
        hd2=max(s2(:,i2));
        if hd2<=0
            pd2=1;
            hdi2=i2;
            hd2=h(i2);
        end
    end
end
i2=find(2*a2<0.1);
end


txt='最新结果';
if isstruct(sol)==1
    txt=sol.t;
    sol=sol.v;
end
pd1=0;
n=ceil(max(t)/365);
if csu<2
b=[n-t0 n];
if xn~=1
    x=x(:,end)';
end
else
 if xn~=1
    x=x(:,1)';
 end
b=[0 t0];
end
j=find(t>=b(1)*365 & t<b(2)*365);
s=sol(j,:);
tm=t(j)-b(1)*365;
tm=tm';
pu=nanmean(s);
k=length(x);
for i=1:k
  %s1=abs(s(:,i)-pu(i));
    %c=trapz(tm,s1);
    %a(i)=c/(365*(b(2)-b(1)))*pi/2;
    a(i)=(max(s(:,i))-min(s(:,i)))/2;
    if pd1==0
        hd=max(s(:,i));
        if hd<=0
            pd1=1;
            hdi=i;
            hd=x(i);
        end
    end
end

if csu1~=0
figure(1);
plot(pu,-x);
title('均温-深度变化');


if nargin==7
    hold on
    plot(pu2,-h);
     legend(txt,txt2)
end     
     
end

i=find(2*a<0.1);

if isempty(i)==0
if csu1~=0
figure(2);
plot(tm,s(:,i(1)));
tx1=[txt,' 序号' ':' num2str(i(1)) '；深度' ':' num2str(x(i(1))) '；温度波幅' ':' num2str(2*a(i(1)))];
end
Pu=nanmean(s(:,i(1)));
Hu=x(i(1));
else
    if csu1~=0
        figure(2);
    tx1=[txt, ' 达不到地温年变化深度' ];
    end
Pu=nan;
Hu=nan;
end

if nargin==7 && csu1~=0
    if isempty(i2)==0
    hold on
    plot(tm2,s2(:,i2(1)),'*');
tx1=[tx1,char(13,10)',txt2,' 序号' ':' num2str(i2(1)) '；深度' ':' num2str(h(i2(1))) '；温度波幅' ':' num2str(2*a2(i2(1)))];
   legend(txt,txt2)
    else 
  tx1=[tx1,char(13,10)',txt2, ' 达不到地温年变化深度'];      

    end
end

try
title(tx1);
catch
end

if pd1==1
if csu1~=0
figure(3);
plot(tm,s(:,hdi));
tx2=[txt,' 序号' ':' num2str(hdi) '；活动层厚度' ':' num2str(hd) ];
end
Hd=hd;
Hdu=nanmean(s(:,hdi));
else
    if csu1~=0
    figure(3);
        tx2=[txt, ' 达不到活动层深度' ];
    end
Hd=nan;
Hdu=nan;
end
if nargin==7 && csu1~=0
    if pd2==1
    hold on
    plot(tm2,s2(:,hdi2),'*');
    tx2=[tx2,char(13,10)',txt2,' 序号' ':' num2str(hdi2) '；活动层厚度' ':' num2str(hd2) ];
legend(txt,txt2)
    else
         tx2=[tx2,char(13,10)',txt2,' 达不到活动层深度' ];
    end
end

try
title(tx2);
catch
end


if csu1~=0
figure(4);
axy1=plot(pu+a,-x,'k',pu-a,-x,'k','LineWidth',2);
if nargin==7
    hold on
   axy2= plot(pu2+a2,-h,'r',pu2-a2,-h,'r','LineWidth',2);
    legend([axy1(1),axy2(1)],txt,txt2)
else

    if isempty(i)==0
line([pu(1)-a(1),pu(1)+a(1)],[-x(i(1)),-x(i(1))],'linestyle',':','LineWidth',2);
text(pu(1)-a(1),-x(i(1))-0.5,['温度年变化深度',num2str(x(i(1))),'m'],'FontSize',20)
    end
    
    if pd1==1
line([pu(1)-a(1),pu(1)+a(1)],[-hd,-hd],'linestyle',':','LineWidth',2);
text(pu(1)-a(1),-hd-0.5,['活动层厚度',num2str(hd),'m'],'FontSize',20)
    end
end
ylabel('h/m ')
xlabel('T/摄氏度')
title('地温包络线');
set(gca,'xaxislocation','top','FontSize',30);
%msgbox(tx)
end
end

