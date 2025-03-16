function [sol,xx] = bxpdesc( x,t,sol0,x0,ys,pd,yc,pq,ku,varargin)
%PDESC 考虑融沉 
% 当前版本缺陷：假定一年期的融沉全部发生在最接近上一年多年冻土上限那一土层中，若那一土层小于设定最小深度1mm，则调整为1mm，剩余融沉分配的下一层土层中。因此若一年期活动层变化跨越多个土层，可能会产生误差。在下一版本pdesc2中得以解决。
%sol0 提供初始第0年冻土上限位置，如果存在则计算在0-1年就发生融沉；如果不存在，则默认以计算的0-1年的冻土上限为基准，即0-1年是不发生融沉
%x0 yc中原始的分层
%fry=@(w)1.2*(w);
%   ry:融沉系数 ry=0.1236;
% fry=@(w)0.5*(w-0.11); %细砾
fry=@(w)0.6*(w-0.14+0.14); %砂土
% fry=@(w)0.7*(w-0.18+0.13); %粉土
% fry=@(w)0.6*(w-0.23); %黏土
clear bxpdes

% global yc
% global ys
t0=1; %时间间隔
n=ceil(max(t)/365);
n0=floor(min(t)/365);

 for i=n0+t0:t0:n

   if i==n0+t0
    tm=t(t>=(i-t0)*365 & t<=i*365);
      
    if isempty(sol0)
     x0=x;
    s0=bxpdes(x,tm,1,0,ys,pd,yc,pq,ku,varargin{:});
    ys=s0(end,:); %更新ys
    [~,~,Hd,~] = pua( s0,tm,x,1,0);%以计算的0-1年的冻土上限为基准
       Hd0=Hd;
        xx=x';
        sol=s0;
        tm0=tm;
        continue;
      else
      s0=sol0(end-365+1:end,:);
      [~,~,Hd0,~] = pua(s0,[1:365],x,1,0); %提取sol0的初始第0年冻土上限位置
      ys=sol0(end,:);
      xx=x';
      tm0=[];
      sol=ys;
    end
    
  else
       tm=t(t>(i-t0)*365 & t<=i*365);
   
  end
   
    tm=[max(tm0),tm];
    % tm不从0开始，但已经补了0位，bxpdes中第35行绝对不能屏蔽&& i>1！！！！
    s0=bxpdes(x,tm,1,0,ys,pd,yc,pq,ku,varargin{:});
    ys=s0(end,:); %更新ys
    tm0=tm;
    sol=[sol;s0(2:end,:)];
    [~,~,Hd,~] = pua( s0,tm,x,1,0);
    
    if ~isnan(Hd)

    [~,~,w]=tyc(yc,x0(x==Hd),0);
     ry=fry(w);
    else
        xx=[xx x'];
        continue;
    end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%
% if nargin==3 
%    if ry<fry(w0)&&Hd<3
%         ry=fry(w0);    %人为控制ry
%     end
% end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ah=ry*(Hd-Hd0);
     
    if ah<=0 || isnan(ah)||ry<0 %只考虑长期融沉
%     if  isnan(ah)||ry<0    %考虑长期冻胀和融沉
         xx=[xx x'];
        continue;
     end
    
    [~,i0]=min(abs(x-Hd0));
   % [~,i1]=min(abs(x-Hd));
   %ahh=ah/(i1-i0);
   %x(i0+1:i1)=x(i0+1:i1)-ahh;
k=1;

   while k
       if  x(i0+k)-ah-x(i0+k-1)<0.001 %最小间距1mm 
          ahh=x(i0+k)-x(i0+k-1)-0.001;
           x(i0+k:end)=x(i0+k:end)-ahh;
          ah=ah-ahh;
          if ah<=0
             break;
          end
       else
           x(i0+k:end)=x(i0+k:end)-ah; 
           break;
       end
       k=k+1;
   end
   Hd0=Hd-ah;
   xx=[xx x'];
   %disp(i);
 end
 clear bxpdes
end

