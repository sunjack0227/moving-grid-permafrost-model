function [sol,xx] = pdesc( x,t,sol0)
%PDESC 考虑融沉
%sol0 提供初始第0年冻土上限位置，如果存在则计算在0-1年就发生融沉；如果不存在，则默认以计算的0-1年的冻土上限为基准，即0-1年是不发生融沉

%   ry:融沉系数

clear pdes

global yc
global ys
t0=1; %时间间隔
n=ceil(max(t)/365);

 for i=t0:t0:n

   if i==t0
    tm=t(t>=(i-t0)*365 & t<=i*365);
      if nargin<3
   
    
    s0=pdes(x,tm,1,0);
    [~,~,Hd,~] = pua( s0,tm,x,1,0);
       Hd0=Hd;
        xx=x';
        sol=s0;
        tm0=tm;
        continue;
      else
      s0=sol0(end-365:end,:);
      [~,~,Hd0,~] = pua( s0,tm,x,1,0);
      ys=sol0(end,:);
      xx=x';
      tm0=[];
      sol=ys;
      end
  else
       tm=t(t>(i-t0)*365 & t<=i*365);
   
  end
   
    tm=[max(tm0),tm];
    % tm不从0开始，但已经补了0位，pdes中第35行绝对不能屏蔽&& i>1！！！！
    s0=pdes(x,tm,1,0);
    tm0=tm;
    sol=[sol;s0(2:end,:)];
    [~,~,Hd,~] = pua( s0,tm,x,1,0);
    
    if ~isnan(Hd)
    [~,~,w]=tyc(yc,Hd,0);
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
     
    if ah<=0 || isnan(ah)||ry<0
        xx=[xx x'];
       continue;
    end
    
    [~,i0]=min(abs(x-Hd0));
   % [~,i1]=min(abs(x-Hd));
   %ahh=ah/(i1-i0);
   %x(i0+1:i1)=x(i0+1:i1)-ahh;
k=1;

   while k
       if  x(i0+k)-ah-x(i0+k-1)<0.001; %最小间距1mm
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
 end
 clear pdes
end

