function[y]=pdeic(x)%建立偏微分方程的初始条件函数
   global ys;
   global xs;
   %yss=roundn(ys,-6);
   y=interp1(xs,ys,x,'PCHIP');
 u0=-1+0.1*rand;
end