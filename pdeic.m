function[y]=pdeic(x)%����ƫ΢�ַ��̵ĳ�ʼ��������
   global ys;
   global xs;
   %yss=roundn(ys,-6);
   y=interp1(xs,ys,x,'PCHIP');
 u0=-1+0.1*rand;
end