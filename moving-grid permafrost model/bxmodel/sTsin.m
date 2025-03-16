function [dwf,ku] = sTsin(tu,sT)
%STSIN 构造地表温度拟合拟合
%   此处显示详细说明

%dwf=@(x)a+b*sin(2*pi/365*x+c)+k*x/365;
%dwf=@(x)a+b*sin(2*pi/365*x+c)+u;

z=find(isnan(sT));
sT(z)=[];
tu(z)=[];
ftype = fittype...
    ('a+b*sin(2*pi/365*x+c)+k*x/365', 'independent', 'x', 'dependent', 'y');
dwf = fit(tu',sT,ftype,...
    'startpoint',[mean(sT),max(abs(sT))-abs(mean(sT)),0,0]);

ku=dwf.k*100;

dwf=@(x)dwf.a+dwf.b*sin(2*pi/365*x+dwf.c);

end

