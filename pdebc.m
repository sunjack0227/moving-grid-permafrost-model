function[pa,qa,pb,qb]=pdebc(xa,ua,xb,ub,t)   %����ƫ΢�ַ��̵ı߽���������
global dwf;
global pd
global u
global tu %modis,8��10��
global sT
global ku
global pq %�����ݶ�
%pd=0; %��������ģʽ��pd=0�������㣬%�ʶ�����
%u=0;

if isempty(dwf)
%�ر��¶�����
T=sT+u;

%�������


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
%pa=ua-dwf(t)-ku/36500*t;qa=0;   %δ������
pa=ua-dwf(t)+ku/2-ku/36500*t;qa=0; %��ʷ����
%pa=ua-dwf(t)+ku/2;qa=0; %50��ǰԤ��
 %pa=ua-(u+12*sin(2*pi/365*t-1.931));qa=0; 

pb=-pq;qb=1;
%pb=ub-u;qb=0;
%pb=0;qb=1;
%pb=-0.0196;qb=1; %�ƹ���
%pb=-0.03;qb=1;
%pb=-0.04;qb=1;
%pb=-0.045;qb=1; %����̲
%pb=-0.05;qb=1;
%pb=-0.065;qb=1;
%pb=-0.076;qb=1;
%pb=-0.035;qb=1; %�����
%pb=-0.02;qb=1; %���ɽ,��´��
%pb=ub+0.8393;qb=0;
%pb=-0.01;qb=1; %����ɽ���
%pb=-0.0056;qb=1;%������

   
    case 1
    pa=ua-dwf(t)+ku/2;qa=0; %50��ǰ
    pb=-pq;qb=1;
    
    case 2
    pa=ua-dwf(t)+ku/2-ku/36500*t;qa=0; %��ʷ����
    pb=-pq;qb=1;
       
    
    case 3
    pa=ua-dwf(t);qa=0; %��ǰ
    pb=-pq;qb=1;
      
    case 4
    pa=ua-dwf(t)-ku/36500*t;qa=0;   %δ������
     pb=-pq;qb=1;
 end
end