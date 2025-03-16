function  bxrunpdes()
%区域模拟
%delete('已存在.mat')
addpath(genpath(pwd));
%Lcs='E:\Matlab\m\bxmodel\';
 t=0:365;

%pi=4;%1：仅预热；2：预热+kun年；3：kun年；4：未来
fbl=1; %空间分辨率（km）
L=load([pwd '\parameter\' num2str(fbl) 'km.mat']);
data=L.data;
% addpath('C:\Users\Administrator.SC-201809111626\Documents\MATLAB\m');
% addpath(genpath('C:\Users\Administrator.SC-201809111626\Documents\MATLAB\gis'));
tu=[];
ku=[];
%kun=[];
%xs=[];
sT=[];
%u=[];
%kun=50;
%pq=0.02;


%ku=2; %升温速率

D=data.v;X=data.lon;Y=data.lat;fD=data.f;
%af=data.af;
%soiltypev=data.soiltypev;
D1=D(:,:,1);
L=load([pwd '\parameter\' 'tu（modis 8天17年).mat']);
tu=L.tu;
L=load([pwd '\parameter\' 'x(30).mat']);

x=L.x;
clear L

 t1=0:365*ceil(max(tu/365));

 YC.yc1=load([pwd '\parameter\' 'yc（1沼泽）.mat']);
 YC.yc2=load([pwd '\parameter\' 'yc（2草原）.mat']);
 YC.yc3=load([pwd '\parameter\' 'yc（3草甸）.mat']);
 YC.yc4=load([pwd '\parameter\' 'yc（4荒漠）.mat']);


% outlj00=[pwd '\结果0\'];
% outlj0=[pwd '\历史升温结果2020\'];
outlj=[pwd '\result\'];

try 
    L=load([pwd '\已存在.mat']); %已存在
    yicunzai=L.yicunzai;
catch
    dirs=dir([pwd '\result']);
    dircell=struct2cell(dirs)';
    yicunzai=dircell(3:end,1);
      %yicunzai=[];
end
F=[1:4];
%F=[1 3 5 7 9 11 13];
%F=[2 4 6 8 10 12];
     yicunzai=[];
   kj=1:length(Y);ki=1:length(X);
     kij=parKij(1,ki,kj);
  for  k = 1:length( kij)
    [i,j]=parKij(2,kij,k);

    T1=D1(j,i);f=fD(j,i);

        %if isnan(T1)||f==0||f>2       %控制计算范围
        if isnan(T1)||~ismember(f,F)%||af(j,i)==0
            continue;
        end
try 
    L=load([pwd '\已存在.mat']); %已存在
    yicunzai=L.yicunzai;
catch
    dirs=dir([pwd '\result']);
    dircell=struct2cell(dirs)';
    yicunzai=dircell(3:end,1);
      %yicunzai=[];
end
       
       outfilename=[num2str(X(i),'%.4f') ',' num2str(Y(j),'%.4f') '.mat'];
       bool_name= ismember(outfilename, yicunzai); %判断是否已经算过了
       
              
       if bool_name==0
          yicunzai{end+1}=outfilename; 
          parsave('已存在.mat','yicunzai',yicunzai); 
%         switch pi
%           case 3
%               try
%                   sol00=load([outlj00 outfilename]');
%                   [~,~,Hd,~] = pua( sol00.sol.solc,t1,x,2,0);
%                   if isnan(Hd)
%                       continue
%                   end
%               catch
%                   msgbox(['警告:未找到',outlj00,outfilename])
%                   continue;
%               end
%           case 4
%               try
%               sol0=load([outlj0 outfilename]);
%               catch 
%                   msgbox(['警告:未找到',outlj0,outfilename])
%                   continue;
%               end
%         end
           
      disp(['正在计算',outfilename]);

%        sol=[];
%        parsave([outlj outfilename],'sol',sol);
%        continue  %测试总计算点数
       
% switch f
%     case 1
%    yc=yc1.yc;pq=yc1.pq;
%     case 2
%    yc=yc2.yc;pq=yc2.pq;
%     case 3
%     yc=yc3.yc;pq=yc3.pq;
%     case 4
%     yc=yc4.yc;pq=yc4.pq;
% end

      YCname=fieldnames(YC);
      PYC=YC.(YCname{f});

sT =gtif(data,X(i),Y(j));     
%MsT=mean(PYC.sT);
yc=PYC.yc;
%ku=PYC.ku;
pq=PYC.pq;
%Mu=PYC.u;

[dwf,ku] = sTsin(tu,sT);
%        try
%            pd0=PYC.pd0;   %确定下边界是第一或第二类条件
%        catch
%            pd0=0;
%        end
% 
%            
%     u=mean(sT)*Mu/MsT;               %修正地表温度平均误差
   
       
 
               

    

%    switch pi
%      case 1 % %kun年前预热
%      pd=1+pd0;
%      ys=x2;
%       T=mean(sT)+u-ku/2;
%       if T<0
%           ys(:)=-1;
%       else
%           ys(:)=1;
%       end  
%      sol=bxpdes(x2,t,500,0,ys,pd,yc,pq,ku,[],[],[],sT,tu,u); 
% %      ys=sol(end,:);
% %      ys(i1:i2)=pq*x2(i1:i2)+ys(i1)-pq*x2(i1);
% %      sol=ys;
%      sol=sol(end-365+1:end,:);
       
      % case 2 %kun年前预热+历史升温
      
     pd=3;
     ys=x;
      T=mean(dwf(0:365));
      if T<0
          ys(:)=-1;
      else
          ys(:)=1;
      end  
     sol=bxpdes(x,t,100,0,ys,pd,yc,pq,ku,dwf); 
     ys=sol(end,:);
   % ys(i1:i2)=pq*x2(i1:i2)+ys(i1)-pq*x2(i1);
     pd=4;
     [solc,xx] =bxpdesc(x,t1,sol,x,ys,pd,yc,pq,ku,dwf);
     sol=struct('solc',solc,'xx',xx);
          

%       case 3            %历史升温
%      pd=2+pd0;
%      ys=sol00.sol(end,:);
%      %ys=sol00.sol.solc(end,:);
%      [solc,xx] =bxpdesc(x2,t1,sol00.sol,x2,ys,pd,yc,pq,ku,[],[],[],sT,tu,u);
%      sol=struct('solc',solc,'xx',xx);
%      
%      
%      
%        case 4           %未来预估
% 
%       sT=sT+ku/100*10;% 从2020年算起，因此先按原来的升温条件升温10年
%            
%       ku=1.7;                   %RCP2.6
%      %ku=6.4;                   %RCP8.5
%      
%      pd=4+pd0;
%       s0= sol0.sol;
%       t2=0:365*80;
%       if isstruct(s0)
%           ys=s0.solc(end,:);
%           x=s0.xx(:,end)';
%           s0=s0.solc;
%       else
%       ys=s0(end,:);
%       x=x2;
%       end
%       [solc,xx] =bxpdesc(x,t2,s0,x2,ys,pd,yc,pq,ku,[],[],[],sT,tu,u);
%       sol=struct('solc',solc,'xx',xx);
%     end
     
     %sol=sol(end-460+1:end,:); %结果截取
    
     parsave([outlj outfilename],'sol',sol);
     disp(['已完成',outfilename]);
     %clear sol;clear sol00;clear sol0
     
    
       end
   end


end

function parsave(fname,x,y)
eval([x,'=y;',])
save(fname,x);
end

function varargout=parKij(csu,varargin)

if csu==1 %制作并行循环kij
    i=varargin{1};
    j=varargin{2};
    kij=[];
   for ki=1:length(i)
    for kj=1:length(j)
        %kij=[kij,str2double([num2str(i(ki)*1e9) , num2str(j(kj))])];
        kij=[kij,[i(ki);j(kj)]];
    end
   end
varargout{1}=kij;

elseif csu==2 %从kij还原ij
    kij=varargin{1};
    k=varargin{2};
    
    %varargout{1}=floor(kij(k)/1e10);
  % varargout{2}=mod(kij(k),1e10);
    varargout{1}=kij(1,k);
    varargout{2}=kij(2,k);
%     varargout{1}=i;
%     varargout{2}=j;
end



end