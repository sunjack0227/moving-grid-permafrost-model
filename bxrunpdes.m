function  bxrunpdes()
%区域模拟
%delete('已存在.mat')
addpath(genpath(pwd));

 t=0:365;



L=load([pwd '\parameter\' 'Inputdata.mat']);
data=L.data;

tu=[];
ku=[];

sT=[];





D=data.v;X=data.lon;Y=data.lat;fD=data.f;

D1=D(:,:,1);
L=load([pwd '\parameter\' 'tu.mat']);
tu=L.tu;
L=load([pwd '\parameter\' 'x.mat']);

x=L.x;
clear L

 t1=0:365*ceil(max(tu/365));

 YC.yc1=load([pwd '\parameter\' 'yc（1沼泽）.mat']);
 YC.yc2=load([pwd '\parameter\' 'yc（2草原）.mat']);
 YC.yc3=load([pwd '\parameter\' 'yc（3草甸）.mat']);
 YC.yc4=load([pwd '\parameter\' 'yc（4荒漠）.mat']);



outlj=[pwd '\result\'];

try 
    L=load([pwd '\已存在.mat']); %已存在
    yicunzai=L.yicunzai;
catch
    dirs=dir([pwd '\result']);
    dircell=struct2cell(dirs)';
    yicunzai=dircell(3:end,1);
     
end
F=[1:4];

     yicunzai=[];
   kj=1:length(Y);ki=1:length(X);
     kij=parKij(1,ki,kj);
  parfor  k = 1:length( kij)
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

           
      disp(['正在计算',outfilename]);

%        sol=[];
%        parsave([outlj outfilename],'sol',sol);
%        continue  %测试总计算点数
       

      YCname=fieldnames(YC);
      PYC=YC.(YCname{f});

sT =gtif(data,X(i),Y(j));     

yc=PYC.yc;

pq=PYC.pq;


[dwf,ku] = sTsin(tu,sT);

   
       
 
               

    


      
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

     pd=4;
     [solc,xx] =bxpdesc(x,t1,sol,x,ys,pd,yc,pq,ku,dwf);
     sol=struct('solc',solc,'xx',xx);
          

    
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