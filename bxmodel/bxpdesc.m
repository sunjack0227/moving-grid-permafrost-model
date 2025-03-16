function [sol,xx] = bxpdesc( x,t,sol0,x0,ys,pd,yc,pq,ku,varargin)
%PDESC �����ڳ� 
% ��ǰ�汾ȱ�ݣ��ٶ�һ���ڵ��ڳ�ȫ����������ӽ���һ����궳��������һ�����У�����һ����С���趨��С���1mm�������Ϊ1mm��ʣ���ڳ��������һ�������С������һ���ڻ��仯��Խ������㣬���ܻ����������һ�汾pdesc2�е��Խ����
%sol0 �ṩ��ʼ��0�궳������λ�ã���������������0-1��ͷ����ڳ�����������ڣ���Ĭ���Լ����0-1��Ķ�������Ϊ��׼����0-1���ǲ������ڳ�
%x0 yc��ԭʼ�ķֲ�
%fry=@(w)1.2*(w);
%   ry:�ڳ�ϵ�� ry=0.1236;
% fry=@(w)0.5*(w-0.11); %ϸ��
fry=@(w)0.6*(w-0.14+0.14); %ɰ��
% fry=@(w)0.7*(w-0.18+0.13); %����
% fry=@(w)0.6*(w-0.23); %���
clear bxpdes

% global yc
% global ys
t0=1; %ʱ����
n=ceil(max(t)/365);
n0=floor(min(t)/365);

 for i=n0+t0:t0:n

   if i==n0+t0
    tm=t(t>=(i-t0)*365 & t<=i*365);
      
    if isempty(sol0)
     x0=x;
    s0=bxpdes(x,tm,1,0,ys,pd,yc,pq,ku,varargin{:});
    ys=s0(end,:); %����ys
    [~,~,Hd,~] = pua( s0,tm,x,1,0);%�Լ����0-1��Ķ�������Ϊ��׼
       Hd0=Hd;
        xx=x';
        sol=s0;
        tm0=tm;
        continue;
      else
      s0=sol0(end-365+1:end,:);
      [~,~,Hd0,~] = pua(s0,[1:365],x,1,0); %��ȡsol0�ĳ�ʼ��0�궳������λ��
      ys=sol0(end,:);
      xx=x';
      tm0=[];
      sol=ys;
    end
    
  else
       tm=t(t>(i-t0)*365 & t<=i*365);
   
  end
   
    tm=[max(tm0),tm];
    % tm����0��ʼ�����Ѿ�����0λ��bxpdes�е�35�о��Բ�������&& i>1��������
    s0=bxpdes(x,tm,1,0,ys,pd,yc,pq,ku,varargin{:});
    ys=s0(end,:); %����ys
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
%         ry=fry(w0);    %��Ϊ����ry
%     end
% end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ah=ry*(Hd-Hd0);
     
    if ah<=0 || isnan(ah)||ry<0 %ֻ���ǳ����ڳ�
%     if  isnan(ah)||ry<0    %���ǳ��ڶ��ͺ��ڳ�
         xx=[xx x'];
        continue;
     end
    
    [~,i0]=min(abs(x-Hd0));
   % [~,i1]=min(abs(x-Hd));
   %ahh=ah/(i1-i0);
   %x(i0+1:i1)=x(i0+1:i1)-ahh;
k=1;

   while k
       if  x(i0+k)-ah-x(i0+k-1)<0.001 %��С���1mm 
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

