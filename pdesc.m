function [sol,xx] = pdesc( x,t,sol0)
%PDESC �����ڳ�
%sol0 �ṩ��ʼ��0�궳������λ�ã���������������0-1��ͷ����ڳ�����������ڣ���Ĭ���Լ����0-1��Ķ�������Ϊ��׼����0-1���ǲ������ڳ�

%   ry:�ڳ�ϵ��

clear pdes

global yc
global ys
t0=1; %ʱ����
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
    % tm����0��ʼ�����Ѿ�����0λ��pdes�е�35�о��Բ�������&& i>1��������
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
%         ry=fry(w0);    %��Ϊ����ry
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
       if  x(i0+k)-ah-x(i0+k-1)<0.001; %��С���1mm
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
