function [T,x0,y0,varargout] = gtif( fp,x,y) %提取数据
%fp tif文件路径或数据结构体,x经度，y纬度

try
[T0,~]=geotiffread(fp);
info=geotiffinfo(fp);
[m,n,~]=size(T0);
r=info.RefMatrix;
varargout{3}=r;
T=T0;


for i=1:m
    [~,y0(i)]=pix2map(r,i,1);
end
    for j=1:n
    [x0(j),~]=pix2map(r,1,j);
    end
    
catch
   try
    T0=fp.v;
   catch
       T0=[];
   end
    x0=fp.lon;
    y0=fp.lat;
    r=[];
%     T1=fp.soiltype;
%     x1=fp.soiltypelon;
%     y1=fp.soiltypelat;
end   

  if nargin==3 %提取特定点的值
   if length(x)==3
    x=dms2degrees(x);
    y=dms2degrees(y);
   end
    [~,i]=min(abs(x0-x));
    [~,j]=min(abs(y0-y));
   if ~isempty(T0)
    T=T0(j,i,:);
    T=squeeze(T);
   else
       T=[];
   end
    
%     [~,i1]=min(abs(x1-x));
%     [~,j1]=min(abs(y1-y));
%     varargout{1}=T1(j1,i1);
      varargout{1}=x0(i);
      varargout{2}=y0(j);
      
  end

    
    

end