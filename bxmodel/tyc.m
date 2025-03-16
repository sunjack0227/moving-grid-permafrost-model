function [y,c,w] = tyc(yc,x,u) %热参数




[i,~]=size(yc);
pw=1000;
l=334.54e3;

% yc(:,14)=yc(:,9).*(yc(:,8)-yc(:,7))+yc(:,11).*(yc(:,10)-yc(:,8))+yc(:,13).*(yc(:,12)-yc(:,10));

%yc(:,13)=-yc(:,14)./yc(:,10);

for k=1:i

if x>=yc(k,1)&&x<yc(k,2)
  try
      w=yc(k,15);   %超限冰含量
  catch
    w=yc(k,9)*(yc(k,8)-yc(k,7))+yc(k,11)*(yc(k,10)-yc(k,8))+...
       yc(k,13)*(yc(k,12)-yc(k,10));   %参与相变总含水含量
  end 
    %if x<1.4
           y=yc(k,4);
      c=yc(k,6)*1e6; 
        %myc21=[yc(k,8),yc(k,11)];
         %myc22=[yc(k,7),yc(k,10)];
      %if t>=77&t<267
         if u<yc(k,7)&&u<yc(k,8)
         y=yc(k,3);
        c=yc(k,5)*1e6;
        
         end
    
         if u>=yc(k,7)&&u<=yc(k,8)
         y=(yc(k,3)+yc(k,4))/2;
        c=(yc(k,5)+yc(k,6))/2*1e6+l*pw*yc(k,9);
         end
        
         if u>=yc(k,8)&&u<yc(k,10)
         y=(yc(k,3)+yc(k,4))/2;
        c=(yc(k,5)+yc(k,6))/2*1e6+l*pw*yc(k,11);
         end
        
         if u>=yc(k,10)&&u<yc(k,12)
         y=(yc(k,3)+yc(k,4))/2;
        c=(yc(k,5)+yc(k,6))/2*1e6+l*pw*yc(k,13);
         end
        
        if u>=yc(k,10)&&u>=yc(k,8)&&u>=yc(k,12)
         y=yc(k,4);
        c=yc(k,6)*1e6;
        end
      return;

end
end
end

