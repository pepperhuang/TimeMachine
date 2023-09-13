function A = GetMeanAngle(angles) 
 n=length(angles);
 x=repmat(NaN,1,n); y=repmat(NaN,1,n);
 for i=1:n
  x(i)=sin(angles(i)*pi/180);
  y(i)=cos(angles(i)*pi/180);
 end;
 meanX=mean(x);
 meanY=mean(y);
 A=atan(meanX/meanY)*180/pi;

 if ((sign(meanX)==1) & (sign(meanY)==1)) 
     A=atan(meanX/meanY)*180/pi; 
 end
 if ((sign(meanX)==1) & (sign(meanY)==-1))
     A=180+atan(meanX/meanY)*180/pi;
 end
 if ((sign(meanX)==-1) & (sign(meanY)==-1))
     A=180+atan(meanX/meanY)*180/pi;
 end
 if ((sign(meanX)==-1) & (sign(meanY)==1))
     A=360+atan(meanX/meanY)*180/pi;
 end

end