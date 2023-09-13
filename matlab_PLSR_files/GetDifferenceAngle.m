function A=GetDifferenceAngle(angle1,angle2)
 n=length(angle1);
 angles=repmat(NaN,1,n);
 for i=1:n
  y1=sin(angle1(i)*pi/180);
  y2=sin(angle2(i)*pi/180);
  x1=cos(angle1(i)*pi/180);
  x2=cos(angle2(i)*pi/180);
  distance12=sqrt((x1-x2)^2+(y1-y2)^2);
  angle12=acos((-(distance12^2)+2)/2)*180/pi;
  angles(i)=angle12;
 end
 A=angles;
end