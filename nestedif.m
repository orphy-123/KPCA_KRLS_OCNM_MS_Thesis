 y=0;
x=40;
 if x>5
        if x>10 
        y=y+1
        break
        end
        y=y-1
 else
    y=y-4
 end

 n=1;
 x=0;
 A=[];
 while x<10000
     x=2^n;
     if x>=10000
         break
     
     end
     A(n)=x;
     n=n+1;
 end
 
 f=10;
 A=5;
 t=0:0.01:1;
 y=A*sin(2*pi*f*t);
 plot(t,y)
 stem(t,y)
 stem(y)
 subplot(2,1,1), plot(y)
 subplot(2,1,2), stem(y)
 stem(y,'r*')
stem(y,'r*')
 stem(y,'r-')
 xlabel('time(ms)');
 ylabel('volt')
 legend('y1','y2')
 
  t=0:0.1:30;
  f1=1/10; A1=2;
  y1=A1*sin(2*pi*f1*t);
  f2=1/15; A2=3;
  y2=A2*sin(2*pi*f2*t);
  plot(t,y1,'r*')
  hold on
  plot(t,y2,'gs')
   xlabel('time(ms)');
 ylabel('volt')
 legend('y1','y2')