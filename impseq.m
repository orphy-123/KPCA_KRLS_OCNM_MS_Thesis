function [x,n]=impseq(n0,n1,n2)
n=n1:n2;
l=length(n);
x=zeros(1,l);
x(n0-n1+1)=1;