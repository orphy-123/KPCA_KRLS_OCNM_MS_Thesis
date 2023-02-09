function [d,x]=f1(r)
A=[5 2*r r; 3 6 2*r-1; 2 r-1 3*r];
b=[2 3 5]';
d=det(A);
x=A\b;