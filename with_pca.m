clc;
clear all;
a=imread('a.jpg');
a1=dwt2(a,'haar');
a2=imresize(a1,0.25,'bil');
b=imread('b.jpg');
b1=dwt2(b,'haar');
b2=imresize(b1,0.25,'bil');
c=imread('c.jpg');
c1=dwt2(c,'haar');
c2=imresize(c1,0.25,'bil');
d=imread('d.jpg');
d1=dwt2(d,'haar');
d2=imresize(d1,0.25,'bil');
e=imread('e.jpg');
e1=dwt2(e,'haar');
e2=imresize(e1,0.25,'bil');
f=imread('f.jpg');
f1=dwt2(f,'haar');
f2=imresize(f1,0.25,'bil');
g=imread('g.jpg');
g1=dwt2(g,'haar');
g2=imresize(g1,0.25,'bil');
h=imread('h.jpg');
h1=dwt2(h,'haar');
h2=imresize(h1,0.25,'bil');
i=imread('i.jpg');
i1=dwt2(i,'haar');
i2=imresize(i1,0.25,'bil');
j=imread('j.jpg');
j1=dwt2(j,'haar');
j2=imresize(j1,0.25,'bil');

[m n]=size(a2);
a2_expand=[];
for r=1:m
a2_expand=[a2_expand a2(r,:)];
end
size(a2_expand);

[m n]=size(b2);
b2_expand=[];
for r=1:m
b2_expand=[b2_expand b2(r,:)];
end
size(b2_expand);

[m n]=size(c2);
c2_expand=[];
for r=1:m
c2_expand=[c2_expand c2(r,:)];
end
size(c2_expand);

[m n]=size(d2);
d2_expand=[];
for r=1:m
d2_expand=[d2_expand d2(r,:)];
end
size(d2_expand);

[m n]=size(e2);
e2_expand=[];
for r=1:m
e2_expand=[e2_expand e2(r,:)];
end
size(e2_expand);

[m n]=size(f2);
f2_expand=[];
for r=1:m
f2_expand=[f2_expand f2(r,:)];
end
size(f2_expand);

[m n]=size(g2);
g2_expand=[];
for r=1:m
g2_expand=[g2_expand g2(r,:)];
end
size(g2_expand);

[m n]=size(h2);
h2_expand=[];
for r=1:m
h2_expand=[h2_expand h2(r,:)];
end
size(h2_expand);

[m n]=size(i2);
i2_expand=[];
for r=1:m
i2_expand=[i2_expand i2(r,:)];
end
size(i2_expand);

[m n]=size(j2);
j2_expand=[];
for r=1:m
j2_expand=[j2_expand j2(r,:)];
end
size(j2_expand);

X=[a2_expand;b2_expand;c2_expand;d2_expand;e2_expand;f2_expand;g2_expand;h2_expand;i2_expand;j2_expand];
%size(X)

[T f] = size(X)

%X = X_P'; %Transpose X, if data matrix is in transposed form
X = X - repmat( mean(X,1) , size(X,1) , 1 ); %Then subtract off each dimension across all timesteps.

C = X'*X;
[V D] = eig(C); %Columns of V are the e-vectors
d = diag(D);
V = fliplr(V);
d = flipud(d);
D = diag(d);

normSquare = sum((X*V).^2);
var = normSquare/sum(normSquare); %Percentage of variances
u = X*V;
u = u ./ repmat( sqrt(normSquare) , size(X,1) , 1 );

r=4; %Number of Principal Components  allocated to normal subspace
R = V(:,1:r);
X_hat = R*R'*X'; %Projections
X_tilde = ( eye(size(R,1)) - R*R' ) * X';
X_state=sum(X_hat.^2,1);
X_residual=sum(X_tilde.^2,1);
subplot(2,1,1);stem(var);
subplot(2,1,2);stem(X_residual);