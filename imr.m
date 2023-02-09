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
size(X)

nu=2; k=5;

X = X./repmat(sqrt(sum(X.*X,2)+eps),1,size(X,2)); %normalize to unit circle (i.e. divide by norm)
n = size(X,1);
h_i = zeros(1,n);

%Get distances to k'th nearest neighbour
clear g g_list;
for i=1:n
    [g(i,:) g_list(i,:)]= M1(X(i,:),X,k); %Sparsity measure M1 is the Euclidean distance
end %for i=1:n

g_sorted = sort(g, 'ascend');

subplot(2,1,1);stem(g);
subplot(2,1,2);stem(g_sorted);

% for i=1:10
%     if g(i)>0.2
%     
%     end
% end
% 
