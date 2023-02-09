clc
clear all

jpgFiles = dir('C:\MATLAB7\work\*.jpeg');
for k = 1:length(jpgFiles)
  filename = jpgFiles(k).name;
  a= double(imread(filename)); 
  a1=dwt2(a,'haar');
%   a1 = fft2(a);
%   a2=imresize(abs(a1),0.1,'bil');

a2=imresize(a1,0.1,'bil');

[m n]=size(a2);
a2_expand=[];
for r=1:m
a2_expand=[a2_expand a2(r,:)];
end
X(k,:)=[a2_expand];
end
% size(X)
% imshow(jpgFiles(8).name)

nu=2; k=3;

X = X./repmat(sqrt(sum(X.*X,2)+eps),1,size(X,2)); %normalize to unit circle (i.e. divide by norm)
n = size(X,1);
h_i = zeros(1,n);

%Get distances to k'th nearest neighbour
clear g g_list;
for i=1:n
    [g(i,:) g_list(i,:)]= M1(X(i,:),X,k); %Sparsity measure M1 is the Euclidean distance
end %for i=1:n

g_sorted = sort(g, 'ascend');

z=1;
for i=1:length(jpgFiles)
    if g(i)>0.03
        figure(i)
        imshow(jpgFiles(i).name)
        if z==1
%         
%         load gong;
%        %load chirp;
%         y1 = y; Fs1 = Fs;
%         
%         wavplay(y1,Fs1,'sync') % The chirp signal finishes before the 
%         wavplay(y,Fs)
%          [y,Fs] = wavread('AVSEQ01');
%             wavplay(y,Fs);
        end
        z=2;
    end
end
figure(i+1)

subplot(2,1,1);stem(g); title('g');
%figure(2)
subplot(2,1,2);stem(g_sorted); title('g_sorted');
