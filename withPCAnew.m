clc
clear all
tic
jpgFiles = dir('C:\Documents and Settings\supriyo.BRACU\Desktop\work\*.jpeg');
for k = 1:length(jpgFiles)
  filename = jpgFiles(k).name;
  a= double(imread(filename)); 
  a1=dwt2(a,'haar');
%    a1 = fft2(a);
%   a2=imresize(abs(a1),0.1,'bil');

a2=imresize(a1,0.1,'bil');

[m n]=size(a2);
a2_expand=[];
for r=1:m
a2_expand=[a2_expand a2(r,:)];
end
X(k,:)=[a2_expand];
end

[T f] = size(X);

trueint=zeros(1,T);

act=[14 20 21 22 30 31 32 34 40 41 44 61 78 105 111 112 113 114 135];
trueint(act)=1;



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

r=3; %Number of Principal Components  allocated to normal subspace
R = V(:,1:r);
X_hat = R*R'*X'; %Projections
X_tilde = ( eye(size(R,1)) - R*R' ) * X';
X_state=sum(X_hat.^2,1);
X_residual=sum(X_tilde.^2,1);


for threshold=.5 : 0.5 : 5
  z=0;
  sumdec=0;

    flagint=zeros(1,T);
        for i=1:length(jpgFiles)
            z=z+1; 
                if X_residual(i)>threshold*10^4
                    sumdec=sumdec+1;
                    flagint(z)=1;
                     flagint(1)=0;
% %         flagint;
                    detected=bitand(flagint,trueint);
                    false=bitxor(flagint,trueint);
                    falsem=false;
                    false(act)=0;
                    missed=bitxor(falsem,false);
                    mis(sumdec)=sum(missed);
                    dec(sumdec)=(sum(detected)/19)*100;
                     fal(sumdec)=(sum(false)/(194-19))*100;
       % scatter(sort(fal),sort(dec));
  %       scatter(fal,dec);
     %   hold on
     
     
%         figure(i)                     % Shows the ANOMALIES
%         imshow(jpgFiles(i).name)
%                if z==1

%          [y,Fs] = wavread('AVSEQ01'); % For Alarm
%             wavplay(y,Fs);
%               end
%       z=2;
                end
       
        end
  %       scatter(sort(fal),sort(dec));
         scatter(fal,dec);
        hold on
end
 




% figure(i+1)
% subplot(3,1,1);stem(var(1:5)); title('Variance');
% subplot(3,1,2);stem(X_state); title('X_state');
% subplot(3,1,3);stem(X_residual); title('X_residual');
toc
















































% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% % % d = d/(T-1);
% % % 
% % % ErrorSub = r+1:size(d,1);
% % % phi1 = sum(d(ErrorSub));
% % % phi2 = sum(d(ErrorSub).^2);
% % % phi3 = sum(d(ErrorSub).^3);
% % % 
% % % h0 = 1 - (2*phi1*phi3)/(3*phi2*phi2);
% % % c=3.090; %p=0.999
% % % 
% % % size(ErrorSub,2);
% % % while (h0<0)
% % %     %%reduce the higher numbered evectors one by one
% % %     %%and recompute h0 with them omitted
% % %     ErrorSub = ErrorSub(1:end-1);
% % %     phi1 = sum(d(ErrorSub));
% % %     phi2 = sum(d(ErrorSub).^2);
% % %     phi3 = sum(d(ErrorSub).^3);
% % %     h0 = 1 - (2*phi1*phi3)/(3*phi2*phi2);
% % %     size(ErrorSub,2);
% % % end %while
% % % 
% % % delta_PCA = phi1 * ( c*sqrt(2*phi2*h0*h0)/phi1 + 1 + phi2*h0*(h0-1)/(phi1*phi1) )^(1/h0);
% % % alarm_PCA = find(X_residual>delta_PCA);
% % % 
% % % if nargout == 1
% % %     alarm_PCA_out = alarm_PCA; 
% % % elseif nargout == 2
% % %     alarm_PCA_out = alarm_PCA; X_residual_out = X_residual;
% % % elseif nargout == 3
% % %     alarm_PCA_out = alarm_PCA; X_residual_out = X_residual; delta_PCA_out = delta_PCA;
% % % end
