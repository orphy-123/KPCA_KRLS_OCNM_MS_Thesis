clc
clear all

jpgFiles = dir('C:\Documents and Settings\supriyo.BRACU\Desktop\work\*.jpeg');
for k = 1:length(jpgFiles)
    filename = jpgFiles(k).name;
    a = imread(filename); 
    a1=dwt2(a,'haar');
    %a1 = fft2(a);
    % a2=imresize(abs(a1),0.1,'bil');
    a2=imresize(a1,0.1,'bil');

    [m n]=size(a2);
    a2_expand=[];
    for r=1:m
        a2_expand=[a2_expand a2(r,:)];
    end
    P(k,:)=[a2_expand];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Training Period
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X = P(1:20,:);
[T f] = size(X);
[T f] = size(X);
[T f] = size(X);

%X = X'; %Transpose X
X = X - repmat( mean(X,1) , size(X,1) , 1 ); %Then subtract of each dimension across 864 timesteps.

C = X'*X;
[V D] = eig(C); %columns of V are the e-vectors
d = diag(D);
V = fliplr(V);
d = flipud(d);
D = diag(d);


normSquare = sum((X*V).^2);
var = normSquare/sum(normSquare); %Percentage of variances
u = X*V;
u = u ./ repmat( sqrt(normSquare) , size(X,1) , 1 );

r=3;
R = V(:,1:r);
X_hat = R*R'*X'; %Projections
X_tilde = ( eye(size(R,1)) - R*R' ) * X';
X_state=sum(X_hat'.^2,1);
X_residual=sum(X_tilde.^2,1);

%%% Threshold calculation %%%

%T = size(X,1) %%%%%%%%%%%
d = d/(T-1); %%%%%%%%

ErrorSub = r+1:size(d,1);
phi1 = sum(d(ErrorSub));
phi2 = sum(d(ErrorSub).^2);
phi3 = sum(d(ErrorSub).^3);

h0 = 1 - (2*phi1*phi3)/(3*phi2*phi2);
c=3.090; %p=0.999

size(ErrorSub,2);
while (h0<0)
    %%reduce the higher numbered evectors one by one
    %%and recompute h0 with them omitted
    ErrorSub = ErrorSub(1:end-1);
    phi1 = sum(d(ErrorSub));
    phi2 = sum(d(ErrorSub).^2);
    phi3 = sum(d(ErrorSub).^3);
    h0 = 1 - (2*phi1*phi3)/(3*phi2*phi2);
    size(ErrorSub,2);
end %while

delta_pca = phi1 * ( c*sqrt(2*phi2*h0*h0)/phi1 + 1 + phi2*h0*(h0-1)/(phi1*phi1) )^(1/h0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run for future timesteps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X = P(21:54,:);
[T f] = size(X);
%X = X'; %Transpose X
X = X - repmat( mean(X,1) , size(X,1) , 1 ); %Then subtract of each dimension across all timesteps.

X_hat = R*R'*X'; %Projections
X_tilde = ( eye(size(R,1)) - R*R' ) * X';
X_state=sum(X_hat'.^2,1);
X_residual=sum(X_tilde.^2,1);


% for i=1:(length(jpgFiles)-20)
%     if X_residual(i)>4*10^4                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
%         figure(i)
%         imshow(jpgFiles(i+19).name)
%     end
% end


figure(i+1)
subplot(3,1,1);stem(var(1:5)); title('Variance');
subplot(3,1,2);stem(X_state); title('X_state');
subplot(3,1,3);stem(X_residual); title('X_residual');














%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear lambda1 lambda2 lambda3 phi1 phi2 phi3 h0 c C D V d normSquare u r R X_hat X_tilde X_state;