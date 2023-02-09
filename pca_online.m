%clear;
load Anukool\AbileneAllTrafficViews-15to21.mat P;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Training Period
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X = P(1:500,:);
[T f] = size(X);
%X = X(1:2*floor(T/2),:); %Ensure even (mod 2) number of timesteps
[T f] = size(X);
%X = [ X(1:2:T-1,:) X(2:2:T,:) ]; %Block 2 timesteps at a time
%X = [ X(1:4:T-3,:) X(2:4:T-2,:) X(3:4:T-1,:) X(4:4:T,:)]; %Block 4 timesteps at a time
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

r=4;
R = V(:,1:r);
X_hat = R*R'*X'; %Projections
X_tilde = ( eye(size(R,1)) - R*R' ) * X';
X_state=sum(X_hat'.^2,1);
X_residual=sum(X_tilde.^2,1);

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

X = P(501:2010,:);
[T f] = size(X);
%X = X'; %Transpose X
X = X - repmat( mean(X,1) , size(X,1) , 1 ); %Then subtract of each dimension across all timesteps.

X_hat = R*R'*X'; %Projections
X_tilde = ( eye(size(R,1)) - R*R' ) * X';
X_state=sum(X_hat'.^2,1);
X_residual=sum(X_tilde.^2,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear lambda1 lambda2 lambda3 phi1 phi2 phi3 h0 c C D V d normSquare u r R X_hat X_tilde X_state;