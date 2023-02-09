%%%%%%%% Online KPCA %%%%%%%%%%

%%%%%%%%% Training Phase %%%%%%%%%%
d = 1; %Kernel parameter

X1 = P(1:400,:);
[T1 F1] = size(X1);
X1 = X1 - repmat( mean(X1,1) , size(X1,1) , 1 ); %Then subtract mean of each dimension across all timesteps.
X1_mean = mean(X1,1); %Store mean value of X_train for each dimension(col)

% Determine Kernel matrix %
for i=1:T1
    for j=1:T1
        K1(i,j) = kernel(X1(i,:), X1(j,:), 3, d);
    end
end


% Find eigenvectors and sort in descending order %
[A1 D1] = eig(K1);
lambda1 = diag(D1);
lambda1(find(lambda1<=0))=eps;
lambda1 = flipud(lambda1);
D1 = diag(lambda1);
A1 = fliplr(A1);

% Normalize eigenvectors to 1/eigenvalue %
A1_norm = A1./sqrt(repmat(lambda1',T1,1)); %Cols of A_norm are the (normalized) eigenvectors
%To check: sum(A_norm.^2).*lambda' must = 1


X1_hat = A1_norm'*K1;
X1_hat = X1_hat'; %Rows are flow vectors for each timestep, cols are projections onto each eigenvector

normSquare1 = sum(X1_hat.^2);
var1 = normSquare1/sum(normSquare1); %Percentage of variances explained by each eigenvector

r = 10;
X1_state = A1_norm(:,1:r)'*K1;
X1_state = X1_state';
X1_tilde = A1_norm(:,r+1:end)'*K1;
X1_tilde = X1_tilde';
X1_residual = sum(X1_tilde.^2,2);

%%%%%%%%% Test Phase %%%%%%%%%%

X2 = P(401:end,:);
%X = P;
[T2 F2] = size(X2);

for t=1:T2
    
    x = X2(t,:); %Arriving input vector
    
    % IMPORTANT: Subtract mean of each dimension of X_train from X_test;
    x = x - X1_mean;
    %Must do this approximation as mean of X_test cannot be known ahead of time
    %in a real-time application
    
    for i=1:T1
        k(i,1) = kernel(X1(i,:), x, 3, d); %Determine kernel of arriving input vector with training pts
    end %for i=1:T1
    
    X2_state = A1_norm(:,1:r)'*k; %Project on to eigenvectors obtained from the training set
    X2_state = X2_state';
    X2_tilde = A1_norm(:,r+1:end)'*k;
    X2_tilde = X2_tilde';
    X2_residual(t+T1) = sum(X2_tilde.^2,2); %Magnitude of proj on to residual space = THE detection statistic
    
end %for t=1:T2

