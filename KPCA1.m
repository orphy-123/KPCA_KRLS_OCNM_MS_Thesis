
X = P(1:400,:);
X = P;
[T F] = size(X);
X = X - repmat( mean(X,1) , size(X,1) , 1 ); %Then subtract off each dimension across all timesteps.


% Determine Kernel matrix %
tic;
for i=1:T
    for j=1:T
        K(i,j) = kernel(X(i,:), X(j,:));
    end
end
toc

% Find eigenvectors and sort in descending order %
[A D] = eig(K);
lambda = diag(D);
lambda(find(lambda<=0))=eps;
lambda = flipud(lambda);
D = diag(lambda);
A = fliplr(A);

% Normalize eigenvectors to 1/eigenvalue %
A_norm = A./sqrt(repmat(lambda',T,1)); %Cols of A_norm are the (normalized) eigenvectors
%To check: sum(A_norm.^2).*lambda' must = 1


X_hat = A_norm'*K;
X_hat = X_hat'; %Rows are flow vectors for each timestep, cols are projections onto each eigenvector

normSquare = sum(X_hat.^2);
var = normSquare/sum(normSquare); %Percentage of variances explained by each eigenvector

r = 10;
X_state = A_norm(:,1:r)'*K;
X_state = X_state';
X_tilde = A_norm(:,r+1:end)'*K;
X_tilde = X_tilde';
X_residual = sum(X_tilde.^2,2);
