% Reset P %
P = R*eye(m);
for i_r=1:r
    k_tilde = zeros(m,1);
    for j=1:m
        k_tilde(j) = kernel(Dictionary(:,j),X(t-i_r,:)'); %Computing k_tilde{t-1}
    end %for j=1:m
    a = K_tilde_inv*k_tilde;
    q = (P*a) / (gamma+a'*P*a);
    P = (1/gamma)*[P - q*a'*P];
    alpha = alpha + K_tilde_inv*q*(Y(t-i_r)-k_tilde'*alpha);
end %for i_r=1:r-1
