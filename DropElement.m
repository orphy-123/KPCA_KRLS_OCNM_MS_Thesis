% DropElement(p) %
%Reorganize K_tilde_p and K_tilde_inv_p, with p'th row/col moved to the end
K_tilde = [ K_tilde(1:p-1,1:p-1) K_tilde(1:p-1,p+1:m) K_tilde(1:p-1,p)  ;  K_tilde(p+1:m,1:p-1) K_tilde(p+1:m,p+1:m) K_tilde(p+1:m,p)  ;  K_tilde(p,1:p-1) K_tilde(p,p+1:m) K_tilde(p,p) ];
K_tilde_inv = [ K_tilde_inv(1:p-1,1:p-1) K_tilde_inv(1:p-1,p+1:m) K_tilde_inv(1:p-1,p)  ;  K_tilde_inv(p+1:m,1:p-1) K_tilde_inv(p+1:m,p+1:m) K_tilde_inv(p+1:m,p)  ;  K_tilde_inv(p,1:p-1) K_tilde_inv(p,p+1:m) K_tilde_inv(p,p) ];
delta_p = 1/(K_tilde_inv(m,m));
a_tilde_p = -delta_p*[K_tilde_inv(1:m-1,m)];
K_tilde_inv = K_tilde_inv(1:m-1,1:m-1)-a_tilde_p*a_tilde_p'/delta_p;
alpha = alpha - (1/delta_p)*[a_tilde_p*a_tilde_p'  -a_tilde_p  ;  -a_tilde_p'  1] *K_tilde*alpha;
alpha = alpha(1:m-1);
K_tilde = K_tilde(1:m-1,1:m-1);
Dictionary(:,p) = [];
drop_index(p) = [];
Lambda(:,p) = [];
dotProd(:,p) = []; %%%
m=m-1;
m_t(t) = m;

ResetP; %Reset P %
