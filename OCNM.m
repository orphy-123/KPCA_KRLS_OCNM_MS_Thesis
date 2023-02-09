function [h_out, g_out g_list_out] = OCNM_2(X,nu,k)
%X is dataset of size n, nu defines percentage of highest density region,
%    and k is nearest-neighbour decider.
%h_i is decision for each i=1:n,
%   g is dist to k'th nearest-neighbour for each pt,
%   g_list is full list of dists (NOT sorted) for every point, to every point

X = X./repmat(sqrt(sum(X.*X,2)+eps),1,size(X,2)); %normalize to unit circle (i.e. divide by norm)
n = size(X,1);
h_i = zeros(1,n);

%Get distances to k'th nearest neighbour
clear g g_list;
for i=1:n
    [g(i,:) g_list{i,:}]= M1(X(i,:),X,k); %Sparsity measure M1 is the Euclidean distance
end %for i=1:n

g_sorted = sort(g, 'ascend');
% Evaluate decision fucntion
if mod(nu*n,1)==0 %nu*n is a Natural number
    rho_star = g_sorted(nu*n);
else %nu*n is NOT a Natural number
    rho_star = g_sorted(floor(nu*n)+1);
end %if mod(nu*n,1)==0
h_i = repmat(rho_star,n,1)-g;
h_i = sign(h_i);

h_out = h_i;
g_out = g;
g_list_out = g_list; %%%Output (optional) distance (NOT sorted) from every pt, in turn, to every point in the sample
