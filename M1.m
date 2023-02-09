function [d_k dist_out] = M1(x, s, k)
if k+1 > size(s,1)
    error('Value of k must be less than the number of pts');
end %if k+1 > size(s,1)

squares = (repmat(x,size(s,1),1) - s).^2;
dist = sqrt(sum(squares,2));
dist_out = dist; %Distance of "x" to every point in "s" (NOT sorted)
[dist, index] = sort(dist, 'ascend');
d_k = dist(k+1); %Note: NOT dist(k), as k=1 then yields dist with itself, equal to 0.