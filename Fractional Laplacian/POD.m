function [V_r,u_mean]= POD(u)

% POD basis of fractional Laplacian

%tolerance
tol = 10^-2;

%set mean to 0
[~,k] = size(u);
u_mean = (1/k)*sum(u,2);

u = u-u_mean;


%SVD
[U,S] = svd(u);
d_s = diag(S);

%Initialize sum and count variables
sum_t = 0;
j = 0;

%Collect Left Singular Vectors based on tolerance
while (sum_t) < tol && j < k
    sum_t = sum_t + d_s(end-j);
    j = j+1;
end

if (sum_t) > tol
    V_r = U(:,1:k-(j-1));
else 
    V_r = U(:,1:k);
end
