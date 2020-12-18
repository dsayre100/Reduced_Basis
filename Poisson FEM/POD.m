function [V_r,samples]= POD(snaps)

% POD basis of fractional Laplacian

%tolerance
tol = 10^(-2);


%Generate solution vectors
samples = 10*rand(1,snaps);
n = size(samples,2);
for i = 1:n
u(:,i) = poissonfem(0,1,0,1,50,50,samples(i));
end

%set mean to 0
[~,k] = size(u);
u_mean = (1/k)*sum(u,2);

u = u-u_mean;


%Use SVD to get "most important" left singular vectors to form basis
[U,S] = svd(u);
d_s = diag(S);

sum_t = 0;
j = 1;
while (sum_t) < tol
    sum_t = sum_t + d_s(end-j);
    j = j+1;
end

V_r = U(:,1:n-(j-1));