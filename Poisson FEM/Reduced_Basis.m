%Reduced Basis Method
clear, clc

%Define tolerance
tol = 10^(-8);
%Generate ROM basis (using POD)
[V_r] = POD(7);
%Generate training set 
theta_samples =  10*rand(1,25);

% Calculate extended basis
j= 1;
e(1)=1;

%Calculate Ay,by,uy for training set of y
for i =1:size(theta_samples,2)
    [u,A,b] = poissonfem(0,1,0,1,50,50,theta_samples(i));
    uy(:,i) = u;
    Ay(:,:,i) = A;
    by(:,i) = b;
end

%Generate error variables
e = zeros(1,size(theta_samples,2));
e_actual = zeros(1,size(theta_samples,2));

for i = 1:size(theta_samples,2)
%Calculate approximate solutions
A_approx(:,:,i) = V_r'*Ay(:,:,i)*V_r;
b_approx = V_r'*by(:,i);
y_approx(:,i) = V_r*(A_approx(:,:,i) \ b_approx);

%Compute error

e(i) = norm(by- Ay(:,:,i)*y_approx(:,i));
e_actual(i) = norm(uy(:,i) - y_approx(:,i));
end


while (max(e) > tol)
    [~,I] = max(e);
    V_r = [V_r,uy(:,I)];
    V_r = orth(V_r);
    
    %clear variables
    clear A_approx
    clear b_approx
    clear y_approx
    e = zeros(1,size(theta_samples,2));
    e_actual = zeros(1,size(theta_samples,2));
    
    
    %Calculate approximate solutions
    for i = 1:size(theta_samples,2)
        A_approx(:,:,i) = V_r'*Ay(:,:,i)*V_r;
        b_approx = V_r'*by(:,i);
        y_approx(:,i) = V_r*(A_approx(:,:,i) \ b_approx);
        
        %Compute error
        
        e(i) = norm(by - Ay(:,:,i)*y_approx(:,i));
        e_actual(i) = norm(uy(:,i) - y_approx(:,i));
    end
end

%Validation Set 
val_samples = 10*rand(1,50);


for i =1:size(val_samples,2)
    [u,A,b] = poissonfem(0,1,0,1,50,50,val_samples(i));
    u_val(:,i) = u;
    A_val(:,:,i) = A;
    b_val(:,i) = b;
    
    
    A_av(:,:,i) = V_r'*A_val(:,:,i)*V_r;
    b_av = V_r'*b_val(:,i);
    y_av(:,i) = V_r*(A_av(:,:,i) \ b_av);
    
    %Compute error
    e(i) = norm(b_val - A_val(:,:,i)*y_av(:,i));
    e_actual(i) = norm(u_val(:,i) - y_av(:,i));
end    
    figure(1)
    subplot(2,1,1)
    plot(val_samples,e_actual,'rs')
    title('Actual Error')
    subplot(2,1,2)
    plot(val_samples,e,'bo')
    title('Calculated Error')
    
    
%create mesh and plot   
h=1/50; 
h2=h^2; 
k=1/50; 
k2=k^2; 
hk=h*k;
x=(0:50)*h; % set mesh values
y=(0:50)*k;


fig2 = figure(2);
w=reshape(u_val(:,1),51,51);
subplot(2,1,1)
mesh(x,y,w')
title('Galerkin Approximation')

w1 = reshape(y_av(:,1),51,51);
subplot(2,1,2)
mesh(x,y,w1')
title('Reduced Basis Approximation')

