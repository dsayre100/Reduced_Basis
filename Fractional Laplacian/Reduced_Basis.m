%Reduced Basis Method
clear, clc

%Define tolerance
tol = 10^(-10);

%Generate ROM basis (using POD)
%Declare number of snapshots to use
snaps = 5;
[u,~,~] = driver_fracLap(snaps);
[V_r,~] = POD(u);

%Determine number of training set samples 
theta_samples =  50;

%Calculate Ay,by,uy for training set of theta_samples
[uy,Ay,by] = driver_fracLap(theta_samples);

%initialize error and count variables
e_calculated(1) = 1; 
k =0;


while max(e_calculated) > tol && k < 2*theta_samples
    k = k+1;
    for i = 1:theta_samples
        
        %Calculate Approximations
        A_bar = V_r'*Ay{i}*V_r;
        b_bar = V_r'*by;
        y_bar = V_r*(A_bar \ b_bar); 
        
        %Record Erros
        e_calculated(i) = norm(by - Ay{i}*y_bar);
        e_actual(i) = norm(uy(:,i) - (y_bar));
    end
    
    %Determine if tolerance has been met
    if max(e_calculated) < tol
        break
    else
        %Add additional vectors to Reduced Basis
        
        %[~,I] = max(e); Found an issue where if we took the greedy
        %approach we can get stuck in an infinite loop due to
        %orthogonalization. Adding a random vector from the training set
        %instead of the vector associated with the max error avoided the 
        %infinite loop
        
        %Pick random solution vector to add to basis
        I = randsample(1:length(e_calculated),1);
        
        %Append solution vector to V_r
        V_r = [V_r,uy(:,I)];
        
        %Orthogonalize the new basis vectors 
        %Note: if this step is skipped the algorithm breaks down due to the
        %condition number growing rapidly
        V_r = orth(V_r);
    end
    
end

%Generate a test set:
%Determine number of test samples
t_samples = 100;
[uy,Ay,by,snapshots] = driver_fracLap(t_samples);

for i = 1:t_samples
        
        %Calculate Approximations
        A_bar = V_r'*Ay{i}*V_r;
        b_bar = V_r'*by;
        y_bar = V_r*(A_bar \ b_bar); 
        
        %Record Erros
        e_calculated(i) = norm(by - Ay{i}*y_bar);
        e_actual(i) = norm(uy(:,i) - (y_bar));
    end


%plot errors for verification
subplot(1,2,1)
semilogy(snapshots,e_actual,'ro')
xlabel('Sample Value')
title('Actual error')


subplot(1,2,2)
semilogy(snapshots,e_calculated, 'bo')
title('Calculated Error')
xlabel('Sample Value')

if max(e_calculated) > tol
    disp('Tolerance not met')
    disp('Maximum Calculated Error is:')
    disp(max(e_calculated))
end