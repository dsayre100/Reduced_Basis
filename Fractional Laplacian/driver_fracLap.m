function [v, A_bar, b,snapshots] = driver_fracLap(samples)

% Solve fractional Poisson problem with non-zero Dirichlet boundary
% conditions (Neumann and other boundary cases are similar, see 
% \cite{}. 
% 
% Authors: Harbir Antil
%          Department of Mathematical Sciences
%          George Mason University
%          Fairfax VA
%
%          Johannes Pfefferer 
%          Chair of Optimal Control, 
%          Center of Mathematical Sciences, 
%          Technical University of Munich, 
%          Garching by Munich, Germany
% 
% Date created:   December 12, 2017
% Last modified:  December 22, 2017
% 
%
% CITATION: IF YOU USE THIS CODE THEN PLEASE CITE THE LIST OF PAPERS GIVEN 
%           AT THE BOTTOM OF THIS FILE. 
% 
% You can find the technical report with implementation details here:
%
%    http://math.gmu.edu/~hantil/Tech_Report/HAntil_JPfefferer_2017a.pdf
% 


set(0, 'defaultaxesfontsize',16,'defaultaxeslinewidth',1,...
      'defaultlinelinewidth',1,'defaultpatchlinewidth',0.7,...
      'defaulttextfontsize',18);


%% generate exact solution data
data = exact_soln;
s = .2;

for ind = 5:5

% Generate the mesh
mesh = rect_grid2(0, 1, 0, 1, 2^(ind+1), 2^(ind+1));

% % display the mesh 
% figure(1); clf
% triplot(mesh.t, mesh.p(:,1), mesh.p(:,2)); drawnow;    
% pause 

% assemble stiffness and mass matrices 
[A,M,FreeNodes,dirichlet] = Stiff_Mass(mesh);

N = size(A,1); 

%% Solve the fractional Poission problem with zero boundary conditions using 
%% A. Bonito and J.E. Pasciak 2015, 2017
rhs = M*data.f(mesh.p);
b = rhs(FreeNodes);
%b = rhs;


u = zeros(N,1);

h = N^(-1/2);    
%c = 2/(pi^2);
k = 1/log(1/h);

Nplus  = ceil(pi^2/(4*s*k^2));
Nminus = ceil(pi^2/(4*(1-s)*k^2));

A_bar = struct([]);
i=1;
snapshots = randsample(-Nminus:1:Nplus,samples);
for ell=snapshots
    y = k*ell;
    %A_bar{i} = A; 
    %A_bar{i}(FreeNodes,FreeNodes) = exp((1-s)*y)*(exp(y)*M(FreeNodes,FreeNodes)+A(FreeNodes,FreeNodes));
    A_bar{i} = exp((1-s)*y)*(exp(y)*M(FreeNodes,FreeNodes)+A(FreeNodes,FreeNodes));
    u(FreeNodes) = A_bar{i} \ b;
    
    %u(FreeNodes) = (sin(s*pi)/pi) * k * u(FreeNodes);
    %v(:,i) = u
    v(:,i) = u(FreeNodes);
    

    i=i+1;
end    



%% Impose Nonzero boundary conditions using 
%% H. Antil, J. Pfefferer and S. Rogovs 2017
% uD = zeros(N,1);
% uD(dirichlet) = data.exactu(mesh.p(dirichlet,:),s); 
% In the above line we are first computing evaluating the boundary data at 
% thenodes which is fine if the boundary data is continuous. In case 
% boundary data is not continuous it is better to do L^2 projection first 
% on the boundary                                                     
                                                                                                                                                            
% v   = uD;
% rhs = -A*uD;
% v(FreeNodes) = A(FreeNodes,FreeNodes) \ rhs(FreeNodes);
% 
% u = u+v;
% 
% % L2-error  
% errL2(ind,1) = sqrt( (u-data.exactu(mesh.p,s))'*M*(u-data.exactu(mesh.p,s)) );
% 
% hmax(ind,1) = h;
end

% trisurf(mesh.t,mesh.p(:,1),mesh.p(:,2),u)
% 
% figure;   
% loglog( hmax, errL2, '-*');
% legend('L^2-error','Location','northwest')
% xlabel('h (meshsize)')
% ylabel('$\|u-u_h\|_{L^2(\Omega)}$','interpreter','latex', 'FontSize', 16)
% axis tight


%% CITATION

%% 1. Please see the technical report for the code

%@TECHREPORT{HAntil_JPfefferer_2017a,
%   author      = {H. Antil and J. Pfefferer},
%   title       = {A short Matlab implementation of fractional Poission equation with nonzero boundary conditions},
%   institution = {Department of Mathematical Sciences},
%   address     = {George Mason University},
%   year        = 2017,
%   note        = "\url{http://math.gmu.edu/~hantil/Tech_Report/HAntil_JPfefferer_2017a.pdf}"
%}


%% 2. Please see the paper H. Antil, J. Pfefferer and S. Rogovs 2017 for fractional PDEs
%     with non-zero boundary conditions

%@article{HAntil_JPfefferer_SRogovs_2017a,
%  title={Fractional Operators with Inhomogeneous Boundary Conditions: Analysis, Control, and Discretization},
%  fauthor={Antil, Harbir and Pfefferer, Johannes and Rogovs, Sergejs},
%  author={Antil, H. and Pfefferer, J. and Rogovs, S.},
%  journal={arXiv preprint arXiv:1703.05256},
%  year={2017}
%}

%% 3. Please see the paper A. Bonito and J.E. Pasciak 2015, 2017 for solving fractional PDEs
%     with zero boundary conditions

%@article {ABonito_JEPasciak_2015a,
%   FAUTHOR = {Bonito, Andrea and Pasciak, Joseph E.},
%    AUTHOR = {Bonito, A. and Pasciak, J.E.},
%     TITLE = {Numerical approximation of fractional powers of elliptic
%              operators},
%   JOURNAL = {Math. Comp.},
%  FJOURNAL = {Mathematics of Computation},
%    VOLUME = {84},
%      YEAR = {2015},
%    NUMBER = {295},
%     PAGES = {2083--2110},
%      ISSN = {0025-5718},
%   MRCLASS = {65N30 (65R20)},
%  MRNUMBER = {3356020},
%MRREVIEWER = {Igor Bock},
%       DOI = {10.1090/S0025-5718-2015-02937-8},
%       URL = {http://dx.doi.org/10.1090/S0025-5718-2015-02937-8},
%}

%@article {ABonito_JEPasciak_2017a,
%   FAUTHOR = {Bonito, Andrea and Pasciak, Joseph E.},
%    AUTHOR = {Bonito, A. and Pasciak, J.E.},    
%     TITLE = {Numerical approximation of fractional powers of regularly
%              accretive operators},
%   JOURNAL = {IMA J. Numer. Anal.},
%  FJOURNAL = {IMA Journal of Numerical Analysis},
%    VOLUME = {37},
%      YEAR = {2017},
%    NUMBER = {3},
%     PAGES = {1245--1273},
%      ISSN = {0272-4979},
%   MRCLASS = {65J10 (26A33 65R10)},
%  MRNUMBER = {3671494},
%MRREVIEWER = {C. Ilioi},
%       URL = {https://doi.org/10.1093/imanum/drw042},
%}



