clear all

addpath('plot/')
addpath('solver/')

p = 0.5; % parameter of constrain
n = 80; % number of nodes

% ------------------------------------------
% CREATE THE PROBLEM
% ------------------------------------------
% discretize the problem
h = 1/(n-1);

A = 2*diag(ones(n,1)) - diag(ones(n-1,1),1) - diag(ones(n-1,1),-1);
%b = -15*ones(n,1);
b = 50*[ones(n/2,1); -ones(n/2,1)];

% apply scaling before dirichlet
A = 1/h*A;
b = h*b;
 
% apply dirichlet - change system so that x(1)=x(N)=0
A(1,1) = 1;
A(n,n) = 1;
A(1,2) = 0;
A(2,1) = 0;
A(n-1,n) = 0;
A(n,n-1) = 0;
b(1) = 0;
b(n) = 0;

A = sparse(A); % bleh!

lb = -ones(n,1);
ub = ones(n,1);

disp(['n = ' num2str(n)]);
    
% ------------------------------------------
% SOLVE THE PROBLEM
% ------------------------------------------

% algorithm settings
x0 = zeros(size(b));
maxit = 1e4;
myeps = 1e-4;
normA = gersgorin(A);

% mprgp
[x1, it1, hess_mult1, gp_norms1] = mprgp(A, b, lb, ub, x0, normA, maxit, myeps);
if true
    % plot algorithm performance
    figure
    hold on
    title('MPRGP - projected gradient norm')
    plot(1:length(gp_norms1),gp_norms1,'r')
    xlabel('$it$','Interpreter','latex')
    ylabel('$\Vert g^P(x_{it}) \Vert_2$','Interpreter','latex')
    
    hold off
end

% spgqp
[x2, it2, hess_mult2, gp_norms2] = spgqp(A, b, lb, ub, x0, normA, maxit, myeps);
if true
    % plot algorithm performance
    figure
    hold on
    title('SPGQP - projected gradient norm')
    plot(1:length(gp_norms2),gp_norms2,'r')
    xlabel('$it$','Interpreter','latex')
    ylabel('$\Vert g^P(x_{it}) \Vert_2$','Interpreter','latex')
    
    hold off
end

% ------------------------------------------
% PLOT THE SOLUTION
% ------------------------------------------
figure
hold on
plot(1:length(x1),x1,'b')
plot(1:length(lb),lb,'r')
plot(1:length(ub),ub,'r')
hold off

