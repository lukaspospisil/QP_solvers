clear all

addpath(genpath('fem/'));
addpath(genpath('plot/'));
addpath(genpath('problem/'));
addpath(genpath('solver/'));

p = 0.5; % parameter of constrain
m = 40; % discretization parameter (size of problem = m^2)

% ------------------------------------------
% CREATE THE PROBLEM
% ------------------------------------------
% discretize the problem
[nodes,edges,idxD,idxN,valuesD,valuesN] = membrane_discretization(1,1,m);
[A,b] = fem2d(nodes,edges,idxD,idxN,valuesD,valuesN);

% get lower bound
l = membrane_get_l( m, p );

n = length(b); % problem dimension
disp(['n = ' num2str(n)]);

% algorithms setting
maxit = 1e4; % max number of iterations
myeps = 1e-4; % precision
u = []; % upper bound
x0 = zeros(size(b)); % initial approximation
normA = gersgorin(A); % estimation of max eigenvalue


% ------------------------------------------
% SOLVE THE PROBLEM
% ------------------------------------------

% matlab solver
%options.Algorithm = 'interior-point-convex';
%options.Display = 'none';
%x = quadprog(A,-b,[],[],B,c,l,[],zeros(size(b)),options);

% mprgp
[x1, it1, hess_mult1, gp_norms1] = mprgp(A, b, l, u, x0, normA, maxit, myeps);
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
[x2, it2, hess_mult2, gp_norms2] = spgqp(A, b, l, u, x0, normA, maxit, myeps);
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
membrane_draw_solution(nodes,edges,x2);


