function [ x, it, hess_mult, gp_norms ] = spgqp(A, b, lb, ub, x0, normA, maxit, myeps)
%SPGQP Spectral Projected Gradient method for convex-constrained QP
%
% More details:
%  E. G. Birgin, J. M. Martinez, and M. M. Raydan. Nonmonotone spectral
%   projected gradient methods on convex sets. SIAM Journal on Optimization,
%   10:1196?1211, 2000.
%  L. Pospisil. Development of Algorithms for Solving Minimizing Problems with Convex Qua-
%   dratic Function on Special Convex Sets and Applications. PhD thesis, VSB-TU Ostrava,
%   <a href="matlab:web('http://am.vsb.cz/export/sites/am/cs/theses/pospisil_phd.pdf')">online</a>.
%

% some magic parameters
m = 12; % the size of generalized Armijo rule buffer
gamma = 0.7; % parameter of Armijo rule
sigma2 = 0.9; % safeguarding parameter
sigma3 = 1e4;

% initialize counters
hess_mult = 0;
it = 0;

x = get_projection_box(x0,lb,ub); % x \in \Omega

g = A*x - b; hess_mult = hess_mult + 1;
f = get_function_value( x, g, b);

% initialize Armijo buffer
fs = f*ones(m,1);

% Barzilai-Borwein step-size
alpha_bar = min(0.95/normA,sigma3); % 0 < alpha_bar <= 2*norm(inv(A));

alpha_bb = alpha_bar; % initial BB step

gp = get_projected_gradient(x,g,alpha_bar,lb,ub);

gp_norms(1) = norm(gp);
fss(1) = f;

while it < maxit && my_stopping_criterion(x,norm(gp),myeps)
    d = get_projection_box(x-alpha_bb*g,lb,ub) - x;

    Ad = A*d; hess_mult = hess_mult + 1;
    dd = dot(d,d);
    dAd = dot(Ad,d); % possibly equal to zero
    
    f_max = max(fs);
    
    xi = (f_max - f)/dAd;
    beta_bar = -dot(g,d)/dAd;
    
    if gamma^2*beta_bar^2 + 2*xi < 0
       disp('error in SPGQP')
    end
    beta_hat = gamma*beta_bar + sqrt(gamma^2*beta_bar^2 + 2*xi);
    
    beta = min([sigma2,beta_hat]);
    
    x = x + beta*d;
    g = g + beta*Ad;
    f = get_function_value( x, g, b);
    
    fs(1:end-1) = fs(2:end);
    fs(end) = f;
    
    alpha_bb = min(max(dd/dAd,0),sigma3);
%    alpha_bb = dd/dAd;
    
    gp = get_projected_gradient(x,g,alpha_bar,lb,ub);
    
    it = it + 1;
    gp_norms(it+1) = norm(gp);
    fss(it+1) = f;
    
end

if false
    figure
    subplot(1,2,1)
    title('projected gradient')
    hold on
    plot(1:length(gp_norms),gp_norms)
    xlabel('$it$','Interpreter','latex')
    ylabel('$\Vert g^P_{it} \Vert_2$','Interpreter','latex')
    set(gca,'yscale','log')
    hold off

    subplot(1,2,2)
    hold on
    title('function value')
    hold on
    plot(1:length(fss),fss)
    xlabel('$it$','Interpreter','latex')
    ylabel('$f(x_it)$','Interpreter','latex')
    set(gca,'yscale','log')
    hold off
end


end

function [result] = my_stopping_criterion(x,norm_gp,myeps)
result = (norm_gp > myeps);
end

function gp = get_projected_gradient(x,g,alpha,lb,ub)
gp = (1/alpha)*(x - get_projection_box(x - alpha*g,lb,ub));
%gp = (x - projection_box(x - alpha*g,lb,ub));
end

function [ fx ] = get_function_value( x, g, b)
% compute value of quadratic function using gradient
fx = 1/2*dot(g-b,x);
end
