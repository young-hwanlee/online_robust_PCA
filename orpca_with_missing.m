% Input:
%    Z = Data matrix with data points as columns
%    N = no. of rows
%    T = no. of columns (iterations)
%    L_prev = initial matrix for the dictionary (basis of the low-dimensional subspace)
%    Omega = missing entries (0: missing, 1: observed)
%    lambda1 = constant for the low-rank component X
%    lambda2 = constant for the sparse component E
%    epsilon = tolerance for convergence
%    r = max. rank
%    rho = step size for SGD for L
%    sgd = use the SGD or not (1st order or 2nd order)
% Output:
%     L = basis of the low-dimensional subspace
%     R = coefficients of the samples w.r.t. the basis L
%     E = sparse noise matrix
%     Delta_Lt = norm of the difference between L at t-1 and L at t
%     fx = instantaneous objective function
%     fx_accumulate = accumulated fx

function [L,R,E] = ...
    orpca_with_missing(Z,N,T,L_prev,Omega,lambda1,lambda2,epsilon,r,rho,sgd)

Omega_c = 1 - Omega;

if sgd == 0
    Lt = L_prev;
    A = zeros(r,r);
    B = zeros(N,r);
end

R = zeros(T,r);
E = zeros(N,T);

Delta_Lt = zeros(1,T);
fx = zeros(1,T);
fx_accumulate = zeros(1,T);

for t = 1:T    
    
    if mod(t,5e3) == 0
        t
    end
    
    zt = Z(:,t);

    [r_traj,e_traj] = data_proj(zt,L_prev,lambda1,lambda2,epsilon,r,Omega(:,t),Omega_c(:,t));
    
    i_last = size(r_traj,2);

    rt = r_traj(:,i_last);
    et = e_traj(:,i_last);
    R(t,:) = rt';
    E(:,t) = et;

%%
    % % update rule for Lt used in [Feng,Xu,Yan,'13] (2nd order)
    if sgd == 0
    % % classic
    % -----------------------------------------------------------------------------------
        A = A + rt*rt';
        B = B + (zt - et)*rt';
    % -----------------------------------------------------------------------------------

    % % scaling the past data (forgetting factor)
    % -----------------------------------------------------------------------------------
    %     rho = 170;
    %     beta = (1-1/t)^rho;
    %     A = beta*A + rt*rt';
    %     B = beta*B + (zt - et)*rt';
    % -----------------------------------------------------------------------------------

        Atilde = A + lambda1*eye(r,r);

        for i = 1:r
            Lt(:,i) = (B(:,i) - Lt*Atilde(:,i))/Atilde(i,i) + Lt(:,i);
        end

    % % update Lt via SGD (1st order)
    elseif sgd == 1
        Lt = L_prev + rho*((zt - L_prev*rt - et)*rt' - lambda1/T*L_prev);
    end
    
    Delta_Lt(t) = norm(Lt-L_prev,'fro');    % change in Lt magnitude
    fx(t) = 0.5*norm(zt-Lt*rt-et)^2 + 0.5*lambda1*(norm(Lt,'fro')^2/T + norm(rt)^2) + lambda2*norm(Omega(:,t).*et,1);
        
    if t == 1
        fx_accumulate(t) = fx(t);
    else
        fx_accumulate(t) = fx_accumulate(t-1) + fx(t);
    end

    L_prev = Lt;

end

L = Lt;
E = E.*Omega;



% min 0.5||zt - L_{t-1} r - e||_2^2 + 0.5*lambda1*||r||_2^2 + lambda2*||e||_1
% using algorithm in Supp. material for [Feng,Xu,Yan,'13]
function [R,E] = data_proj(z,L,lambda1,lambda2,epsilon,r,omega,omega_c)

N = length(z);
e = zeros(N,1);
conv_metric = epsilon + 999;

r_prev = zeros(r,1);
e_prev = zeros(N,1);
k = 0;

R = [];
E = [];

while conv_metric >= epsilon
    
    k = k + 1;

    rk = (L'*L + lambda1*eye(r,r))\L'*(z-e_prev);
    
    tmp = z - L*rk;
%     e = soft_thresh(tmp, lambda2);       # original method in [Feng,Xu,Yan,'13]
    e = soft_thresh(tmp,lambda2).*omega - L*rk.*omega_c;
    
    conv_metric = max([norm(rk - r_prev), norm(e - e_prev)])/norm(z);     % original method in [Feng,Xu,Yan,'13]
%     conv_metric = max([norm(rk - r_prev), norm(e - e_prev)]);       % alternative method
    
    r_prev = rk;
    e_prev = e;
    
    R(:,k) = rk;
    E(:,k) = e;
end



function y = soft_thresh(x,lambda)

y = max(abs(x) - lambda,0).*sign(x);
