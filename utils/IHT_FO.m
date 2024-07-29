%% Author : TAO ZHANG  * zt1996nic@gmail.com *
% Created Time : 2023-07-20 00:00
% Last Revised : TAO ZHANG ,2023-08-01 00:00
% Remark : Iterative hard thresholding to solve the Fractional Chaotic System

function [Xi,q,outputs,outoptions] = IHT_FO(X, h, t0, tfinal, Theta, options)

if nargin < 3 || isempty(options) 
    options.maxit = 100;
    options.del = 100;
    options.dq = 0.001;
    options.q0 = 0.9; 
    options.lambda = 0.05;   
end

options = setDefaultField(options, 'maxit', 100);
options = setDefaultField(options, 'del', 10);
options = setDefaultField(options, 'dq', 0.01);
options = setDefaultField(options, 'q0', 0.9);
options = setDefaultField(options, 'lambda', 0.05);
maxit = options.maxit;
del = options.del;
dq = options.dq;
q0 = options.q0;
lambda = options.lambda;

% Main Algorithm
[N_X,Dim] = size(X);
iter = 0;
q=q0;
Dalpha1 = ban(q,N_X,h);
for i=1:Dim
    D_qx(:,i) = Dalpha1 * X(:,i);
end
Xi = Theta(del:end-del,:)\D_qx(del:end-del,:);
while ((iter <= maxit) && (q <= 1))
    iter = iter+1;
    
    Xi = Theta(del:end-del,:)\D_qx(del:end-del,:);
    for iter_inner = 1:10 
        smallinds = (abs(Xi)<lambda); 
        Xi(smallinds) = 0; 
        for ind = 1:size(D_qx,2)
            biginds = ~smallinds(:,ind);
            Xi(biginds,ind) = Theta(del:end-del,biginds)\D_qx(del:end-del,ind);
        end
    end
    
    RXi = Xi(:);
    [~, SY] = fde12(q,@(t, y)sparse_model(t, y, RXi),t0,tfinal,X(1,:)',h);  
    
    error_x = norm(X - SY', 'fro')/norm(X , 'fro');
    Ob2(iter)=error_x;
    Rerror = norm(D_qx(del:end-del,:) - Theta(del:end-del,:)*Xi, 'fro')/norm(D_qx(del:end-del,:), 'fro'); 
%     Rerror = norm(D_qx(del:end-del,:) - Theta(del:end-del,:)*Xi, 'fro');  
    rel_err(iter) = Rerror;
    
    Xi_P{iter} = Xi;
    q_i(iter)=q;
    Ob1(iter)= nnz(Xi);
    

    q=q+dq;
    Dalpha1 = ban(q,N_X,h);
    for i=1:Dim
        D_qx(:,i) = Dalpha1 * X(:,i);
    end
end

%  Outputs and options
outputs.iter = iter;
outputs.rerror_model = rel_err;
outputs.Xi = Xi_P;
outputs.q = q_i;
outputs.ObFu1 = Ob1;
outputs.ObFu2 = Ob2;
outputs.ObFu3 = rel_err;

outoptions.maxit = maxit;
outoptions.coef_thres = lambda;
outoptions.q0 = q0;
outoptions.dq = dq;
outoptions.del = del;

end

function s = setDefaultField(s, field, defaultValue)
    if ~isfield(s, field)
        s.(field) = defaultValue;
    end
end