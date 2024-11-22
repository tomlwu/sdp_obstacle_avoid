function solution = rank_min_mip_obstacle_avoid(A_ieq,A_eq,b_ieq,b_eq,objective_fun,n_obs,d,options)
%RANK_MIN_MIP_OBSTACLE_AVOID Summary of this function goes here
%   Detailed explanation goes here
%% settings
maxIter = 100;
tol_sigma = 1e-5;
tol_e = -1e-5;
tol_U = 1e-7;
trace_fixed_Y = 2;
flag_alter_rc = false; % a flag that is true if we want to alter steps between rank and cost
if isfield(options,'alterrankcost')
    flag_alter_rc = options.alterrankcost;
end

sqstd = @(Y) (kron(eye(4),[0,0,1,0])*vec(Y)+1/2)'*(kron(eye(4),[0,0,1,0])*vec(Y)+1/2); % squared standard deviation of theta
%% initialize
cvx_begin sdp
cvx_precision high
cvx_solver mosek_3

variable x(n_obs*d,1)
variable Y(2,2,n_obs) symmetric
% minimize(objective_fun(x,Y));
% minimize(sum_max_eig_3d(Y));
subject to
for i = 1:n_obs
    Y(:,:,i) >= 0;
    Y(1,1,i) == 1;
    Y(2,2,i) == 1;
end
A_ieq*[x;vec(Y)]<=b_ieq;
A_eq*[x;vec(Y)]==b_eq;

cvx_end
x = full(x);
Y = full(Y);
solution.xInit = x;
solution.YInit = Y;
flag_solver_fail = false;
%% rank minimization
max_eigen_all = zeros(maxIter,n_obs);
max_eigen_all(1,:) = max_eigen(Y);
obj_cost = zeros(maxIter,1);
obj_cost(1) = objective_fun(x,Y);
standard_diviation_all = zeros(maxIter,1);
grad_sigma_1_Y=grad_singularvalue(Y,ones(n_obs,1),'psd');
index_large_rank_Y = abs(max_eigen(Y)-trace_fixed_Y)>=tol_sigma;
Ux = 100*ones(n_obs*d,1);
Uy = 100*repmat(eye(2),[1,1,n_obs]);
k= 1;
if flag_alter_rc
    alter_status = 0; % status flag for alternation: 0 for rank, 1 for cost minimization
end

while k<maxIter && (any(index_large_rank_Y) && (norm(Ux,"fro")>=tol_U || norm(Uy,"fro")>=tol_U))

    cvx_begin sdp
    %     cvx_precision high
    cvx_solver mosek_3

    variable Ux(n_obs*d,1)
    variable Uy(2,2,n_obs) symmetric
    variable c1
    if flag_alter_rc
        if alter_status == 0
            minimize(c1);%objective_fun(M+U)+
        elseif alter_status == 1
            minimize(objective_fun(x+Ux,Y+Uy));
        end
    else
        minimize(objective_fun(x+Ux,Y+Uy)+n_obs*c1);%
    end
    subject to
    vec(Uy)'*vec(grad_sigma_1_Y)>=(c1-1).*(sum(max_eigen_all(k,:))-trace_fixed_Y*n_obs)+tol_e;
%     obj*[Ux;vec(Uy)]<=(c2-1)*objective_fun(x,Y);
    for i = 1:n_obs
        Y(:,:,i)+Uy(:,:,i) >= 0;
        Uy(1,1,i) == 0;
        Uy(2,2,i) == 0;
    end
    c1 >= 0;
    if flag_alter_rc
        if alter_status == 0
            c1 <= 1;
        elseif alter_status == 1
            c1 <= 2;
        end
    else
        c1 <= 1;
    end
    A_ieq*[x+Ux;vec(Y+Uy)]<=b_ieq;
    A_eq*[Ux;vec(Uy)]==zeros(size(b_eq));

    cvx_end
    if strcmp(cvx_status,'Failed') || strcmp(cvx_status,'Infeasible') || strcmp(cvx_status,'Unbounded') || strcmp(cvx_status,'Overdetermined')
        flag_solver_fail = true;
        break
    end
    k = k+1;
    x = x+Ux;
    Y = Y+Uy;
    if flag_alter_rc
        alter_status = mod(alter_status+1,2);
    end
    grad_sigma_1_Y=grad_singularvalue(Y,ones(n_obs,1),'psd');
    index_large_rank_Y = abs(max_eigen(Y)-trace_fixed_Y)>=tol_sigma;
    max_eigen_all(k,:) = max_eigen(Y);
    obj_cost(k) = objective_fun(x,Y);
    standard_diviation_all(k) = sqstd(Y);
end
if flag_solver_fail
    disp('solver failed')
    exit_flag = 4;
else
    if k>=maxIter
        disp('max iteration reached!')
        exit_flag = 3;
    elseif ~any(index_large_rank_Y) %&& obj_cost(k)<tol_cost
        disp('Low rank solution found!')
        exit_flag = 1;
    elseif norm(Ux,"fro")<tol_U || norm(Uy,"fro")<tol_U
        disp('Stopped with small improvement!')
        exit_flag = 2;
    end
end
solution.x = full(x);
solution.Y = full(Y);
solution.max_eigen_all = max_eigen_all(1:k,:);
solution.obj_cost = obj_cost(1:k,:);
solution.exitFlag = exit_flag;
solution.standard_diviation = standard_diviation_all(1:k);
end

function sum_eig_max = sum_max_eig_3d(M)
% sum of largest eigenvalues
sum_eig_max = 0;
for i=1:size(M,3)
    sum_eig_max = sum_eig_max+lambda_max(M(:,:,i));
end
end
