clear
clc
close all

%% build the world
% build convex polygons
pCenters = [-2, 0, 2, 0;
            0, 2, 0, -2];
r = 1.75;
r_b = 0.2;
nc = 8; % number of sides for each obstacle
Aps = zeros(8,2,4);
bps = zeros(8,4);
A_norms = zeros(8,4);
desired_direction = [-1;2];
for i = 1:4
    V = circle(pCenters(:,i),r,'n',nc);
%     Aps(:,:,i) = V';
%     bps(:,i) = zeros(size(V,2),1);%r^2.*ones(size(V,2),1);
    [AV,bV] = polygonConstraints(V);
    Aps(:,:,i) = AV;
    bps(:,i) = bV;
    for j = 1:nc
        A_norms(j,i) = norm(Aps(j,:,i),1);
    end
end

A_ieq = [blkdiagfrom3Dmat(Aps), blkdiagfrom3Dmat(permute(A_norms.*r_b-bps,[1,3,2]))];
b_ieq = zeros(size(A_ieq,1),1);
A_eq = [zeros(1,2*4),ones(1,4)];
b_eq = 1;

% Q = blkdiag(eye(2*4),zeros(4));
% obj = zeros(1,2*4+4);
obj = -desired_direction'*[repmat(eye(2),[1,4]),zeros(2,4)];%[ones(1,4*2),zeros(1,4)];

%% Solve with a mip solver
model.A = sparse([A_ieq;A_eq]);
% model.Q = sparse(Q);
model.obj = obj;
model.rhs = [b_ieq;b_eq];
model.sense = [repmat('<',size(A_ieq,1),1);repmat('=',size(A_eq,1),1)];
model.vtype = [repmat('C',[2*4,1]);repmat('B',[4,1])];
model.modelsense = 'min';
model.lb = -10*ones(size(obj,2),1);
model.ub = 10*ones(size(obj,2),1);

% gurobi_write(model, 'mip1.lp');

params.outputflag = 1;

result = gurobi(model, params);

disp(result);

yVals = zeros(2,4);
for i = 1:4
    yVals(:,i) = result.x((i-1)*2+1:i*2);
end
p_mip = sum(yVals,2);
%% Solve with relaxed SDP
M_trans_i = [0 0 1/2 0];
N_trans_i = 1/2;
obj2 = -desired_direction'*[repmat(eye(2),[1,4]),zeros(2,4*4)];%[ones(1,4*2),zeros(1,4*4)];
objective_fun = @(x,Y) obj2*[x;vec(Y)];%x'*x;

A_ieq2 = [blkdiagfrom3Dmat(Aps), blkdiagfrom3Dmat(matmult3D(permute(A_norms.*r_b-bps,[1,3,2]),repmat(M_trans_i,[1,1,4])))];
b_ieq2 = vec(matmult3D(permute(A_norms.*r_b-bps,[1,3,2]),repmat(-N_trans_i,[1,1,4])));
A_eq2 = [zeros(1,2*4),kron(ones(1,4),M_trans_i)];
b_eq2 = 1-4/2;
options.alterrankcost = false;

solution = rank_min_sdp_obstacle_avoid(A_ieq2,A_eq2,b_ieq2,b_eq2,objective_fun,4,2,options);
p_sdp = sum(reshape(solution.x,2,[]),2);

Y_one = [1,1;1,1];
Y_zero = [1,-1;-1,1];
test_variable = [0;0;0;1;0;0;0;0;vec(Y_zero);vec(Y_one);vec(Y_zero);vec(Y_zero)];
%% plot
figure
for i = 1:4
    draw2dcircle(pCenters(:,i),r,'n',nc);
    hold on
end

h1 = draw2dcircle(p_mip,r_b,'r');
plot(p_mip(1),p_mip(2),'rx')
h2 = draw2dcircle(p_sdp,r_b,'b--');
plot(p_sdp(1),p_sdp(2),'b+')
plotArrows([0;0],desired_direction,'g-')
hold off
axis equal
grid on
legend([h1,h2],{'MIP','SDP'})

figure
plot(1:size(solution.max_eigen_all,1),solution.max_eigen_all)
figure
plot(1:size(solution.obj_cost,1),solution.obj_cost)
% figure
% plot(1:size(solution.standard_diviation,1),solution.standard_diviation)