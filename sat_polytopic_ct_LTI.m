clear
close all
clc

warning off

% system
A = [-1, 0, 1;
    0, 2, -1;
    2, 0, -3];

B = [1;
    1;
    0];

C = eye(size(A));

D = zeros(size(C,1),size(B,2));

% initial condition
x0 = [2; 0.65; 0];

X = sdpvar(size(A,1), size(A,2), 'symmetric', 'real');
Y = sdpvar(size(B,2), size(A,2), 'full', 'real');
Z = sdpvar(size(B,2), size(A,2), 'full', 'real');

LMI1 = X - 1e-5*eye(size(X)) >= 0;

% polytope
LMIpoly = [];
m = size(B,2);

D1 = 1;
D1_ = 1 - D1;
D2 = 0;
D2_ = 1 - D2;

Dp = [D1, D2];
Dp_ = [D1_, D2_];

for i = 1:2^m
    LMIi = (X*A' + A*X + B*Dp(i)*Y + + B*Dp_(i)*Z + Y'*Dp(i)'*B' + Z'*Dp_(i)'*B') + 1e-5*eye(size(A)) <= 0;
    LMIpoly = [LMIpoly, LMIi];
end

% ellipsoid
LMIellip = [];
for i = 1:m
    LMIi = [X, Z(i,:)';
        Z(i,:), 1] >= 0;
    LMIellip = [LMIellip, LMIi];
end

LMIs = [LMI1, LMIpoly, LMIellip];

options = sdpsettings('solver','sedumi','verbose',1);
sol = optimize(LMIs, [], options);

if sol.problem == 0
    [primal,~] = check(LMIs);
    if (min(primal) >= 0 && all(primal(1) > 0))
        disp('Sucessfully solved LMIs without problems'); 
    else
        disp('LMIs not solved');
    end
else
    [primal,~] = check(LMIs);
    if (min(primal) >= 0 && all(primal(1) > 0))
        disp(['Sucessfully solved LMIs, but solver acused ' yalmiperror(sol.problem)]);
    else
        disp(['LMIs not solved. Solver acused ' yalmiperror(sol.problem)]);
    end
end

X = value(X); % X = inv(P)
Y = value(Y); % Y = F*inv(P)
Z = value(Z); % Z = H*inv(P)

H = Z/X;

% feedback gain
K = Y/X;

Tsim = 15;
sim = sim('sim_ssf_sat');

figure
plot(sim.x.time, sim.x.signals.values)
axis tight

figure
plot(sim.u.time, sim.u.signals.values)
axis tight

return