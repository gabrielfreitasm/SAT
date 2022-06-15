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
Y = sdpvar(size(B,2), size(B,2), 'diag', 'real');
Z = sdpvar(size(B,2), size(A,2), 'full', 'real');
W = sdpvar(size(B,2), size(A,2), 'full', 'real');

LMI1 = X - 1e-5*eye(size(X)) >= 0;

LMI2 = Y - 1e-5*eye(size(Y)) >= 0;

LMI3 = [X*A' + W'*B' + A*X + B*W, B*Y - Z';
    Y*B' - Z, -2*Y] + 1e-5*eye(size(A,1) + size(Y,1)) <= 0;

% ellipsoid
m = size(B,2);

LMIellip = [];
for i = 1:m
    LMIi = [X, (W(i,:) - Z(i,:))';
        W(i,:) - Z(i,:), 1] >= 0;
    LMIellip = [LMIellip, LMIi];
end

LMIs = [LMI1, LMI2, LMI3, LMIellip];

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
Y = value(Y); % Y = inv(R)
Z = value(Z); % Z = G*inv(P)
W = value(W); % W = K*inv(P)

G = Z/X;

% feedback gain
K = W/X;

Tsim = 15;
sim = sim('sim_ssf_sat');

figure
plot(sim.x.time, sim.x.signals.values)
axis tight

figure
plot(sim.u.time, sim.u.signals.values)
axis tight

return