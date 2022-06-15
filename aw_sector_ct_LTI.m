clear
close all
clc

warning off

% system
A = [0.1, - 0.1;
    0.1, - 3];

B = [5, 0;
    0, 1];

C = eye(size(A));

D = zeros(size(C,1),size(B,2));

% stabilizing controller
Ac = [- 171.2, 27.2;
    - 68, - 626.8];

Bc = [- 598.2, 5.539;
    - 4.567, 149.8];

Cc = [0.146, 0.088;
    - 6.821, - 5.67];

Dc = zeros(size(B,2), size(Ac,2));

% closed-loop
Acl = [A + B*Dc*C, B*Cc;
    Bc*C, Ac];

F = [Dc*C, Cc];

% initial condition
x0 = [0; 0; 0];

W = sdpvar(size(Acl,1), size(Acl,2), 'symmetric', 'real');
R = sdpvar(size(B,2), size(B,2), 'diag', 'real');
Y = sdpvar(size(B,2), size(Acl,2), 'full', 'real');
Z = sdpvar(size(Ac,1), size(B,2), 'full', 'real');

LMI1 = W - 1e-5*eye(size(W)) >= 0;

LMI2 = R - 1e-5*eye(size(R)) >= 0;

LMI3 = [W*Acl' + Acl*W, [B; zeros(size(B))]*R + [zeros(size(R)); eye(size(R))]*Z - Y';
    R*[B; zeros(size(B))]' + Z'*[zeros(size(R)); eye(size(R))]' - Y, - 2*R] + 1e-5*eye(size(Acl,1) + size(R,1)) <= 0;

% ellipsoid
m = size(B,2);

LMIellip = [];
for i = 1:m
    LMIi = [W, W*F(i,:)' - Y(i,:)';
        F(i,:)*W - Y(i,:), 1] >= 0;
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

W = value(W); % W = inv(P)
R = value(R); % R = inv(T)
Y = value(Y); % Y = G*inv(P)
Z = value(Z); % Z = Ec*inv(T)

G = Y/W;

% aw gain
Ec = Z/R;

return