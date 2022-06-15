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

Bcl = [B;
    zeros(size(B))];

F = [Dc*C, Cc];

% initial condition
x0 = [0; 0; 0];

W = sdpvar(size(Acl,1), size(Acl,2), 'symmetric', 'real');
Z = sdpvar(size(B,2), size(Acl,2), 'full', 'real');

LMI1 = W - 1e-5*eye(size(W)) >= 0;

% polytope
LMIpoly = [];
m = size(B,2);

D1 = eye(m);
D1_ = eye(m) - D1;
D2 = [1, 0;
    0, 0];
D2_ = eye(m) - D2;
D3 = [0, 0;
    0, 1];
D3_ = eye(m) - D3;
D4 = zeros(m,m);
D4_ = eye(m) - D4;

Dp = [D1, D2, D3, D4];
Dp_ = [D1_, D2_, D3_, D4_];

for i = 1:m:size(Dp,2)
    LMIi = (Acl - Bcl*F + Bcl*Dp(:, i:i + 1)*F)*W + Bcl*Dp_(:, i:i + 1)*Z + ((Acl - Bcl*F + Bcl*Dp(:, i:i + 1)*F)*W + Bcl*Dp_(:, i:i + 1)*Z)' + 1e-5*eye(size(Acl)) <= 0;
    LMIpoly = [LMIpoly, LMIi];
end

% ellipsoid
LMIellip = [];
for i = 1:m
    LMIi = [W, Z(i,:)';
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

W = value(W); % W = inv(P)
Z = value(Z); % Z = H*inv(P)

H = Z/W;

return