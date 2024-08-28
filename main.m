%% --------------- System Model -----------------------------------------------
% parameters
syms l1 l2 l3 m1 m2 m3 g 'real'

param_vars = [  l1   l2   l3 m1 m2 m3    g];
param_vals = [0.50 0.43 0.35 30 25 20 9.81];

% state variables
syms theta1 theta2 theta3 dtheta1 dtheta2 dtheta3 'real'

% inputs
T = sym('T',[3 1],'real');

% Lagrange formulation
q  = [ theta1  theta2  theta3]';
dq = [dtheta1 dtheta2 dtheta3]';
u  = T;

% Kinetic and Potential energy
p1 =  0 + l1 * [cos(              theta1); sin(              theta1)];
p2 = p1 + l2 * [cos(       theta1+theta2); sin(       theta1+theta2)];
p3 = p2 + l3 * [cos(theta1+theta2+theta3); sin(theta1+theta2+theta3)];

v1 = jacobian(p1, q) * dq;
v2 = jacobian(p2, q) * dq;
v3 = jacobian(p3, q) * dq;

% Kinetic energy
K = 0.5 * m1 * v1'*v1 + 0.5 * m2 * v2'*v2 + 0.5 * m3 * v3'*v3;

% Potential energy, center off mass in the middle of each link.
pm1 = 0.5 * p1;
pm2 = 0.5 * (p1 + p2);
pm3 = 0.5 * (p2 + p3);
U = g * (m1 * pm1(2) + m2 * pm2(2) + m3 * pm3(2) );
% Hamiltonian
H = K + U;
Hfn = simplify(subs(H, param_vars, param_vals));
Hfn = matlabFunction(Hfn, 'Vars', {[q;dq]});

DK_Ddq     = simplify( jacobian(K, dq)' );
DK_Dq      = simplify( jacobian(K, q)' );
D2K_DqDdq  = simplify( jacobian(DK_Ddq, q)' );
D2K_D2dq   = simplify( jacobian(DK_Ddq, dq)' );
DU_Dq      = simplify( jacobian(U, q)' );

f = (D2K_D2dq) \ (u - D2K_DqDdq' * dq - DU_Dq + DK_Dq);
f = [dq; simplify(f)];
x = [q; dq];

f = simplify(subs(f, param_vars, param_vals));

f = matlabFunction(f, 'Vars', {x, u});
