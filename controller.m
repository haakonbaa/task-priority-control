tic
% A(theta) and B(theta, dtheta) eq 8, Nakamura 1987
%   tau = A(theta) * ddtheta + B(theta,dtheta);
A = simplify(D2K_D2dq);
B = simplify(D2K_DqDdq' * dq + DU_Dq + DK_Dq);

Afn = matlabFunction(subs(A, param_vars, param_vals), 'Vars', {q});
Bfn = matlabFunction(subs(B, param_vars, param_vals), 'Vars', {q, dq});

% taskts
r1 = p3;
r2 = cos(theta1 + theta2 + theta3);

dr1 = simplify(jacobian(r1, q) * dq);
dr2 = simplify(jacobian(r2, q) * dq);

t = sym('t', 'real');

% reference tasks
r1o = [0.5*(1+sin(t-pi/2)); 0.8];

T = 20;

sigmoid = @(x) 1./(1+exp(-x));
alpha = @(x) 1-sigmoid(10*(x-0.5*T));

alpha2 = 1-0.5*(1+sin((pi/(0.5*T))*t-pi/2));
r1o1 = (1-alpha2)*[1; 0.8] + alpha2*[  0;  0.8];
r1o2 = (1-alpha2)*[1; 0.8] + alpha2*[0.5;-0.15];
r1o = r1o1*alpha(t) + r1o2*(1-alpha(t));


r2o = cos(pi/2);
dr1o  = diff( r1o, t);
ddr1o = diff(dr1o, t);
dr2o  = diff( r2o, t);
ddr2o = diff(dr2o, t);

J1 = simplify(jacobian(r1, q));
J2 = simplify(jacobian(r2, q));

dJ1 = matrix_dot('dJ1', J1, q, dq);
dJ2 = matrix_dot('dJ2', J2, q, dq);

J1pinv = simplify(pinv(J1));
J2pinv = simplify(pinv(J2));

G11 =  10;
G21 =  10;
G12 =  20;
G22 = 100;
h1 = simplify(- dJ1 * dq + ddr1o + G11*(dr1o - dr1) + G21*(r1o - r1));
h2 = simplify(- dJ2 * dq + ddr2o + G12*(dr2o - dr2) + G22*(r2o - r2));

r1ofn = matlabFunction(subs(r1o, param_vars, param_vals), 'Vars', {t});
h1fn = matlabFunction(subs(h1, param_vars, param_vals), 'Vars', {q, dq, t});
h2fn = matlabFunction(subs(h2, param_vars, param_vals), 'Vars', {q, dq, t});
J1fn = matlabFunction(subs(J1, param_vars, param_vals), 'Vars', {q});
J2fn = matlabFunction(subs(J2, param_vars, param_vals), 'Vars', {q});

ddtheta_o = J1pinv * h1 + (eye(3) - J1pinv*J1)*J2pinv * h2;
ddtheta_ofn = matlabFunction(subs(ddtheta_o, param_vars, param_vals), 'Vars', {q, dq, t});

toc
% matrix_dot differentiates a symbolic matrix and returns a new symbolic
% variable that is the time derivative of that matrix. 
% inputs
%   name: name of new symbolic variable
%      M: matrix to be differentiated
%   vars: M is a function of 'vars', 'vars' are functions of time
%  dvars: the time derivative of 'vars'
function M_dot = matrix_dot(name, M, vars, dvars)
    M_dot = sym(name, size(M), 'real');
    for i = 1:size(M, 1)
        for j = 1:size(M, 2)
            M_dot(i, j) = jacobian(M(i, j), vars) * dvars;
        end
    end
    M_dot = simplify(M_dot);
end
