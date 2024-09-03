%% --------------- Run Simulation ---------------------------------------------

ddtheta_o = @(theta, dtheta, t) pinv(J1fn(theta)) * h1fn(theta,dtheta,t) + ... 
        (eye(3) - pinv(J1fn(theta)) *J1fn(theta)) * pinv(J2fn(theta)) * h2fn(theta,dtheta,t);
%ddtheta_o = @(theta, dtheta, t)pinv(J1fn(theta)) * h1fn(theta,dtheta,t)

%ddtheta_o = @(theta, dtheta, t) pinv(J2fn(theta)) * h2fn(theta,dtheta,t);

%u = @(x,t) Bfn(x(1:3),x(4:6)) + Afn(x(1:3))*[0;0;0];
u = @(x,t) max(-1000*ones(3,1), min(1000*ones(3,1), Bfn(x(1:3),x(4:6)) + Afn(x(1:3))*ddtheta_o(x(1:3),x(4:6), t)));
%u = @(x,t) zeros(3,1);

%eqn = [p2(1); p3(1); p3(2)] == [0; 0; 0.8];
%theta_init = solve(subs(eqn, param_vars, param_vals), q);
%theta_init = double(subs(q, theta_init(1)))

theta_init = [0.6369;2.1418;-1.2078];

theta_0  = theta_init;
dtheta_0 = zeros(3,1);
x0 = [theta_0; dtheta_0];

f_closed_loop = @(t, x)f(x,u(x,t));
fprintf("Running simulation\n");
opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
[tout, yout] = ode45(f_closed_loop, [0:0.01:20], x0, opts);

%% --------------- Plot Results -----------------------------------------------

angles = yout(:,1:3);

l1 = 0.5;
l2 = 0.43;
l3 = 0.35;
lt = l1 + l2 + l3;
pts = zeros(2,4);
figure(1); grid;
set(gcf, 'Position', [100 -200 1000 1000]);

idx = 1;
tic;

% Calculate Desired Path
path_des = zeros(2, length(tout));
inputs = zeros(3, length(tout));
for i = 1:length(tout)
    path_des(:,i) = r1ofn(tout(i));
    inputs(:,i) = u(yout(i,:)',tout(i));
end

for i = 1:length(tout)

end

% Plot system in real time
while true
    t = toc;
    while tout(idx) < t
        idx = idx + 1;
        if idx > size(tout)
            tic;
            idx = 1;
            break;
        end
    end
    t1 = angles(idx, 1);
    t2 = angles(idx, 2);
    t3 = angles(idx, 3);

    pts(:,2) = l1 * [cos(t1); sin(t1)];
    pts(:,3) = pts(:,2) + l2 * [cos(t1+t2); sin(t1+t2)];
    pts(:,4) = pts(:,3) + l3 * [cos(t1+t2+t3); sin(t1+t2+t3)];
    ref = r1ofn(t);
    clf;

    subplot(3, 2, [1:4])
    plot(ref(1), ref(2),'o'); hold on;
    plot(pts(1,:), pts(2,:)); hold on;
    plot(pts(1,:), pts(2,:), 'x'); hold on;
    plot(path_des(1,:), path_des(2,:),'--'); hold on;

    title(sprintf('Time = %.2f, Energy = %.2e',t, Hfn(yout(idx,:)')));
    xlabel('x');
    ylabel('y');
    xlim([-0.2 1.4]);
    ylim([-0.6 1]);
    grid;

    subplot(3, 2, [5 6]);
    plot(tout, inputs(1,:),'r'); hold on;
    plot(tout, inputs(2,:),'g'); hold on;
    plot(tout, inputs(3,:),'b'); hold on;

    in_now =u(yout(idx,:)',tout(idx));
    plot(tout(idx), in_now(1),'ro');
    plot(tout(idx), in_now(2),'go');
    plot(tout(idx), in_now(3),'bo');
    xline(t);

    legend('\tau_1','\tau_2','\tau_3');
    xlabel('t [s]');
    ylabel('\tau [Nm]')

    title('Torques');
    grid on;

    pause(0.05);
end

function pts = forward_kinematics(theta, l1, l2, l3)
    t1 = theta(1);
    t2 = theta(2);
    t3 = theta(3);
    pts = zeros(2,3);
    pts(:,1) = l1 * [cos(t1); sin(t1)];
    pts(:,2) = pts(:,1) + l2 * [cos(t1+t2); sin(t1+t2)];
    pts(:,3) = pts(:,2) + l3 * [cos(t1+t2+t3); sin(t1+t2+t3)];
end
