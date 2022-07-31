% Taylor Howell
% Research-Manchester
% June 19, 2018
% iLQR for control of a pendulum (simple) 
    % Midpoint discretization
    % "Iterative Linear Quadratic Regulator Design for Nonlinear Biological
    % Movement Systems" derivation

clear; clc; close all

% system 

L = [1,2];
r = L./2;
m = [1,2]; %[3,2,1]; %[1,2,3]; %
MI = m.*L.^3./12;
g = 9.81;
n = 4; % dimensions of system
%m = 1 % mass
%l = 1; % length
%g = 9.8; % gravity
%mu = 0.01; % friction coefficient

fc = @(x,u) [omg1;omg2; 
               
     



            (g/l*sin(x(1)) - mu/m/(l^2)*x(2) + 1/m/(l^2)*u)];
        
dynamics_midpoint = @(x,u,dt) x + fc(x + fc(x,u)*dt/2,u)*dt;

A_midpoint = @(x,dt) [(1 + g*cos(x(1))*(dt^2)/(2*l)) (dt - mu*(dt^2)/(2*m*l^2));
                      (g*cos(x(1) + x(2)*dt/2)*dt/l - mu*g*cos(x(1))*(dt^2)/(2*m*l^3)) (1 + g*cos(x(1) + x(2)*dt/2)*(dt^2)/(2*l) - mu*dt/(m*l^2) + (mu^2)*(dt^2)/(2*(m^2)*l^4))];

B_midpoint = @(x,dt) [(dt^2)/(2*m*l^2); 
                      (-mu*(dt^2)/(2*(m^2)*l^4) + dt/(m*l^2))];

% initial conditions
x0 = [0; 0];

% goal
xf = [pi; 0]; % (ie, swing up)

% costs
Q = 1e-5*eye(2);
Qf = 25*eye(2);
R = 1e-5*eye(1);

e_dJ = 1e-12;

% simulation
dt = 0.01;
tf = 1;
N = floor(tf/dt);
t = linspace(0,tf,N);
iterations = 100;

% initialization
u = zeros(1,N-1);
x = zeros(n,N);
x_prev = zeros(n,N);
x(:,1) = x0;

% first roll-out
for k = 2:N
    x(:,k) = dynamics_midpoint(x(:,k-1),u(k-1),dt);
end

% original cost
J = 0;
for k = 1:N-1
    J = J + (x(:,k) - xf)'*Q*(x(:,k) - xf) + u(k)'*R*u(k);
end
disp('Original cost:')
J = 0.5*(J + (x(:,N) - xf)'*Qf*(x(:,N) - xf))

%% iterations of iLQR using Todorov derivation
for i = 1:iterations
    %% Backward pass
    S = zeros(n,n,N);
    v = zeros(n,1,N);
    S(:,:,N) = Qf;
    v(:,N) = Qf*(x(:,N) - xf);
    K = zeros(1,n,N-1);
    Kv = zeros(1,n,N-1);
    Ku = zeros(1,N-1);
    for k = N-1:-1:1
        K(1,:,k) = (B_midpoint(x(:,k),dt)'*S(:,:,k+1)*B_midpoint(x(:,k),dt) + R)\B_midpoint(x(:,k),dt)'*S(:,:,k+1)*A_midpoint(x(:,k),dt);
        Kv(1,:,k) = (B_midpoint(x(:,k),dt)'*S(:,:,k+1)*B_midpoint(x(:,k),dt) + R)\B_midpoint(x(:,k),dt)';
        Ku(k) = (B_midpoint(x(:,k),dt)'*S(:,:,k+1)*B_midpoint(x(:,k),dt) + R)\R;
        S(:,:,k) = A_midpoint(x(:,k),dt)'*S(:,:,k+1)*(A_midpoint(x(:,k),dt) - B_midpoint(x(:,k),dt)*K(1,:,k)) + Q;
        v(:,k) = (A_midpoint(x(:,k),dt) - B_midpoint(x(:,k),dt)*K(1,:,k))'*v(:,k+1) - K(1,:,k)'*R*u(k) + Q*x(:,k);
    end
    
    % update control, roll out new policy, calculate new cost
    x_prev = x;
    J_prev = J;
    J = Inf;
    alpha = 1;
    iter = 0;
    while J > J_prev 
        x = zeros(n,N);
        x(:,1) = x0;
        for k = 2:N
            u_(k-1) = u(k-1) -K(1,:,k-1)*(x(:,k-1) - x_prev(:,k-1)) + alpha*(-Kv(1,:,k-1)*v(:,k) - Ku(k-1)*u(k-1));
            x(:,k) = dynamics_midpoint(x(:,k-1),u_(k-1),dt);
        end

        J = 0;
        for k = 1:N-1
            J = J + (x(:,k) - xf)'*Q*(x(:,k) - xf) + u_(k)'*R*u_(k);
        end
        J = 0.5*(J + (x(:,N) - xf)'*Qf*(x(:,N) - xf));
        alpha = alpha/2;
        iter = iter + 1;
    end
    disp('New cost:')
    J
    u = u_;
    
    if abs(J - J_prev) < e_dJ
        disp(strcat('eps criteria met at iteration: ',num2str(i)))
        break
    end
end


%% Results

% Animation
r = 1;
figure
U = [u(1) u];
for i = 1:N
    p1 = subplot(1,2,1);
    X = r*cos(x(1,i) - pi/2);
    Y = r*sin(x(1,i) - pi/2);
    plot([0 X],[0 Y],'k-')
    hold on
    plot(X,Y,'ko','MarkerFaceColor', 'k')
    xlabel('pendulum (simple)')
    axis([-1.5*r 1.5*r -1.5*r 1.5*r])
    axis square
    
    p2 = subplot(1,2,2);
    stairs(t(1:i),U(1:i))
    xlabel('t')
    ylabel('u(t)')
    axis([0 tf min(u) max(u)])
    axis square

    drawnow
    %pause(dt)

    if i ~= N
        cla(p1);
        cla(p2);
    end
     
end

figure
hold on
plot(linspace(0,tf,N),x(1,:))
plot(linspace(0,tf,N),x(2,:))
legend('\theta','\theta_{d}')
