%% Inverted Pendulum Problem
% Behnam Moradi
% Student ID: 200 433 555
% Camputer Aided Process Engineering
clc; clear all; close all;
% CasADi v3.4.5
%% Symbol Definition
import casadi.* % Importing Casadi Library
T = 0.1; % Sampling time [s]
N = 14; % Prediction horizon
F_max = 12; F_min = -F_max; % Saturation Constraints
x1 = SX.sym('x1'); % State X1
x2 = SX.sym('x2'); % State X2
x3 = SX.sym('x3'); % State X3
x4 = SX.sym('x4'); % State X4
states = [x1;x2;x3;x4]; % States Vestor
n_states = length(states);
F = SX.sym('F'); % Input Variavble or Decision Variable
controls = [F]; % Input Vectors
n_controls = length(controls);
rhs = [x2;(-0.1431*x2 + 2.095*sin(x3)*cos(x3) + 0.0856*sin(x3)*x4^2 +1.3373*F)/(1 -
0.2143*(cos(x3))^2); x4; (-24.5*sin(x3) + 0.3579*x2*cos(x3) -
0.2146*sin(x3)*cos(x3)*x4^2 - 3.32*cos(x3)*F)/(1-0.2143*(cos(x3))^2)];
f = Function('f',{states,controls},{rhs}); % nonlinear mapping function f(x,u)
U = SX.sym('U',n_controls,N); % Decision variables (controls)
P = SX.sym('P',n_states + n_states); % parameters (which include the initial and the
reference state of the robot)
% P = SX.sym('P',N*(n_states + n_controls));
X = SX.sym('X',n_states,(N+1)); % A Matrix that represents the states over the
optimization problem.
%% compute solution symbolically
X(:,1) = P(1:4); % initial state
for k = 1:N
 st = X(:,k);
 con = U(:,k);
 f_value = f(st,con);
 st_next = st+ (T*f_value);
 X(:,k+1) = st_next;
end
% this function to get the optimal trajectory knowing the optimal solution
ff=Function('ff',{U,P},{X});
obj = 0; % Objective function
g = []; % constraints vector
Q = zeros(4,4);
Q(1,1) = 100;
Q(2,2) = 0;
Q(3,3) = 50;
Q(4,4) = 0;% weighing matrices (states)
R = 1; % weighing matrices (controls)

%% Compute Objective
for k=1:N
 st = X(:,k);
 con = U(:,k);
 obj = obj+(st-P(5:8))'*Q*(st-P(5:8)) + con'*R*con; % calculate obj
% obj = obj+(st-P(5*(k-1)+1:5*k-1))'*Q*(st-P(5*(k-1)+1:5*k-1)) + (con -
P(5*k))'*R*(con - P(5*k));
end
%% Compute Constraints
for k = 1:N+1
 g = [g ; X(1,k)]; %state x1
 g = [g ; X(2,k)]; %state x3
 g = [g ; X(3,k)]; %state x1
 g = [g ; X(4,k)]; %state x3
end
% make the decision variables one column vector
OPT_variables = reshape(U,N,1);
nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', P);
opts = struct;
opts.ipopt.max_iter = 100;
opts.ipopt.print_level =0;%0,3
opts.print_time = 0;
opts.ipopt.acceptable_tol =1e-8;
opts.ipopt.acceptable_obj_change_tol = 1e-6;
solver = nlpsol('solver', 'ipopt', nlp_prob,opts);
args = struct;
args.lbg(1:4:4*(N+1),1) = -0.5; % lower bound of the states x1
args.ubg(1:4:4*(N+1),1) = 0.5; % upper bound of the states x1
args.lbg(2:4:4*(N+1),1) = -inf; % lower bound of the states x2
args.ubg(2:4:4*(N+1),1) = inf; % upper bound of the states x2
args.lbg(3:4:4*(N+1),1) = -pi; % lower bound of the states x3
args.ubg(3:4:4*(N+1),1) = 2*pi; % upper bound of the states x3
args.lbg(4:4:4*(N+1),1) = -inf; % lower bound of the states x4
args.ubg(4:4:4*(N+1),1) = inf; % upper bound of the states x4
% input constraints
args.lbx(1:1:N,1) = F_min;
args.ubx(1:1:N,1) = F_max;
%% THE SIMULATION LOOP SHOULD START FROM HERE
%-------------------------------------------
t0 = 0;
x0 = [-0.25 ; 0 ;pi-pi/20; 0]; % initial condition.
xs = [0 ; 0 ; pi; 0]; % Reference posture.
xref1 = -0.25;
xs(1) = xref1;
xx(:,1) = x0; % xx contains the history of states
t(1) = t0;
u0 = zeros(N,1); % two control inputs
sim_tim = 20; % Maximum simulation time

%% Start MPC
mpciter = 0;
xx1 = [];
u_cl=[];
main_loop = tic;
xref1 = -0.25;
uref = 0;
j = 1;
for i = 1:5
 xs(1) = xs(1) + 0.1;
while norm((x0 - xs),3) > 1e-4
 args.p = [x0;xs]; % set the values of the parameters vector
 args.x0 = reshape(u0',N,1); % initial value of the optimization variables
 sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
 'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);

 u = reshape(full(sol.x)',1,N)';
 ff_value = ff(u',args.p); % compute OPTIMAL solution TRAJECTORY
 xx1(:,1:4,mpciter+1)= full(ff_value)';
 u_cl= [u_cl ; u(1,:)];
 t(mpciter+1) = t0;
 [t0, x0, u0] = shift(T, t0, x0, u,f); % get the initialization of the next
optimization step
 xx(:,mpciter+4) = x0;
 mpciter
 mpciter = mpciter + 1;
 setP(mpciter) = xs(1);
 itr(mpciter) = mpciter;
end;
end
main_loop_time = toc(main_loop);
ss_error = norm((x0-xs),3)
average_mpc_time = main_loop_time/(mpciter+1);
u_cl(1) = 0;
xx(1,2) = xx(1,1);
xx(1,3) = xx(1,4);
xx(3,2) = xx(3,1);
xx(3,3) = xx(3,4);
figure;
plot(xx(1,:), 'LineWidth',2);
hold on
plot(itr,setP, 'LineWidth',2);
title('The Cart Position and Setpoint Tracking');
legend('The Cart Position','Setpoints to be Tracked');
xlabel('Number of Iterations');
ylabel('The Cart Position(m)');
figure;
plot(xx(3,:)*180/pi, 'LineWidth',2);
title('The Penulum Angle Response');
legend('The Pendulum Angle');
xlabel('Number of Iterations');
ylabel('The Pendulum Angle(degree)');
figure;
stairs(u_cl, 'Linewidth', 3);
title('The Inpute Voltage applied to the DC Motor');
legend('Applied Voltage');
xlabel('Number of Iterations');
ylabel('Applied Voltage(Volts)');