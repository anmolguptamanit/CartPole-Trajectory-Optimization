clear all
clc
m = 1;
M = 1;
g = 10;
dt = 0.05;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Linear Equality Variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Aeq = zeros(9, 200);
Aeq(1,1) = 1; % d = 0 for t = 0
Aeq(2,2) = 1; % theta = 0 for t = 0
Aeq(3,3) = 1; % d_dot = 0 for t = 0
Aeq(4,4) = 1; % theta_dot = 0 for t = 0

Aeq(5,158) = 1; % theta = pi at t = 40
Aeq(6,159) = 1; % d_dot = 0 at t = 40
Aeq(7,160) = 1; % theta_dot = 0 at t = 40
Aeq(8,161) = 1; % initially u =0
Aeq(9,200) = 1;

beq = zeros(9, 1);
beq(5,1) = 3.14;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Linear Inequality Variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = [];
b = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Lower Bound %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lb = -100 * ones(200, 1);

for i=1:4:157
  lb(i) = -3;
end

for i=2:4:158
  lb(i) = -3.14;
end

for i=161:200
  lb(i) = -20;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Upper Bonds %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ub = 100 * ones(200, 1);

for i=1:4:157
  ub(i) = 3;
end

for i=2:4:158
  ub(i) = 3.14;
end

for i=161:200
  ub(i) = 20;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Non Linear Constraint %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nonlcon = @dynamic_constraint;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initial Guess %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0 = zeros(200,1);
for i=1:4:157
    val = floor(i/4);
    val = val/39;
    x0(i:i+3) = val*[1;3.14;0;0];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

solution = fmincon(@objective, x0, A, b, Aeq, beq, lb, ub, nonlcon)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d_traj = [];
for i=1:4:157
    idx = ceil(i/4);
    d_traj(idx) = solution(i);
end

theta_traj = [];
for i=2:4:158
    idx = ceil(i/4);
    theta_traj(idx) = solution(i);
end

u_traj = [];
for i=161:200
    u_traj(i-160) = solution(i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting the Trajectories %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = [0:0.05:1.95];
figure
plot(t, d_traj, t, theta_traj, '.-');
xlabel('time'), ylabel('d/theta'), title('Trajectory'), legend('d', 'theta')

figure
plot(t, u_traj);
xlabel('time'), ylabel('u'), title('Trajectory'), legend('u')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% System Dynamics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d_ddot = lin_acc(x_i, m, M, l, g, u_i)
  
  numerator = (l*m*sin(x_i(2))*x_i(4)^2) + u_i + (m*g*cos(x_i(2))*sin(x_i(2)));
  denomenator = M + m*(1 - cos(x_i(2))^2);
  
  d_ddot = numerator/denomenator;
  
end

function theta_ddot = ang_acc(x_i, m, M, l, g, u_i)
  
  numerator = -1*((l*m*cos(x_i(2))*sin(x_i(2))*x_i(4)^2) + u_i*cos(x_i(2)) + ((m+M)*g*sin(x_i(2))));
  
  denomenator = l*M + l*m*(1-cos(x_i(2))^2);
  
  theta_ddot = numerator/denomenator;
end


function q_dot = dynamics(q_prev, m, M, l, g, u_prev, u_now, dt)
  
  acc_prev = lin_acc(q_prev, m, M, l, g, u_prev);
  ang_acc_prev = ang_acc(q_prev, m, M, l, g, u_prev);
  
  d_dot = q_prev(3) + acc_prev*dt;
  theta_dot = q_prev(4) + ang_acc_prev*dt;
  
  d = q_prev(1) + q_prev(3)*dt + 0.5*acc_prev*dt*dt;
  theta = q_prev(2) + q_prev(4)*dt + 0.5*ang_acc_prev*dt^2;
  
  q_now = [d theta d_dot theta_dot];
  
  acc_now = lin_acc(q_now, m, M, l, g, u_now);
  ang_acc_now = ang_acc(q_now, m, M, l, g, u_now);
  
  q_dot = [d_dot theta_dot acc_now ang_acc_now];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Optimization Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f = objective(x)
  cost = 0;
  for i=161:200
    cost = cost + x(i)^2;
  end
  f = cost;
end



function [c, ceq] = dynamic_constraint(x)
  f = zeros(160,1);
  cons = zeros(160, 1);
  for i = 5:4:160
    u_idx = ceil(i/4) + 160;
    f(i:i+3) = dynamics(x(i-4:i-1), 0.5, 1, 0.5, 10, x(u_idx-1), x(u_idx), 0.05);
    cons(i:i+3) = x(i:i+3) - x(i-4:i-1) - 0.5*0.05*(f(i-4:i-1) + f(i:i+3));
  end
  c = [];
  ceq = cons;
end

