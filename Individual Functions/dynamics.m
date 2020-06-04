function q_dot = dynamics(q_prev, m, M, l, g, u_prev, u_now, dt)
  
  acc_prev = lin_acc(q_prev, m, M, l, g, u_prev);
  ang_acc_prev = ang_acc(q_prev, m, M, l, g, u_prev);
  
  d_dot = q_prev(3) + acc_prev*dt;
  theta_dot = q_prev(3) + ang_acc_prev*dt;
  
  d = q_prev(1) + q_prev(3)*dt + 0.5*acc_prev*dt*dt;
  theta = q_prev(2) + q_prev(4)*dt + 0.5*ang_acc_prev*dt^2;
  
  q_now = [d theta d_dot theta_dot];
  
  acc_now = lin_acc(q_now, m, M, l, g, u_now);
  ang_acc_now = ang_acc(q_now, m, M, l, g, u_now);
  
  q_dot = [d_dot theta_dot acc_now ang_acc_now];
end
