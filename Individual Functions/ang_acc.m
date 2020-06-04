function theta_ddot = ang_acc(x_i, m, M, l, g, u_i)
  
  numerator = -1*((l*m*cos(x_i(2))*sin(x_i(2))*x_i(4)^2) + u_i*cos(x_i(2)) + ((m+M)*g*sin(x_i(2))));
  
  denomenator = l*M + l*m*(1-cos(x_i(2))^2);
  
  theta_ddot = numerator/denomenator;
end
