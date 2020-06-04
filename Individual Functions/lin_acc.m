function d_ddot = lin_acc(x_i, m, M, l, g, u_i)
  
  numerator = (l*m*sin(x_i(2))*x_i(4)^2) + u_i + (m*g*cos(x_i(2))*sin(x_i(2)));
  denomenator = M + m*(1 - cos(x_i(2))^2);
  
  d_ddot = numerator/denomenator;
  
end
