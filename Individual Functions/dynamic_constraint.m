function [c, ceq] = dynamic_constraint(x)
  f = zeros(40,4);
  constraint = zeros(40,4);
  
  for i=2:40
    f(i, :) = dynamics(x(i,:), 1, 1, 1, 10, x(41+i-1,1), x(41+i, 1), 0.05);
    
    constraint(i, :) = x(i, :) - x(i-1,:) - 0.5*0.05*(f(i-1,:)+f(i,:));
    
  end
  
  c=[];
  ceq=constraint;
end
