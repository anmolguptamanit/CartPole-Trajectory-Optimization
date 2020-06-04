function f = objective(x)
  cost = 0;
  for i = 41:80
    cost = cost + x(i,1)^2;
  end
  f = cost;
end

    