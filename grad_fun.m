function eval = grad_fun(x, Q)
    eval = Q*x/(1 + x'*Q*x);
end