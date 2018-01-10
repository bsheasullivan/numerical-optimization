function p = cauchy_point(x, delt, Q)
    g = grad_fun(x, Q);
    b = hess_fun(x, Q);
    check = g'*b*g;
    if (check <= 0)
        alpha = delt/(norm(g));
    else
        alpha = min(delt/(norm(g)), (norm(g)).^2/check);
    end
    p = -alpha*g;

end