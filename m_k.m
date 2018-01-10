function m = m_k(p, x_k, Q)
    f = fun(x_k, Q);
    g = grad_fun(x_k, Q);
    b = hess_fun(x_k, Q);
    m = f + g'*p + .5*p'*b*p;
end