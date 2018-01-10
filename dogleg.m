function p = dogleg(x, delt, Q)
    g = grad_fun(x, Q);
    b = hess_fun(x, Q);
    p_b = -inv(b)*g;
    p_u = (-(norm(g).^2)/(g'*b*g))*g;
    if (norm(p_b) <= delt)
        p = p_b;
    else
        if (norm(p_u) >= delt)
            p = (-delt/norm(g))*g;
        else
            a = (norm(p_b - p_u)).^2;
            b = 2*p_u'*(p_b - p_u);
            c = (norm(p_u)).^2 - delt.^2;
            tau_wink = (-b + sqrt(b.^2 - 4*a*c))/2*a;
            p = p_u + tau_wink*(p_b - p_u);
        end
    end


end