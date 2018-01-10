%Trust region algorithim
n = 100;
[A, R] = qr(rand(n, n));
v = n*(rand(n, 1));
D = diag(v);
Q = A'*D*A;

delt_hi = 1;
delt_0 = delt_hi/2;
eta = .2;

x_0 = 2*rand(n, 1);
x_k = x_0;

IT = [];
RE = [];

delt_k = delt_0;

maxiter = 1000;
tol = .00000000000001;
err = 1000000000;
k = 0;

while (err > tol && k < maxiter)
    %p_k = cauchy_point(x_k, delt_k, Q);
    p_k = dogleg(x_k, delt_k, Q);
    rho_k = (fun(x_k, Q) - fun((x_k + p_k), Q))/(fun(x_k, Q) - m_k(p_k, x_k, Q));
    if rho_k < .25
        delt_k = .25*delt_k;
    else
        if (rho_k > .75 && norm(p_k) == delt_k)
            delt_k = min(2*delt_k, delt_hi);
        else
            delt_k = delt_k;
        end
    end
    if (rho_k > eta)
        x_k = x_k + p_k;
    else 
        x_k = x_k;
    end
    err = fun(x_k, Q);
    IT(end + 1) = k;
    RE(end + 1) = err;
    k = k + 1;
    
end

plot(IT, log10(RE));
