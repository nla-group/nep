% Newton_trace
F  = @(z) [exp(1i*z.^2) 1; 1 1];
Fp = @(z) [2i*z*exp(1i*z.^2) 0; 0 0];
tol = 1e-8; maxit = 20; lam = 2.2 + 1e-4i;
for k = 0:maxit
    [L,U] = lu(F(lam));
    if abs(prod(diag(U)))/norm(F(lam),'fro')<tol, break, end
    corr = trace(U\(L\Fp(lam)));
    lam = lam - 1/corr;
end
if k < maxit, nbr_iter = k, lambda = lam, end
