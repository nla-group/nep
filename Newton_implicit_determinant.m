% Newton_implicit_determinant
n = 2; F = @(z) [exp(1i*z.^2) 1; 1 1];
Fp = @(z) [2i*z*exp(1i*z.^2) 0; 0 0];
tol = 1e-8; maxit = 20; b = [0; 1]; c = b; lam = 2.2 + 1e-4i;
for k = 0:maxit
    [L,U] = lu([F(lam) b; c.' 0]);
    xf = U\(L\[zeros(n,1); 1]);
    if abs(xf(n+1))/norm(F(lam),'fro') < tol, break, end
    xfp = U\(L\[-Fp(lam)*xf(1:n); 0]);
    lam = lam - xf(n+1)/xfp(n+1);
end
if k < maxit, nbr_iter = k, lambda = lam, end
