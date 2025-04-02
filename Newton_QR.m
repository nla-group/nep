% Newton_QR
n = 2; F = @(z) [exp(1i*z.^2) 1; 1 1];
Fp = @(z) [2i*z*exp(1i*z.^2) 0; 0 0];
tol = 1e-8; maxit = 20; lam = 2.2 + 1e-4i;
for k = 0:maxit
    [Q,R,P] = qr(F(lam));
    if abs(R(n,n))/norm(F(lam),'fro') < tol, break, end
    p = R(1:n-1,1:n-1)\R(1:n-1,n);
    lam = lam - R(n,n)/(Q(:,n)'*Fp(lam)*P*[-p; 1]);
end
if k < maxit, nbr_iter = k, lambda = lam, end
