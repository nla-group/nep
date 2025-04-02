% inverse_iteration
n = 2; F = @(z) [exp(1i*z.^2) 1; 1 1];
Fp = @(z) [2i*z*exp(1i*z.^2) 0; 0 0];
tol = 1e-8; maxit = 20; lam = 2.2 + 1e-4i;
v = [1; 1]; v = v/norm(v); u = [1; 0];
for k = 0:maxit-1
	if norm(F(lam)*v) < tol, break; end
	vt = F(lam)\Fp(lam)*v;
	lam = lam - u'*v/(u'*vt)
	v = vt/norm(vt);
end
if k < maxit, nbr_iter = k, lambda = lam, end