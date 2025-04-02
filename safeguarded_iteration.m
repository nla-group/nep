% safeguarded_iteration
n = 100; C1 = n*gallery('tridiag',n); C1(end) = C1(end)/2;
C2 = (abs(gallery('tridiag',n)) + 2*speye(n))/(6*n);
C2(end) = C2(end)/2; C3 = sparse(n,n); C3(n,n) = 1;
F = @(z) C1 - z*C2 + C3*z/(z-1); tol = n*eps;
lam = 1.1; maxit = 5; nb_evs = 5;
fprintf('eigenvalue  #iter   residual\n')
for j = 1:nb_evs
    for k = 0:maxit
        [X,E] = eigs(F(lam),nb_evs+0,'sa'); v = X(:,j);
        res = norm(F(lam)*v)/norm(F(lam),'fro');
        if res < tol, break, end
        c1 = v'*C1*v; c2 = v'*C2*v; c3 = v'*C3*v;
        f = @(z) c1-c2*z+c3*z/(z-1); lam = fzero(f,lam);
    end
    fprintf('%9.2e %5.0f %12.2e\n',lam,k,res)
end
