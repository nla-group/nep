% Newton_QR_banded
n = 100; lambda0 = 4.0;
ind(1,:) = [1 (1:n-1)]; ind(2,:) = [(2:n) n]; ind(3,:) = ind(2,:);
opts.tol = 5.0e-10; opts.maxitn = 20; opts.rr = 0; opts.nevs = 0;
fun = @(z) f_loaded_string(z,n);
for k = 1:5
    [lam,v,w,h] = NQR4UB(n,fun,ind,lambda0,opts);
    opts.evs(k) = lam; opts.nevs = opts.nevs+1;
end
evs = opts.evs

function [F Fp] = f_loaded_string(z,n)
% Return F(z) and derivative F'(z) for the loaded_string problem.
% Both F(z) and F'(z) are compactly stored as required by NQR4UB.
F(1,:) = [2*n-2*z/3/n (-n-z/6/n)*ones(1,n-1)];
F(2,:) = [-n-z/(6*n)  (2*n-2*z/3/n)*ones(1,n-2) n-z/3/n+z/(z-1)];
F(3,:) = [0 (-n-z/6/n)*ones(1,n-2) 0];
Fp(1,:) = [-2/3/n (-1/6/n)*ones(1,n-1)];
Fp(2,:) = [-1/(6*n)  (-2/3/n)*ones(1,n-2) -1/3/n-1/((z-1)^2)];
Fp(3,:) = [0 (-1/6/n)*ones(1,n-2) 0];
F = sparse(F); Fp = sparse(Fp);
