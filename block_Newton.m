% block_Newton
A(:,:,1) = [0 1; 1 1]; A(:,:,2) = [1 0; 0 0];
n = 2; k = 4; ell = k; maxit = 30; tol = n*eps;
% construct initial pair (V,J)
a = sqrt(2*pi); % eigenvalues are [a, -a, 1i*a, -1i*a]
d = [a -a a*1i -a*1i]+1e-2*(randn(1,k)+1i*randn(1,k)); J = diag(d);
V = diag([1 -1])*ones(n,k) + 1e-2*(randn(n,k)+1i*randn(n,k));
for iter = 0:maxit
  Z = zeros(n*ell,k); Z(1:n,:) = V;
  for j = 2:ell, Z((j-1)*n+1:j*n,:) = Z((j-2)*n+1:(j-1)*n,:)*J; end
  [Q,R] = qr(Z); R = R(1:k,1:k); V = V/R; J = R*(J/R);
  W(:,:,1) = V; for j = 2:ell, W(:,:,j) = W(:,:,j-1)*J; end
  Res = A(:,:,1)*V*f(1,J) + A(:,:,2)*V*f(2,J);
  if norm(Res,'fro') < tol, break, end
  [DV,DJ] = nlevp_newtonstep(A,@f,V,J,W,Res,zeros(k));
  V = V - DV; J = J - DJ;
end
if k < maxit, nbr_iter = k, evs = eig(J), end

function X = f(j,M)
if j == 1, X = eye(size(M)); end
if j == 2, X = expm(1i*M*M); end
end