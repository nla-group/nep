% just a test on Marc Van Barel's remark on contour based method
n = 2; F = @(z) [exp(1i*z.^2), 1; 1, 1];

m = 60; % nr of quadrature points
b = randn(n,1);
z = exp(2i*pi*(0:m-1)/m); plot([z,z(1)],'k:')
V = zeros(N,m);
for j = 1:m,
    V(:,j) = F(z(j))\b;
end

c = randn(n,1);
f = c'*V;

nc = 35; 
g = [];
for k = 0:2*nc-1, % could use FFT for this!
    g(k+1) = sum((z.^(k+1)).*f)/m; 
end
H1 = hankel(g(1:nc),g(nc:2*nc-1));
H2 = hankel(g(2:nc+1),g(n+1:2*nc));
[VV,DD] = eig(H2,H1);
plot(diag(DD)+1i*eps,'g*')
