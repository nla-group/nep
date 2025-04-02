% contour_eigensolver for the string problem
n = 100; C1 = n*gallery('tridiag',n); C1(end) = C1(end)/2;
C2 = (abs(gallery('tridiag',n))+2*speye(n))/(6*n); C2(end) = C2(end)/2;
C3 = sparse(n,n); C3(n,n) = 1; F = @(z) C1 - z*C2 + C3*z/(z-1);
ell = 2; r = 2; L = rand(n,ell); R = rand(n,r); % probing matrices
gam = 10; rad = 10; nc = 500; % circle center & radius, nbr of nodes
w = exp(2i*pi*(1:nc)/nc); z = gam+rad*w; % unit roots and quad pts
pbar = 1; A = zeros(ell,r,2*pbar);       % matrices of moments
for k = 1:nc
    Fz = L'*(F(z(k))\R);
    for j = 0:2*pbar-1
        A(:,:,j+1) = A(:,:,j+1) + ((exp(w(k)-3).^1).*w(k)^j*rad*w(k)/nc)*Fz;
    end
end
A = A(:,:); B0 = zeros(pbar*ell,pbar*r); B1 = B0;
for j = 0:pbar-1
    B0(1+j*ell:(j+1)*ell,:) = A(:,1+j*r:pbar*r+j*r);
    B1(1+j*ell:(j+1)*ell,:) = A(:,1+(j+1)*r:pbar*r+(j+1)*r);
end
[V,Sig,W] = svd(B0); mbar = find(diag(Sig)/Sig(1)>1e-12,1,'last');
V0 = V(:,1:mbar); Sig0 = Sig(1:mbar,1:mbar); W0 = W(:,1:mbar);
M = (V0'*B1*W0)/Sig0; evs = eig(M); evs = gam+rad*evs(abs(evs)<1)
