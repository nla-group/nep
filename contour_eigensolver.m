% contour_eigensolver
n = 2; F = @(z) [exp(1i*z.^2), 1; 1, 1];
ell = 1; r = 1; L = rand(n,ell); R = rand(n,r); % probing matrices
gam = 0; rad = 3; nc = 32; % circle center & radius, nbr of nodes
w = exp(2i*pi*(1:nc)/nc); z = gam+rad*w; % unit roots and quad pts
pbar = 6; A = zeros(ell,r,2*pbar);       % matrices of moments
for k = 1:nc
    Fz = L'*(F(z(k))\R);
    for j = 0:2*pbar-1
        A(:,:,j+1) = A(:,:,j+1) + (w(k)^j*rad*w(k)/nc)*Fz;
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
