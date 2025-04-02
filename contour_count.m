% contour_count
F = @(z) [exp(1i*z.^2) 1; 1 1];
Fp = @(z) [2i*z*exp(1i*z.^2) 0; 0 0];
gam = 10; rad = 1.25; nc = 900; % center, radius, nr of nodes on circle
w = exp(2i*pi*(1:nc)/nc); z = gam+rad*w; % unit roots and quad pts
nr = 0;
for k = 1:nc
    nr = nr + rad*w(k)*trace(F(z(k))\Fp(z(k)));
end
nbr_evs = nr/nc % number of eigenvalues
