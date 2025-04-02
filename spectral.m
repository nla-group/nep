% spectral
n = 2; F = @(z) [ exp(1i*z.^2) 1; 1 1];
npts = 150; a = 5; ax = [-a a -a a];
levels = [];
nptsx = npts; nptsy = npts;
x = linspace(-a,a,npts);
[xx, yy] = meshgrid(x, x);
z = complex(xx, yy);
Smin = zeros(nptsy, nptsx);
for j=1:nptsx
    for i=1:nptsy
        Smin(i,j) = min(svd(F(z(i,j))));
    end
end
figure
z = log10(Smin + eps);
pcolor(x,x,z);
map = hot; l = length(map); 
ind = [(1:1:24) (26:4:l)]; 
colormap(map(ind,:)), colorbar
shading interp, axis('square')
hold off, caxis([-1.5,0.2]), shg
set(gca,'XTick',-5:5,'YTick',-5:5)
