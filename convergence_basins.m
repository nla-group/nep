% convergence_basins
% This code produces the plot in the upper left of Figure 4.6. 
% It can easily be modified to produce the other plots in that figure.

n = 2; F = @(z) [ exp(1i*z.^2) 1; 1 1 ];
Fp = @(z) [ 2i*z*exp(1i*z.^2) 0; 0 0 ];
tol =  1e-8; maxit = 20;
a = -3; b = 3; nbpt = 200; % in paper we used nbpt = 1000;
x = [a+(0:nbpt-2)*(b-a)/(nbpt-1), b];
[xx, yy] = meshgrid(x,x);
z = complex(xx,yy);
u = [0 sqrt(2*pi) 1i*sqrt(2*pi) -sqrt(2*pi) -1i*sqrt(2*pi)]'; % exact evs

S = zeros(size(z));
warning off, hb = waitbar(0,'Please wait...');
for i1 = 1:nbpt
    waitbar(i1/nbpt,hb)
    for i2 = 2:nbpt
        % The following is Newton-trace iteration with starting guess lam.
        % Change it to Newton-QR, or any of the other methods, to produce 
        % the other plots in Figure 4.6.
        lam = z(i1,i2);
        for k = 0:maxit-1
            [L,U] = lu(F(lam));
            if abs(prod(diag(U)))/norm(F(lam),'fro') < tol, break, end
            corr = trace(U\(L\Fp(lam)));
            lam = lam - 1/corr;
            if isnan(lam), break; end
        end
        it = k;
        if isnan(lam) || abs(lam) > 3 || norm(F(lam),'fro') > (1/eps)
            it = maxit;
        end
        if it < maxit
            [tmp,S(i1,i2)] = min(abs(lam-u));
        end
    end
end
close(hb), warning on

% plotting
colspec = parula(5); figure
for i = 1:5
    k = find(S == i);
    plot(real(z((k))),imag(z((k))),'.','MarkerSize',5,'Color',.9*colspec(i,:))
    hold on
end
axis square; axis([-3,3,-3,3])
title('Newton-trace method','Fontweight','Normal','Interpreter','latex')
hold on; plot(u,'p','Color','none','MarkerFaceColor',[1,1,1],'MarkerSize',14)
