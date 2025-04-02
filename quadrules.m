% quadrules

% This code requires the Schwarz-Christoffel Toolbox by Toby Driscoll
% http://www.math.udel.edu/~driscoll/SC/
% and the Chebfun system
% http://www.chebfun.org/
% to be in the MATLAB path.

figure(1)
p = [ 1+1i , -1+1i , -1-1i , 1-1i , 1+1i  ];
pp = polygon(p);
f = diskmap(pp);
f = center(f,0);
fv = feval(inv(f),1);
r = 1.1; s = 1.08;
Psi = @(w) s*f(w/r);
Psid = @(w) s/r*feval(diff(f),w/r);
z = exp(1i*linspace(0,2*pi,1e4+1));
h1 = fill(real(Psi(r*z)),imag(Psi(r*z)),'y','FaceColor',[.9,.9,.9],'LineStyle','none');
hold on
fill(real(Psi(z/r)),imag(Psi(z/r)),'w','LineStyle','none')
h2 = plot(Psi(z),'b-');
h3 = plot(p,'r--');
evs = [0.3-.4i,-.354,0.2922+0.904i];
h4 = plot(evs,'kx');
tt = linspace(0,2*pi);
rt = 1.105;
bern = (rt*exp(1i*tt)+exp(-1i*tt)/rt)/2;
a = -1+1i; b = 1+1i;
bern = (a-b)/2*bern + (a+b)/2;
h5 = plot(bern,'g:','Color', [0 .5 0]);
legend([h4,h1,h2,h3,h5],'eigenvalues','domain \Omega','conformal','square','Bernstein ellipse','Location',[.5,.5,0.2,.2])

% apply rule
F = @(z) diag(evs)-z*eye(3);
R = [1;1;1];
mom = 0; ex = [ -1;-1;-1];
Nc = 8:8:160;
for k = 1:length(Nc),
    nc = Nc(k);
    % circular trapezoid
    zl = exp(2i*pi*(0.5:nc)/nc);
    nodes = Psi(zl);
    if nc == 16,
        figure(1)
        plot(nodes,'bo')
    end
    omega = (nodes).^mom.*Psid(zl).*zl/nc;
    Ap = 0*R;
    for j = 1:nc,
        Ap = Ap + omega(j)*(F(nodes(j))\R);
    end
    err1(k) = norm(Ap-ex);
    
    % trapezoid on square with endpoints
    pts = 2*(0:nc/4)/(nc/4)-1; % nc/4+1 pts on [-1,1]
    nodes = [ pts-1i , 1+1i*pts , 1i-pts , -1-1i*pts ];
    if nc == 16,
        figure(1)
        plot(nodes,'r+')
    end
    omega1 = nodes(1:nc/4+1).^mom*2/(nc/4)/(2i*pi);
    omega1(1) = omega1(1)/2; omega1(end)=omega1(end)/2;
    omega2 = nodes(nc/4+2:2*nc/4+2).^mom*2/(nc/4)/(2i*pi);
    omega2(1) = omega2(1)/2; omega2(end)=omega2(end)/2;
    omega3 = nodes(2*nc/4+3:3*nc/4+3).^mom*2/(nc/4)/(2i*pi);
    omega3(1) = omega3(1)/2; omega3(end)=omega3(end)/2;
    omega4 = nodes(3*nc/4+4:4*nc/4+4).^mom*2/(nc/4)/(2i*pi);
    omega4(1) = omega4(1)/2; omega4(end)=omega4(end)/2;
    omega = [ omega1, 1i*omega2, -omega3, -1i*omega4 ];
    Ap2 = 0*R;
    for j = 1:length(nodes),
        Ap2 = Ap2 + omega(j)*(F(nodes(j))\R);
    end
    err2(k) = norm(Ap2-ex);
    
    % Gauss on square
    [pts, wgt] = legpts(nc/4); % pts and weights on (-1,1)
    pts = pts(:).';
    nodes = [ pts-1i , 1+1i*pts , 1i-pts , -1-1i*pts ];
    omega = nodes.^mom.*[wgt,wgt,wgt,wgt]/(2i*pi);
    omega(nc/4+1:nc/2) = 1i*omega(nc/4+1:nc/2);
    omega(2*nc/4+1:3*nc/4) = -1*omega(2*nc/4+1:3*nc/4);
    omega(3*nc/4+1:nc) = -1i*omega(3*nc/4+1:nc);
    Ap3 = 0*R;
    for j = 1:nc,
        Ap3 = Ap3 + omega(j)*(F(nodes(j))\R);
    end
    err3(k) = norm(Ap3-ex);
end

figure(1)
axis([-1.2,1.2,-1.2,1.2])
figure(2)
h1 = semilogy(Nc,err1,'b-o');
hold on
semilogy(Nc,.5*r.^(-Nc),'b--')
h2 = semilogy(Nc,err2,'r-+');
semilogy(Nc,4./(Nc.^2),'r--')
h3 = semilogy(Nc,err3,'g-','Color', [0 .5 0]);
semilogy(Nc,.5*rt.^(-Nc/2),'g:','Color', [0 .5 0])
legend([h1,h2,h3],'conformal','trapezoid','Gauss-Legendre','Location','SouthWest')
axis([0,nc+1,1e-5,2])
grid on
