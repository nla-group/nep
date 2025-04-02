% sampling1 - Sampling example using the loaded_string problem by Sovolev.

clear all, close all

if ~exist('util_nleigs','file'), % download and install RKToolbox
    disp('Press any key to download and install RKToolbox.')
    disp('Press Ctrl+C to abort.'), pause
    unzip('http://guettel.com/rktoolbox/rktoolbox.zip'); 
    cd('rktoolbox'); addpath(fullfile(cd)); 
    addpath(fullfile(cd,'utils')); savepath; cd('..')
end

%% problem definition
n = 100; C1 = n*gallery('tridiag',n); C1(end) = C1(end)/2;
C2 = (abs(gallery('tridiag',n))+2*speye(n))/(6*n); C2(end) = C2(end)/2;
C3 = sparse(n,n); C3(n,n) = 1; F = @(z) C1 - z*C2 + C3*z/(z-1);
Nmax  = 50; zz = linspace(4,296,50);
 
%% polynomial Leja on circle
Sigma = 150 + 146*exp(linspace(0,2i*pi,400));  Xi = inf;
tol  = 0; R = util_nleigs(F, Sigma, Xi, tol, Nmax); 
disp('POLYNOMIAL LEJA SAMPLING ON CIRCLE')
disp('  Please wait, evaluating interpolant for all degrees m=0,...,100.')
disp('  This is slow and only done for producing a convergence plot.')
for newN = 0:Nmax,
    for j = 1:length(zz),
        Q = R(zz(j),newN);
        err(j) = norm(Q - F(zz(j)),'fro');
    end
    errN(newN+1) = max(err);
end
figure
semilogy(0:(length(errN)-1),errN,'k-+','Color',.6*[1,1,1])
hold on

%% polynomial Leja on interval
Sigma = [4,296];  Xi = inf;
tol  = 1e-12; R = util_nleigs(F, Sigma, Xi, tol, Nmax); 
disp('POLYNOMIAL LEJA SAMPLING ON INTERVAL')
disp('  Please wait, evaluating interpolant for all degrees m=0,...,100.')
disp('  This is slow and only done for producing a convergence plot.')
errN = [];
for newN = 0:Nmax,
    for j = 1:length(zz),
        Q = R(zz(j),newN);
        err(j) = norm(Q - F(zz(j)),'fro');
    end
    errN(newN+1) = max(err);
end
semilogy(0:(length(errN)-1),errN,'b-.')

%% rational Leja-Bagby on interval
Sigma = [4,296];  Xi = [1,inf]; Nmax = 2;
tol  = 1e-15; R = util_nleigs(F, Sigma, Xi, tol, Nmax); 
disp('RATIONAL LEJA-BAGBY SAMPLING ON INTERVAL')
disp('  Evaluating interpolant for all degrees m=0,...,2.')
errN = [];
for newN = 0:2, 
    for j = 1:length(zz),
        Q = R(zz(j),newN);
        err(j) = norm(Q - F(zz(j)),'fro');
    end
    errN(newN+1) = max(err);
end
figure(1)
semilogy(0:(length(errN)-1),errN,'r-o')
legend('(a) polynomial Leja on disk', '(b) polynomial Leja on interval', ...
    '(c) rational Leja-Bagby','Location',[.505,.138,.2,.17])
grid on
semilogy(0:50,(146/149).^(0:50),'k--','Color',.6*[1,1,1])
rt = 1/1.2243;
semilogy(0:50,rt.^(0:50),'b--')
title('\texttt{loaded\_string:} approximation error',...
    'Fontweight','Normal','Interpreter','latex')
axis([0,50,1e-11,100])

%%
figure(2)
AB = linearize(R);
[Am,Bm] = AB.get_matrices();
ee = eig(full(Am),full(Bm));
plot(real(ee),imag(ee),'ro')
axis([-50,350,-200,200])
legend('eigenvalues of "exact" linearization','Location','SouthEast')
title('\texttt{loaded\_string:} eigenvalues',...
    'Fontweight','Normal','Interpreter','latex')

disp('DONE.')
