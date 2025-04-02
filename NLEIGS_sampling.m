% NLEIGS_sampling

if ~exist('util_nleigs','file'), % download and install RKToolbox
    disp('Press any key to download and install RKToolbox.')
    disp('Press Ctrl+C to abort.'), pause
    unzip('http://guettel.com/rktoolbox/rktoolbox.zip'); 
    cd('rktoolbox'); addpath(fullfile(cd)); 
    addpath(fullfile(cd,'utils')); savepath; cd('..')
end

n = 100; C1 = n*gallery('tridiag',n); C1(end) = C1(end)/2;
C2 = (abs(gallery('tridiag',n)) + 2*speye(n))/(6*n);
C2(end) = C2(end)/2; C3 = sparse(n,n); C3(n,n) = 1;
F = @(z) C1 - z*C2 + C3*z/(z-1); tol = 0; mmax = 50;
Sigma = 150+146*exp(2i*pi*(0:99)/100); Xi = inf; % Leja on circle
Sigma = [4,296]; Xi = inf;                       % Leja on interval
Sigma = [4,296]; Xi = [1,inf]; mmax = 2;         % Leja-Bagby
Rm = util_nleigs(F, Sigma, Xi, tol, mmax)        % NLEIGS sampling
Lm = linearize(Rm); [Am,Bm] = Lm.get_matrices(); % linearization

figure % plot linearization pencil
subplot(1,2,1), spy(Am), title('A_m')
subplot(1,2,2), spy(Bm), title('B_m')