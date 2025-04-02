% sampling2 - Sampling example for the Gun problem taken from
% http://guettel.com/rktoolbox/examples/html/example_nlep.html

clear all, close all

if ~exist('nlevp','file'), % download and install NLEVP collection
    disp('Press any key to download and install NLEVP collection.')
    disp('Press Ctrl+C to abort.'), pause
    mkdir('nlevp'); cd('nlevp'); 
    unzip('http://www.cl.eps.manchester.ac.uk/medialand/maths/software/NLEVP/nlevp_111222.zip');
    addpath(fullfile(cd)); savepath; cd('..'); 
end

if ~exist('util_nleigs','file'), % download and install RKToolbox
    disp('Press any key to download and install RKToolbox.')
    disp('Press Ctrl+C to abort.'), pause
    unzip('http://guettel.com/rktoolbox/rktoolbox.zip'); 
    cd('rktoolbox'); addpath(fullfile(cd)); 
    addpath(fullfile(cd,'utils')); savepath; cd('..')
end

%% The Gun problem
% We demonstrate NLEIGS linearization and the rational Krylov solution of
% the resulting eigenvalue problem on the gun problem
% 
% $$ \displaystyle A(\lambda)\mathbf x = [K-\lambda M + i\sqrt{\lambda} W_1 
%   + i\sqrt{\lambda-108.8774^2} W_2]\mathbf x = \mathbf 0, $$
% 
% where $K, M, W_1, W_2$ are sparse $9956\times 9956$ matrices. Let us 
% define a function handle to the NEP:

[coeffs, fun] = nlevp('gun');
n = size(coeffs{1}, 1);
A = @(lam) 1*coeffs{1} - lam*coeffs{2} + ...
    1i*sqrt(lam)*coeffs{3} + ...
    1i*sqrt(lam-108.8774^2)*coeffs{4};

%%
% The target set $\Sigma$ in this example is an upper-half disk with
% center $62500$ and radius $50000$.  
% Note that the definition of $A(\lambda)$ involves two branch cuts
% $(-\infty,0]$ and $(-\infty,108.8774^2]$ caused by the square
% roots, and the union of these two is a good choice for the
% singularity set $\Xi$. 
%
% We can now use the utility function |util_nleigs| to sample the
% NLEP on the target set $\Sigma$ using poles from the singularity set
% $\Xi$. The function requires as inputs a function handle to $A(\lambda)$,
% the vertices of $\Sigma$ and $\Xi$ represented by polygons, a tolerance
% for the sampling procedure, and the maximal number of terms. 

Nmax  = 50;
Sigma = 62500 + 50000*exp(1i*pi*[1, linspace(0, 1)]);
Xi   = [-inf, 108.8774^2];
tol  = 1e-15;
QN   = util_nleigs(A, Sigma, Xi, tol, Nmax); 
disp(QN)

%%
% As we can see, the output of |util_nleigs| is an RKFUNM object |QN| 
% representing $Q_N(\lambda)$, a rational matrix-valued function 
% which interpolates $A(\lambda)$ at the nodes $\sigma_j$ (i.e., 
% $Q_N(\sigma_j) = A(\sigma_j)$ for all $j=0,1,\ldots,N$). We can evaluate 
% $Q_N$ at any point $z$ in the complex plane by typing |QN(z)|.
% The |linearize| function can be used to convert $Q_N$ into an
% equivalent linear matrix pencil structure $L_N(z) = A_N - z B_N$ with the 
% same eigenvalues as $Q_N$. Via $L_N$ we can also access the norms $\|D_j\|_F$ 
% of the matrices $D_j$ in the expansion of $Q_N$. Apparently, a degree of 
% $N=32$ was sufficient to represent $A(\lambda)$ to accuracy 
% $\mathrm{tol} = 10^{-15}$:

LN   = linearize(QN); 
figure(1), semilogy(0:LN.N, LN.nrmD/LN.nrmD(1), 'r-'), grid on
legend('relative Frobenius norm of D_j'); xlabel('j')
axis([0, LN.N, 1e-16, 1])

%% 
% Luckily this $N=32$ is quite small due to the Leja-Bagby sampling strategy 
% employed by |util_nleigs|. However, the full linearisation matrices
% $(A_N,B_N)$ are of size $Nn\times Nn$ and hence quite large. Here is
% a spy plot of $(A_N,B_N)$. 

[AN,BN] = LN.get_matrices();
figure(2)
subplot(1,2,1), spy(AN), title('A_N')
subplot(1,2,2), spy(BN), title('B_N')

%%
% The |LN| structure  provides two function handles |multiply| and
% |solve|, which can be used by the |rat_krylov| function to compute a
% rational Krylov basis for $(A_N,B_N)$ without forming these matrices
% explicitly [7, 8]. For the rational Arnoldi algorithm we choose, rather 
% arbitrarily, $5$ cyclically repeated shifts in the interior of $\Sigma$.

shifts = [9.6e+4, 7.9e+4+1.7e+4i, 6.3e+4, 4.6e+4+1.7e+4i, 2.9e+4];

%%
% Let us first plot the target set, the sampling points $\sigma_j$, the 
% poles $\xi_j$, and the shifts of the rational Krylov space:

figure(3)
fill(real(sqrt(Sigma)), imag(sqrt(Sigma)), [1 1 .6])
hold on
plot(sqrt(LN.sigma), 'gx', 'Color', [0 .5 0])
plot(sqrt(LN.xi(LN.xi>0))+1i*eps, 'r.', 'MarkerSize', 14)
plot(sqrt(shifts), 'mo')
xlabel('Re sqrt(lambda)'), ylabel('Im sqrt(lambda)')
legend('target set \Sigma', 'interpolation nodes \sigma_j', ...
       'poles \xi_j', 'RK shifts', 'Location', 'NorthWest')
axis([0,350,-10,110])

%%
% The computation of the rational Arnoldi decomposition $A_N V_{m+1}
% \underline{K_m} = B_N V_{m+1} \underline{H_m}$ is conveniently performed
% by providing the pencil structure |LN| as the first input argument 
% to the |rat_krylov| function. Here, $m=70$ and the starting vector is 
% chosen at random. 
%
% *Note:* _The shifts in this example are cyclically repeated and the
% |solve| function provided in the |LN| structure attempts to 
% reuse LU factors of $n\times n$ matrices whenever possible. 
% The five LU factors required for this example are stored automatically as
% |persistent| variables within the |solve| function._

v = randn(LN.N*n, 1);
shifts = repmat(shifts, 1, 14);
[V, K, H] = rat_krylov(LN, v, shifts, struct('waitbar', 1));

%%
% From the rational Arnoldi decomposition we can easily compute the Ritz
% pairs for the linearisation $(A_N,B_N)$. In the following we extract 
% the Ritz values in the interior of $\Sigma$ and find 
% that there are 21 Ritz values. The leading $n$ elements of the
% corresponding Ritz vectors (normalised to unit norm) are then
% approximations to the eigenvectors $\mathbf{x}$ of $A(\lambda)$.

[X, D] = eig(H(1:end-1, :), K(1:end-1, :));
ritzval = diag(D);
ind = inpolygon(real(ritzval), imag(ritzval), ...
                real(Sigma), imag(Sigma));
ritzval = ritzval(ind); 
ritzvec = V(1:n, 1:end)*(H*X(:, ind));
ritzvec = ritzvec/diag(sqrt(sum(abs(ritzvec).^2)));
disp(length(ritzval))

%%
% Let us compute the nonlinear residual norm $\| A(\lambda)\mathbf{x}\|_2$
% for all 21 Ritz pairs $(\lambda,\mathbf{x})$:

res = arrayfun(@(j) norm(A(ritzval(j))*ritzvec(:, j), 'fro'), ...
               1:length(ritzval));

figure(4)
semilogy(res, 'b-o'), xlim([1, length(res)])
legend('residual norm of Ritz pairs')
xlabel('index of Ritz pair'), hold on

%% 
% We find that all but $4$ Ritz pairs are good approximations to the
% eigenpairs of the nonlinear problem. Let us run five more rational
% Arnoldi iterations by using as shift the mean of the four nonconverged Ritz 
% values. This can be done by simply extending the existing rational 
% Arnoldi decomposition.

% Use mean of Ritz values.
shifts = repmat(mean(ritzval(res>1e-8)), 1, 5);
% Extend the decomposition.
[V, K, H] = rat_krylov(LN, V, K, H, shifts, struct('waitbar',1));

%%
% Now let us compute the improved Ritz pairs and the corresponding nonlinear 
% residuals exactly as above:

[X, D] = eig(H(1:end-1, :), K(1:end-1, :));
ritzval = diag(D);
ind = inpolygon(real(ritzval), imag(ritzval), ...
                real(Sigma),   imag(Sigma));
ritzval = ritzval(ind); 
ritzvec = V(1:n, 1:end)*(H*X(:, ind));
ritzvec = ritzvec/diag(sqrt(sum(abs(ritzvec).^2)));

res = arrayfun(@(j) norm(A(ritzval(j))*ritzvec(:, j), 'fro'), ...
               1:length(ritzval));

figure(4)
semilogy(res, 'r-o')
legend('residual norm of Ritz pairs','residual norm (extended)')

%%
% All wanted Ritz pair are now of sufficiently high accuracy and we are 
% done. Finally, here is a plot of the Ritz values:

figure(3)
plot(sqrt(ritzval), 'b+')
legend('target set \Sigma', 'interpolation nodes \sigma_j', ...
       'poles \xi_j', 'RK shifts', '21 Ritz values', ...
       'Location', 'NorthWest')

   