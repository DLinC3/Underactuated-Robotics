%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Aerodynamic residual learning (offline):
%   - Baseline: flat-plate lift coefficient C_{L,fp} (loaded as Clfp)
%   - Residual: ΔC_L(α, δ_e) ≈ θ^T φ([α; δ_e]) with RBF features + bias
%   - Fit θ by ridge regression: θ = (Φ'Φ + γI)^{-1} Φ'y
%
% notes for my future self:
%   1) loads data-driven lift residuals (Cl) and flat-plate baseline (Clfp)
%   2) constructs a tiled grid of RBF centers μ_i over (α, δ_e)
%   3) fits θ via regularized least squares
%   4) compares: (data + baseline) vs (baseline) vs (fit + baseline)
%      and also compares against a saved model in model_coeffs.mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fit_coeffs
close all
clear all

% Load aerodynamic coefficients:
%   - aero_coeffs: learned/identified residual lift (Cl) from data
%   - aero_coeffs_fp: flat-plate baseline lift (Clfp)
addpath("data/")
aero_coeffs=load('aero_coeffs.mat');
aero_coeffs_fp=load('aero_coeffs_fp.mat');

Clfp=aero_coeffs_fp.Cltfp;    % baseline C_{L,fp}(α, δ_e)
Cl=aero_coeffs.Clt;           % residual ΔC_L(α, δ_e) from data
alpha=aero_coeffs.alphat;     % α samples
el=aero_coeffs.elt;           % δ_e samples (elevator)

% -------------------------------------------------------------------------
% Build an RBF center grid over the (α, δ_e) domain
% -------------------------------------------------------------------------
numalphas=4;
numels=4;

maxalpha=pi;
minalpha=0;
maxel=pi/4;
minel=-pi/4;

dalpha=(maxalpha-minalpha)/numalphas;
del=(maxel-minel)/numels;

alphas=minalpha:dalpha:maxalpha;
els=minel:del:maxel;

muCl=make_grid(alphas,els);   % μ_i tiling the domain

% RBF spread Σ (shared across centers here)
sigmaCl=[.1 0;0 .1];

paramf.sigma=sigmaCl;
paramf.mu=muCl;

display('Fitting Aerodynamic Coefficients Started...')

% -------------------------------------------------------------------------
% Ridge regression on the linear-in-parameters model:
%   ΔC_L(α,δ_e) ≈ θ^T φ([α;δ_e])
% where φ stacks Gaussian RBFs evaluated at each center, plus a bias term.
% -------------------------------------------------------------------------
thetaCl=regularized_least_squares(@gaussian,[alpha;el],Cl',paramf);

% Evaluate learned residual on the training samples
for k=1:length(alpha)
   phi= gaussian([alpha(k);el(k)],paramf);   % φ([α;δ_e])
   Clf(k)=thetaCl'*phi;                      % θ^T φ
end

display('Lift Coefficients Fit...')

% Compare:
%   - Data total:     (Cl + Clfp)
%   - Baseline:       Clfp
%   - Fit total:      (Clf + Clfp)
figure()
plot(alpha,Cl+Clfp,'*b',alpha,Clfp,'*r',alpha,Clf+Clfp,'og')
xlabel('Alpha')
ylabel('Lift Coefficient')
legend('Data', 'Flat Plat Theory', 'Model')

% Print fitted weights θ
thetaCl

% -------------------------------------------------------------------------
% Compare against a saved model (model_coeffs.mat)
% -------------------------------------------------------------------------
coeffs = load('model_coeffs.mat');
coeffs.thetaCl

coeffs = load('model_coeffs.mat');

Cl_off = zeros(size(Cl));
for k = 1:length(alpha)
    x = [alpha(k); el(k)];
    phi_off = gaussian(x, paramf);           % same feature map φ
    Cl_off(k) = coeffs.thetaCl' * phi_off;   % saved θ applied to φ
end

figure()
plot(alpha, Cl + Clfp, '*b', alpha, Clfp, '*r', ...
     alpha, Clf + Clfp, 'og', alpha, Cl_off + Clfp, 'oy');
xlabel('Alpha');
ylabel('Lift Coefficient');
legend('Data', 'Flat Plat Theory', 'My Model', 'modelcoeffs.mat Model');
title('Lift Coefficient Compare');

rmpath("data/")
end



%==========================================================================
% regularized_least_squares
%
% Fit θ for a linear model y ≈ Φ θ using ridge regression:
%   θ = argmin ||Φθ - y||_2^2 + γ||θ||_2^2
%     = (Φ'Φ + γI)^{-1} Φ'y
%
% Inputs:
%   BASIS: function handle for feature map φ(x) ∈ R^{Np+1}
%   x:     input samples stacked as columns (d×N)
%   y:     targets (N×1 or 1×N)
%   param: carries RBF centers μ and spread Σ
%==========================================================================
function theta = regularized_least_squares(BASIS, x, y, param)

    % N  = number of samples
    % Np = number of RBF centers (excluding bias)
    N  = size(x,2);
    Np = size(param.mu,2);
    
    % Design matrix Φ ∈ R^{N × (Np+1)} where row k is φ(x_k)^T
    Phi = zeros(N, Np+1);
    for k = 1:N
        phi = BASIS(x(:,k), param);   % (Np+1)×1
        Phi(k,:) = phi';              % 1×(Np+1)
    end
   
    gamma = 0.1; % ridge regularization strength (matches your note)
    theta = (Phi' * Phi + gamma * eye(Np+1)) \ (Phi' * y(:))
end

%==========================================================================
% gaussian
%
% Gaussian RBF feature map with shared Σ:
%   φ_i(x) = exp( -1/2 (x-μ_i)^T Σ^{-1} (x-μ_i) ),  i=1..Np
% plus a constant bias feature appended at the end.
%
% Output:
%   phi ∈ R^{Np+1}
%==========================================================================
function phi = gaussian(x, param)
    mu    = param.mu;     % 2×Np centers
    sigma = param.sigma;  % 2×2 spread matrix Σ
    
    [~, Np] = size(mu);
    phi = zeros(Np+1,1);
    for i = 1:Np
        diff = x - mu(:,i);
        % Using (sigma \ diff) = Σ^{-1} diff (numerically nicer than inv)
        phi(i) = exp( -0.5 * (diff' * (sigma \ diff)) );
    end
    phi(Np+1) = 1; % bias term
end

%==========================================================================
% make_grid
%
% Create a Cartesian grid of RBF centers μ_i over two axes x and y.
% Output mu is 2×(mx*my) with columns [x_i; y_j].
%==========================================================================
function mu=make_grid(x,y)
    k=1;
    mx=length(x);
    my=length(y);
    for i=1:mx
        for j=1:my
            mu(:,k)=[x(i);y(j)];
            k=k+1;
        end
    end
end
