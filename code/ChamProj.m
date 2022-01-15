function res = ChamProj(f, lambda, params, PrevData)
%% Solves E(u) = J_TV(u) + Lambda/2 * ||u - f||^2
% This function evolves an image u[n] to u[n+1] based on a simple
% gradient descent method, using the functional:
%% input: f - input image; lambda - ratio between energies;
%% output: u - processed image; 
%% example: u = AHMOD(f,lambda)

% parse inputs and defaults
if ~exist('params','var'), params=struct(); end
if ~isfield(params,'SmoothParam')
    params.SmoothParam.Derv_HWinSize = 3;
    params.SmoothParam.Derv_Sigma = 1.5;
    params.SmoothParam.ST_HWinSize = 3;
    params.SmoothParam.ST_Sigma = 1.5;
end
if ~isfield(params,'print'), params.print = 0; end
if ~isfield(params,'numIterations'), params.numIterations = 1000; end
if ~isfield(params,'TransformType'), params.TransformType = 'TV'; end
if ~isfield(params,'k'), params.k = 1; end
if ~isfield(params,'tau'), params.tau = 0.2; end
if ~isfield(params,'A') && ~strcmp(params.TransformType,'TV'), params.A = CalcA(f, params.k, params.SmoothParam); end

% code inits
[ny,nx]=size(f); 
u = zeros(ny,nx); 
if exist('PrevData') && isfield(PrevData,'px') && isfield(PrevData,'py')
	q = cat(3,PrevData.px,PrevData.py);
else
    q = zeros(ny,nx,2);
end

% algo inits
absVector = @(v) sqrt(sum(v.^2,3));
if strcmp(params.TransformType,'TV')
    grad_func = @ (u) grad(u);
    div_func = @ (q) div(q);
else
    grad_func = @ (u) gradA(params.A, u);
    div_func = @ (q) divA(params.A, q);
end
tau = params.tau;
lami=1/lambda;  % here the algorithm works with inverse of lambda

% epsilon1 = 1e-5;
% if (max(f(:))-min(f(:)))< epsilon1, params.numIterations = 100; end

for i=1:params.numIterations 
    divP = div_func(q);   
    G = grad_func(divP - f/lami);
    aG = absVector(G);
    q(:,:,1) = (q(:,:,1) + tau*G(:,:,1))./(1+tau*aG);
    q(:,:,2) = (q(:,:,2) + tau*G(:,:,2))./(1+tau*aG);
end 

u_prev = f-lami*divP; 
u=f-lami*div_func(q);

res.u = u;
res.lastErr = norm(u-u_prev)/norm(u);
res.PrevData.px = q(:,:,1);
res.PrevData.py = q(:,:,2);
end