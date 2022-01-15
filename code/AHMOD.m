function res = AHMOD(f, lambda, params, PrevData)
%% solves min F(Kx)+G(x)
%         x
% for:
%   G(u) = (lambda/2)*||u-f||^2_2
%   F(Du) = ||Du||_1 (=TV(u) or =A2TV(u))
% using the AHMOD variant described in Chambolle-Pock
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
ProxFS = @ (q,sigma) q./max(1,absVector(q));
ProxG = @ (u,tau) (u+tau*lambda*f)/(1+tau*lambda);
if strcmp(params.TransformType,'TV')
    K = @ (u) grad(u);
    KS = @ (q) - div(q);
else
    K = @ (u) gradA(params.A, u);
    KS = @ (q) - divA(params.A, q);
end
L = 8; 
gamma = 0.7*lambda;
tau = params.tau;
sigma = 1/(tau*L^2);%should satisfy sigma*tau*norm(K)^2<1

% epsilon1 = 1e-5;
% if (max(f(:))-min(f(:)))< epsilon1, params.numIterations = 100; end

for i=1:params.numIterations
    u_prev = u;
    q = ProxFS(q + sigma*K(u),sigma);
    u = ProxG(u - tau*KS(q),tau);       
    theta = 1/sqrt(1 + 2*gamma*tau);
    tau = theta*tau;
    sigma = sigma/theta;
end

res.u = u;
res.lastErr = norm(u-u_prev)/norm(u);
res.PrevData.px = q(:,:,1);
res.PrevData.py = q(:,:,2);
end