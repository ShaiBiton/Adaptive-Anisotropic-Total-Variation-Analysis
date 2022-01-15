function XTV = spec2D_evolve(f, Max_time, dt, params, PrevData)
% Computes TV spectrum of an image, returns S, time interval T, Phi(T) 
% and residual image f_r (to be added in the reconstruction).
% params is an optional struct specifying numerical method and params
% Example: XTV = spec2D_evolve(f, Max_time, dt)

% parse inputs and defaults
if ~exist('params','var'), params=struct(); end
if ~isfield(params,'TransformType'), params.TransformType = 'TV'; end
if ~isfield(params,'NumericalMethod'), params.NumericalMethod = 'ChambolleProjection'; end
if ~isfield(params,'tau'), params.tau = 0.1; end
if ~isfield(params,'numIterations'), params.numIterations = 5000; end
if ~isfield(params,'UsePrevData'), params.UsePrevData = 1; end
if ~isfield(params,'maxError'), params.maxError = 0; end
if ~isfield(params,'CalcAInside'), params.CalcAInside = 0; end
if ~isfield(params,'scale'), params.scale = 1; end
if ~isfield(params,'k'), params.k = 1; end
if ~isfield(params,'verbose'), params.verbose = 1; end
if ~isfield(params,'SmoothParam')
    params.SmoothParam.Derv_HWinSize = 3;
    params.SmoothParam.Derv_Sigma = 1.5;
    params.SmoothParam.ST_HWinSize = 5;
    params.SmoothParam.ST_Sigma = 2;
end

% code inits
A=[];
Phi = zeros(size(f,1),size(f,2));
Zx = zeros(size(f,1),size(f,2));
Zy = zeros(size(f,1),size(f,2));

% algo inits
scale = params.scale;
mu = 1/(2*dt);
NumIter = round(Max_time/dt);
if strcmp(params.NumericalMethod,'ChambolleProjection')
    NumericalFunction = @ChamProj;
else
    NumericalFunction = @AHMOD;
end
if ~exist('PrevData','var'), PrevData=struct(); end
if ~isfield(PrevData,'Alast') && ~strcmp(params.TransformType,'TV')
    params.A = CalcA( f , params.k, params.SmoothParam);
elseif isfield(PrevData,'Alast')
    params.A = PrevData.Alast;
end
if ~isfield(PrevData,'ulast')
    u0 = f*scale;
else
    u0 = ulast;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

res_first = NumericalFunction(f, 2*mu, params);
if params.UsePrevData
    res = NumericalFunction(res_first.u, 2*mu, params, res_first.PrevData);
else
    res = NumericalFunction(res_first.u, 2*mu, params);
end
u1 = res_first.u;
u2 = res.u;
Zx(:,:,1) = res_first.PrevData.px;
Zy(:,:,1) = res_first.PrevData.py;
Zx(:,:,2) = res.PrevData.px;
Zy(:,:,2) = res.PrevData.py;
u(:,:,1) = u0;
u(:,:,2) = u1;
u(:,:,3) = u2;

% setting waitbar
h = waitbar(0,'Please wait...','Name','Approximating XTV Flow...',...
            'CreateCancelBtn',...
            'setappdata(gcbf,''canceling'',1)');
setappdata(h,'canceling',0)

for i=1:NumIter
    startCyc = tic;
    ddu = (u0+u2-2*u1)/(dt*dt);
    t = i*dt;
    phi = ddu*t;
    Phi(:,:,i) = phi;
    S(i) = sum(abs(phi(:)));
    if (i<NumIter) % not last iteration
        u0=u1;
        u1=u2;
        if params.UsePrevData
            res = NumericalFunction(u2, 2*mu, params, res.PrevData);
        else
            res = NumericalFunction(u2, 2*mu, params);
        end
        u2 = res.u;
        Zx(:,:,i+2) = res.PrevData.px;
        Zy(:,:,i+2) = res.PrevData.py;
        lastErr(i+2) = res.lastErr;
        u(:,:,i+3) = u2;
        if params.CalcAInside == 1
            params.A = CalcA(u2 , params.k, params.SmoothParam);
        end
    end
    
    if getappdata(h,'canceling')
        NumIter = i;
        break
    end
    TCyc = toc(startCyc);
    rmnTotalTime = (NumIter - i)*TCyc;
    Hours = fix(rmnTotalTime/3600);
    Min = fix(mod(rmnTotalTime,3600)/60);
    Sec = fix(mod(rmnTotalTime,60));
    remainingTimeStr = [num2str(Hours,'%02d') ':' num2str(Min,'%02d') ':' num2str(Sec,'%02d')];
    if params.verbose > 1
        figure(100); plot((1:length(S))*dt,S);title('S(t)');
        figure(101); imagesc(u2);colorbar;title('u(t)');
    end
    waitbar(i/NumIter,h,['Iter ' num2str(i) ' Out of ' num2str(NumIter) '. Remaning Time: ' remainingTimeStr])
end
delete(h) 

if params.verbose > 1
    figure(200); plot(lastErr);title('lastErr');
end

f_r = (NumIter+1)*u1-NumIter*u2;  % residual image
T = (1:NumIter)*dt;
% rescale
S = S/scale; Phi = Phi/scale; f_r = f_r/scale;u = u/scale;

for j=1:length(S)
    p(:,:,j) = (u(:,:,j)-u(:,:,j+1))./dt;
end

XTV.S = S;
XTV.T = T;
XTV.Phi = Phi;
XTV.u = u;
XTV.f_r = f_r;
XTV.lastErr = lastErr;
XTV.Zx = Zx;
XTV.Zy = Zy;
XTV.p = p;
XTV.params = params;
XTV.ScaleParam.dt = dt;
XTV.ScaleParam.Max_time = Max_time;
XTV.ScaleParam.NumIter = NumIter;
if exist('A')
    XTV.A = A;
end
end