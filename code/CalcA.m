function A = CalcA( f , k, params)
%% Calculate the matrix A using the structure tensor
% parse inputs and defaults
if ~exist('params','var'), params=struct(); end
if ~isfield(params,'verbose'), params.verbose = 0; end
if ~isfield(params,'interflag'), params.interflag = false; end
if ~isfield(params,'interX'), params.interX = 3; end
if ~isfield(params,'binary_image'), params.binary_image = false; end
if ~isfield(params,'Derv_HWinSize'), params.Derv_HWinSize = 3; end
if ~isfield(params,'Derv_Sigma'), params.Derv_Sigma = 1.5; end
if ~isfield(params,'ST_HWinSize'), params.ST_HWinSize = 3; end
if ~isfield(params,'ST_Sigma'), params.ST_Sigma = 1.5; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
interX = params.interX;

if params.interflag
    f = imresize(f,interX);
end

hSmooth = fspecial('gaussian', params.Derv_HWinSize,params.Derv_Sigma);
df = grad(f);
fx = df(:,:,1);
fy = df(:,:,2);

sfx = imfilter(fx, hSmooth);
sfy = imfilter(fy, hSmooth);

hSmooth = fspecial('gaussian', params.ST_HWinSize,params.ST_Sigma);
fx2 = imfilter(sfx.*sfx, hSmooth);
fxfy = imfilter(sfx.*sfy, hSmooth);
fy2 = imfilter(sfy.*sfy, hSmooth);

[szy, szx]=size(f);
ST = [fx2(:) fxfy(:) fxfy(:) fy2(:)]';

% Create Stats for choosing K 
parfor n = 1 : size(ST,2)
    structureTensor = reshape(ST(:,n),2,2);
    [~, D] = eig(structureTensor);
    maxEigs(n) = max(D(:));
end
maxImEigs = reshape(maxEigs,szy,[]);
K = k * mean(maxEigs);

% Create A
maxImD = zeros(1,size(ST,2));
NewST = zeros(size(ST));
parfor n = 1 : size(ST,2)
    structureTensor = reshape(ST(:,n),2,2);
    [V, D] = eig(structureTensor);
    if (D(1,1) > D(2,2))
        D(1,1) = WeickertG(D(1,1),K).^(0.5);
        maxImD(n) = D(1,1) ;
        D(2,2) = 1;
    else
        D(2,2) = WeickertG(D(2,2),K).^(0.5);
        maxImD(n) = D(2,2) ;
        D(1,1) = 1;
    end
    
    ST_temp = V*D/V;
    NewST(:,n) = ST_temp(:);
end

maxImD = reshape(maxImD,szy,[]);
A{1,1} = reshape(NewST(1,:),szy,[]);
A{1,2} = reshape(NewST(2,:),szy,[]);
A{2,1} = reshape(NewST(3,:),szy,[]);
A{2,2} = reshape(NewST(4,:),szy,[]);

if params.interflag
    A{1,1}= imresize( A{1,1},1/interX, 'nearest');
    A{1,2}= imresize( A{1,2},1/interX, 'nearest');
    A{2,1}= imresize( A{2,1},1/interX, 'nearest');
    A{2,2}= imresize( A{2,2},1/interX, 'nearest');
    maxImEigs = imresize( maxImEigs,1/interX, 'nearest');
    maxImD = imresize( maxImD,1/interX, 'nearest');
end
if params.verbose > 1
    figure(201); imshow(maxImEigs,[]);
    figure(202); imshow(maxImD,[]);
    drawnow;pause(0.1);
end
end


