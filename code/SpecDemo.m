%% Calculate the Spectrum Transform for the "FlowerBox" image
close('all','force')
clc
% Choose Transform Type
TransformType = 'TV';
% TransformType = 'A2TV';
NumericalMethod = 'ChambolleProjection';
% NumericalMethod = 'ChambollePock';

%% Perform the chosen Spectral Transform
image = normalize_image(rgb2gray(imread('FlowerBox.jpg')));
figure(1);imshow(image,[]);

if strcmp(TransformType, 'TV')
    Max_time = 2;
elseif strcmp(TransformType, 'A2TV')
    Max_time = 10;
end
Num_of_bands = 100;
dt = Max_time/Num_of_bands;

% define params
params.TransformType = TransformType;
params.NumericalMethod = NumericalMethod;
params.numIterations = 1000;

tic; XTV = spec2D_evolve(image, Max_time, dt, params);toc;  % evolve image

figure(2); plot(XTV.T,XTV.S);title('S(t)');
%% Show the Phi Bands
Phi = XTV.Phi;
figure(3); montage(permute(Phi,[1 2 4 3]),'size',[5 20],'DisplayRange',[-1 1]); title('All the Phi Bands');
%% Filter Bands
H = ones(size(XTV.T));
f_H = [];
if strcmp(TransformType, 'TV')
    freq = [ 44 + [-2:2], 70 + [-2:2]]; %Small & Big Flowers freq
elseif strcmp(TransformType, 'A2TV')
    freq = [ 27 + [-2:2], 58 + [-2:2]]; %Small & Big Flowers freq
end
if length(H) >= max(freq)
    for ii=1:12
        H1 = H; 
        H1(freq)=(ii-1)/2;
        f_H(:,:,ii) = spec2D_filter( XTV.Phi, H1, XTV.f_r, XTV.ScaleParam.dt );
    end
    figure(4); montage(permute(f_H,[1 2 4 3]),'size',[4 3],'DisplayRange',[-0.5 2]); title('Filtered Bands');
end
