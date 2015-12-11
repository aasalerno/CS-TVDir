% this is a script to demonstrate the original experiment by Candes, Romberg and Tao
%
% (c) Michael Lustig 2007
close all
clearvars
rand('twister',2000);
addpath(strcat(pwd,'/utils'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L1 Recon Parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = [256,256]; 		% image Size
DN = [256,256]; 	% data Size
%N = [128 128];
%DN = [128 128];
pctg = [0.25];  	% undersampling factor
P = 5;			% Variable density polymonial degree
TVWeight = 0.001; 	% Weight for TV penalty
xfmWeight = 0.;	% Weight for Transform L1 penalty
Itnlim = 8;		% Number of iterations


% generate variable density random sampling
pdf = genPDF(DN,P,pctg , 2 ,0.1,0);	% generates the sampling PDF
k = genSampling(pdf,10,60);		% generates a sampling pattern

%generate image
%im = (phantom(N(1)))  + randn(N)*0.01 + i*randn(N)*0.01;
load brain6.1-zpad-ksp.mat

%generate Fourier sampling operator
ph = phCalc(im,0,0);
FT{1} = p2DFT(k, N, ph, 2);
data = FT{1}*im;

%generate transform operator

XFM = Wavelet('Daubechies',6,4);	% Wavelet
%XFM = TIDCT(8,4);			% DCT
%XFM = 1;				% Identity transform 	

% initialize Parameters for reconstruction
param = init;
param.FT = FT;
param.XFM = XFM;
param.TV = TVOP;
param.data = data;
param.TVWeight =TVWeight;     % TV penalty 
param.xfmWeight = xfmWeight;  % L1 wavelet penalty
param.Itnlim = Itnlim;

%im_dc = FT{1}'*(data./pdf);	% init with zf-w/dc (zero-fill with density compensation)
im_dc = imsharpen(im,'radius',0.5);
h1 = figure('Position',[0 0 1920 1080]);
imshow(abs(im_dc),[]);drawnow;

res = XFM*im_dc;

% do iterations
h2 = figure('Position',[0 0 1920 1080]);
tic
for n=1:8
	res = fnlCg(res,param);
	im_res = XFM'*res;
	%figure(100), imshow(abs(im_res),[]), drawnow
    
    subplot(2,4,n)
    imshow(abs(im_res),[]);
end
toc

saveas(h1,'~/Desktop/12.11.15/Sharpen-Radius0.5-TV_Only_0.001-StartingPoint.tif');
saveas(h2,'~/Desktop/12.11.15/Sharpen-Radius0.5-TV_Only_0.001.tif');


