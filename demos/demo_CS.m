function outs = demo_CS(filename)
rand('twister',2000);
addpath(strcat(pwd,'/utils'));


% This one is for reality checking
% load brain.6.01-zpad.mat
load(filename)
N = size(im);
if length(N) == 2
    N = [N 1];
end
DN = [N(1) N(2)];
data = zeros(N);
im_dc = zeros(N);
res = zeros(N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Direction Recon Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = '/micehome/asalerno/Documents/CompressedSensing/GradientVectorMag.txt'; % Vector file
thresh = 0.9; % Minimum dot product we'll accept
sigma = 2; % Standard deviation of the gaussian to control thickness (this is the weight)
dirWeight = 0.00;   % Weight for directionally similar penalty

% Check to make sure that the number of directions is the same as number of
% slices (one per direction) in our stack!
% Comment this out when doing reality checks
dimcheck = load(filename);
if (isempty(find(size(dimcheck) == N(3),1))) && dirWeight ~= 0
    error('The data does not comply with the number of directions')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L1 Recon Parameters -- in normal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%N = [256,256]; 		% image Size
%DN = [256,256]; 	% data Size
% N = [128 128];
% DN = [128 128];
pctg = [0.25];  	% undersampling factor
P = 5;			% Variable density polymonial degree
TVWeight = 0.01; 	% Weight for TV penalty
xfmWeight = 0.1;	% Weight for Transform L1 penalty
Itnlim = 8;		% Number of iterations

% generate variable density random sampling
pdf = genPDF(DN,P,pctg , 2 ,0.1,0);	% generates the sampling PDF
k = genSampling(pdf,10,60);		% generates a sampling pattern

%generate transform operator
XFM = Wavelet('Daubechies',6,4);	% Wavelet
%XFM = TIDCT(8,4);			% DCT
%XFM = 1;				% Identity transform

%  ------------- END


% AS
for kk=1:N(3)
    % calculate the phase:
    ph = phCalc(squeeze(im(:,:,kk)),0,0);
    % FT
    trans.FT{kk} = p2DFT(k, [N(1) N(2)], ph, 2);
    FT = trans.FT{kk};
    data(:,:,kk) = reshape(FT*squeeze(im(:,:,kk)),[N(1) N(2) 1]);
    im_dc(:,:,kk) = reshape(FT'*(squeeze(data(:,:,kk))./pdf),[N(1) N(2) 1]);
    res(:,:,kk) = reshape(XFM*(squeeze(im_dc(:,:,kk))),[N(1) N(2) 1]);
end
% ---------------------------


% initialize Parameters for reconstruction
param = init;
param.FT = trans.FT;
param.XFM = XFM;
param.TV = TVOP;
param.data = data;
param.TVWeight =TVWeight;     % TV penalty
param.xfmWeight = xfmWeight;  % L1 wavelet penalty
param.dirWeight = dirWeight;  % directional weight
param.Itnlim = Itnlim;

if param.dirWeight
    [param.dirPair, param.dirPairWeight] = dotThresh(filename,thresh,sigma); % AS
end

tic
for n=1:8
    res = fnlCg(res,param);
	im_res = XFM'*res(:,:,1);
	%figure(100), imshow(abs(im_res),[]), drawnow
    figure(3)
    subplot(2,4,n)
    imshow(abs(im_res),[])
end
toc

% % create a low-res mask
% mask_lr = genLRSampling_pctg(DN,pctg,1,0);
% im_lr = ifft2c(zpad(fft2c(im).*mask_lr,N(1),N(2)));
% 
% im_full = ifft2c(zpad(fft2c(im),N(1),N(2)));
% figure
% imshow(abs(cat(2,im_full,im_lr,im_dc,im_res)),[]);
% c = caxis;
% clf
% subplot(221)
% imshow(abs(im_full),c);
% title('Original (Fully Sampled)')
% subplot(222)
% imshow(abs(im_lr),c);
% title('Low Res')
% subplot(223)
% imshow(abs(im_dc),c);
% title('Zero Filled')
% subplot(224)
% imshow(abs(im_res),c);
% title(['CS with ',num2str(round(100*pctg)),'% pts']);
% 
% figure, plot(1:N(1), abs(im_full(end/2,:)),1:N(1), abs(im_lr(end/2,:)), 1:N(2), abs(im_dc(end/2,:)), 1:N(2), abs(im_res(end/2,:)),'LineWidth',2);
% legend('original', 'LR', 'zf-w/dc', 'TV');


w = whos;
for a = 1:length(w)
outs.(w(a).name) = eval(w(a).name);
end
