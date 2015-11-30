function [im_res,diffRMS] = demo_TVxfmDir(TVWeight,xfmWeight,dirWeight,threshAcc)

% Load in the data
load brain.6-zpad-ksp.mat
dat_share = 'data_added_nearest.mat';
load(dat_share);
load('sampPattern.mat');
origDataSize = [180 180];


N = size(im);
if length(N) == 2
    N = [N 1];
end
for i = 1:N(3)
     im(:,:,i) = im(:,:,i)/max(max(im(:,:,i)));
end

DN = [N(1) N(2)];
Itnlim = 8;		% Number of iterations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Direction Recon Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = '/micehome/asalerno/Documents/CompressedSensing/GradientVectorMag.txt'; % Vector file
if nargin < 4
    threshAcc = 4;
end

% Check to make sure that the number of directions is the same as number of
% slices (one per direction) in our stack!
% Comment this out when doing reality checks
dimcheck = load(filename);
if (isempty(find(size(dimcheck) == N(3),1))) && dirWeight ~= 0
    error('The data does not comply with the number of directions')
end

param = init;
if dirWeight
    dirInfo = lsqA(filename,threshAcc);
    param.dirInfo = dirInfo;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L1 Recon Parameters -- in normal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add in a method to create a set of sampling data if the kspace data
% doesn't exist already
if ~exist('samp','var')
    pctg = 0.25;  	% undersampling factor
    radius = 0.2;
    disp('Creating sampling pattern');
    samp = genSampToAdd(DN,pctg,radius,threshAcc+1,filename);
    disp('Sampling pattern created.');
end

%generate transform operator
XFM = Wavelet('Daubechies',20,4);	% Wavelet



% Make sure we are working with data that is 2^n
if 2.^nextpow2(DN(1)) ~= DN(1) || (2.^nextpow2(DN(2)) ~= DN(2))
    im_hold = zeros(N);
    for i = 1:N(3)
        im_hold(:,:,i) = fft2c(zpad(ifft2c(im(:,:,i)),2.^nextpow2(DN(1)),2.^nextpow2(DN(2))));
    end
    im = im_hold;
end


data = zeros(N);
for kk=1:N(3)
    % calculate the phase:
    ph = phCalc(squeeze(im(:,:,kk)),0,0);
    % FT
    trans.FT{kk} = p2DFT(samp(:,:,kk), DN, ph, 2);
    FT = trans.FT{kk};
    data(:,:,kk) = reshape(FT*squeeze(im(:,:,kk)),[N(1) N(2) 1]);
end

% Now we add the k-space data together to create im_dc
data_dc = ksp_add(samp,data,filename,origDataSize,threshAcc+1,1);
im_dc = zeros(N);
res = zeros(N);
for kk = 1:N(3)
    im_dc(:,:,kk) = reshape(ifft2c(squeeze(data_dc(:,:,kk))),[N(1) N(2) 1]);
    res(:,:,kk) = reshape(XFM*(squeeze(im_dc(:,:,kk))),[N(1) N(2) 1]);
end
    
% % Then have the  term be our data that's NOT added together.



% initialize Parameters for reconstruction
param.FT = trans.FT;
param.XFM = XFM;
param.TV = TVOP;
param.data = data;
param.TVWeight =TVWeight;     % TV penalty
param.xfmWeight = xfmWeight;  % L1 wavelet penalty
param.dirWeight = dirWeight;  % directional weight
param.Itnlim = Itnlim;

steps = zeros([size(res) 8]);
tic
for n=1:8
    res = fnlCg(res,param);
    n
    steps(:,:,:,n) = res;
% 	im_hold = XFM'*res(:,:,1);
% % 	%figure(100), imshow(abs(im_res),[]), drawnow
%     figure(3)
%     subplot(2,4,n)
%     imshow(abs(im_hold),[])
end
toc
for i=N(3):-1:1
    im_res(:,:,i) = XFM'*res(:,:,i);
end


diffRMS = rms(im(:)-im_res(:));


mat2mnc(abs(im_res),['/projects/egerek/asalerno/CS-TVDir/TV' num2str(TVWeight) 'DIR' num2str(dirWeight) 'XFM' num2str(xfmWeight) '.mnc'])
