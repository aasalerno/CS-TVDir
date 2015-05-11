function ph = phCalc(data,rl)
% This code will calculate the phase for the data -- it will use different
% styles of calculation depending on if the data is real or complex.
%
% data - this is the dataset
% rl   - Logical. Is the data real?

if nargin < 2; rl = isreal(data); end
data = ifftshift(ifft2(data));
N = size(data);

if rl
    disp('Lustig''s way')
    phmask = zpad(hamming(6)*hamming(6)',N(1),N(2)); %mask to grab center frequency
    phmask = phmask/max(phmask(:));
    ph = exp(1i*angle((ifft2c(data.*phmask))));
else
    disp('Brian''s Way')
%     F = fspecial('gaussian',[5 5],2);
%     filtdata = imfilter(data,F,'same');
%     imshow(filtdata);
%    ph = conj(filtdata)./abs(filtdata);
    ph = conj(data)./abs(data);
end