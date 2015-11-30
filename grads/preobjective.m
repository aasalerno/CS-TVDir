function [FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx, XFMtx, XFMtdx] = preobjective(x, dx, params)

% precalculates transforms to make line search cheap

% AS - preallocation and for looping
FTXFMtx = zeros(size(x));
FTXFMtdx = zeros(size(x));

DXFMtx = zeros([size(x,1) size(x,2) 2 size(x,3)]);
DXFMtdx = DXFMtx;

XFMtx = zeros(size(x));
XFMtdx = zeros(size(x));


for kk = 1:size(x,3)
    %x1 = squeeze(x(:,:,kk));
    %dx1 = squeeze(dx(:,:,kk));
    FTXFMtx(:,:,kk) = params.FT{kk}*(params.XFM'*x(:,:,kk));
    FTXFMtdx(:,:,kk) = params.FT{kk}*(params.XFM'*dx(:,:,kk));
end


if params.TVWeight
    for kk = 1:size(x,3)
        %x1 = squeeze(x(:,:,kk));
        %dx1 = squeeze(dx(:,:,kk));
        DXFMtx(:,:,:,kk) = params.TV*(params.XFM'*x(:,:,kk));
        DXFMtdx(:,:,:,kk) = params.TV*(params.XFM'*dx(:,:,kk));
        XFMtx(:,:,kk) = params.XFM'*(x(:,:,kk));
    end
    % else
    %     DXFMtx = 0;
    %     DXFMtdx = 0;
end

% if isfield(params,'dirWeight') && params.dirWeight ~= 0
%     for kk = 1:size(x,3)
% %         x1 = squeeze(x(:,:,kk));
% %         dx1 = squeeze(dx(:,:,kk));
%         XFMtx(:,:,kk) = (params.XFM'*x(:,:,kk));
%         XFMtdx(:,:,kk) = (params.XFM'*dx(:,:,kk));
%     end
%     
%     % This part is what will be used in the new version of the TV operator
%     % Since we need to build the vectors in a very specific way, for a
%     % single direction and pixel, we can get the information we need from this matrix
%     % as squeeze(dDirx(i,j,k,:))
%     dDirx = zeros([size(x,1) size(x,2) size(params.dirInfo.inds)]);
%     
%     for i = 1:size(params.dirInfo.inds,1)
%         for j = 1:size(params.dirInfo.inds,2)
%             dDirx(:,:,i,j) = params.XFM'*(x(:,:,params.dirInfo.inds(i,j))-x(:,:,i));
%         end
%     end
%     
%         
%     
% end