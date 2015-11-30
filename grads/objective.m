function [res, obj, RMS] = objective(FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx, XFMtx, XFMtdx, x,dx,t, param)
%calculated the objective function

p = param.pNorm;
N = size(FTXFMtx);

if length(N) == 2
    N = [N 1];
end

obj = zeros(N(3),1);
objPre = FTXFMtx + t*FTXFMtdx - param.data;


for kk = 1:N(3)
    mul = objPre(:,:,kk);
    obj(kk) = (mul(:)'*mul(:))';
end

if param.TVWeight
    %     TV = zeros(N(1)*N(2)*2,N(3));
    %     for kk = 1:N(3)
    %         DXFMtx1 = DXFMtx(:,:,:,kk);
    %         DXFMtdx1 = DXFMtdx(:,:,:,kk);
    %         w = DXFMtx1(:) + t*DXFMtdx1(:);
    %         TV(:,kk) = (w.*conj(w)+params.l1Smooth).^(p/2);
    %     end
    w = reshape(DXFMtx + t*DXFMtdx,[N(1)*N(2)*2,N(3)]);
    dir = reshape(XFMtx + t*param.gDir,[N(1)*N(2),N(3)]);
    if param.dirWeight
        TV = (sum(w(:).*conj(w(:))+param.l1Smooth) + param.dirWeight.*(sum(dir(:).*conj(dir(:))))/(sum(cellfun('prodofsize', param.dirInfo.Ause)))).^(p/2);
    else
        TV = (w.*conj(w)+param.l1Smooth).^(p/2);
    end
else
    TV = 0;
end

if param.xfmWeight
    %   XFM = zeros(N(1)*N(2),N(3));
    %     for kk = 1:N(3)
    %         x1 = x(:,:,kk);
    %         dx1 = dx(:,:,kk);
    %         w = x1(:) + t*dx1(:);
    %         XFM(:,kk) = (w.*conj(w)+params.l1Smooth).^(p/2);
    %     end
    w = reshape(x + t*dx,[N(1)*N(2),N(3)]);
    XFM = (w.*conj(w)+param.l1Smooth).^(p/2);
else
    XFM=0;
end

% if isfield(params,'dirWeight') && params.dirWeight~=0
%     wgt = params.dirPairWeight;
%     dirPair = params.dirPair;
%     n = length(wgt);
%     dir = zeros(1,n);
%     for kk = 1:n
%         Xi = XFMtx(:,:,dirPair(kk,1));
%         DXi = XFMtdx(:,:,dirPair(kk,1));
%         Xj = XFMtx(:,:,dirPair(kk,2));
%         DXj = XFMtdx(:,:,dirPair(kk,2));
%         val = (Xi(:) + t*DXi(:) - Xj(:) - t*DXj(:));
%         dir(kk) = wgt(kk).*(val'*val);
%     end
%     
%     % Separate by diffusion direction?
%     if 1 == 0
%         n = size(dirPair,1);
%         nums = (1:n)';
%         dirDiff = zeros(1,N(3));
%         for kk = 1:max(dirPair(:))
%             inRow = nums(any(kk == dirPair,2).*nums ~= 0);
%             dirDiff(kk) = sum(dir(inRow),2);
%         end
%         
%         % How shall I sum them?
%         % Each diff gets it's own
%         res = sum(obj,1) + sum(params.xfmWeight(:).*XFM,1) + sum(params.TVWeight(:).*TV,1)...
%             + params.dirWeight(:).*dirDiff;
%         
%         RMS = sqrt(obj/sum(abs(params.data(:))>0));
%     else
        
%         TV = sum(TV.*params.TVWeight(:));
%         XFM = sum(XFM.*params.xfmWeight(:));
%         dirDiff = sum(dir.*params.dirWeight(:));
%         RMS = sqrt(obj/sum(abs(params.data(:))>0));
%         
%         
%         res = sum(obj(:) + TV(:) + XFM(:) + dirDiff(:)/n);
% %         fprintf('Objective: %1.3e \n',(sum(obj(:))))
% %         fprintf('XFM: %1.3e \n',(sum(XFM(:))))
% %         fprintf('TV: %1.3e \n',(sum(TV(:))))
% %         fprintf('Directions: %1.3e \n',(sum(dirDiff(:))))
%         obj = sum(obj);
%         RMS = sum(RMS);
%     end
%     
    
% else % If dirWeight doesn't exist.
    
    TV = sum(TV.*param.TVWeight(:));
    XFM = sum(XFM.*param.xfmWeight(:));
    RMS = sum(sqrt(obj/sum(abs(param.data(:))>0)));
    fprintf('Objective: %1.3e \n',(sum(obj(:))))
    fprintf('XFM: %1.3e \n',(sum(XFM(:))))
    fprintf('TV: %1.3e \n',(sum(TV(:))))
    
    res = sum(obj(:) + TV(:) + XFM(:));
% end
