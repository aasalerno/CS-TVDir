function [grad,param] = wGradient(x,param)

gradXFM = 0;
gradTV = 0;

gradObj = gOBJ(x,param);

if param.xfmWeight
    gradXFM = gXFM(x,param);
end

if isfield(param,'dirWeight') && param.dirWeight ~= 0
    dDirx = zeros([size(x,1) size(x,2) size(param.dirInfo.inds)]);
    
    for i = 1:size(param.dirInfo.inds,1)
        for j = 1:size(param.dirInfo.inds,2)
            dDirx(:,:,i,j) = param.XFM'*(x(:,:,param.dirInfo.inds(i,j))-x(:,:,i));
        end
    end
    
    param.dDirx = permute(dDirx,[4 3 1 2]); % Put it such that it will produce a column matrix for :,columnWeNeed,i,j
    
    for i = 1:size(x,3)
        
    end
end

if param.TVWeight
    if param.dirWeight
        [gradTV, param] = gTV(x,param);
    else
        gradTV = gTV(x,param);
    end
end

% AAS
% if isfield(params,'dirWeight') && params.dirWeight
%    gradDir = gDir(x,params);
%    grad = (gradObj + params.xfmWeight.*gradXFM + params.TVWeight.*gradTV ...
%                 + params.dirWeight.*gradDir);
% %    gdir = gradDir(:,:,1);
% %    save(['/projects/muisjes/asalerno/CS/data/dirArtefactData/dx_grad.' num2str(cnt) '.mat'],'gdir');
% %    cnt = cnt + 1;
% else
    %NORMAL
   grad = (gradObj + param.xfmWeight.*gradXFM + param.TVWeight.*gradTV);
% end

