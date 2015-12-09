function [grad,param] = wGradient(x,param)

gradXFM = 0;
gradTV = 0;

gradObj = gOBJ(x,param);

if param.xfmWeight
    gradXFM = gXFM(x,param);
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

