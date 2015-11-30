function [grad,param] = gTV(x,param)
% compute gradient of TV operator
%
% Edited from the original to also include the information from the
% directional term as well
%
% The directional term's gradient will be multiplied by a factor that
% compares the gradient term and the directional term -- it will be
% dirWeight/TVWeight

p = param.pNorm;


grad = zeros(size(x));

for kk = 1:size(x,3)
    % TV
    x1 = squeeze(x(:,:,kk));
    Dx = param.TV*(param.XFM'*x1);
    G = p*Dx.*(Dx.*conj(Dx) + param.l1Smooth).^(p/2-1);
    grad(:,:,kk) = (param.TV'*G);
end

if param.dirWeight
    %     Ahat = params.dirInfo.Ahat;
    M = param.dirInfo.M;
    Ause = param.dirInfo.Ause;
    dI = param.dirInfo.dI;
    dDirx = param.dDirx;
    %inds = params.dirInfo.inds;
    %indsPos = params.dirInfo.indsPos;
    %indsNeg = params.dirInfo.indsNeg;
    
    % Preallocate for speed and having recursive addition
    gradHold = zeros(size(x));
    %     betaHold = zeros([size(x) 3]);
    
    for kk = 1:size(x,3)
       
        for d = 1:length(Ause{kk})
            colUse = Ause{kk}(d); % Which column are we looking at 
            dIM = dI(:,colUse,kk)'*M(:,:,colUse); % dI(:,d,kk) gives the d'th column in the kk'th set -- that is it tells us the column of where the image from direction kk shows up
            for i = 1:size(x,1)
                for j = 1:size(x,2)
                    
                    gradHold(i,j,kk) = gradHold(i,j,kk) + 2*dIM*dDirx(:, colUse, i,j); % Pixel i,j in our dataset, all of the differences (column), for column d which is told to us by the above
                    
                end
            end
        end
       
        
    end
    grad = grad + gradHold;
    param.gDir = gradHold;
end

for kk = 1:size(x,3)
    grad(:,:,kk) = param.XFM*grad(:,:,kk);
end