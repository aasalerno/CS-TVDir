function grad = gTV(x,params,dDirx)
% compute gradient of TV operator
%
% Edited from the original to also include the information from the
% directional term as well
%
% The directional term's gradient will be multiplied by a factor that
% compares the gradient term and the directional term -- it will be
% dirWeight/TVWeight

p = params.pNorm;


grad = zeros(size(x));

for kk = 1:size(x,3)
    % TV
    x1 = squeeze(x(:,:,kk));
    Dx = params.TV*(params.XFM'*x1);
    G = p*Dx.*(Dx.*conj(Dx) + params.l1Smooth).^(p/2-1);
    grad(:,:,kk) = (params.TV'*G);
end

if params.dirWeight
    Ahat = params.dirInfo.Ahat;
    inds = params.dirInfo.inds;
    indsPos = params.dirInfo.indsPos;
    %indsNeg = params.dirInfo.indsNeg;
    % Preallocate for speed and having recursive addition
    gradHold = zeros(size(x));
    betaHold = zeros([size(x) 3])
    
    for kk = 1:size(x,3)

        
        % Have the work done with respect to addition (when it's the
        % positive term in the calculation)
%         for i = 1:numel(indsPos{kk})
%             gradHold(:,:,i) = gradHold(:,:,i) + x(:,:,kk) - x(:,:,inds(indsPos{kk}(i)));
%         end
        
        % Have the calculations done for the negative sets. My only worry
        % with this part is that this will lead to a total cancellation for
        % the ones that are closest together -- does this part need to
        % exist?
%         for i = 1:numel(indsNeg{kk})
%             gradHold(:,:,i) = gradHold(:,:,i) +x(:,:,inds(indsNeg{kk}(i))) - x(:,:,kk);
%         end
        [rows,cols] = find(inds == kk);

        for i = 1:size(x,1)
            for j = 1:size(x,2)
                
                gradHold(i,j,kk,:) = gradHold(i,j,kk,:) - reshape(Ahat(:,:,kk)*squeeze())
                
            end
        end
        
        
        
    end
    grad = grad + gradHold;
end

for kk = 1:size(x,3)
    grad(:,:,kk) = params.XFM*grad(:,:,kk);
end