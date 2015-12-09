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

for dir = 1:size(x,3)
    % TV
    x1 = squeeze(x(:,:,dir));
    Dx = param.TV*(param.XFM'*x1);
    G = p*Dx.*(Dx.*conj(Dx) + param.l1Smooth).^(p/2-1);
    grad(:,:,dir) = (param.TV'*G);
end

if param.dirWeight
    
    % Here, we build a matrix of differences for each of the combinations that
    % we get in the indicies.
    %
    % It's a bit slow, as we need to do all 30 directions amd there are ~4 sets
    % each. So it overall takes about a second after some cleaning up.
    if isfield(param,'dirWeight') && param.dirWeight ~= 0
        dDirx = zeros([size(x,1) size(x,2) size(param.dirInfo.inds)]);
        
        for ind_q = 1:size(param.dirInfo.inds,1)
            for ind_r = 1:size(param.dirInfo.inds,2)
                dDirx(:,:,ind_q,ind_r) = param.XFM'*(x(:,:,param.dirInfo.inds(ind_q,ind_r))-x(:,:,ind_q));
            end
        end
        
%         dDirx = permute(dDirx,[4 3 1 2]); % Put it such that it will produce a column matrix for :,columnWeNeed,i,j
        
    end
    %     Ahat = params.dirInfo.Ahat;
    %M = param.dirInfo.M; % Ahat'*Ahat
    %dI = param.dirInfo.dI; % Either -1, 0, or 1.
    %inds = params.dirInfo.inds;
    %indsPos = params.dirInfo.indsPos;
    %indsNeg = params.dirInfo.indsNeg;
    
    Ause = param.dirInfo.Ause; % Which 'A's do we need (theres a different one for each grad direction)
    dIM = param.dirInfo.dIM;
    
    % Preallocate for speed and having recursive addition
    gradHold = zeros(size(x));
        
    for dir = 1:size(x,3)
        
        % Comb is which "combination" we're looking at -- that is the
        % combinations that are closest to dir
%         for comb = 1:length(Ause{dir})
%             
%             columnOfData = Ause{dir}(comb); % Which column are we looking at
%             
%             % Run through each pixel 
%             for pixel_i = 1:size(x,1)
%                 for pixel_j = 1:size(x,2)
%                     
%                     gradHold(pixel_i,pixel_j,dir) = gradHold(pixel_i,pixel_j,dir)...
%                                                    + 2*dIM(columnOfData,:,dir)*dDirx(:, columnOfData, pixel_i,pixel_j); % Pixel i,j in our dataset, all of the differences (column), for column d which is told to us by the above
%                     
%                 end
%             end
%         end
        
        % OR %
        
        for comb = 1:length(Ause{dir})
            
            columnOfData = Ause{dir}(comb); % Which column are we looking at
            
            for qr = 1:size(param.dirInfo.inds,2)
            
                gradHold(:,:,dir) = dIM(dir,qr,columnOfData)*dDirx(:,:,dir,qr) + gradHold(:,:,dir);
            end
            
        end
    end
    
    % Add the original gTV and the gTV from our direction term
    grad = grad + gradHold;
    
    % Save it for future looking if we want
    %param.gDir = gradHold;
    
    % Convert it to XFM space
    grad(:,:,dir) = param.XFM*grad(:,:,dir);
end