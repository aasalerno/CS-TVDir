function samp = genSampToAdd(imSize,pctg,radius,nmins,filename)
%
% samp = genSampToAdd(imSize,pctg,radius,nmins,filename)
%
% samp is out sampling pattern that we wish to obtain
% imSize - the image size that we are working with
% pctg - how much data we want to collect for each dataset
% radius - How much of the centre we want fully sampled
% nmins - How many directions we want to be grouped together
% 
% Figure out how the vectors are going to work -- what combinations should
% be used
%
%
% This is the function that cdoes all of the combinations from 1:nDirs in
% our datatset

if nargin < 5
    nmins = 5;
end

dirs = load(filename);

% Make sure all of the data is on one side -- this may not need to happen,
% since we're doing the next step

n = length(dirs);
r = zeros(n);

% Calculate radius -- this is being done using an electrostatics approach
% We figure out which distance is shorter, between the two directions, or
% the direction and the negative of the other
for i = 1:n
    for j = 1:n
        r(i,j) = sqrt(sum((dirs(i,:)-dirs(j,:)).^2));
        r(i,j) = min(sqrt(sum((-dirs(i,:)-dirs(j,:)).^2)),r(i,j));
    end
end

invR = 1./(r+eps);

% Find all possible combos
k = floor(n*pctg); % How many vectors can be used for the given p
combs = combnk(1:n,k); % Different k-vectored combinations
vecs = combnk(1:k,2); % What are all the vector combos
engStart = zeros(1,length(combs));

% Figure out energies for said combos -- that is, the energy if the point
% is in each one of the datasets for the vectors as stated in combs
for i = 1:length(combs)
    for j = 1:length(vecs)
        engStart(i) = engStart(i) + invR(combs(i,vecs(j,1)),combs(i,vecs(j,2)));
    end
end
%save(['eng-' num2str(n) 'dirs-' num2str(nmins) 'mins.mat'],'engStart');


% Ensure that we don't have major differences in how many counts for each
% vector there are
[eng,ind] = sort(engStart);
vecsind = combs(ind,:)';
locs = zeros(n,nmins);
vecsMin = zeros(k,n*nmins);


% Figure out how many times each number was chosen
for i = 1:n
    [~,locs(i,:)] = find(vecsind==i,nmins);
    vecsMin(:,nmins*(i-1)+1:nmins*(i)) = vecsind(:,locs(i,:));
end
vecsMin = unique(vecsMin','rows');
amts = zeros(1,n);
for i = 1:n
    amts(i) = numel(vecsMin(vecsMin==i));
end
[cts,srt] = sort(amts);
qEng = quantile(eng,4);

% if theres a huge difference, tack more of the lower counts on, but make
% sure that we aren't hitting too high energy sets
while cts(end)/cts(1) >= 1.1
    
    [~,indx] = ismember(sort(srt(1:k)),vecsind','rows');
    
    if eng(indx) < qEng(2)
        vecsMin(end+1,:) = sort(srt(1:k));
    else
        %         [~,indx] = ismember(sort(srt([1:6,8])),vecsind','rows');
        %         if eng(indx) < qEng(2)
        %             vecsMin(end+1,:) = sort(srt([1:6,8]));
        %         else
        while eng(indx) >= qEng(2)
            arr = [];
            while length(arr) ~= k
                st = ceil(0.7*n*rand(1));
                if ~any(arr==st)
                    arr(end+1) = st;
                end
            end
            [~,indx] = ismember(sort(srt(arr)),vecsind','rows');
        end
        vecsMin(end+1,:) = sort(srt(arr));
        %         end
    end
    
    for i = 1:30
        amts(i) = numel(vecsMin(vecsMin==i));
    end
    
    [cts,srt] = sort(amts);
    %srthold = srthold(srt);
end
% for i = 1:30
%     amts(i) = numel(vecsMin(vecsMin==i));
% end


% -------------------------------------------------------------------
% Now we generate the sampling pattern!
%
[x,y] = meshgrid(linspace(-1,1,imSize(2)),linspace(-1,1,imSize(1)));
r = sqrt(x.^2+y.^2);
if 1 == 0
    r = r/max(abs(r(:)));
end

%r(r > 1 | r < radius) = 0;
[rows, cols] = find(r <= 1 & r > radius);
[rx,ry] = find(r <= radius);

samp = zeros([imSize,n]);
% samp(rx,ry,:) = 1;
nSets = size(vecsMin,1);

for i = 1:numel(rows)
    val(i) = ceil(nSets*rand(1));
    choice = vecsMin(val(i),:);
    
    samp(rows(i),cols(i),choice) = 1;
end

for i = 1:numel(rx)
    samp(rx(i),ry(i),:) = 1;
end
% 
% if nargin > 6 | endSize ~= imSize
%     samp_final = zeros([endSize,n]);
%     for i = 1:n
%         samp_final(:,:,i) = zpad(samp(:,:,i),endSize(1),endSize(2));
%     end
%     samp = samp_final;
% end
