function signalGrid = interpolateDeadElectrodes(signalGrid, deadElectrodes)
% Interpolate some positions to be the median value of surrounding nodes.
% Input signalGrid can be a 2D matrix with dimension [X Y Z] or a 4D matrix
% with dimension [X Y Z T], where t gives time. badElectrodes is
% an array containing the indices of the dead electrodes for the current
% data set.
%

sg = size(signalGrid);
signalGrid = signalGrid(:,:,:,:);

[nr, nc, nz, nt] = size(signalGrid);

% If no dead electrodes are supplied, interpolate over corners and any
% channels with NaNs
if nargin == 1
    nanChans = find(any(isnan(signalGrid), 4));
    deadElectrodes = nanChans;
end

signalGrid = reshape(signalGrid, [], nt);

ii = 1;
ndead = length(deadElectrodes(:));
fprintf('channel %d\n', ndead);

deadElectrodes = deadElectrodes(:)';
tic;

while ii <= length(deadElectrodes(:))
    id = deadElectrodes(ii);
    ii = ii+1;
  
    
    % Define surrounding nodes in the form [down up left right]
    surrounding = [id-nr-1, id-nr, id-nr+1, id-1, id+1 id+nr-1 id+nr,id+nr+1];
    surrounding = [surrounding, surrounding-nr*nc, surrounding+nr*nc];
    
    % Remove nodes that wrap around the edges
%     if mod(id,nr) == 0
%         surrounding(1) = [];
%     elseif mod(id,nr) == 1
%         surrounding(2) = [];
%     end

    % Remove out of bounds nodes
    outOfBounds = surrounding<1 | surrounding>nr*nc*nz;
    surrounding(outOfBounds) = [];
    
    % Remove surrounding nodes that are corners or other bad nodes unless
    % these are the only available nodes
    badSurround = ismember(surrounding,deadElectrodes);
    if ii<=ndead || sum(badSurround) < length(surrounding)
        surrounding(badSurround) = [];
    end
    
    % If there are no valid surroudning electrodes, push it to the back of
    % the queue
    if nargin == 2
        if isempty(surrounding)
            deadElectrodes = [deadElectrodes, id];
            continue
        end
    end
    
    % Set to median value of surrounding electrodes
    if ~isempty(surrounding)
        signalGrid(id, :) = nanmean(signalGrid(surrounding, :), 1);
    end
end

% Reshape to original dimensions
signalGrid = reshape(signalGrid, sg);
toc;

end