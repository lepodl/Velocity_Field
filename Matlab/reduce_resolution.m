function re_signal = reduce_resolution(signalGrid,reduce_size)

signalGrid(signalGrid==0)=nan;
sg = size(signalGrid);
if nargin < 2
    reduce_size = ceil(sg(1:end-1) / 2);
end

[nr, nc, nz, nt] = size(signalGrid);
re_sg = [reduce_size, nt];
re_signal = zeros(re_sg);
re_signal = reshape(re_signal, [], nt);
signalGrid = reshape(signalGrid, [], nt);

sphere_sz = sg(1:end-1);
fprintf('reduce size %d\n', [reduce_size, nt]);
tic;

ll_id = 0;
for k=1:2:nz
    for j=1:2:nc
        for i=1:2:nr
            ll_id = ll_id + 1;
            id = sub2ind(sphere_sz, i, j, k);
            surrounding = [id-nr-1, id-nr, id-nr+1, id-1, id,  id+1 id+nr-1 id+nr,id+nr+1];
            surrounding = [surrounding, surrounding-nr*nc, surrounding+nr*nc];
            outOfBounds = surrounding<1 | surrounding>nr*nc*nz;
            surrounding(outOfBounds) = [];
            re_signal(ll_id, :) = nanmean(signalGrid(surrounding, :), 1);
        end
    end
end
if ll_id ~= size(re_signal, 1)
    error(['shape not correct' num2str(ll_id) '!=' size(re_signal, 1)]);
end
                     
re_signal = reshape(re_signal, re_sg);

end

