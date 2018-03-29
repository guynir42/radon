function R = partial2final(P)
% extracts the final Radon image from a cell array of partial transforms. 

    if nargin==0, help('radon.partial2final'); return; end
    
    if isempty(P)
        R = [];
    elseif iscell(P)
        R = permute(P{end}, [1,3,2]);
    else
        R = permute(P, [1,3,2]);
    end

end
