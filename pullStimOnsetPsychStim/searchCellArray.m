function targetIndex = searchCellArray(cellArray,target)
% searches a cell array for a string and returns the index, or 0 if not
% found.


%%%%%% code written by Adam Ranson 2014

targetIndex  = find(ismember(cellArray,target));
if isempty(targetIndex)
    targetIndex = 0;
end
end

