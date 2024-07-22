% Remove points that are not close to the mask (processes in the blue chanel - blueIm)
% Input:
% pos - position of molecules
% thDist - threshold for the distance to processes
% blueIm - mask for processes
function newpos = RemoveDistantBlue(pos, blueIm, thDist, removeclose)
    if ~isempty(pos) 
    d = bwdist(blueIm);
     id = find(d(sub2ind(size(blueIm),pos(:,2),pos(:,1)))>thDist);

     if removeclose
         newpos = pos(id,:);
     else
        pos(id,:) = [];
        newpos = pos;   
     end
    else
        newpos = [];
    end
end 