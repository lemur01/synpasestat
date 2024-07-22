% Remove points that are not close to reference points
% Input:
% pos - position of molecules
% poref - position of reference molecules
% thDist - threshold for the distance to processes
  
function newpos = RemoveDistant(pos, posRef, thDist, removeclose)
    
     [posP distP] = knnsearch(posRef, pos);
     id = find(distP>thDist);
     if removeclose
         newpos = pos(id,:);
     else
        pos(id,:) = [];
        newpos = pos;   
     end
end % Remove points that are not close to any point of a specific color
