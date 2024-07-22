% creates lists of positions from local maxima in the images
function [idx idy] = FindInitPos(data, ch)
     
        if ch~=1

            pmax = data.mx{ch};     
            M = pmax-bwareaopen(pmax,5);

        else
  
            data.mx{ch} = data.bIm;
            M =  data.mx{ch} ;

        end
        
        [idy idx] =  find((M>0));
    end
