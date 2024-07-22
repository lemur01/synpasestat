
% Read in series from lif file
% if not found or not 4 channels per FOV returns error
% Output: data
function data = OpenImage( nm, pt)

        data = struct;
        
        try      
            fn = fullfile(pt,nm)
            %imTIF = Tiff(strcat(pt,nm),'r');
            dlif = bfopen(fn);
            data.nrseries = size(dlif,1);
            nrZ = 1; 
            nrch = 4;
            data.A = cell(nrch,data.nrseries);
            data.im = zeros(size(dlif{1,1}{1,1},1),size(dlif{1,1}{1,1},2),4);
            data.mask = zeros(size(data.im(:,:,1)));
            ct = 1;
            for j = 1:data.nrseries
            for i = 1:size(dlif{j,1},1)
                z = str2num(dlif{j,1}{i,2}(findstr('Z=', dlif{j,1}{i,2})+2));
                if isempty(z)
                    break;
                end
                ch = str2num(dlif{j,1}{i,2}(findstr('C=', dlif{j,1}{i,2})+2));            
                data.A{ch,j}(:,:,z) = dlif{j,1}{i,1};
                ct = ct+1;
            end
        end
        catch
            disp('File selection ore read failed...')
        end

end
