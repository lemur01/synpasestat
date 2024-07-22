function data = PreprocessIm(d, nrs, thRemove, TVgap, gapcorr)
        data = d;
        data.bIm = [];
        data.im = [];
        
        try
        for ch = 1:4
            data.im(:,:,ch) = max(double(data.A{ch, nrs}),[],3);
            if ch ==1
                aux = imgaussfilt(data.im(:,:,ch), 3);
                aux = aux/max(aux(:));
                data.nuclei = imdilate(imopen(im2bw(aux, graythresh(aux)),strel('disk',11)), strel('disk',2));                
            end         
        end
        data.D = cell(4,1);
        data.mx = cell(4,1);
        version = 0;
        if version == 1
            for ch = 1:4
                [data.D{ch} data.Det B1 n1 w1] = FindPeakWav(data.im(:,:,ch), 0.01, 0, 3, 1, 1); 
                %[data.D{ch} data.Det Backgr1 nm1 wm1] = FindPeakModif(data.im(:,:,ch), data.im(:,:,1)>0, 0.01, 0, 3, 2, 0);
                data.bIm = [];
                data.mx{ch} = imregionalmax(data.D{ch}*data.Det);
                data.upper(ch) = max(data.D{ch}(:));
                data.thr(ch) = 0;               
            end            
        else
            for ch = 1:4
                if ch == 1
                    num_iter = 15;
                    delta_t = 1/7;
                    kappa = 200;
                    option = 2;
                   
                    data.D{ch} = fibermetric(data.im(:,:,1)/4096 .*(~imdilate(data.nuclei, strel('disk',3))));
                    % normalize all blue                    
                    [data.blueIm RR Backgr1 nm1 wm1] = FindPeakWav(data.D{ch}, 0.01, 0, 3, 1, 0);
                    RR = bwareaopen(RR,thRemove);
                    data.bIm = RR;
                    data.blueIm = RR;
                    data.mx{1}= RR;
                else %Line 65,66,67
                     [data.D{ch} data.Det(:,:,ch) Backgr1 nm1 wm1] = FindPeakWav(data.im(:,:,ch), 0.01, 0, 3, 1, 0);
                    
                    % support condition
                    data.mx{ch} = imregionalmax(data.D{ch}.*bwareaopen(data.Det(:,:,ch),3));
                    data.upper(ch) = max(data.D{ch}(:));
                    data.thr(ch) = 0;
                end
               
            end
        end
        data.mask = zeros(size(data.im(:,:,1)));
     
    end 
end
function  re = TV(img, TVgap) % not used
    % Run the tensor voting framework        
    T = find_features(img,TVgap); %18
    % Threshold un-important data that may create noise in the
    % output.
    [e1,e2,l1,l2] = convert_tensor_ev(T);
    z = l1-l2;
    l1(z<0.3) = 0;
    l2(z<0.3) = 0;
    T = convert_tensor_ev(e1,e2,l1,l2);
   
    diffL = l1-l2;
    % Run a local maxima algorithm on it to extract curves
    re = calc_ortho_extreme(T,15,pi/8);     
 end
