 
addpath ..\..\bfmatlab\

% parameters
thDist = 2;         % allowed distance to processes (pixels)
thRemove = 3*50;    % process fragments smaller than the threshold are removed
 
TVgap = 9;          % tensor voting parameter - covers gaps smaller than this threshold. Not used
gapcorr = 7;
%%
% input data folder
pt = 'D:\Projects\Helene_eva\testData'
% result folder
respt = 'testResult'
mkdir(respt)
allf = dir(strcat(pt, '\*.lif'));   %all lif files
%% processing
for f = 1:length(allf)
    
    nm = allf(f).name;
    fn = fullfile(respt, nm)
    clear data
    data = OpenImage( nm, pt)

    for nrs = 1:size(data.A,2)
        close all
        clear data.mx data.ch data.bIm data.blueIm
        data = PreprocessIm(data, nrs, removelength, TVgap, gapcorr);
        if ~isempty(data.im)
        
            [idxB idyB] = FindInitPos(data,1);      
            [sky skx] = find(data.bIm>0);            
                        
            %green
            [idxG idyG] = FindInitPos(data,3);    
            posG = RemoveDistantBlue([idxG idyG], data.blueIm, thDist,0);
            
            if (~isempty(posG)) &(size(posG,1)>10)
                % white
                      
                [idxW idyW] = FindInitPos(data,2);        
                posW = RemoveDistantBlue([idxW idyW], data.blueIm, thDist,0);
                
                %red                
                [idxR idyR] = FindInitPos(data,4); 
                posR = RemoveDistantBlue([idxR idyR], data.blueIm, thDist,0);
                
                posRtoW = RemoveDistant(posR, posW, thDist,0); % reds close to white
                posWtoR = RemoveDistant(posW, posR, thDist,0); % whites close to red
                posWtoG = RemoveDistant(posW, posG, thDist,0);
                posRtoG = RemoveDistant(posR, posG, thDist,0);
                posRtoWG = RemoveDistant(posRtoW, posG, thDist,0); 
                posWtoRG = RemoveDistant(posWtoR, posG, thDist,0);
                posGtoR = RemoveDistant(posG, posR, thDist,0);
                posGtoW = RemoveDistant(posG, posW, thDist,0);
                posGtoRW = RemoveDistant(posGtoR, posW, thDist,0);
                
                % onlies
                posGalone = RemoveDistant(posG, posW, thDist,1);
                posGalone = RemoveDistant(posGalone, posR, thDist,1);
                posRalone = RemoveDistant(posR, posW, thDist,1);
                posRalone = RemoveDistant(posRalone, posG, thDist,1);
                posWalone = RemoveDistant(posW, posG, thDist,1);
                posWalone = RemoveDistant(posWalone, posR, thDist,1);
                
                
                
                % save results
                colNames = {'Total_Green','Total_Yellow','Total_Red', 'Only_Green', 'Only_Yellow', 'Only_Red', ...
                'Green_to_Red', 'Green_to_Yellow', 'Yellow_to_Red', 'Yellow_to_Green', 'Red_to_Green', 'Red_to_Yellow',...
                'Green_to_RedYellow','Red_to_GreenYellow','Yellow_to_RedGreen','Thresh_Blue','Thresh_Yellow',...
                'Thresh_Green','Thresh_Red', 'Total_Blue_pixels','Skeleton_cyan_pixels'};
                T = table(length(posG), length(posW),  length(posR), length(posGalone), length(posWalone),  length(posRalone), ...
                length(posGtoR) , length(posGtoW), length(posWtoR), length(posWtoG), length(posRtoG), length(posRtoW), ...
                length(posGtoRW), length(posRtoWG), length(posWtoRG), data.thr(1), data.thr(2), data.thr(3), data.thr(4), sum(data.blueIm(:)), sum(data.bIm(:)),'VariableNames',colNames);
                
                fname = strcat( strtok(fn, '.'),'_', 'Series', num2str(nrs));
                writetable(T, strcat( strtok(fname, '.'), 'Res.txt'), 'WriteRowNames',true, 'Delimiter', ';');
                
                h = figure('units','normalized','outerposition',[0 0 1 1]);  
                clf
                plot(idxB, idyB, 'b.'); hold on % switched y and x
                plot(posG(:,1)-0.5, posG(:,2)-0.5, 'gx', 'LineWidth', 2.5); hold on        
                plot(posW(:,1)+0.5, posW(:,2)+0.5, '+', 'LineWidth', 2.5, 'Color', [1 0.5 0]); hold on
                plot(posR(:,1), posR(:,2), 'ro', 'LineWidth', 2.5); %title(strcat('Pairs: ', num2str(length(g))));
                title(strcat('Measure blue: total pix.' ,num2str(sum(data.blueIm(:))),'- skeleton', num2str(sum(data.bIm(:)))))
                axis equal tight;
                plot(skx, sky, 'c.')  
                set(gca, 'YDir','reverse')
                F = getframe(h);
                [Final, Map] = frame2im(F);      
                imwrite(Final, strcat(strtok(fname, '.'), 'Final.png'), 'png');
                                
                clf
                imshow(data.im(:,:,2),[]); hold on
                plot(posW(:,1), posW(:,2), 'o', 'LineWidth', 2.5, 'Color', [1 0.5 0]); 
                title(strcat('Yellow: ', num2str(length(posW))));
                c = caxis;
                caxis([c(1) c(2)/2])
                F = getframe(gcf);
                [WhiteIm, Map] = frame2im(F);        
                imwrite(WhiteIm, strcat(strtok(fname, '.'), 'Yellow.png'), 'png');        
                
                clf
                imshow(data.im(:,:,4),[]); hold on
                plot(posR(:,1), posR(:,2), 'ro', 'LineWidth', 5); 
                title(strcat('Red: ', num2str(length(posR))));
                c = caxis;
                caxis([c(1) c(2)/2])
                F = getframe(gcf);
                [RedIm, Map] = frame2im(F);
                imwrite(RedIm, strcat(strtok(fname, '.'), 'Red.png'), 'png');     
                                
                clf
                imshow(data.im(:,:,3),[]); hold on
                plot(posG(:,1), posG(:,2), 'gx', 'LineWidth',5); 
                title(strcat('Green: ', num2str(length(posG))));
                c = caxis;
                caxis([c(1) c(2)/2])
                F = getframe(gcf);
                [GreenIm, Map] = frame2im(F);
                imwrite(GreenIm, strcat(strtok(fname, '.'), 'Green.png'), 'png');
                
                clf
                imshow(data.im(:,:,1),[]); hold on
                plot(skx, sky, 'c.','LineWidth', 2.5);      
                title(strcat('Cyan: ', num2str(length(skx))));
                c = caxis;
                caxis([c(1) c(2)/2])
                F = getframe(gcf);
                [CyanIm, Map] = frame2im(F);
                imwrite(CyanIm, strcat(strtok(fname, '.'), 'Cyan.png'), 'png');
            end                      
            
            RGB2(:,:,1) = max(max(0,data.im(:,:,4)/4095)*2, data.im(:,:,2)/4095);
            RGB2(:,:,2) = max(data.im(:,:,3)/4095*2, data.im(:,:,2)/4095);
            RGB2(:,:,3) =max(0, data.im(:,:,1)/4095*2);
            figure(1); 
            clf
            imshow(RGB2,[]); hold on            
            plot(posRtoG(:,1), posRtoG(:,2), 'ro', 'LineWidth', 2);
            %plot(posG(:,1), posG(:,2), 'go', 'LineWidth',2);
            plot(posWtoG(:,1), posWtoG(:,2), 'o', 'LineWidth', 2, 'Color', [1 0.5 0]);
        end
    end
end
