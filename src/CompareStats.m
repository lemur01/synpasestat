
% creates summary statistics and figures.


clear all
close all

% folder with outputs of Synapse Stats - has to contain lif files with
% 100Pa and 10kPa in their name (the two conditions compared)
folder = 'testResult'

folderBAD = ''
suffix = '.png';
fn = dir(strcat(folder, '\*Res.txt'))
[match,noMatch] = regexp(folder,'\','match','split')
%%

if isempty(folderBAD)
    bad = [];
else
    bad = dir(strcat(folderBAD, '\*.png'))
end

compfname =@(x,y)(strcmp(x(1:strfind(x, 'DIV')),y.name(1:strfind(y.name, 'DIV')) ) && ...
    (str2num(x(strfind(x, '_Series')+7:strfind(x, 'Res')-1) ) == str2num(y.name(strfind(y.name, 'DIV')+3:strfind(y.name, 'DIV')+5)) ));


%%
allT = {}
kk = 1;
for k = 1:length(fn)
   x = fn(k).name
   isbad = 0;
   for ii = 1:size(bad,1)
       if strcmp(x(1:strfind(x, 'DIV')),bad(ii).name(1:strfind(bad(ii).name, 'DIV')) ) && ...
             (str2num(x(strfind(x, '_Series')+7:strfind(x, 'Res')-1) ) == str2num(bad(ii).name(strfind(bad(ii).name, 'DIV')+3:strfind(bad(ii).name, 'DIV')+5)) )

           isbad = 1
       end
   end

   if ~isbad
       k
       T =  readtable(fullfile(fn(k).folder, fn(k).name))
       fnkk{kk} = fn(k).name;
       allT = [allT; T];
       G(kk) = T.Total_Green;
       Sk(kk) = T.Skeleton_cyan_pixels*180/10^6;
       B(kk) = T.Total_Blue_pixels*180/10^6;
       YtoG(kk) =  T.Yellow_to_Green;
       YtoRG(kk) = T.Yellow_to_RedGreen;
       RtoG(kk) = T.Red_to_Green;
       RtoYG(kk) = T.Red_to_GreenYellow;
       Y1(kk) = (YtoG(kk)-YtoRG(kk))/ Sk(kk);
       R1(kk) = (RtoG(kk)-RtoYG(kk))/ Sk(kk);
       is10k(kk) = contains(fn(k).name,'0kPa');
       is100k(kk) = contains(fn(k).name,'150911');
       kk = kk+1;
   end
end
%%

h = figure; % ('units','normalized','outerposition',[0 0 1 1]);  
boxplot(G./B, is10k,'Labels',{'100Pa','10kPa'})
title('Total Green')
F = getframe(h);
[Final, Map] = frame2im(F);      
imwrite(Final, strcat(noMatch{length(noMatch)}, '_GCorr10kPa',suffix), 'png');

h = figure; % ('units','normalized','outerposition',[0 0 1 1]);  
boxplot(B, is10k,'Labels',{'100Pa','10kPa'})
ylabel('mm')
title('Total Cyan')
F = getframe(h);
[Final, Map] = frame2im(F);      
imwrite(Final, strcat(noMatch{length(noMatch)}, '_CCorr10kPa',suffix), 'png');

%%

figure 
fY = (YtoG-YtoRG)./B;
boxplot( fY, is10k,'Labels',{'100Pa','10kPa'})
[h p] = ttest2(fY(is10k==0),fY(is10k==1) )
title(strcat('Total Y (close to G, no R) per mm (p= ', num2str(p),')'))

figure
fR = (RtoG-RtoYG)./B;
boxplot(fR,  is10k,'Labels',{'100Pa','10kPa'})
[h p] = ttest2(fR(is10k==0),fR(is10k==1) )
title(strcat('Total R (close to G, no Y) per mm (p= ', num2str(p),')'))
%%
hh = figure; % ('units','normalized','outerposition',[0 0 1 1]);  
fR = RtoG./ B;
boxplot(fR, is10k,'Labels',{'100Pa','10kPa'})
[h p] = ttest2(fR(is10k==0),fR(is10k==1) )
title(strcat('Total R (close to G) per mm (p= ', num2str(p),')'))
F = getframe(hh);
[Final, Map] = frame2im(F);      

imwrite(Final, strcat(noMatch{length(noMatch)}, '_RtoGCorr10kPa',suffix), 'png');


hh = figure; % ('units','normalized','outerposition',[0 0 1 1]);  
fY = YtoG./ B;
boxplot(fY,  is10k,'Labels',{'100Pa','10kPa'})
[h p] = ttest2(fY(is10k==0),fY(is10k==1) )
title(strcat('Total Y (close to G) per mm (p= ', num2str(p),')'))
F = getframe(hh);
[Final, Map] = frame2im(F);      
imwrite(Final, strcat(noMatch{length(noMatch)}, '_YtoGCorr10kPa', suffix), 'png');
%%
T = table(fnkk',B', (G./B)', fR', fY', 'VariableNames',{'Sample','Processes_mm', 'Neurolignin_per_mm','Gaba_per_mm','Glutamate_per_mm'});
writetable(T,strcat(noMatch{length(noMatch)}, '_Corr10kPa',strtok(suffix,'.'),'.csv'));
%%

figure;
scatter(Sk, fR,12, is10k)
xlabel('Sk')
ylabel('R')
hold on
scatter(Sk(is100k==1), fR(is100k==1),12, 'rx')
colormap winter
title('RED')
figure;
scatter(Sk, fY,12, is10k)
xlabel('Sk')
ylabel('Y');
hold on
scatter(Sk(is100k==1), fY(is100k==1),12, 'rx')
colormap winter
title('YELLOW')
%%
figure; scatter( allT.Skeleton_cyan_pixels, allT.Yellow_to_Green, [], is10k); title('yellow');hold on
scatter(allT.Skeleton_cyan_pixels(is100k==1), allT.Yellow_to_Green(is100k==1),12, 'rx')
figure; scatter( allT.Skeleton_cyan_pixels, allT.Red_to_Green, [], is10k); title('red'); hold on
scatter(allT.Skeleton_cyan_pixels(is100k==1), allT.Red_to_Green(is100k==1),12, 'rx')


%%
id = find(Sk>0.01);
fYid = fY(id);
fRid = fR(id);
[h p] = ttest2(fYid(is10k(id)==0),fYid(is10k(id)==1) )

[h p] = ttest2(fRid(is10k(id)==0),fRid(is10k(id)==1) )
%%

RGB(:,:,1) = min(min(bg/4095,1-data.im(:,:,4)/4095),1-data.im(:,:,2)/4095);
RGB(:,:,2) = min(min(bg/4095,1-data.im(:,:,2)/4095),1-data.im(:,:,3)/4095);
RGB(:,:,3) = min(bg/4095,1-data.im(:,:,4)/4095);
figure; imshow(max(bg,1-RGB),[])

hold on; scatter(posWtoG(:,1),posWtoG(:,2), 15,'y')
hold on; scatter(posRtoG(:,1),posRtoG(:,2), 15, 'm')
hold on; scatter(posG(:,1),posG(:,2), 15, 'g')

RGB(:,:,1) = max(max(0,data.im(:,:,4))/4095*4,data.im(:,:,2)/4095*2);
RGB(:,:,2) = max(max(0/4095,data.im(:,:,2)/4095*2),data.im(:,:,3)/4095*4);
RGB(:,:,3) = max(max(0/4095,data.im(:,:,4)/4095*4), data.im(:,:,1)/4095);
figure; imshow(RGB,[])


hold on; plot(posWtoG(:,1),posWtoG(:,2),'yo', 'MarkerSize',15)
hold on; plot(posRtoG(:,1),posRtoG(:,2), 'mo', 'MarkerSize',15)
hold on; plot(posG(:,1),posG(:,2), 'go', 'MarkerSize',15)


hold on; scatter(posWtoG(:,1),posWtoG(:,2), 15,'y')
hold on; scatter(posRtoG(:,1),posRtoG(:,2), 15, 'm')
hold on; scatter(posG(:,1),posG(:,2), 15, 'g')



hold on; plot(posWtoG(:,1),posWtoG(:,2),'yo', 'MarkerSize',15)
hold on; plot(posRtoG(:,1),posRtoG(:,2), 'mo', 'MarkerSize',15)
hold on; plot(posG(:,1),posG(:,2), 'go', 'MarkerSize',15)
%%
bg = 4095-data.im(:,:,1);

RGB(:,:,1) = max(max(0,data.im(:,:,4)/4095), data.im(:,:,2)/4095);
RGB(:,:,2) = max(0/4095,data.im(:,:,3)/4095, data.im(:,:,2)/4095);
RGB(:,:,3) =max(0/4095,data.im(:,:,4)/4095);

figure; imshow(1-max(bg,RGB),[])
