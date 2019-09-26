% this goes round a folder and runs detectBugs4 on it. The results are then
% ordered by cell size across all cells and the average cell intensity
% distributions are plotted.

function [allAverageIntensity,rounded] = detectBugs4Folder(folder,saveFolder)

% if no folder input ask for folder
if nargin < 2
    folder = uigetdir();
    saveFolder = uigetdir();
end

plotting = 0; % this controls whether individual images and their intensity profiles are saved

% read in all the files with the desired ending from the folder
%fileList = dir(fullfile(folder,'*.ome.tif'));
fileList = dir(fullfile(folder,'*.tif'));

% find out how many files there are
numFiles = numel(fileList);

lengths = zeros(1);
unsortedImageList = {};
unsortedImageMaskList = {};
unsortedImageNameList = {};

index = 1;

% loop through each image
for i=1:numFiles
    disp(['looping through images, this is image ',num2str(i)]);
    % run the detect bugs code to get the mean profile data for each cell
    [miniLengths,miniImageList,miniImageMaskList,miniImageNameList] = detectBugs4([folder,'/',fileList(i).name]);
    
%     % make sure no non-identified cells are included
%     deleteIndices = numel(miniImageMaskList) == 0;
%     miniLengths(deleteIndices) = [];
%     miniImageList(deleteIndices) = [];
%     miniImageNameList(deleteIndices) = [];
%     miniImageMaskList(deleteIndices) = [];

    numNewCells = numel(miniLengths);

    lengths(index:index + numNewCells - 1,1) = miniLengths;
    unsortedImageList(index:index + numNewCells - 1,1) = miniImageList;
    unsortedImageMaskList(index:index + numNewCells - 1,1) = miniImageMaskList;
    unsortedImageNameList(index:index + numNewCells - 1,1) = miniImageNameList;
    
    index = index + numNewCells;   
end

altLengths = zeros(numel(unsortedImageList),1);
for i=1:numel(altLengths)
    altLengths(i) = size(unsortedImageList{i},1);
end

[~,sortedLengthIndices] = sort(altLengths);
imageList = unsortedImageList(sortedLengthIndices);
imageMaskList = unsortedImageMaskList(sortedLengthIndices);
imageNameList = unsortedImageNameList(sortedLengthIndices);

allAverageIntensity = NaN(size(imageList{end},1));

close all;

for i=1:numel(imageList)
    
    % convert to double as NaN not defined for uint32
    noZeroIm = double(imageList{i});
    noZeroIm(noZeroIm == 0) = NaN;

    clusters = bwconncomp(imageMaskList{i});

    %assert(clusters.NumObjects == 1);

    maskedIm = noZeroIm;
    maskedIm(clusters.PixelIdxList{1} == 0) = NaN;

    averageIntensity = nanmean(maskedIm,2);
    
    allAverageIntensity(1:numel(averageIntensity),i) = averageIntensity;

    if plotting == 1
        found = strfind(imageNameList{i},'/');
        indexBegin = found(end);
        shortNameFile = imageNameList{i}(indexBegin+1:end);
        fig = figure;
        imshow(maskedIm,[min(maskedIm(:)),max(maskedIm(:))]);  

        saveas(fig,[saveFolder,'/cell-',num2str(i),'-name-',shortNameFile,'.png']);       
        close(fig)

        fig = figure;
        plot(averageIntensity,'linewidth',2)

        saveas(fig,[saveFolder,'/cell-',num2str(i),'-averageIntensity.png']);        
        close(fig)
    end

    
end

% save the results so the whole image analysis does not have to be redone
% every time to change the plots
save([saveFolder,'/results.mat'])

% here the plotting occurs

% double check that there are no sketchy zeros
allAverageIntensity(allAverageIntensity == 0) = NaN;

meanIntensity = mean(allAverageIntensity(~isnan(allAverageIntensity)));
stdIntensity = std(allAverageIntensity(~isnan(allAverageIntensity)));

bottomLimit = meanIntensity - stdIntensity;
topLimit = meanIntensity + stdIntensity;

figure;
imagesc(allAverageIntensity)
caxis([bottomLimit,topLimit]);
colorbar()

rounded = allAverageIntensity;
for j=1:size(rounded,2)
    rounded(:,j) = (rounded(:,j)-min(rounded(:,j)))./(max(rounded(:,j))-min(rounded(:,j)));
end

roundedCentred = NaN(size(rounded));
for j=1:size(roundedCentred,2)
    column = rounded(:,j);
    miniColumn = column(~isnan(column));
    
    halfway = floor(numel(miniColumn)/2);
    maxStart = mean(miniColumn(1:halfway));
    maxEnd = mean(miniColumn(halfway+1:end));
    if maxStart > maxEnd
        miniColumn = miniColumn(end:-1:1);
    end
    
    numExtra = numel(column(isnan(column)));
    halfExtra = floor(numExtra/2);
    roundedCentred(halfExtra+1:halfExtra+numel(miniColumn),j) = miniColumn;
end

figure;
imagesc(roundedCentred)
c = colorbar();
c.Label.String = 'Normalised Fluorescence Intensity';
set(gca,'Fontsize',16)
xticklabels({''})
yticklabels({''})

% peak finding code

countPeaks = zeros(numel(imageList),1);

for k=3:numel(imageList)

    image = imageList{k};
    
    % this finds the peaks
    BW = imregionalmax(image);
    
    linIndexes = find(BW==1);
    
    % convert indexes to subscripts   
    %[x,y] = ind2sub(size(image),linIndexes);
    
    % mask out sketchy bits
    %image(imageMaskList{k} == 0) = NaN;
    
    % remove peaks that are below the average
    
    thresh = double(nanmean(image(:)));
    
    chosenPeaksLin = linIndexes(image(linIndexes)>thresh);
    
    [x,y] = ind2sub(size(image),chosenPeaksLin);
    countPeaks(k) = numel(chosenPeaksLin);
    
    newX = x;
    newY = y;
    
    toBeDeleted = [];

    for i = 1:numel(newX)
        for j=1:numel(newX)
            if i~=j
                if sqrt((newX(i)-newX(j))^2 + (newY(i)-newY(j))^2) < 3
                    %disp('found to delete')
                    if image(newX(j),newY(j)) > image(newX(i),newY(i))
                        %delete i
                        toBeDeleted(numel(toBeDeleted)+1) = i;
                    else
                        toBeDeleted(numel(toBeDeleted)+1) = j;
                    end
                end
            end
        end
    end
    
    newX(toBeDeleted) = [];
    newY(toBeDeleted) = [];
% 
%     figure;
%     imshow(image,[nanmin(image(:)),nanmax(image(:))]);
%     hold on;
%     plot(newY,newX,'ro')
%    
end

figure;
plot(lengths,countPeaks,'o')

% plot the number of peaks as a function of cell length

minLength = min(round(lengths));
maxLength = max(round(lengths));

averagePeaksAtLength = zeros(maxLength-minLength+1,1);
stdPeaksAtLength = zeros(maxLength-minLength+1,1);

for i=minLength:maxLength
    index = i-minLength+1;
    
    averagePeaksAtLength(index) = mean(countPeaks(round(lengths)==i));
    stdPeaksAtLength(index) = std(countPeaks(round(lengths)==i));
end

figure;
errorbar((minLength:maxLength)*0.097,averagePeaksAtLength,stdPeaksAtLength,'x-')
