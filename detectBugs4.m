% updated version of detectBugs developed by Stephen Cross.

% main file. Binarises each image, identifies cells based on an intensity
% and pixel size threshold and then crops each cell out and returns the
% cropped cell image and the length of the cell

% INPUTS
% filename; a full filepath to a tiff image file

% OUTPUTS
% lengths; 1xC vector, length of each cell
% images; 1xC cell array, cropped image of each cell
% imageMasks; 1xC cell array, masks of each cell region in the cropped image
% imageNames; 1xC cell array, the name of the image file for each cell

function [lengths,images,imageMasks,imageNames] = detectBugs4(filename)

min_px_in_bug = 300; %The minimum number of pixels per detected cell, which
%should remove the majority of background noise
max_px_in_bug = 800; %The maximum area a detected cell can take, it
%should remove any large "blobs" of dirt on the sample, or groups of cells

disp (filename); %Displaying the current image file being worked on

%Loading the current image to "loaded_im"
image = tiffread(filename);
loaded_im = image.data();

% show the original image
figure;imshow(loaded_im,[min(loaded_im(:)) max(loaded_im(:))]);title('original')
original_im = loaded_im; % save the original image for later

% we want to threshold the image. The images are close to background, and
% the default threshold fails. So we calculate our own threshold as the avg
% + N*s.d. for some N

avg = mean(loaded_im(:)); % calculate avg intensity

sd = std(double(loaded_im(:))); % SD of intensity (double cause type issues)

%thresh = (double(avg)+2.5*double(sd))/double(2^32); % threshold
%thresh = (double(avg)+2*double(sd))/double(2^32); % threshold

thresh = (double(avg)+4*double(sd))/double(2^16); % threshold

% binarise the image
filtered_im = imbinarize(loaded_im,thresh);
filtered_im = imfill(filtered_im,'holes'); % deal a bit with speckly cells

figure;imshow(filtered_im)

% identify individual cells by finding connected components in the
% filtered_im
clusters = bwconncomp(filtered_im);

% remove clusters that are too big or too small to be bugs
numPixels = cellfun(@numel,clusters.PixelIdxList);

for i=1:length(clusters.PixelIdxList)
    if numPixels(i) > max_px_in_bug || numPixels(i) < min_px_in_bug
        filtered_im(clusters.PixelIdxList{i})=0;
    end
end

% find clusters again now we may have removed some
clusters = bwconncomp(filtered_im);

% identify information about each region
s = regionprops(clusters,'BoundingBox','MajorAxisLength','MinorAxisLength','Orientation','Area','PixelList','Extrema','Centroid');

% extrema [top-left top-right right-top right-bottom bottom-right bottom-left left-bottom left-top]

numAnalysedBugs = numel(s);

fig = figure;
hold on;
imshow(original_im,[min(loaded_im(:)) max(loaded_im(:))]);title('chosen bugs')

for k = 1:numAnalysedBugs
    rectangle('Position', s(k).BoundingBox,'EdgeColor','r', 'LineWidth', 1)
    text(s(k).BoundingBox(1)-5,s(k).BoundingBox(2)-5,num2str(k),'Color','g');
    disp(['Num bug: ',num2str(k),' Long axis: ',num2str(s(k).MajorAxisLength),' Short axis: ',num2str(s(k).MinorAxisLength)])
end

%This cade saves the analysed bug image if desired
%saving it with the suffix "analysedBugs"
%set(disp_fig,'position',screen_size);

% print('-dtiff', '-r200', strcat(filename, num2str(curr_num),...
%     'analysedBugs','.','tif'));

%close(fig)

% we do not preallocate size as some cells may be ignored and there are
% too few cells per image to make it worthwhile memory-wise
images = {};
imageMasks = {};
%points = zeros(100,100);
lengths = [];
%peaks = [];

imageNames = cell(1);

% do the analysis of each cell. The peak finding code is currently
% commented out.

index = 1;
for k=1:numAnalysedBugs
    
    % NEW APPROACH
    % first rotate the image by the theta that will make this cell vertical
    theta = s(k).Orientation;
    
    rotIm = imrotate(loaded_im,90-theta);
    % rotate the masked image too so we can find connected components
    % again
    rotImMask = imrotate(filtered_im,90-theta);
    
    % also rotate the centroid of the current cell, so that we can identify
    % which component in the rotated image corresponds to our cell
    % we have to turn the rotation angle (90-theta) into radians to do a
    % matrix rotation
    thetaRad = ((90-theta)*pi)/180;
    % NOTE: this form of the rotation matrix takes into account that the y
    % axis is flipped in images
    rotMatrix = [[cos(thetaRad),sin(thetaRad)];[-sin(thetaRad),cos(thetaRad)]];

    % imrotate rotates around the image centre. So first we must find the
    % image centre
   
    imageCentre = [round(size(loaded_im,2)/2),round(size(loaded_im,1)/2)]';
    
    % make the co-ordinates relative to the image centre
    relativeCentroid = s(k).Centroid' - imageCentre;
    
    % rotate the co-ordinates
    newRelativeCentroid = rotMatrix*relativeCentroid; % note column vector
    
    % find the new image centre
    newCentre = [round(size(rotIm,2)/2),round(size(rotIm,1)/2)]';
    
    % add this back on to make the co-ordinates no longer relative
    newCentroid = newRelativeCentroid + newCentre;
    
    % now we find the connected components
    sTemp = findBugsInImage(rotImMask,min_px_in_bug,max_px_in_bug);
    
    rightBug = 0;
    
    for bug = 1:numel(sTemp)
        if abs(sTemp(bug).Centroid(1) - newCentroid(1)) < 2 && abs(sTemp(bug).Centroid(2) - newCentroid(2)) < 2
            rightBug = bug;
            break
        end
    end
    try
        
        assert(rightBug ~= 0,'Could not find bug in rotated image')
        
        % is optional ability to test whether the long axis is a lot longer
        % than the short axis (more rod shaped) to get rid of very small
        % round cells on end/junk/not stuck down proper cells
        
        %if sTemp(rightBug).MajorAxisLength > 1.5*sTemp(rightBug).MinorAxisLength
        box = sTemp(rightBug).BoundingBox;
        %rectangle('Position', sTemp(rightBug).BoundingBox,'EdgeColor','r', 'LineWidth', 1)

        lengths(index) = sTemp(rightBug).MajorAxisLength;

        % crop out the cell
        % have to rotate first or get bad edges
        images{index} = imcrop(rotIm,box);
        imageMasks{index} = imcrop(rotImMask,box);
        imageNames{index} = filename(1:end-4);
        index = index + 1;
        %end

    catch
        disp('Could not find bug in rotated image')
        
        figure; hold on;
        imshow(rotIm,[min(rotIm(:)),max(rotIm(:))])
        plot(newCentroid(1),size(rotIm,1)-newCentroid(2),'ro')
        title('failed')
    end
    
    % PEAK FINDING CODE
    
    %figure;imshow(images{k},[min(images{k}(:)),max(images{k}(:))]);
    
    % apply a median filter
    
    %images{k} = medfilt2(images{k},[2 2]);
    
    %find the peaks, overlay them and save this image
    
%     BW = imregionalmax(images{k});
%     
%     linIndexes = find(BW==1);
%     
%     % convert indexes to subscripts
%     
%     [x,y] = ind2sub(size(images{k}),linIndexes);
%     
% %     RGB = cat(3,BW*255,zeros(size(BW)),zeros(size(BW)));
% %     imshow(RGB)
% %     
% %     alpha(.25)
%     
%     plot(y,x,'o','MarkerSize',0.5);
%     
%     % remove peaks that are below the average
%     
%     thresh = double(mean(images{k}(:)));
%     
%     chosenPeaksLin = linIndexes(images{k}(linIndexes)>thresh);
%     
%     [x,y] = ind2sub(size(images{k}),chosenPeaksLin);
%     
%     plot(y,x,'ro','MarkerSize',0.5)
%     
%     % save this figure
%     
%     print('-dtiff', '-r400', strcat(filename,'-cell',num2str(k),'.','tif'));
%     
%     % count how many peaks there are
%     
%     numPeaks = length(x);
%     
%     % save this
%     
%     lengths(k) = s(oldNum).MajorAxisLength;
%     peaks(k) = numPeaks;
%     
%     % calculate the normalised peak locations
%     % this should be done by dividing x and y by the size of the rectangle
%     % *100 (for percent)
%     
%     normX = (x/size(images{k},1))*100;
%     normY = (y/size(images{k},2))*100;
%     
%     normXr = round(normX);
%     normYr = round(normY);
%     
%     for i=1:numPeaks
%         points(normXr(i),normYr(i)) = points(normXr(i),normYr(i))+1;
%     end
%     
%     %figure(f);hold on;plot(normYr,100-normXr,'o');  
end

end

function s = findBugsInImage(filtered_im,min_px_in_bug,max_px_in_bug)
    clusters = bwconncomp(filtered_im);

    % remove clusters that are too big or too small to be bugs
    numPixels = cellfun(@numel,clusters.PixelIdxList);

    for i=1:length(clusters.PixelIdxList)
        if numPixels(i) > max_px_in_bug || numPixels(i) < min_px_in_bug
            filtered_im(clusters.PixelIdxList{i})=0;
        end
    end

    % find clusters again now we may have removed some

    clusters = bwconncomp(filtered_im);

    % test doing a bounding box

    s = regionprops(clusters,'BoundingBox','MajorAxisLength','MinorAxisLength','Orientation','Area','PixelList','Extrema','Centroid');

end
