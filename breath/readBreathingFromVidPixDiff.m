[file,directory] = uigetfile()
vidObj = VideoReader([directory file]);

i=0;
while hasFrame(vidObj)
    vidFrame = readFrame(vidObj);
    i = i+1;
    if i == 1
        
       mask = roipoly(vidFrame)

    end

    meanPix(i) = mean(vidFrame(mask));

end

pixDiff = diff(meanPix);

