function imgs = fileIO_lvtiff2mat(fileImgIn, imgNum, repNum, cropSize, bgValue)
% convert image and zernike files from microscope to matlab variables
% % Output
% imgs: diversity images for phase reconstruction;
% % Input
% fileImgIn: strings of the input path/image name;
% imgNum: number of total phases (including initial aberrated phase);
% repNum: number of images for each phase;
% cropSize: image size for calculation;
% bgValue: background of the input image;
 
% by Min Guo
% Jul. 30, 2020;

%% image preprocessing

flagSmoothEdge = 1; % smooth the edges of the image;
imgsRaw = single(ReadTifStack(fileImgIn));
[Sx1, Sy1, rawNum] = size(imgsRaw);

if(rawNum<imgNum*repNum)
    error('image number in TIFF is less than image number');
end

%% average if acquistion repeated for each phase
if(repNum>1)
    imgsAve = zeros(Sx1,Sy1,imgNum, 'single');
    for i = 1:imgNum
        imgAve = zeros(Sx1,Sy1, 'single');
        for j = 1:repNum
            iSlice = (i-1)*repNum + j;
            if(iSlice <=rawNum)
                imgAve = imgAve + imgsRaw(:,:,iSlice);
            end
        end
        imgsAve(:,:,i) = imgAve;
    end
    imgsAve = imgsAve/repNum;
else
    imgsAve = imgsRaw;
end

%% rotation and crop
imgs = zeros(cropSize,cropSize,imgNum, 'single');
for i = 1:imgNum
    img = imgsAve(:,:,i);
    Soxy = round((size(img)-cropSize)/2+1);
    Sox = Soxy(1);
    Soy = Soxy(2);
    imgs(:,:,i) = img(Sox:Sox+cropSize-1,Soy:Soy+cropSize-1);
end
imgs = imgs - bgValue;
imgs = max(imgs,0.01);

%% smooth the boundary of images
if(flagSmoothEdge==1)   
    sKernel = fspecial('gaussian',[10 10],3); 
    for i = 1:imgNum
        img = imgs(:,:,i);
        img = edgetaper(img,sKernel);
        imgs(:,:,i) = img;
    end
end