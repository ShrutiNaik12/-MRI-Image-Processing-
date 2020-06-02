%%  Image Processing with MATLAB: applications in Medicine and biology project
%
%%
clear;
close all;
clc;
imtool close all;
%% extracting the files from the folder
filefolder = fullfile(pwd,'11'); % files in this variable
files = dir(fullfile(filefolder, '*.dcm')); % extracting every dicom file
filenames = {files.name};
A = filenames.';
%% Examine data from the images of the file
for i = 1:48
info = dicominfo (fullfile(filefolder,filenames{i}));
voxel_size = [info.PixelSpacing; info.SliceThickness];
% read one file to get size
I = dicomread(fullfile(filefolder,filenames{i}));
class_I = class(I);
size_I = size(I);
numImages = length(filenames);
 slice_loc(i)=info.SliceLocation;   %slice location of each MRI slice
 slice_thk(i)=info.SliceThickness;  %slice thickness of each MRI slice    
end
%% Read slice Images - 3D matrix
%creating array to read and store the DICOM images
MRI = zeros(size_I(1), size_I(2), numImages, class_I);  %store images in new array
for i = length(filenames):-1:1
fname = fullfile(filefolder,filenames{i});
MRI(:,:,i) = uint16(dicomread(fname));
end
%% Display in Montage
figure
montage(reshape(uint16(MRI), [size(MRI,1),size(MRI,2),1,size(MRI,3)]));title('Brain MRI slices in Montage View');
set(gca,'clim', [0 300]); % for better resolution
drawnow;
shg;
%% Anatomically organized montage view
sort_slice_loc=sort(slice_loc,'descend');
for i=1:length(sort_slice_loc)
    for j=1:length(slice_loc)
        if sort_slice_loc(i)== slice_loc(j)
            d=(string(fullfile(filefolder,(A(j)))));
            sorted_mri(:,:,i)=int16(dicomread(d));
        end
    end
end
figure;
montage(reshape((sorted_mri), [size(MRI,1),size(MRI,2),1,size(MRI,3)]));title('Anatomically Arranged MRI slices in Montage View');
set(gca,'clim', [0 300]); % for better resolution  
drawnow;
shg;
%% exploring image data
im = sorted_mri(:,:,15); % pick convinent slice with tumor
max_level = double(max(im(:)));
imt = imtool(im, [0, max_level]);
%% Brain tissue thresholding by creating a binary mask
threshhold_value = 70;
binary_sorted_mri = im > threshhold_value;
figure
imshow(binary_sorted_mri), title('binary image')
binaryImage = bwareafilt(binary_sorted_mri,1);
binaryImage = imopen(binaryImage, true(5));
binaryImage = bwareafilt(binaryImage,1);
binaryImage = imfill(binaryImage, 'holes');
binaryImage = imdilate(binaryImage,true(5));
figure
imshow(binaryImage, []); title('binary image after removing skull')
%% segment brain tissue
brain_noskull = im;
brain_noskull(~binaryImage) = 0; % Remove skull from image
figure
imshow(brain_noskull, []), title('Image after removing skull')
%% region growing
img = brain_noskull; %sorted_mri(:,:,15);
seedmask = zeros(256,256);
seedmask(167,176)=64;
seedintensity=img(167,176);
seedrangemin=seedintensity-40;
seedrangemax=seedintensity+40; 
oldseed=1;
newseed=0;
while newseed~=oldseed
oldseed=newseed;
newseed=0;
for i=2:255
for j=2:255
if seedmask(i,j)>0
intens = img((i-1),j);
if (intens>=seedrangemin) && (intens<=seedrangemax)
newseed = newseed+1;
seedmask((i-1),j) = 64;
end
intens = img((i+1),j);
if (intens>=seedrangemin) && (intens<=seedrangemax)
newseed = newseed+1;
seedmask((i+1),j) = 64;
end
intens = img(i,(j-1));
if (intens>=seedrangemin) && (intens<=seedrangemax)
newseed = newseed+1;
seedmask(i,(j-1)) = 64;
end
intens = img(i,(j+1));
if (intens>=seedrangemin) && (intens<=seedrangemax)
newseed = newseed+1;
seedmask(i,(j+1)) = 64;
end
end
end
end
end
figure
imagesc(seedmask), colormap(gray); title('seedmask')
A1=double(brain_noskull)+seedmask; % to display image with tumor prominent
figure;imagesc(A1),colormap(gray),title('Tumor and Brain can be Differentiated');
%% opening of the image
% erosion followed by dilation
A=A1; %bw(:,:,10);
modimg = zeros(256,256);
for i=2:255
for j=2:255
if A(i,j)>0
modimg(i,j)=1;
if (A(i-1,j)==0) || (A(i+1,j)==0) || (A(i,j-1)==0) || (A(i,j+1)==0)
modimg(i,j)=0;
end
end
end
end
modimg2 = zeros(256,256);
for i=2:255
for j=2:255
if modimg(i,j)>0
modimg2(i,j)=1;
end
if (modimg(i-1,j)>0)||(modimg(i+1,j)>0)||(modimg(i,j-1)>0)||(modimg(i,j+1)>0)
modimg2(i,j)=1;
end
end
end
figure
image(256*modimg2),colormap(gray),title('After opening of the image')
%% Display tumor and healthy tissue with different colors
% identify blobs and display them
L = bwlabeln(modimg2);
stats = regionprops(L,'Area','Centroid'); % gives detail of image
LL=L +1; %AA=L(:,:,15)+logical(seedmask);
cmap = hsv(length(stats));cmap=[0 0 0;cmap];
LL=cmap(LL,:);
LL=reshape(LL,[256,256,3]);
AA=LL+logical(seedmask); %stats = regionprops(AA,'Area','Centroid');
figure
imshow(AA),title('Tumor in the selected slice') 
%different colors in above figure shows different whole regions
%% contrast map of brain image
level=thresh_tool((A1+double(seedmask)),'gray');
mriBrainPartition = uint8(zeros(size(sorted_mri)));
mriBrainPartition(sorted_mri>0 & A1+ double(sorted_mri)<level) = 2;
mriBrainPartition(sorted_mri>=level)=3;
figure
imshow(uint16(mriBrainPartition(:,:,15)),[0 0 0;0 0 0;0.25,0.25,0.25]);
% selecting 286 in thresh tool gives us better differentiation
% so set 286 in intensity distribution and click done
%% contour of the image
mriAdjust(:,:,15)=sorted_mri(:,:,15)+int16(seedmask);
image_num=15;
cm = brighten(jet(60),-0.5);
figure('Colormap',cm),title('contour of the image')
contourslice(mriAdjust,[],[],image_num)
axis ij tight
daspect([1,1,1])
%% create a 3D display
Ds = imresize(mriBrainPartition,0.25,'nearest');
Ds=flip(Ds,1);
Ds=flip(Ds,2);
Ds=permute(Ds,[3 2 1]);
voxel_size2 = voxel_size([1 3 2]).*[4 1 4];
white_vol = isosurface(Ds,2.5);
gray_vol = isosurface(Ds,1.5);
%% the output of isosurface as 3D display
h=figure('visible','off','outerposition',[0 0 800 600],'renderer','openGL','windowstyle','normal');
patch(white_vol,'Facecolor','b','Edgecolor','none');
patch(gray_vol,'Facecolor','y','Edgecolor','none','FaceAlpha',0.5);
view(45,15);
daspect(1./voxel_size2(1:3));
daspect([1 1 1]);
axis tight;
axis off;
camlight;
camlight(-80,-10);
lighting phong;
movegui(h,'center'),title('3D image of the slice');
set(h,'visible','on');