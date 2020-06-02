# -MRI-Image-Processing-

1. Load MRI data and generate montage view of MRI data- First with extract the dicom files from the folder with the help of dir function. dir gives us the names and numbers of the files. Now examine data of the images and collect the information of slice location and slice thickness of each MRI slice. Then create array to store and read dicom images. With the help of montage function, displayed full resolution view of the entire stack of the dicom images.
2. Replot the montage view with the image slices in an anatomically organized order- The displayed montage is not anatomically organized. Therefore, organize it in anatomical order with the help of sort and for loop. The slice location which we collected earlier plays a key roll in anatomically arranging the montage.
3. Segment the brain tissue. Eliminate the skull and scalp- After that select one slice from stack and view image with imtool. With the help of thresholding convert brain image into binary image. Then using morphological operations, remove the skull from binary image. Then convert binary image into gray scale image which shows only brain tissue without skull.
4. Perform the contrast map of the segmented brain image- Then perform region growing to differentiate tumor and healthy brain tissue and displayed tumor and healthy brain tissue in different colors with HSV function. Using thresh_tool function of matlab we create contrast map of brain image. In which we can select different threshold value in histogram of input image. 286 value will give us better differentiation.
5. Perform the 2D contour view of above image- To get contour of image use contour slice function which gives X,Y,Z coordinate data in single surface.
6. Generate a 3-D image of both the damaged tissue and normal tissue- Now create a 3D display of the brain slices in which different color shows different regions, this is done by patch command. Isosurface allows us to form a 3D stack of slices. movegui function helps to move figure to specified location on screen.