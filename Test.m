% Kamran Siddiqui
% Department of Mechanical and Materials Engineering
% University of Western Ontario
%====================================================================
% M-file to detect the object in PIV images and save it as binary image
%=====================================================================
% Dated: January 12, 2008
% Last Modified: March 10, 2017
%=====================================================================

clearvars

fname00='_Camera_'; %input(' Input the file prefix, fname= ');
is=1; %input(' Enter the start file number, is= ');
%ie=input(' Enter the ending file number, ie= ');

foldernum = input(' Enter the folder number to process: ');

pthr=['F:\Paraffin Wax\18Deg - Test 1',num2str(foldernum),'\Input\'];
pthw=['F:\Paraffin Wax\18Deg - Test 1',num2str(foldernum),'\Input - Cropped'];
home=['F:\Paraffin Wax\18Deg - Test 1',num2str(foldernum),'\'];

D = dir([pthr, '\*.tif']);
ie = length(D(not([D.isdir])));
ie = ie-1;

cd(home)

mkdir(home,'Input - Cropped')

h = waitbar(0,'Please wait...');

for I=is:ie
    
   waitbar(I / ie)
    
    % add the No. of the files
    ss=int2str(I);
    if I>=0 & I<=9
        ss=strcat('00',ss);
    end
    if I>=10 & I<=99
        ss=strcat('0',ss);
    end
    if I>=100 & I<=999
        ss=strcat(ss);
    end
        


    %loading the image file
    imfile1=[pthr num2str(foldernum) fname00 ss '.tif'];
    
    im=imread(imfile1);
    
    
    im1 = rot90(rot90(rot90(im)));
    
      
    % output file name
    imfile1=[pthw '\' num2str(foldernum) fname00 ss '.tif'];
    
    % saving the images on the tiff format
    imwrite(im1,imfile1,'tif','compression','none');
    
    
end

delete(h)
    
    