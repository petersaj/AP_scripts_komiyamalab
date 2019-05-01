%% concatenate selected tiffs to make one large file

[tiff_filename,tiff_path]=uigetfile('*.tif','Choose TIFF files','Multiselect','on');

savename = [tiff_path tiff_filename{1}(1:11) '_Concat.tif'];
w = waitbar(0,'Combining images');

for i = 1:length(tiff_filename);
    img_filename = [tiff_path tiff_filename{i}];
    
    imageinfo=imfinfo(img_filename,'tiff');
    numframes=length(imageinfo);
    M=imageinfo(1).Width;
    N=imageinfo(1).Height;
    
    %     im_tensor = TIFFStack(img_filename);
    
    disp('Loading file (native format)....')
    clear im
    for loadframe = 1:numframes
        im(:,:,loadframe) = imread(img_filename,'tiff',loadframe);
    end
    
    for u = 1:numframes
            if u == 1 && i == 1;
                for windows_pause = 1:10
                    try
                        imwrite(im(:,:,u),savename,'tif','Compression','none','WriteMode','overwrite', ...
                    'Description', imageinfo(1).ImageDescription);
                        waitbar(u/numframes,w,['Combining images: ' num2str(i) '/' num2str(length(tiff_filename))])
                        break;
                    catch me
                        pause(0.2);
                    end
                    disp('Windows Locked! Didn''t Write')
                    
                end
            else
                for windows_pause = 1:10
                    try
                        imwrite(im(:,:,u),savename,'tif','Compression','none','WriteMode','append', ...
                    'Description', imageinfo(1).ImageDescription);
                        waitbar(u/numframes,w,['Combining images: ' num2str(i) '/' num2str(length(tiff_filename))])
                        break;
                    catch me
                        pause(0.2);
                    end
                    disp('Windows Locked! Didn''t Write')
                end
            end
    end
end
close(w);
