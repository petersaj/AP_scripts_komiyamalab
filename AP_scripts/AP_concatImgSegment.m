%% make small segment of image for whole duration

img = {1 2 3 4 6 7 8 9 10 11};% 12 13 14 15 16 17 18 19 20 21 22};
img_str = cellfun(@num2str, img, 'UniformOutput', false);
img_path = 'C:\Users\Andy\Documents\KomiyamaLab\Data\TempFilesCorrected\110728\';
savename = sprintf('%sAP_110728_18_Turboreg_StitchWhole.tif',img_path);
w = waitbar(0,'Combining images');

for i = 1:length(img);
    % get img filename
    curr_img_str = sprintf('00%s',img_str{i});
    curr_img_str = curr_img_str(end-2:end);
    img_filename = sprintf('%sAP_110728_18_%s_Turboreg.tif',img_path,curr_img_str);
    
    x_borders = [200 320];
    y_borders = [60 100];
    size_sm = [diff(y_borders)+1 diff(x_borders)+1];
    
    imageinfo=imfinfo(img_filename,'tiff');
    numframes=length(imageinfo);
    M=imageinfo(1).Width;
    N=imageinfo(1).Height;
    
    im_lg=zeros(N,M);
    im_sm = im_lg;
    %im_sm = zeros(size_sm(1),size_sm(2));
    
    
    for u = 1:numframes
            im_lg(:,:)=imread(img_filename,'tiff',u);
            im_sm(:,:) = im_lg;%(y_borders(1):y_borders(2), ...
                %x_borders(1):x_borders(2));
            im_sm=uint16(round(im_sm));
            if u == 1 && i == 1;
                for windows_pause = 1:10
                    try
                        imwrite(im_sm,savename,'tif','Compression','none','WriteMode','overwrite');%, 'Description', ImageDescriptio]);
                        waitbar(u/numframes,w,['Combining images: ' num2str(i) '/' num2str(length(img))])
                        break;
                    catch me
                        pause(0.2);
                    end
                end
            else
                for windows_pause = 1:10
                    try
                        imwrite(im_sm,savename,'tif','Compression','none','WriteMode','append');
                        waitbar(u/numframes,w,['Combining images: ' num2str(i) '/' num2str(length(img))])
                        break;
                    catch me
                        pause(0.2);
                    end
                end
            end
    end
end
close(w);




