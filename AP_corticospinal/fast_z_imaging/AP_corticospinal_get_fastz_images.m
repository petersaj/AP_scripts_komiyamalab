% Get the dendrite and soma plane images for all cells (run this locally)

% Copy over summed movies and roi files for dendrite/soma
data_path = '/usr/local/lab/People/Andy/Data/corticospinal_fastz_imaging';
animal_paths = {'151212_AP162','151212_AP163'};

cell_idx = 1;

soma_fig = figure('Name','Somas');
dendrite_fig = figure('Name','Dendrites');

for curr_animal = 1:length(animal_paths)
    
    curr_path = [data_path filesep animal_paths{curr_animal}];
    curr_dir = dir(curr_path);
    curr_regions = {curr_dir(3:end).name};
    
    for curr_region = 1:length(curr_regions)
        
        curr_region_path = [curr_path filesep curr_regions{curr_region}];
        region_dir = dir(curr_region_path);
        region_dir_folders = {region_dir([region_dir.isdir]).name};
        region_subregions = region_dir_folders(3:end);
        
        for curr_subregion = 1:length(region_subregions);
            
            curr_subregion_path = [curr_region_path filesep region_subregions{curr_subregion} filesep ...
                'motion_corrected' filesep 'summed'];
            curr_summed_filename = dir([curr_subregion_path filesep '*.tif']);
            im = AP_load_tiff([curr_subregion_path filesep curr_summed_filename.name]);
            
            curr_dendrite_im = nanmean(im(:,:,1:2:end),3);
            curr_soma_im = nanmean(im(:,:,2:2:end),3);
            
            roi_dir = dir([curr_subregion_path filesep '*.roi']);
            roi_filenames = sort({roi_dir.name});
            
            curr_dendrite_rois = load([curr_subregion_path filesep roi_filenames{1}],'-MAT');
            curr_soma_rois = load([curr_subregion_path filesep roi_filenames{2}],'-MAT');
            
            if length(curr_soma_rois.polygon.ROI) > 1;
                roilabel_filenames = dir([curr_subregion_path filesep '*.roilabel']);
                curr_roilabels = load([curr_subregion_path filesep roilabel_filenames.name],'-MAT');
                
            end
            
            for curr_cell = 1:length(curr_soma_rois.polygon.ROI)
                
                curr_soma_roi = curr_soma_rois.polygon.ROI{curr_cell};
                
                if length(curr_soma_rois.polygon.ROI) > 1;
                    curr_dendrite_roi_idx = ...
                        cellfun(@(x) strcmp(num2str(curr_cell),x),curr_roilabels.roi_labels);
                else
                    curr_dendrite_roi_idx = true(size(curr_dendrite_rois.polygon.ROI));
                end
                
                curr_dendrite_rois_cell = curr_dendrite_rois.polygon.ROI(curr_dendrite_roi_idx);
                
                % Plot plane image and ROI
                figure(dendrite_fig);
                subplot(7,7,cell_idx); hold on;
                imagesc(curr_dendrite_im);colormap(gray);
                for curr_dendrite = 1:length(curr_dendrite_rois_cell);
                    plot(curr_dendrite_rois_cell{curr_dendrite}(:,1), ...
                        curr_dendrite_rois_cell{curr_dendrite}(:,2),'b');
                end
                axis square; axis off;
                px_leeway = 20;
                curr_dendrite_cat = vertcat(curr_dendrite_rois_cell{:});
                axis_dim = ceil(max(max(curr_dendrite_cat) - ...
                    min(curr_dendrite_cat)));
                min_x = floor(min(curr_dendrite_cat(:,1)));
                min_y = floor(min(curr_dendrite_cat(:,2)));
                xlim([min_x-px_leeway,min_x+axis_dim+px_leeway]);
                ylim([min_y-px_leeway,min_y+axis_dim+px_leeway]);
                
                figure(soma_fig);
                subplot(7,7,cell_idx); hold on;
                imagesc(curr_soma_im);colormap(gray);
                plot(curr_soma_roi(:,1), ...
                    curr_soma_roi(:,2),'b');
                axis square; axis off;
                px_leeway = 20;
                axis_dim = ceil(max(max(curr_soma_roi) - ...
                    min(curr_soma_roi)));
                min_x = floor(min(curr_soma_roi(:,1)));
                min_y = floor(min(curr_soma_roi(:,2)));
                xlim([min_x-px_leeway,min_x+axis_dim+px_leeway]);
                ylim([min_y-px_leeway,min_y+axis_dim+px_leeway]);
                
                cell_idx = cell_idx + 1;
            end           
        end
    end
end 




