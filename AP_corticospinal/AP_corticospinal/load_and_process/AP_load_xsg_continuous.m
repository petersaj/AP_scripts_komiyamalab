% Written by AK, something or another changed by AP

function data = AP_load_xsg_continuous(dirname)
    if(nargin==0)
        lastwarn('Directory was not given.');
        dirname = uigetdir;
    end
    
    files = dir(fullfile(dirname, '*.xsglog'));
    
    name_first='';
    epoch_first='';
    
    for i_file = 1:numel(files)
        fn = files(i_file).name;
        [name, epoch] = parse_filename(fn);
        if(i_file==1)
            name_first = name;
            epoch_first = epoch;
        else
            if(~strcmp(name_first,name)||~strcmp(epoch_first,epoch))
                error('Files does not match. Only one epoch is allowed. Delete/move unused files.');
            end
        end
    end
    
    data.name = name_first;
    data.epoch = epoch_first;
    data.file_infos = files;
    
    fn = [name_first '.txt'];
    fid = fopen(fullfile(dirname,fn),'r');
    if(fid<0)
        warning(['Could not open file ''' fullfile(dirname,fn) '''']);
        data.txt = '';
    else
    	data.txt = fread(fid,inf,'*char');
        fclose(fid);
    end
    
    data.channels = zeros(0,numel(files));
    data.channel_names = cell(numel(files),1);
    for i_file = 1:numel(files)
        fn = files(i_file).name;
        [~, ~, channel_name] = parse_filename(fn);
        
        fid = fopen(fullfile(dirname,fn),'r');
        if(fid<0)
            error(['Could not open file ''' fullfile(dirname,fn) '''']);
        end
        fdata =  fread(fid,inf,'double');
        fclose(fid);
%         m = memmapfile(fullfile(dirname,fn),'Format','double');
%         fdata = m.Data;
        
        if(numel(fdata)~=size(data.channels,1))
            data.channels(end+1:numel(fdata),:)=NaN; % changed from data to data.channels
        end
        data.channels(1:numel(fdata),i_file)=fdata;
            
        data.channel_names{i_file} = channel_name; % fixed the indexing parenthesis
    end
    
    if(any(any(isnan(data.channels))))
%     if(any(diff(cellfun(@numel,data.channels.values))))
        lastwarn('Mismatched data size.'); 
    end
    
end

function [name, epoch, channel] = parse_filename(fn) 
    matches = regexp(fn,'([A-z]{2}\d{4})([A-z]{4}\d{4})_(.*)\.xsglog','tokens');
    if(isempty(matches))
        error(['Wrong file name. Correct or delete ''' fn '''']);
    end
    name = matches{1}{1};
    epoch = matches{1}{2};
    channel = matches{1}{3};
end
