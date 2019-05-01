function SI4_frames = associate_frame_tags_to_triggers(SI4_frame_triggers, SI4_frame_tags)

    delay_limit = 30*1E-3;
    
    SI4_frames = SI4_frame_triggers;
    assert(numel(SI4_frames)>0);
    
    SI4_frames(1).SI4_repeats_done = [];
    SI4_frames(1).SI4_frame_tag = [];
    
    first_repeats_done = [];
    first_offset = [];
    
    for i_tag = 1:length(SI4_frame_tags)
        index_of_associated_trigger = find(SI4_frame_tags(i_tag).xsg_sec - delay_limit < [SI4_frame_triggers.xsg_sec_end] ...
            & [SI4_frame_triggers.xsg_sec_end] < SI4_frame_tags(i_tag).xsg_sec);
        
        if(numel(index_of_associated_trigger)~=1)
            lastwarn(['could not find an associated trigger to the ' num2str(i_tag) 'th frame tag.']);
        else
            SI4_frames(index_of_associated_trigger).SI4_repeats_done = SI4_frame_tags(i_tag).SI4_repeats_done;
            SI4_frames(index_of_associated_trigger).SI4_frame_tag = SI4_frame_tags(i_tag).SI4_frame_tag;
            if(isempty(first_repeats_done))
                first_repeats_done = SI4_frame_tags(i_tag).SI4_repeats_done;
                first_offset = SI4_frame_tags(i_tag).SI4_frame_tag - index_of_associated_trigger;
            end
        end
    end
    
    current_repeats_done = first_repeats_done;
    current_offset = first_offset;
    
    for i_frame = 1:length(SI4_frames)
        if(isempty(SI4_frames(i_frame).SI4_repeats_done))
            SI4_frames(i_frame).SI4_repeats_done = current_repeats_done;
            SI4_frames(i_frame).SI4_frame_tag = i_frame+current_offset;
        else
            if(current_repeats_done == SI4_frames(i_frame).SI4_repeats_done)
                if(current_offset ~= SI4_frames(i_frame).SI4_frame_tag - index_of_associated_trigger)
                    lastwarn('some frames are not properly aligned. check the data.');
                    current_offset = SI4_frames(i_frame).SI4_frame_tag - index_of_associated_trigger;
                end                
            else 
                current_repeats_done = SI4_frames(i_frame).SI4_repeats_done;
                current_offset = SI4_frames(i_frame).SI4_frame_tag - index_of_associated_trigger;
            end
        end
    end
end
