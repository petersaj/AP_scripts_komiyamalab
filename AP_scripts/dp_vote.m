% %% DP vote
% disp('voting...')
% underscore = 1;
% curr_vote = 0;
% while 1;
%     
%     if underscore == 1
%         urlread('http://www.kprifm.com/pages/battle_of_the_fans','post', ...
%             {'poll_selection','53293','Vote','Vote','pollID','14675', ...
%             'pollDiv','14675','columns','','module','poll','data_only','1', ...
%             'user_function[poll_place_vote]','1'});
%             curr_vote = curr_vote+1;
%     end
%     
%     % check current percentages
%     ranks = urlread('http://www.kprifm.com/pages/battle_of_the_fans','get', ...
%         {'poll_selection','53293','Vote','Vote','pollID','14675', ...
%         'pollDiv','14675','columns','','module','poll','data_only','1', ...
%         'user_function[poll_place_vote]','1'});
%     
%     idx_s = strfind(ranks,'</tr></table>')+13;
%     idx_e = strfind(ranks,'%')-1;
%     curr_ranks = zeros(4,1);
%     for j = 1:4
%         curr_ranks(j) = str2num(ranks(idx_s(j):idx_e(j)));
%     end
%     rank_diff = curr_ranks(4) - curr_ranks(1:3);
%     underscore = curr_ranks(4) < 70; 
%     
% end

%% DP vote cont

disp('voting continuously...')
underscore = 1;
curr_vote = 0;
while 1;
    
    urlread('http://www.kprifm.com/pages/battle_of_the_fans','post', ...
        {'poll_selection','53293','Vote','Vote','pollID','14675', ...
        'pollDiv','14675','columns','','module','poll','data_only','1', ...
        'user_function[poll_place_vote]','1'});
    curr_vote = curr_vote+1;
    
end
