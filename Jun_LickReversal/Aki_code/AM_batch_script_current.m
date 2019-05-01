% mousename = 'JL062';
% try 
%     AM_batch_get_dFoF_from_img(mousename);
%     sessions = AM_batch_behaviordispatcher(mousename);
%     sessions = AM_batch_load_behavior(sessions);
%     sessions = AM_batch_load_analysis(sessions);
%     save([mousename '_data.mat'], 'sessions');
% catch err 
%     save([mousename '_err.mat'], 'err');
% end
% 
% 
