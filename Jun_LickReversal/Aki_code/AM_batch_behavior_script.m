%%
mousename = 'JL041';

    load([mousename '_dspt.mat']);
    sessions = AM_batch_load_behavior(sessions);
    save([mousename '_behavior.mat'], 'sessions');




%%
mousename = 'JL042';
    load([mousename '_dspt.mat']);
    sessions = AM_batch_load_behavior(sessions);
    save([mousename '_behavior.mat'], 'sessions');




%%
mousename = 'JL043';
    load([mousename '_dspt.mat']);
    sessions = AM_batch_load_behavior(sessions);
    save([mousename '_behavior.mat'], 'sessions');




%%
mousename = 'JL048';
    load([mousename '_dspt.mat']);
    sessions = AM_batch_load_behavior(sessions);
    save([mousename '_behavior.mat'], 'sessions');



%%
mousename = 'JL051';
    load([mousename '_dspt.mat']);
    sessions = AM_batch_load_behavior(sessions);
    save([mousename '_behavior.mat'], 'sessions');

%%
mousename = 'JL053';
    load([mousename '_dspt.mat']);
    sessions = AM_batch_load_behavior(sessions);
    save([mousename '_behavior.mat'], 'sessions');


%%
mousename = 'JL055';
    load([mousename '_dspt.mat']);
    sessions = AM_batch_load_behavior(sessions);
    save([mousename '_behavior.mat'], 'sessions');


%%
mousename = 'JL056';     % bad behavior
% try 
% 
% %     AM_batch_get_dFoF_from_img(mousename);
% %     sessions = AM_batch_behaviordispatcher(mousename);
% %     sessions = AM_batch_load_behavior(sessions);
% %     sessions = AM_batch_load_analysis(sessions);
% %     save([mousename '_data.mat'], 'sessions');
% %     save([mousename '_aligned.mat'], 'data');
% catch err 
%     save([mousename '_err.mat'], 'err');
% end










