% AP: debugging oopsi

framerate = 28.4098; % Hz

foopsi = {};

% set options/parameters
clear V
V.dt = 1/framerate;
V.est_sig = 1;
V.est_lam = 1;
V.est_gam = 1;
V.est_b = 1;
V.est_a = 1;

clear P
% P.a  = 1;
% P.b = 0;
% P.sig = 0.04;
% tau = 1;
% P.gam = 1-(V.dt/tau);
% P.lam = 1;

% set of possible taus
taus = 1; %[0.1 0.3 0.7 1];

% run fast oopsi over params
for curr_cell = 1:length(oopsi_test_cells)
    
    F = oopsi_test_cells{curr_cell};
    F = F - min(F);
    percentiles = quantile(F,[0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75]);
    
    lambdas = [0.1 0.5 1];
    
    for curr_param = 1:length(lambdas)
        tau = 1;
        P.gam = 1-(V.dt/tau);
        P.lam = lambdas(curr_param);
        %P.a = percentiles(end-1);
        [n_out P_out V_out] = fast_oopsi(F,V,P);
        foopsi{curr_cell}{curr_param} = n_out/max(n_out);      
    end
end

% run fast oopsi with good params
for curr_cell = 1:length(oopsi_test_cells)
    F = oopsi_test_cells{curr_cell};
    F = F - min(F);
    tau = 1;
    P.gam = 1-(V.dt/tau);
    P.gam = 0.99
    if curr_cell < 5
        P.lam = 0.5;
    else
        P.lam = 0.5;
    end
    [n_out P_out V_out] = fast_oopsi(F,V,P);
    foopsi{curr_cell}{1} = n_out/max(n_out);
end

%% plot results: fluor on top, deconvolved on bottom
for curr_cell = 1:length(oopsi_test_cells)
    p = [];
    figure
    p(1) = subplot(1+length(lambdas),1,1);
    plot(oopsi_test_cells{curr_cell},'k')
    for i = 1:length(lambdas)
        p(i+1) = subplot(1+length(lambdas),1,i+1);
        plot([foopsi{curr_cell}{i}],'k');  
    end
    linkaxes(p,'x');
end

   
   
   
   
   
   