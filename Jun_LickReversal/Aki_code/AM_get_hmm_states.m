function states = AM_get_hmm_states(result)
    p_transition = 0.3;
    p_variability = 0.1;

    initial_transition_probabilities = ones(4)*p_transition/3 + diag([1 1 1 1])*(1-p_transition);
    clear p_trans

    inital_emission_probabilities = [1-p_variability p_variability 1-p_variability p_variability; 1-p_variability 1-p_variability p_variability p_variability; p_variability p_variability 1-p_variability 1-p_variability; p_variability 1-p_variability p_variability 1-p_variability];
    clear p 
    
    [estimated_transition_probabilities, estimated_emission_probabilities] = ...
        hmmtrain(result,initial_transition_probabilities,inital_emission_probabilities);
    
    states = hmmviterbi(result,estimated_transition_probabilities,estimated_emission_probabilities);
end