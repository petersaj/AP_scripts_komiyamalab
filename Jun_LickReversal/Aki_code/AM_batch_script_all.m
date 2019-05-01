%%
m = mousedata();

for n_mouse = 1:length(m)
    try
        AM_batch_script_single(m(n_mouse));
    catch err 
        disp(err);
        save([mousename '_err.mat'], 'err');
    end
end












