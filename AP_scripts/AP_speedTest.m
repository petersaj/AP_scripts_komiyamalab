clear all

temp = magic(2000);

matrixmath.for = zeros(100,1);
for i = 1:100
    tic;
    temp*temp; 
    matrixmath.for(i) = toc;
end

matrixmath.single = 0;
tic;
temp*temp;
matrixmath.single = toc;

clear a
memalloc.for = zeros(100,1);
for i = 1:100;
    tic;
    a = zeros(2000);
    memalloc.for(i) = toc;
end

memalloc.single = 0;
tic;
a = zeros(2000);
memalloc.single = toc;

bench_temp = bench(5);
close all;
matlab_bench = median(bench_temp(:,1:4));

clearvars -except matrixmath memalloc matlab_bench 
