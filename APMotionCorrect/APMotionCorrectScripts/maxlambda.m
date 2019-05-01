
%before running this file, load up the _PI.mat file that was created by
%split_createPI

%powers of 10 to scan over for lambda
%sampling uniform in log space determines lambda to a uniform percentage
%precision, i have intentionally scanned over a larger range of lambda than
%is neccesary in order to illustrate what happens at small and large values of lambda.
h = waitbar(0.0,'Scanning Lambda Values...');
% set(h,'Position',[50 0 360 72]);
set(h,'Name','Expectation Maximization');

n=-3:.1:.5;

%values of lambda to sample
lambdas=10.^n;

%initialize index for lambda values
jj=0;

%clear previous results if any
clear saveprobs

%loop over values of lambda
for lambda=lambdas
  %increment index
  jj=jj+1;
  dl=1/length(lambdas);
  waitbar((jj-1)/length(lambdas),h,['Lambda=' num2str(lambda)]);
  %run an abbreviated version of the HMM on the first 20 frames
  %totprob will be the overall probability for that run
  %offsets predicted for that value are plotted at the completion of each
  %iteration.
  short_markov;
  
  %save the total probablity for that value of lambda
  saveprobs(jj)=totprob;
end

%plot the curve of overall probabilities


%pick out the value of lambda which is the overall most probable
[maximumvalue,maximumindex]=max(saveprobs);
lambda=lambdas(maximumindex);

figure(99);
clf;
set(gcf,'Position',[800 200 1280-800 500]);
hold on;
plot(lambdas,saveprobs);
plot(lambda,maximumvalue,'rx');
xlabel('\lambda');
ylabel('total probability');


waitbar(0.0,h,['Run with Lambda=' num2str(lambda)]);
set(h,'Name',['HMM Lambda=' num2str(lambda)]);

%run the HMM algorithm, now that lambda is set to its most probable value,
%over all frames of the movie. this will result in a _L file containing the
%offsets
markov_on_PIsave;
if exist('h') delete(h); end
%save the probabilities into the _PI file for posterity
save([savedir filename(find(filename == '\',1,'last'):end-4) '_PI.mat'],'saveprobs','lambdas','-append');