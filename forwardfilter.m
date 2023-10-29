function  [p, xhat, sigsq, xhatold, sigsqold] = forwardfilter(I, sigE, xguess, sigsqguess, mu)
% forwardfilter is a helper function that implements the forward recursive
% filtering algorithm to estimate the learning state (hidden process) at
% trial k as the Gaussian random variable with mean x{k|k} (xhat) and
% SIG^2{k|k} (sigsq).
%
% variables:
%   xhatold      x{k|k-1}, one-step prediction (equation A.6)*
%   sigsqold     SIG^2{k|k-1}, one-step prediction variance (equation A.7)*
%   xhat         x{k|k}, posterior mode (equation A.8)*
%   sigsq        SIG^2{k|k}, posterior variance (equation A.9)*
%   p            p{k|k}, observation model probability (equation 2.2)*
%   N            vector of number correct at each trial
%   Nmax         total number that could be correct at each trial
%   K            total number of trials
%   number_fail  saves the time steps if Newton's Method fails

K = size(I,2);
N = I(1,:);
Nmax = I(2,:);

xhatold = nan(1,K); sigsqold = nan(1,K); sigsq = nan(1,K); xhat = nan(1,K); flagfail = nan(1,K);

% Initial conditions: use values from previous iteration
xhat(1)   = xguess;
sigsq(1)  = sigsqguess;
% number_fail = [];

% for each trial, compute estimates of the one-step prediction, the
% posterior mode, and the posterior variance

for loopk = 1:K,
    k = loopk+1;
    
    % Compute the one-step prediction estimate of mean and variance
    xhatold(k)  = xhat(loopk);
    sigsqold(k) = sigsq(loopk) + sigE^2;
    
    xhat(k) = xhatold(k) + sigsqold(k)*(N(loopk) - Nmax(loopk)*(exp(xhatold(k)+mu)/(1+exp(xhatold(k)+mu))));
    denom       = -1/sigsqold(k) - Nmax(loopk)*exp(mu)*exp(xhat(k))/(1+exp(mu)*exp(xhat(k)))^2;
    sigsq(k)    = -1/denom;
end

number_fail = find(flagfail==1);

if ~isempty(number_fail)
    fprintf(2,'Newton convergence failed at times % d \n', number_fail)
%     pause
end

% Compute the observation model probability estimate
% xhat
% keyboard
p = exp(mu)*exp(xhat)./(1+exp(mu)*exp(xhat));

