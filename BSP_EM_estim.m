function [Responses, ptile,xnew,signewsq,A,p, x, s] = BSP_EM_estim(bsr_binary, Fs, number_of_samples)


% % Input
% bsr_binary: binary series
% Fs: sampling frequency
% number_of_samples: binomial size, =1 for bernoulli
% % Output
% Responses = analyzed series (binary if number_of_samples =1)
% ptile: struct containing bsp and confidence intervals
% xnew,signewsq,A : needed for trialtotrial_comparisons

%----------------------Initialization
Responses = doSum(bsr_binary,Fs,number_of_samples/Fs); MaxResponse=number_of_samples; BackgroundProb=0.5; SigE = 1e-4;

% make it a row vector
Responses = Responses(:)';
% Make the second row be the max number of correct responses
I = [Responses; MaxResponse*ones(1,length(Responses))];

SigsqGuess  = SigE^2; mu = log(BackgroundProb/(1-BackgroundProb));
% Criterion for convergence
CvgceCrit = 1e-8;

xguess = 0; NumberSteps = 20000; NumberMove = 100;
newsigsq = nan(1,NumberSteps); xnew1save = newsigsq;
% Parameters for the gamma regularization
alpha = 50000;
beta = 2;

%----------------------EM
for loopstep = 1:NumberSteps,
    for loopmove = 1:NumberMove,
        loopstepmove = (loopstep-1)*100+loopmove;
        [p, x, s, xold, sold] = forwardfilter(I, SigE, xguess, SigsqGuess, mu);
         
        [xnew, signewsq, A]   = backwardfilter(x, xold, s, sold);
        xnew(1) = xnew(2); signewsq(1) = signewsq(2);
        [newsigsq(loopstepmove)] = em_bino(I, xnew, signewsq, A, 2,alpha,beta);
        xnew1save(loopstepmove) = xnew(1);
        SigE = sqrt(newsigsq(loopstepmove)); xguess = xnew(1); SigsqGuess = signewsq(1);
    end

% Check for convergence of EM
    a1 = abs(newsigsq(loopstepmove) - newsigsq(loopstepmove-1));
    a2 = abs(xnew1save(loopstepmove) -xnew1save(loopstepmove-1));
    if (a1 < CvgceCrit && a2 < CvgceCrit)
        fprintf(2, 'EM estimates of learning state process variance and start point converged after %d steps   \n',  loopstepmove)
        break
    end
end
newsigsq = newsigsq(~isnan(newsigsq)); xnew1save = xnew1save(~isnan(xnew1save));

if (loopstep == NumberSteps)
    fprintf(2,'failed to converge after %d steps; convergence criterion was %f \n', loopstepmove, CvgceCrit)
end
% keyboard
% -------------------------------------------
%integrate and do change of variables to get confidence limits
%%
    ptile = pdistn(xnew, signewsq, mu, BackgroundProb);

% -------------------------------------------------------------------------------------
