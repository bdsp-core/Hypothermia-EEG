function varargout = pdistn(x, s, mu, background_prob)
% pdist is a helper function that calculates the confidence limits of the EM
% estimate of the learning state process.  For each trial, the function
% constructs the probability density for a correct response.  It then builds
% the cumulative density function from this and computes the p values of
% the confidence limits
%
% variables:
%   xx(ov)   EM estimate of learning state process
%   ss(ov)   EM estimate of learning state process variance
%   pmatrix  vector of the level of certainty the ideal observer has that performance is better than chance at each trial
%   dels     bin size of the probability density p values
%   pr       bins of the probability density distribution
%   fp       p{k|j}, probability density of the probability of a correct response at trial k     (equation B.3)*
%   pdf      probability density function
%   sumpdf   cumulative density function of the pdf
%   lowlimit index of the p value that gives the lower 95% confidence
%            bound
%   highlimit   index of the p value that gives the upper 95% confidence
%               bound
%   middlimit   index of the p value that gives the
%   p05      the p value that gives the lower 95% confidence bound
%   p95      the p value that gives the upper 95% confidence bound
%   p50     the p value that gives the 50% confidence bound
%   pmode    the p value that gives the highest probability density

dels = 1e-4; pr = dels:dels:1-dels; pmatrix = nan(numel(x),numel(pr));
p050 = nan(size(x)); p950 = p050; p500 = p050; pmode = p050;
p025 = p050; p975 = p050;

for ov = 1:size(x,2),
    xx = x(ov); ss = s(ov);
    term1 = 1./(sqrt(2*pi*ss) * (pr.*(1-pr)));
    term2 = exp(-1/(2*ss) * (log (pr./((1-pr)*exp(mu))) - xx).^2);
    pdf = term1 .* term2;
    pdf = dels * pdf;
    
    % Integrate the pdf
    sumpdf = cumtrapz(pdf);  % sumpdf = cumsum(pdf);
    
    q025  = find(sumpdf>0.025,1);
    if isempty(q025), q025 = 1; end
    q050  = find(sumpdf>0.05,1);
    if isempty(q050), q050 = 1; end
    q500 = find(sumpdf>0.5,1);
    if isempty(q500), q500 = length(pr); end
    q950 = find(sumpdf<0.95,1,'last');
    if isempty(q950), q950 = length(pr); end
    q975 = find(sumpdf<0.975,1,'last');
    if isempty(q975), q975 = length(pr); end
    
    p025(ov)   = pr(q025);
    p050(ov)   = pr(q050);
    p950(ov)   = pr(q950);
    p500(ov)   = pr(q500);
    p975(ov)   = pr(q975);
    [maxwhy,maxeye]     = max(pdf);
    pmode(ov) = pr(maxeye);
    
    pmatrix(ov,:) = sumpdf;
end
inte = fix(background_prob/dels); pmatrix = pmatrix(:, inte);

if nargout==1,
    varargout = {struct('p025',p025,'p050',p050,'p500',p500,'p950',p950,'p975',p975,'pmode',pmode,'pmatrix',pmatrix)};
else
    varargout = {p050, p950, p500, pmode, pmatrix};
end

% [EOF] pdistn.m