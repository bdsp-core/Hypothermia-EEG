function sum_sig = doSum(signal,Fs,sec)

samples_per_batch = floor(Fs*sec);
n = length(signal);
% length of the vector must be a multiple of samples_per_batch
samples_to_consider = n - rem(n,samples_per_batch);
sum_sig = sum(reshape(signal(1:samples_to_consider),samples_per_batch,samples_to_consider/samples_per_batch),1) ;

end

