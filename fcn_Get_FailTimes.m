function t=fcn_Get_FailTimes(bT,b,L); 

ct=0; 
ind=find(abs(bT-b)<=.05); 
if isempty(ind); 
    ind=find(abs(bT-b)<=.1); 
end
t=L(ind); 