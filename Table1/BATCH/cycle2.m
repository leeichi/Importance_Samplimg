clear
n0 = [1000];
p = [0.9];
k = [10];
nb = [300];
delta=1;
va=100;
trial=500;

resultVAR=[];
result = [];
VAR=[]; %batch var while it ends.
count=1;
for z = 1 : length(k)
    for i = 1:length(n0)
        for j = 1:length(nb)
                [ANS,PCS,CPU_TIME VAR] = fsp_SampleVar(n0(i),p,k(z),nb(j),delta,va,trial);
               
                result(count,1) = k(z);
                result(count,2) = n0(i);
                result(count,3) = nb(j);
                result(count,4) = PCS;
                result(count,5) = ANS;
                result(count,6) = CPU_TIME;
                
                resultVAR(count,:)=mean(VAR,2);
                count = count + 1 ;
                
                
        end
    end
end
