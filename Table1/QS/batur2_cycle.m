clear


n0 = [200 400 1000];
p = [0.9];
K = [3 10 50];
nb = [100];
result = [];
delta=1;
count=1;
replication = 500;
for z = 1 : length(K)
for i = 1:length(n0)
    for j = 1:length(nb)
            [ANS,PCS,CPU_TIME] = baturcodeori(n0(i),p,delta,K(z),replication);
            
     %       result(i+3*(j-1),1+3*(z-1)) = PCS;
    %        result(i+3*(j-1),2+3*(z-1)) = ANS;
   %         result(i+3*(j-1),3+3*(z-1)) = CPU_TIME;
            result(count,1) = K(z);
            result(count,2) = n0(i);
            result(count,3) = nb(j);
            result(count,4) = PCS;
            result(count,5) = ANS;
            result(count,6) = CPU_TIME;   
            count = count + 1 ;



    end
end
end
