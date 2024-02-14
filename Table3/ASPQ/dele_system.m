function [delesys,estimate,var_sample,temp_quantile_TIME,temp_dele_TIME] = dele_system(index,obser,k,r,delta,p,b,w,quantile_TIME,dele_TIME)
    para_a = -log(2*0.05/(k-1));%%%step 0的a
    estimate = [];%%%各個系統的分位數
    var_sample = [];%%%各個系統的變異數
    

timer2=tic;
    for i = 1:k
        estimate(i) = quantile(obser(i,:),p);
        var_sample(i) = asym_var(obser(i,:),b,w,r,p,estimate(i));
    end
temp_quantile_TIME = quantile_TIME + toc(timer2);

timer3=tic;
    temp1 = estimate';%%%%%step 3 刪除判斷式
    temp1(1:k,2) = 1;%%%% i
    temp2(2,1:k) = -estimate;%%%% h
    temp2(1,1:k) = 1;
    estimate_diff = temp1*temp2;%%%內積Yi - Yh，對角線因為先減變為0
    clear temp1 temp2;

    temp1 = var_sample';%%%%%step 3 刪除判斷式
    temp1(1:k,2) = 1;
    temp2(2,1:k) = var_sample;
    temp2(1,1:k) = 1;
    indiff = temp1*temp2;%%%內積Si2 + Sh2，對角線因為先加不為0
    clear temp1 temp2;
    
    indiff= delta/2-indiff*para_a/delta;%%%step 3 min 裡面的東西
    indiff = indiff-indiff.*eye(size(indiff)); %對角線歸零para_a
    indiff(indiff>0)=0; %min{0,複雜的式子}
    delesys=[estimate_diff < indiff]; %判斷刪除式子(step3)，logic 1是 or 0否

    for i =1:k
        if index(i)>=1
            delesys(i,1:k)=false;%%%%%%%%%0
        end
    end
temp_dele_TIME = dele_TIME + toc(timer3);
end
