function [delesys estimate,var_sample,temp_quantile_TIME,temp_dele_TIME] = dele_system(index,obser,k,r,delta,s,n0,L,p,MU_Y,VAR_Y,b,w,quantile_TIME,dele_TIME)
    para_a = -log(2*0.05/(k-1));
       
    estimate = [];
    var_sample = [];

timer3=tic;
    for i = 1:k
        estimate(i) = calQuantile(obser(i,:,:),s,L,n0,r,p);
        var_sample(i) = asym_var(obser(i,:,:),b,w,s,L,n0,r,p,estimate(i),MU_Y(i),VAR_Y(i));
    end
temp_quantile_TIME = quantile_TIME + toc(timer3);


timer4=tic;
    temp1 = estimate';
    temp1(1:k,2) = 1;
    temp2(2,1:k) = -estimate;
    temp2(1,1:k) = 1;
    estimate_diff = temp1*temp2;
    clear temp1 temp2;
    
    
    temp1 = var_sample';
    temp1(1:k,2) = 1;
    temp2(2,1:k) = var_sample;
    temp2(1,1:k) = 1;
    indiff = temp1*temp2;
    clear temp1 temp2;
    
    indiff= delta/2-indiff*para_a/delta;
    indiff = indiff-indiff.*eye(size(indiff)); %對角線歸零para_a
    indiff(indiff>0)=0; %min{0,複雜的式子}
    delesys=[estimate_diff < indiff]; %判斷刪除式子(step3)，成立是1

    for i =1:k %%%%已經被刪的甭再算了
        if index(i)>=1
            delesys(i,1:k)=false;
        end
    end
temp_dele_TIME = dele_TIME + toc(timer4);
end
