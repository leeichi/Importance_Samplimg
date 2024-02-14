function [delesys,var_sample,overall_estimate] = dele_system(b_size,index,obser,k,r,delta,p,estimate)
    para_a = -log(2*0.05/(k-1));
    overall_estimate = [];
    b_num = r/b_size;

    for i = 1:k
        overall_estimate(i) = quantile(obser(i,:),p);
        var_sample(i) = 1/(b_num-1)*sum((estimate(i,:)-overall_estimate(i)).^2);
 
    end
    
    temp1 = overall_estimate';
    temp1(1:k,2) = 1;
    temp2(2,1:k) = -overall_estimate;
    temp2(1,1:k) = 1;
    estimate_diff = temp1*temp2;
    clear temp1 temp2;
    
    
    temp1 = var_sample';
    temp1(1:k,2) = 1;
    temp2(2,1:k) = var_sample;
    temp2(1,1:k) = 1;
    indiff = temp1*temp2;
    clear temp1 temp2;
    
    indiff= delta/2-indiff/b_num*para_a/delta;
    indiff = indiff-indiff.*eye(size(indiff)); %�﨤�u�k�spara_a
    indiff(indiff>0)=0; %min{0,���������l}
    delesys=[estimate_diff < indiff]; %�P�_�R�����l(step3)

    for i =1:k
        if index(i)>=1
            delesys(i,1:k)=false;
        end
    end
end
