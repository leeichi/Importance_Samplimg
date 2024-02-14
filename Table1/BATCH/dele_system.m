function [delesys,var_sample,sample_mean] = dele_system(index,estimate,k,r,delta,nb)
    para_a = -log(2*0.05/(k-1));
    
    sample_mean = mean(estimate,2)'; 
    batch_num = r/nb;

    temp1=sample_mean';
    temp1(1:k,2)=1;
    temp2(2,1:k)=-sample_mean;
    temp2(1,1:k)=1;
    mean_diff=temp1*temp2;%不同系統theta_bar的差
    
    
    var_sample=var(estimate,0,2)';
    temp1=var_sample';
    temp1(1:k,2)=1;
    temp2(2,1:k)=var_sample;
    temp2(1,1:k)=1;
    indiff=temp1*temp2;
    clear temp1 temp2;
    
    indiff= delta/2-indiff/batch_num*para_a/delta;
    indiff = indiff-indiff.*eye(size(indiff));%對角線歸零
    indiff(indiff>0)=0;%min{0,複雜的式子}
    delesys=[mean_diff<indiff];%判斷刪除式子(step3)
    for i =1:k
        if index(i)>=1
            delesys(i,1:k)=false;
        end
    end


end