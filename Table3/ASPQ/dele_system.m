function [delesys,estimate,var_sample,temp_quantile_TIME,temp_dele_TIME] = dele_system(index,obser,k,r,delta,p,b,w,quantile_TIME,dele_TIME)
    para_a = -log(2*0.05/(k-1));%%%step 0��a
    estimate = [];%%%�U�Өt�Ϊ������
    var_sample = [];%%%�U�Өt�Ϊ��ܲ���
    

timer2=tic;
    for i = 1:k
        estimate(i) = quantile(obser(i,:),p);
        var_sample(i) = asym_var(obser(i,:),b,w,r,p,estimate(i));
    end
temp_quantile_TIME = quantile_TIME + toc(timer2);

timer3=tic;
    temp1 = estimate';%%%%%step 3 �R���P�_��
    temp1(1:k,2) = 1;%%%% i
    temp2(2,1:k) = -estimate;%%%% h
    temp2(1,1:k) = 1;
    estimate_diff = temp1*temp2;%%%���nYi - Yh�A�﨤�u�]�������ܬ�0
    clear temp1 temp2;

    temp1 = var_sample';%%%%%step 3 �R���P�_��
    temp1(1:k,2) = 1;
    temp2(2,1:k) = var_sample;
    temp2(1,1:k) = 1;
    indiff = temp1*temp2;%%%���nSi2 + Sh2�A�﨤�u�]�����[����0
    clear temp1 temp2;
    
    indiff= delta/2-indiff*para_a/delta;%%%step 3 min �̭����F��
    indiff = indiff-indiff.*eye(size(indiff)); %�﨤�u�k�spara_a
    indiff(indiff>0)=0; %min{0,���������l}
    delesys=[estimate_diff < indiff]; %�P�_�R�����l(step3)�Alogic 1�O or 0�_

    for i =1:k
        if index(i)>=1
            delesys(i,1:k)=false;%%%%%%%%%0
        end
    end
temp_dele_TIME = dele_TIME + toc(timer3);
end
