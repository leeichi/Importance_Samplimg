function [delesys,estimate,estimate_bw,var_sample,psi_all_stage,temp_quantile_TIME,temp_dele_TIME] = dele_system_IS(index,obser,allobser,k,r,delta,p,b,w,quantile_initial,quantile_initial_bw,ak,stage,psi_all_previous,n0,C,alpha,quantile_TIME,dele_TIME)
    para_a = -log(2*0.05/(k-1));%%%step 0��a
    estimate = [];%%%�U�Өt�Ϊ������
    estimate_bw=[];
    var_sample = [];%%%�U�Өt�Ϊ��ܲ���
    psi_all_stage =0;
    
%     b = [0.1 / sqrt(r) , 0.2 / sqrt(r) ] ;
%     w = [4/3 , -1/3];
timer3=tic;
    for i = 1:k
        if (index(i))==0
            estimate(i) = CalquantileSA_IS(obser(i,:,:),p,stage,quantile_initial(i),ak,C,alpha);

            if(p+b > 1)%  �p��_�l�� bandwidth
              estimate_bw(i,1) = quantile(allobser(i,:,1),1-(1-p)/30);% ��l���q�p�G�ܨu�� �o��| = 0
              estimate_bw(i,2) = quantile(allobser(i,:,1),2*p-1+(1-p)/30);
            else
              estimate_bw(i,1) = quantile(allobser(i,:,1),p+b);   % �p��_�l�� bandwidth_1 (p+b)
              estimate_bw(i,2) = quantile(allobser(i,:,1),p-b);   % �p��_�l�� bandwidth_2 (p-b)    
            end

            psi_all_stage(i) = CalPsi(obser(i,:,:),p,stage,estimate(i),psi_all_previous(i),n0);
            %%%%�[�Wn0�Ppsi_all_previous
            
            var_sample(i) = asym_varSA_IS(p,b,w,r,estimate_bw(i,:),psi_all_stage(i),stage);
        else
            estimate(i) = 0;
            estimate_bw(i,1:2) = 0;
            var_sample(i) = 0;
        end
    end
temp_quantile_TIME = quantile_TIME + toc(timer3);

timer4=tic;
    
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
    
temp_dele_TIME = dele_TIME + toc(timer4);
end
