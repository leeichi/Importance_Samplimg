function [PCS, ANS,CPU_TIME,quantile_TIME,optmize_TIME,dele_TIME,sample_TIME]=ASPQ_IS_OptionSA(p,k,n0,nb,delta,mu,va,trial,c,v,ak,tk,C,alpha,input_S0,input_K,sigma,T,rate,t1)%%��Jn0 p k c v nb delta va
% n0 : �_�l��˼�
% p : quantile
% k : �t�μ�
% c,v : CFD bandwidth parameter b = c / r ^(v)
% nb : ���ᶥ�q�n�⪺�˥���
% delta : ���t�Τ��L�t�Φhdelta�AIZ parameter
% va : �ܲ���

%------�ۭq-----
% clear;
rng(2)
timer1=tic;
optmize_TIME = 0;
quantile_TIME = 0;
sample_TIME = 0;
dele_TIME = 0;
LR_TIME = 0;
VAR_Y = []; % Original state parameter
MU_Y = []; % Original state parameter
MU_C = []; %�̷s�@������˰Ѽ�
final_quantile = [];

p = p;
k = k; %�t�έӼ�

n0 = n0;%�Ψ�ĴX�Ӽ˥��}�lsetp2
nb = nb;%�h�⪺�ƶq
t1=t1;
delta = delta;
c = c; %fd�Ѽ�
v = v;
ak = ak;

tk = tk;
w = 1;


for i =1:k
    MU_Y(i) = mu;%%������
    VAR_Y(i) = va;%%�ܲ���
end


%����a�Ѽ�,�bfunction�̭�,�H�ߤ���=0.05
correct=zeros(1,k);%���T��ܦ���
numwsample=[];%#wholesample
trial=trial;%�`�����զ���
FileName =  ['IS_OptionSA_N(0,',num2str(VAR_Y(1)^(1/2)),') , p = ',num2str(p),' , n0 = ',num2str(n0),' , nb = ',num2str(nb),' , k = ',num2str(k),' , trial = ',num2str(trial),' , ak = ',num2str(ak),' , tk = ',num2str(tk),'.mat'];

%-----�`������------
for t = 1:trial
    clear obser count_num stage MU_C k_count allobser
    MU_C = MU_Y ;
    count_num =0;
    stage = 1;
    psi_all_previous = zeros(k,1);
    index(1:k) = 0;%�ĴX�Өt�άO�_�w�簣
    k_count = 0;%�w�R���t�έӼ�;
    r = n0; % sample counter
    
%-------------------------------���ͪ�l���q�˥��P�p��_�l��-------------------------obser(system_code,sample,LR)
    b = c / r^(v);

    for i =1:k %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%�_�l�ѩ�� �ζ��ǲέp��
        timer5=tic;%%sample_path(�t��,�˥���,path)
        
        [obser(i,1:n0,:),sample_path(i,1:n0),temp_LR_TIME] = OptionSampling_IS(input_S0(i),input_K(i),n0,sigma,T,rate,MU_Y(i),VAR_Y(i),MU_C(i),VAR_Y(i),LR_TIME,t1);
        LR_TIME = temp_LR_TIME;
        
        count_num = count_num + n0;
        sample_TIME = sample_TIME + toc(timer5);
        
        timer3=tic;
        quantile_first(i) = quantile(obser(i,:,1),p);   %%%%�p��_�l��
        quantile_initial(i) = quantile_first(i);  

        if(p+b > 1)%%%%%%%%%% �p��_�l�� bandwidth
            quantile_first_bw(i,1) = quantile(obser(i,:,1),1-(1-p)/30);%%��l���q�p�G�ܨu�� �o��| = 0
            quantile_first_bw(i,2) = quantile(obser(i,:,1),2*p-1+(1-p)/30);
            quantile_initial_bw(i,1) = quantile_first_bw(i,1);
            quantile_initial_bw(i,2) = quantile_first_bw(i,2);
        else
            quantile_first_bw(i,1) = quantile(obser(i,:,1),p+b);   %%%%�p��_�l�� bandwidth_1 (p+b)
            quantile_first_bw(i,2) = quantile(obser(i,:,1),p-b);   %%%%�p��_�l�� bandwidth_2 (p-b)
            quantile_initial_bw(i,1) = quantile_first_bw(i,1);  
            quantile_initial_bw(i,2) = quantile_first_bw(i,2);
        end
    quantile_TIME = quantile_TIME + toc(timer3); 
    allobser(i,:)=obser(i,:,1);
    end

%-------------------------------���ͪ�l���q�˥��P�p��_�l��-------------------------obser(system_code,sample,LR)
%-------�R�t��-------
while true
    
    
    
%------------------------------���ͤU�@���q����ˤ��t-------------------------------- Algorithm 3.2
 timer2=tic;%%%%�p��̨ΤƮɶ�   
    for i = 1 : k    

        [temp_MU_gradient] = update_Optiontheta_gradient(obser(i,:,:),p,stage,MU_C(i),quantile_initial(i),tk,MU_Y(i),VAR_Y(i),MU_C(i),C,alpha,sample_path(i,:));
        MU_C(i) = temp_MU_gradient;
        clear temp_MU_gradient
        

    end
 optmize_TIME = optmize_TIME + toc(timer2);   
%------------------------------���ͤU�@���q����ˤ��t-------------------------------- Algorithm 3.2 
    
%------------------------------��U�@���q�˥�---------------------------------------
    obser = [];
    sample_path = [];
    timer5=tic;
    for i = 1 : k
        if (index(i))==0
            [obser(i,1:nb,1:2),sample_path(i,1:nb),temp_LR_TIME] = OptionSampling_IS(input_S0(i),input_K(i),nb,sigma,T,rate,MU_Y(i),VAR_Y(i),MU_C(i),VAR_Y(i),LR_TIME,t1);
            count_num = count_num + nb;%%%%�p��ANS
            
            LR_TIME = temp_LR_TIME;
            allobser(i,r+1:r+nb)=obser(i,1:nb,1);

        else
            obser(i,1:nb,1:2) = 0;
            sample_path(i,1:nb) = 0;
        end
    end
    sample_TIME = sample_TIME + toc(timer5);
    r = r + nb;
    b = c / r^(v);
%------------------------------��U�@���q�˥�---------------------------------------    

%-------------------------------------�R�t��---------------------------------------
    [delesys,temp_quan,temp_quanbw,variance,temp_all_stage,temp_quantile_TIME,temp_dele_TIME] = dele_system_IS(index,obser,allobser,k,r,delta,p,b,w,quantile_initial,quantile_initial_bw,ak,stage,psi_all_previous,n0,C,alpha,quantile_TIME,dele_TIME);
    
         for i=1:k
                 if (index(i))==0
                    if  temp_quan(i)>quantile_first(i)*1.3
                        temp_quan(i)=quantile_first(i)*1.3 ;
                    end
                    if  temp_quan(i)<quantile_first(i)*0.7
                        temp_quan(i)=quantile_first(i)*0.7 ;
                    end
                end
         end


    quantile_initial = temp_quan;
    quantile_initial_bw = temp_quanbw;
    psi_all_previous = temp_all_stage;
    quantile_TIME = temp_quantile_TIME;
    dele_TIME = temp_dele_TIME;%%%%�z��ɶ�
    clear temp_quan temp_quanbw temp_all_stage temp_quantile_TIME temp_dele_TIME
    
    %%%��X�U�@���q���_�l��
    stage = stage+1;
    
    timer4=tic;
    if(any(any(delesys)) == true) %determine whether there's a system to be deleted�A1
        [a,b] = find(delesys == true);%%a�����Өt��
        a = unique(a); %%%���ƥX�{���R��
        
            for i =1:length(a)
                index(a(i)) = index(a(i)) + 1;%�����Q�R�����t��
                    if a(i) == 1%%test
                        check_qwe(t) = stage;
                    end
            end
        k_count = k_count+length(a);%�����Q�R�t�μ�
        
        if(k_count >= k-1)%%%%�Ѥ@�Өt�Ϋh����
            break;
        end
    end
    clear a b delesys;
    dele_TIME = dele_TIME + toc(timer4);
%-------------------------------------�R�t��---------------------------------------

end
%-------------------------------------�p�⥿�T���---------------------------------------
if (index(1))==0
    correct(1) = correct(1)+1;
elseif(index(2))==0
    correct(2) = correct(2)+1;
elseif(index(3))==0
    correct(3) = correct(3)+1;
elseif(index(4))==0
    correct(4) = correct(4)+1;
elseif(index(5))==0
    correct(5) = correct(5)+1;
elseif(index(6))==0
    correct(6) = correct(6)+1;
elseif(index(7))==0
    correct(7) = correct(7)+1;
elseif(index(8))==0
    correct(8) = correct(8)+1;
elseif(index(9))==0
    correct(9) = correct(9)+1;
elseif(index(10))==0
    correct(10) = correct(10)+1;    
end

numwsample(t)=count_num;
% numwsample(t)=length(find(obser));%%%�䤣��0������
end

PCS=correct(1)/trial;
ANS=sum(numwsample)/trial/k;
CPU_TIME=toc(timer1);
save(FileName)
end
