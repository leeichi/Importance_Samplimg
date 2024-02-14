function [ANS, PCS,CPU_TIME,quantile_TIME,optmize_TIME,dele_TIME,sample_TIME]=fsp_IS_OptionSAA(n0,p,k,nb,delta,mu,va,trial,input_S0,input_K,sigma,T,rate,t1,amount)%%%fsp_CMC_General(n0,p,k,c,v,nb,delta,va)�A�֤Fc v va
%------�ۭq-----
% clear;

%%%%%freq fix �C�����|�]��
%%%%%bandwidth code 134-167
%%%%%%%%%code 156 ���O�nL1*L2�� temp1(:,2).^2
%%%%%%%%%code 160 psi

n0 = n0;%�Ψ�ĴX�Ӽ˥��}�lsetp2VAR_Y = []; % Original state parameter
MU_Y = []; % Original state parameter
MU_C = []; % Current sampling parameter
final = [];
p = p;
delta = delta ;
k = k; %�t�έӼ�
L = nb; %��s�W�v
t1=t1;

for i =1:k
    MU_Y(i) = mu;
    VAR_Y(i) = va;
end


timer1=tic;
optmize_TIME = 0; %%%�̨ΤƮɶ�
quantile_TIME = 0; %%%����Ʈɶ�
dele_TIME = 0; %%%�z��ɶ�
sample_TIME = 0;
LR_TIME = 0;


%����a�Ѽ�,�bfunction�̭�,�H�ߤ���=0.05
correct=0;%���T��ܦ���
numwsample=0;%#wholesample
trial=trial;%�`�����զ���
options = optimset('disp', 'off','TolX',1e-10);
c = 0.5;
v = 1/2;

FileName = ['IS_OptionSAA_N(0,',num2str(VAR_Y(1)^(1/2)),') , p = ',num2str(p),' , n0 = ',num2str(n0),' , nb = ',num2str(L),' , k = ',num2str(k),' , trial = ',num2str(trial),'.mat'];

%-----�`������------
for t = 1:trial
clear obser MU_C ini quantil temp1 temp2 delesys
MU_C = MU_Y    ;%%%�{�b�����󤧫e��

index(1:k) = 0;%�ĴX�Өt�άO�_�w�簣
k_count = 0;%�w�R���t�έӼ�;
r = n0; % sample counter
r_last = r;
s = 1; % update stage
s_last = 1;
calculate = 1;%�Ψӧ���᭱�ռƪ�obser&sample
amount = amount;
sample_path = [];
band = c / r^(v); %CFD
w = 1;
%---Generate sample batch#1----obser(system_code,sample,LR,batch)
timer3=tic;
    for i =1:k
        
        [obser(i,1:n0,:),sample_path(i,:),temp_LR_TIME] = OptionSampling_IS(input_S0(i),input_K(i),n0,sigma,T,rate,MU_Y(i),VAR_Y(i),MU_C(i),VAR_Y(i),LR_TIME,t1);
        LR_TIME =  temp_LR_TIME;
        clear temp_LR_TIME
        
    end
obser_last = obser;
sample_path_last = sample_path;
sample_TIME = sample_TIME + toc(timer3); 
    
    
%-------�R�t��-------
v = 1/2;
while true
    
    b = c / r_last^(v); %CFD
    w = 1;
    
    [delesys,quantil,variance,temp_quantile_TIME,temp_dele_TIME] = dele_system(index,obser_last,k,r_last,delta,s_last,n0,L,p,MU_Y,VAR_Y, band,w,quantile_TIME,dele_TIME,s,amount,r,obser);
    quantile_TIME = temp_quantile_TIME;
    dele_TIME = temp_dele_TIME;
    clear temp_quantile_TIME temp_dele_TIME
    
timer5=tic; %%%% �z��ɶ�    
    if(any(any(delesys)) == true) %determine whether there's a system to be deleted
        [a,b] = find(delesys == true);
        a = unique(a); 
        
            for i =1:length(a)
                index(a(i)) = index(a(i)) + 1;%�����Q�R�����t��
            end
        k_count = k_count+length(a);%�����Q�R�t�μ�
        if(k_count >= k-1)
            break;
        end
    end
    clear a b delesys;
dele_TIME = dele_TIME + toc(timer5);%%%% �z��ɶ�


% check whether the parameters need to be updated 
timer2=tic;%%%% �̨ΤƮɶ�
    freq = (r-n0)/L;
    if fix (freq) == freq && freq >= 0 % update parameter�Afix ���ܾ��
        for i = 1:k

            if index(i) == 0
                temp1 = squeeze(obser(i,:,:));%%%dimension length�O1���R���A�]�N�O�Ĥ@���A�t�μơA�u�d�Usample LR
                fun = @(MU_New) varfun(MU_New,temp1,s,L,n0,r,p,MU_Y(i),VAR_Y(i),quantil(i),band,w,sample_path(i,:));%%%MU_New���ܼơA�n�D��
                [x,var1] = fmincon(fun,[MU_C(i)],[],[],[],[],[],[],[],options);
                MU_C(i) = x;%%%%�U�@�����|����W�@����output
                
            end
        end
        s = s + 1;
        s_last = s_last + 1;
   end    
optmize_TIME = optmize_TIME + toc(timer2);

%-------�h��Ӽ˥��^�ĤG�B--------
timer3=tic;
    for i =1 : k
       if index(i) == 0
            [obser(i,r+1:r+L,:),sample_path(i,r+1:r+L),temp_LR_TIME] = OptionSampling_IS(input_S0(i),input_K(i),nb,sigma,T,rate,MU_Y(i),VAR_Y(i),MU_C(i),VAR_Y(i),LR_TIME,t1);
            LR_TIME =  temp_LR_TIME;
            clear temp_LR_TIME
       end
    end


    if s > amount
        obser_last = obser(:,n0+1+L*(calculate-1):n0+L*(calculate-1+amount),:);
        sample_path_last = sample_path(:,n0+1+L*(calculate-1):n0+L*(calculate-1+amount));
        calculate = calculate + 1;
        
        r = r + L;%�h�⪺�˥�
        r_last = amount*L;
        s_last = amount;
        
    else
        obser_last=obser;
        sample_path_last = sample_path;
        r = r + L;%�h�⪺�˥�
        r_last = r;

    end
sample_TIME = sample_TIME + toc(timer3);

    v = 1/2;
end
if (index(1))==0%SC
%if(index(1)*index(2)==0)%MDM
    correct = correct+1;
end
numwsample(t)=length(find(obser))/2;
final = [final;[quantil,variance]];
% save(FileName)
end

PCS=correct/trial;
ANS=sum(numwsample)/trial/k;
CPU_TIME=toc(timer1);
save(FileName)
end

function [variance] = varfun(MU_New,obs,s,L,n0,r,p,MU_Y,VAR_Y,quantil,b,w,sample_path)
    obs = squeeze(obs);
    sample_path = squeeze(sample_path);
    sample_path = sample_path';
    phi = zeros(1,length(b));
    new_LR = normpdf(sample_path(:,:),MU_Y,VAR_Y^(1/2)) ./ normpdf(sample_path(:,:),MU_New,VAR_Y^(1/2));
    new_LR = prod(new_LR,2);
    obs_new = obs;
    obs_new(:,2) = obs_new(:,2) .* new_LR ;%%%% L1*L2
   for i = 1 : length(phi)
        if(p+b(i)>1)
            %phi(i) = ((calQuantile(obs,s,L,n0,r,1-(1-p)/30)) - calQuantile(obs,s,L,n0,r,2*p-1+(1-p)/30))/(29*(1-p)/15); 
            phi(i) = ((quantile(obs(:,1),1-(1-p)/30)) - quantile(obs(:,1),2*p-1+(1-p)/30))/(29*(1-p)/15);
        else
            %phi(i) = (calQuantile(obs,s,L,n0,r,p+b(i)) - calQuantile(obs,s,L,n0,r,p-b(i)))/(2*b(i)); 
            phi(i) = (quantile(obs(:,1),p+b(i)) - quantile(obs(:,1),p-b(i)))/(2*b(i));
        end
        
   end
   Phi = phi';
   for i = 0: s-1   
        if(i == 0)%%%original
            temp1 = obs_new((1:n0),:);
        else
            temp1 = obs_new((n0+L*(i-1)+1:min(n0 +L*i,r)),:);%%%�C�����}�B�z�A�@L��
        end
        temp2(i+1) = 1/length(temp1)*((temp1(:,1)>quantil)'*temp1(:,2)) - (1-p)^2; %%%psi���l%%%%%%%%%%%%%%%%%%%%%%%%%��W��temp1���讳���C
   end
   
  
   Psi = r * sum(temp2) / ((s+1)^2);%%%p26(3.18)
   
    
    variance = (Psi*Phi^2)/r;
    clear temp1 temp2 obs_new
% end

end