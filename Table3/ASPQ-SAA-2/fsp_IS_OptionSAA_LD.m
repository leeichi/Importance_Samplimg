function [ANS, PCS,CPU_TIME,quantile_TIME,optmize_TIME,dele_TIME,sample_TIME]=fsp_IS_OptionSAA_LD(n0,p,k,nb,delta,mu,va,trial,input_S0,input_K,sigma,T,rate,t1,amount)%%%fsp_CMC_General(n0,p,k,c,v,nb,delta,va)�A�֤Fc v va
%------�ۭq-----
% clear;

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

%%%%%%%%%%%LD�ѼƳ]�w
q=0;S=input_S0;K=input_K;deltaT=1;
tau=T-t1;

FileName = ['IS_OptionSAA_N(0,',num2str(VAR_Y(1)^(1/2)),') , p = ',num2str(p),' , n0 = ',num2str(n0),' , nb = ',num2str(L),' , k = ',num2str(k),' , trial = ',num2str(trial),'.mat'];

%-----�`������------
for t = 1:trial
clear obser MU_C ini quantile temp1 temp2 delesys
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
    
    [delesys,quantile,variance,temp_quantile_TIME,temp_dele_TIME] = dele_system(index,obser_last,k,r_last,delta,s_last,n0,L,p,MU_Y,VAR_Y, band,w,quantile_TIME,dele_TIME,s,amount,r,obser);
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
                %%%%%%%LD%%%%%%%%%          
                [theta,LDdelta,gamma]=theta_delta_gamma(rate,q,sigma,S(i),K(i),tau);
                a0=-(-theta.c)*1*deltaT;
                a=-(-LDdelta.c)*ones(1,1);
                A=-0.5*(diag((-gamma.c)*ones(1,1)));
                Sigma=diag((S(i)*sigma)^2*ones(1,1)*deltaT);
                [V,D]=eig(Sigma);
                Ctil=V*sqrt(D);
                [U,Lambda]=eig(Ctil'*A*Ctil);
                C=Ctil*U;
                LDb=a'*C;
                LDb=LDb';
                lambda=diag(Lambda);

                temp=fzero(@(the)phi_prime(the,lambda,LDb,quantile(i),a0),0);%%LD get theta
                theta_IS=temp;
                LDphi=phi_f(theta_IS,lambda,LDb);
                sigma2_new=(1./(1-2*theta_IS*lambda));
                mu_new=theta_IS*LDb.*sigma2_new;

        MU_C(s+1,i)=mu_new(1);
        VAR_C(i)=sigma2_new(1);
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
            [obser(i,r+1:r+L,:),sample_path(i,r+1:r+L),temp_LR_TIME] = OptionSampling_IS(input_S0(i),input_K(i),nb,sigma,T,rate,MU_Y(i),VAR_Y(i),MU_C(s,i),VAR_C(i),LR_TIME,t1);%%%%���t�令 LD �᪺���t%%
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
    correct = correct+1;
end
numwsample(t)=length(find(obser))/2;
final = [final;[quantile,variance]];
% save(FileName)
end

PCS=correct/trial;
ANS=sum(numwsample)/trial/k;
CPU_TIME=toc(timer1);
save(FileName)
end