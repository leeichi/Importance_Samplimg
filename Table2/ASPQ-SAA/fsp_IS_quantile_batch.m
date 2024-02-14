function [ANS, PCS,CPU_TIME,quantile_TIME,optmize_TIME,dele_TIME,sample_TIME]=fsp_IS_quantile_batch(n0,p,k,nb,delta,mu,va,trial)%%%fsp_CMC_General(n0,p,k,c,v,nb,delta,va)�A�֤Fc v va
%------�ۭq-----
% clear;

n0 = n0; %��l���q�˥�
MU_Y = []; % ��l��������
MU_C = []; % ��l���q��������
final = [];
p = p;
delta = delta ; %�L�t�O�϶�
k = k; %�t�έӼ�
L = nb; %��s�W�v

for i =1:k %%�]�w�U�t�Υ����ƻP�ܲ���
    MU_Y(i) = mu;
    VAR_Y(i) = va;
end

MU_Y(1) = MU_Y(1) + delta;  %%�t��1���L�t�Υ����Ƥj


timer1=tic; %%%���������ɶ�
optmize_TIME = 0; %%%�̨ΤƮɶ�
quantile_TIME = 0; %%%����Ʈɶ�
dele_TIME = 0; %%%�z��ɶ�
sample_TIME = 0;
LR_TIME = 0;

%����a�Ѽ�,�bfunction�̭�,�H�ߤ���=0.05
correct=0;%���T��ܦ���
numwsample=0;%%�������˥���
trial=trial;%�`�����զ���
options = optimset('disp', 'off','TolX',1e-10);%%%matlab���س̨Τƪ��Ѽ�
c = 0.5; %%%bandwidth�Ѽ�
v = 1/2; %%%bandwidth�Ѽ�

%%%�x�s��
FileName = ['IS_N(0,',num2str(VAR_Y(1)^(1/2)),') , p = ',num2str(p),' , n0 = ',num2str(n0),' , nb = ',num2str(L),' , k = ',num2str(k),' , trial = ',num2str(trial),'.mat'];


%-----�`������------
for t = 1:trial
disp(t)
clear obser MU_C ini quantile temp1 temp2 delesys
MU_C = MU_Y    ;%%%�Ĥ@���q�P��l���q�ۦP

index(1:k) = 0;%�P�_�ĴX�Өt�άO�_�w�簣
k_count = 0;%�w�R���t�έӼ�;
r = n0; % sample counter
s = 1; % update stage

band = c / r^(v); %CFD
w = 1;


%---Generate sample batch#1----obser(system_code,sample,LR,batch)

    for i =1:k
        timer3=tic;
        obser(i,1:n0,1) = normrnd(MU_C(i), VAR_Y(i)^(1/2), [1 n0]);%%% ��˭�
        sample_TIME = sample_TIME + toc(timer3);  
        
        timer6=tic;
        obser(i,1:n0,2) = normpdf(obser(i,1:n0,1),MU_Y(i),VAR_Y(i)^(1/2)) ./ normpdf(obser(i,1:n0,1),MU_C(i),VAR_Y(i)^(1/2));%%%% ������
        LR_TIME = LR_TIME + toc(timer6); 
    end
  
%-------�R�t��-------
v = 1/2;
while true
    
    b = c / r^(v); %CFD bandwidth
    w = 1; %%% bandwidth�v��
    
    %%% �R�t��
    [delesys,quantile,variance,temp_quantile_TIME,temp_dele_TIME] = dele_system(index,obser,k,r,delta,s,n0,L,p,MU_Y,VAR_Y, band,w,quantile_TIME,dele_TIME);
    %%% ����Ʈɶ��P�z��ɶ�
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

 
timer2=tic;%%%% �̨ΤƮɶ�
    freq = (r-n0)/L;%%%% ��s�W�v
    if fix (freq) == freq && freq >= 0 % update parameter�Afix ���ܾ��
        for i = 1:k
            if index(i) == 0
                temp1 = squeeze(obser(i,:,:));%%%dimension length�O1���R���A�]�N�O�Ĥ@���A�t�μơA�u�d�Usample LR
                fun = @(MU_New) varfun(MU_New,temp1,s,L,n0,r,p,MU_Y(i),VAR_Y(i),quantile(i),band,w);%%%MU_New���ܼơA�n�D��
                [x,var1] = fmincon(fun,[MU_C(s,i)],[],[],[],[],[],[],[],options);
                MU_C(s+1,i) = x;%%%%�U�@�����|����W�@����output
            end
        end
        s = s + 1;
    end 
optmize_TIME = optmize_TIME + toc(timer2);

%-------�h��Ӽ˥��^�ĤG�B--------

    for i =1 : k
       if index(i) == 0
            timer3=tic;
            obser(i,r+1:r+L,1) = normrnd(MU_C(s,i), VAR_Y(i)^(1/2),[1,L]);%%%���
            sample_TIME = sample_TIME + toc(timer3);
            
            timer6=tic;
            obser(i,r+1:r+L,2) = normpdf(obser(i,r+1:r+L,1),MU_Y(i),VAR_Y(i)^(1/2)) ./ normpdf(obser(i,r+1:r+L,1),MU_C(s,i),VAR_Y(i)^(1/2));%%%%LR
            LR_TIME = LR_TIME + toc(timer6); 
       end
    end


    r = r + L;%�h�⪺�˥�
    v = 1/2;
end
if (index(1))==0%SC
%if(index(1)*index(2)==0)%MDM
    correct = correct+1;
end
numwsample(t)=length(find(obser))/2;%%% �p��һݼ˥���
final = [final;[quantile,variance]];
% save(FileName)
end

PCS=correct/trial;
ANS=sum(numwsample)/trial/k;
CPU_TIME=toc(timer1);
save(FileName)%%% �x�s
end

function [variance] = varfun(MU_New,obs,s,L,n0,r,p,MU_Y,VAR_Y,quantile1,b,w)
    obs = squeeze(obs);
    phi = zeros(1,length(b));
    aphi = zeros(1,length(b));
    new_LR = normpdf(obs(:,1),MU_Y,VAR_Y^(1/2)) ./ normpdf(obs(:,1),MU_New,VAR_Y^(1/2));
    obs_new = obs;
    obs_new(:,2) = obs_new(:,2) .* new_LR ;%%%% ��ⶵ L1*L2    
    
    for i = 1 : length(phi)
        if(p+b(i)>1)
            phi(i) = ((quantile(obs(:,1),1-(1-p)/30)) - quantile(obs(:,1),2*p-1+(1-p)/30))/(29*(1-p)/15);
        else
            phi(i) = (quantile(obs(:,1),p+b(i)) - quantile(obs(:,1),p-b(i)))/(2*b(i));
        end    
    end
    Phi = w * phi';%%% �p��phi�Aw�O1�A�Ǻ�
    
    for i = 0: s-1 
        if(i == 0)%%%original
            temp1 = obs_new((1:n0),:);
        else
            temp1 = obs_new((n0+L*(i-1)+1:min(n0 +L*i,r)),:);%%%�C�����}�B�z�A�@L��
        end
        temp2(i+1) = 1/length(temp1)*((temp1(:,1)>quantile1)'*temp1(:,2)) - (1-p)^2; %%%psi���l
        
%%%%%%%%%%%%%%%%%%%%%%%%%��W��temp1���讳���C
    end
    Psi = r * sum(temp2) / (s^2);%%% �p��psi
    
    variance = (Psi*Phi^2)/r;
    clear temp1 temp2 obs_new
% end

end