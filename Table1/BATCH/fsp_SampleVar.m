function [ ANS,PCS CPU_TIME VAR]=fsp_SampleVar(n0,p,k,nb,delta,va,trial)
% %------�ۭq-----
% clear;
timer1=tic;
n0 = n0;%�Ψ�ĴX�Ӽ˥��}�lsetp2
VAR_Y = []; % Original state parameter
MU_Y = []; % Original state parameter
p = p;
nb = nb;
delta = delta ;
k = k; %�t�έӼ�

for i =1:k
    MU_Y(i) = 0;
    VAR_Y(i) = va;
end

MU_Y(1) = MU_Y(1)+delta;

%����a�Ѽ�,�bfunction�̭�,�H�ߤ���=0.05
correct=0;%���T��ܦ���
numwsample=[];%#wholesample

trial=trial;%�`�����զ���

FileName = ['(nb=',num2str(nb),'p=',num2str(p),')SC_FSPSampleVar_n0=',num2str(n0),'k=',num2str(k),'var=',num2str(VAR_Y(1)),'del=',num2str(delta),'.mat'];

%-----�`������------
for t = 1:trial
    disp(t)
    clear obser delesys index estimate
    estimate = [];
    index(1:k) = 0;%�ĴX�Өt�άO�_�w�簣
    k_count = 0;%�w�R���t�έӼ�;
    r = n0; % sample counter
    %---Generate samples
    for i =1:k
        obser(i,1:n0,1) = normrnd(MU_Y(i), VAR_Y(i)^(1/2), [1 n0]);
        
        for j =1:r/nb
            temp = obser(i,1+nb*(j-1):nb*(j));
            temp = sort(temp);
            estimate(i,j) = temp(ceil(nb*p));
            clear temp
        end
        
    end
    
    %-------�R�t��-------
    while true
        
        [delesys,variance,quan] = dele_system(index,estimate,k,r,delta,nb);
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
        
        
        %-------�h��Ӽ˥��^�ĤG�B--------
        for i =1 : k
            if index(i) == 0
                obser(i,r+1:r+nb) = normrnd(MU_Y(i), VAR_Y(i)^(1/2),[1 nb]);
                temp = obser(i,r+1:r+nb);
                temp = sort(temp);
                estimate(i,(r/nb)+1) = temp(ceil(nb*p));
                clear temp
            end
        end
        r = r + nb;%�h��Ӽ˥�
    end
    if (index(1))==0%SC
        %if(index(1)*index(2)==0)%MDM
        correct = correct+1;
    end
    VARR(t,:)=variance;
    numwsample(t)=length(find(obser));
end
VAR=(mean(VARR,1));
PCS=correct/trial;
ANS=sum(numwsample)/trial/k;
CPU_TIME=toc(timer1);
save(FileName)
end
