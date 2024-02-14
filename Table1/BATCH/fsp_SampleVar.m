function [ ANS,PCS CPU_TIME VAR]=fsp_SampleVar(n0,p,k,nb,delta,va,trial)
% %------自訂-----
% clear;
timer1=tic;
n0 = n0;%用到第幾個樣本開始setp2
VAR_Y = []; % Original state parameter
MU_Y = []; % Original state parameter
p = p;
nb = nb;
delta = delta ;
k = k; %系統個數

for i =1:k
    MU_Y(i) = 0;
    VAR_Y(i) = va;
end

MU_Y(1) = MU_Y(1)+delta;

%有個a參數,在function裡面,信心水準=0.05
correct=0;%正確選擇次數
numwsample=[];%#wholesample

trial=trial;%循環測試次數

FileName = ['(nb=',num2str(nb),'p=',num2str(p),')SC_FSPSampleVar_n0=',num2str(n0),'k=',num2str(k),'var=',num2str(VAR_Y(1)),'del=',num2str(delta),'.mat'];

%-----循環測試------
for t = 1:trial
    disp(t)
    clear obser delesys index estimate
    estimate = [];
    index(1:k) = 0;%第幾個系統是否已剔除
    k_count = 0;%已刪除系統個數;
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
    
    %-------刪系統-------
    while true
        
        [delesys,variance,quan] = dele_system(index,estimate,k,r,delta,nb);
        if(any(any(delesys)) == true) %determine whether there's a system to be deleted
            [a,b] = find(delesys == true);
            a = unique(a);
            
            for i =1:length(a)
                index(a(i)) = index(a(i)) + 1;%紀錄被刪除的系統
            end
            k_count = k_count+length(a);%紀錄被刪系統數
            
            if(k_count >= k-1)
                break;
            end
        end
        clear a b delesys;
        
        
        %-------多抽個樣本回第二步--------
        for i =1 : k
            if index(i) == 0
                obser(i,r+1:r+nb) = normrnd(MU_Y(i), VAR_Y(i)^(1/2),[1 nb]);
                temp = obser(i,r+1:r+nb);
                temp = sort(temp);
                estimate(i,(r/nb)+1) = temp(ceil(nb*p));
                clear temp
            end
        end
        r = r + nb;%多抽個樣本
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
