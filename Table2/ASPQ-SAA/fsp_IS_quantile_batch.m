function [ANS, PCS,CPU_TIME,quantile_TIME,optmize_TIME,dele_TIME,sample_TIME]=fsp_IS_quantile_batch(n0,p,k,nb,delta,mu,va,trial)%%%fsp_CMC_General(n0,p,k,c,v,nb,delta,va)，少了c v va
%------自訂-----
% clear;

n0 = n0; %初始階段樣本
MU_Y = []; % 初始的平均數
MU_C = []; % 其餘階段的平均數
final = [];
p = p;
delta = delta ; %無差別區間
k = k; %系統個數
L = nb; %更新頻率

for i =1:k %%設定各系統平均數與變異數
    MU_Y(i) = mu;
    VAR_Y(i) = va;
end

MU_Y(1) = MU_Y(1) + delta;  %%系統1比其他系統平均數大


timer1=tic; %%%紀錄全部時間
optmize_TIME = 0; %%%最佳化時間
quantile_TIME = 0; %%%分位數時間
dele_TIME = 0; %%%篩選時間
sample_TIME = 0;
LR_TIME = 0;

%有個a參數,在function裡面,信心水準=0.05
correct=0;%正確選擇次數
numwsample=0;%%全部的樣本數
trial=trial;%循環測試次數
options = optimset('disp', 'off','TolX',1e-10);%%%matlab內建最佳化的參數
c = 0.5; %%%bandwidth參數
v = 1/2; %%%bandwidth參數

%%%儲存檔
FileName = ['IS_N(0,',num2str(VAR_Y(1)^(1/2)),') , p = ',num2str(p),' , n0 = ',num2str(n0),' , nb = ',num2str(L),' , k = ',num2str(k),' , trial = ',num2str(trial),'.mat'];


%-----循環測試------
for t = 1:trial
disp(t)
clear obser MU_C ini quantile temp1 temp2 delesys
MU_C = MU_Y    ;%%%第一階段與初始間段相同

index(1:k) = 0;%判斷第幾個系統是否已剔除
k_count = 0;%已刪除系統個數;
r = n0; % sample counter
s = 1; % update stage

band = c / r^(v); %CFD
w = 1;


%---Generate sample batch#1----obser(system_code,sample,LR,batch)

    for i =1:k
        timer3=tic;
        obser(i,1:n0,1) = normrnd(MU_C(i), VAR_Y(i)^(1/2), [1 n0]);%%% 抽樣值
        sample_TIME = sample_TIME + toc(timer3);  
        
        timer6=tic;
        obser(i,1:n0,2) = normpdf(obser(i,1:n0,1),MU_Y(i),VAR_Y(i)^(1/2)) ./ normpdf(obser(i,1:n0,1),MU_C(i),VAR_Y(i)^(1/2));%%%% 概似比
        LR_TIME = LR_TIME + toc(timer6); 
    end
  
%-------刪系統-------
v = 1/2;
while true
    
    b = c / r^(v); %CFD bandwidth
    w = 1; %%% bandwidth權重
    
    %%% 刪系統
    [delesys,quantile,variance,temp_quantile_TIME,temp_dele_TIME] = dele_system(index,obser,k,r,delta,s,n0,L,p,MU_Y,VAR_Y, band,w,quantile_TIME,dele_TIME);
    %%% 分位數時間與篩選時間
    quantile_TIME = temp_quantile_TIME;
    dele_TIME = temp_dele_TIME;
    clear temp_quantile_TIME temp_dele_TIME
    
timer5=tic; %%%% 篩選時間
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
dele_TIME = dele_TIME + toc(timer5);%%%% 篩選時間

 
timer2=tic;%%%% 最佳化時間
    freq = (r-n0)/L;%%%% 更新頻率
    if fix (freq) == freq && freq >= 0 % update parameter，fix 取至整數
        for i = 1:k
            if index(i) == 0
                temp1 = squeeze(obser(i,:,:));%%%dimension length是1的刪掉，也就是第一項，系統數，只留下sample LR
                fun = @(MU_New) varfun(MU_New,temp1,s,L,n0,r,p,MU_Y(i),VAR_Y(i),quantile(i),band,w);%%%MU_New為變數，要求解
                [x,var1] = fmincon(fun,[MU_C(s,i)],[],[],[],[],[],[],[],options);
                MU_C(s+1,i) = x;%%%%下一階的會等於上一階的output
            end
        end
        s = s + 1;
    end 
optmize_TIME = optmize_TIME + toc(timer2);

%-------多抽個樣本回第二步--------

    for i =1 : k
       if index(i) == 0
            timer3=tic;
            obser(i,r+1:r+L,1) = normrnd(MU_C(s,i), VAR_Y(i)^(1/2),[1,L]);%%%抽樣
            sample_TIME = sample_TIME + toc(timer3);
            
            timer6=tic;
            obser(i,r+1:r+L,2) = normpdf(obser(i,r+1:r+L,1),MU_Y(i),VAR_Y(i)^(1/2)) ./ normpdf(obser(i,r+1:r+L,1),MU_C(s,i),VAR_Y(i)^(1/2));%%%%LR
            LR_TIME = LR_TIME + toc(timer6); 
       end
    end


    r = r + L;%多抽的樣本
    v = 1/2;
end
if (index(1))==0%SC
%if(index(1)*index(2)==0)%MDM
    correct = correct+1;
end
numwsample(t)=length(find(obser))/2;%%% 計算所需樣本數
final = [final;[quantile,variance]];
% save(FileName)
end

PCS=correct/trial;
ANS=sum(numwsample)/trial/k;
CPU_TIME=toc(timer1);
save(FileName)%%% 儲存
end

function [variance] = varfun(MU_New,obs,s,L,n0,r,p,MU_Y,VAR_Y,quantile1,b,w)
    obs = squeeze(obs);
    phi = zeros(1,length(b));
    aphi = zeros(1,length(b));
    new_LR = normpdf(obs(:,1),MU_Y,VAR_Y^(1/2)) ./ normpdf(obs(:,1),MU_New,VAR_Y^(1/2));
    obs_new = obs;
    obs_new(:,2) = obs_new(:,2) .* new_LR ;%%%% 拆成兩項 L1*L2    
    
    for i = 1 : length(phi)
        if(p+b(i)>1)
            phi(i) = ((quantile(obs(:,1),1-(1-p)/30)) - quantile(obs(:,1),2*p-1+(1-p)/30))/(29*(1-p)/15);
        else
            phi(i) = (quantile(obs(:,1),p+b(i)) - quantile(obs(:,1),p-b(i)))/(2*b(i));
        end    
    end
    Phi = w * phi';%%% 計算phi，w是1，甭管
    
    for i = 0: s-1 
        if(i == 0)%%%original
            temp1 = obs_new((1:n0),:);
        else
            temp1 = obs_new((n0+L*(i-1)+1:min(n0 +L*i,r)),:);%%%每階分開處理，共L個
        end
        temp2(i+1) = 1/length(temp1)*((temp1(:,1)>quantile1)'*temp1(:,2)) - (1-p)^2; %%%psi式子
        
%%%%%%%%%%%%%%%%%%%%%%%%%把上面temp1平方拿掉。
    end
    Psi = r * sum(temp2) / (s^2);%%% 計算psi
    
    variance = (Psi*Phi^2)/r;
    clear temp1 temp2 obs_new
% end

end