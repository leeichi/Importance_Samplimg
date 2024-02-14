function [PCS, ANS,CPU_TIME,quantile_TIME,optmize_TIME,dele_TIME,sample_TIME]=ASPQ_IS_OptionSA(p,k,n0,nb,delta,mu,va,trial,c,v,ak,tk,C,alpha,input_S0,input_K,sigma,T,rate,t1)%%輸入n0 p k c v nb delta va
% n0 : 起始抽樣數
% p : quantile
% k : 系統數
% c,v : CFD bandwidth parameter b = c / r ^(v)
% nb : 之後階段要抽的樣本數
% delta : 此系統比其他系統多delta，IZ parameter
% va : 變異數

%------自訂-----
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
MU_C = []; %最新一期的抽樣參數
final_quantile = [];

p = p;
k = k; %系統個數

n0 = n0;%用到第幾個樣本開始setp2
nb = nb;%多抽的數量
t1=t1;
delta = delta;
c = c; %fd參數
v = v;
ak = ak;

tk = tk;
w = 1;


for i =1:k
    MU_Y(i) = mu;%%平均數
    VAR_Y(i) = va;%%變異數
end


%有個a參數,在function裡面,信心水準=0.05
correct=zeros(1,k);%正確選擇次數
numwsample=[];%#wholesample
trial=trial;%循環測試次數
FileName =  ['IS_OptionSA_N(0,',num2str(VAR_Y(1)^(1/2)),') , p = ',num2str(p),' , n0 = ',num2str(n0),' , nb = ',num2str(nb),' , k = ',num2str(k),' , trial = ',num2str(trial),' , ak = ',num2str(ak),' , tk = ',num2str(tk),'.mat'];

%-----循環測試------
for t = 1:trial
    clear obser count_num stage MU_C k_count allobser
    MU_C = MU_Y ;
    count_num =0;
    stage = 1;
    psi_all_previous = zeros(k,1);
    index(1:k) = 0;%第幾個系統是否已剔除
    k_count = 0;%已刪除系統個數;
    r = n0; % sample counter
    
%-------------------------------產生初始階段樣本與計算起始解-------------------------obser(system_code,sample,LR)
    b = c / r^(v);

    for i =1:k %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%起始解抽樣 用順序統計算
        timer5=tic;%%sample_path(系統,樣本數,path)
        
        [obser(i,1:n0,:),sample_path(i,1:n0),temp_LR_TIME] = OptionSampling_IS(input_S0(i),input_K(i),n0,sigma,T,rate,MU_Y(i),VAR_Y(i),MU_C(i),VAR_Y(i),LR_TIME,t1);
        LR_TIME = temp_LR_TIME;
        
        count_num = count_num + n0;
        sample_TIME = sample_TIME + toc(timer5);
        
        timer3=tic;
        quantile_first(i) = quantile(obser(i,:,1),p);   %%%%計算起始解
        quantile_initial(i) = quantile_first(i);  

        if(p+b > 1)%%%%%%%%%% 計算起始解 bandwidth
            quantile_first_bw(i,1) = quantile(obser(i,:,1),1-(1-p)/30);%%初始階段如果很罕見 這邊會 = 0
            quantile_first_bw(i,2) = quantile(obser(i,:,1),2*p-1+(1-p)/30);
            quantile_initial_bw(i,1) = quantile_first_bw(i,1);
            quantile_initial_bw(i,2) = quantile_first_bw(i,2);
        else
            quantile_first_bw(i,1) = quantile(obser(i,:,1),p+b);   %%%%計算起始解 bandwidth_1 (p+b)
            quantile_first_bw(i,2) = quantile(obser(i,:,1),p-b);   %%%%計算起始解 bandwidth_2 (p-b)
            quantile_initial_bw(i,1) = quantile_first_bw(i,1);  
            quantile_initial_bw(i,2) = quantile_first_bw(i,2);
        end
    quantile_TIME = quantile_TIME + toc(timer3); 
    allobser(i,:)=obser(i,:,1);
    end

%-------------------------------產生初始階段樣本與計算起始解-------------------------obser(system_code,sample,LR)
%-------刪系統-------
while true
    
    
    
%------------------------------產生下一階段的抽樣分配-------------------------------- Algorithm 3.2
 timer2=tic;%%%%計算最佳化時間   
    for i = 1 : k    

        [temp_MU_gradient] = update_Optiontheta_gradient(obser(i,:,:),p,stage,MU_C(i),quantile_initial(i),tk,MU_Y(i),VAR_Y(i),MU_C(i),C,alpha,sample_path(i,:));
        MU_C(i) = temp_MU_gradient;
        clear temp_MU_gradient
        

    end
 optmize_TIME = optmize_TIME + toc(timer2);   
%------------------------------產生下一階段的抽樣分配-------------------------------- Algorithm 3.2 
    
%------------------------------抽下一階段樣本---------------------------------------
    obser = [];
    sample_path = [];
    timer5=tic;
    for i = 1 : k
        if (index(i))==0
            [obser(i,1:nb,1:2),sample_path(i,1:nb),temp_LR_TIME] = OptionSampling_IS(input_S0(i),input_K(i),nb,sigma,T,rate,MU_Y(i),VAR_Y(i),MU_C(i),VAR_Y(i),LR_TIME,t1);
            count_num = count_num + nb;%%%%計算ANS
            
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
%------------------------------抽下一階段樣本---------------------------------------    

%-------------------------------------刪系統---------------------------------------
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
    dele_TIME = temp_dele_TIME;%%%%篩選時間
    clear temp_quan temp_quanbw temp_all_stage temp_quantile_TIME temp_dele_TIME
    
    %%%輸出下一階段的起始解
    stage = stage+1;
    
    timer4=tic;
    if(any(any(delesys)) == true) %determine whether there's a system to be deleted，1
        [a,b] = find(delesys == true);%%a為哪個系統
        a = unique(a); %%%重複出現的刪掉
        
            for i =1:length(a)
                index(a(i)) = index(a(i)) + 1;%紀錄被刪除的系統
                    if a(i) == 1%%test
                        check_qwe(t) = stage;
                    end
            end
        k_count = k_count+length(a);%紀錄被刪系統數
        
        if(k_count >= k-1)%%%%剩一個系統則停止
            break;
        end
    end
    clear a b delesys;
    dele_TIME = dele_TIME + toc(timer4);
%-------------------------------------刪系統---------------------------------------

end
%-------------------------------------計算正確選擇---------------------------------------
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
% numwsample(t)=length(find(obser));%%%找不為0的長度
end

PCS=correct(1)/trial;
ANS=sum(numwsample)/trial/k;
CPU_TIME=toc(timer1);
save(FileName)
end
