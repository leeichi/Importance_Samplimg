function [PCS, ANS,CPU_TIME,quantile_TIME,optmize_TIME,dele_TIME,sample_TIME]=ASPQ_SA_IS(p,k,n0,nb,delta,mu0,va,trial,c,v,ak,C,alpha)%%輸入n0 p k c v nb delta va
% n0 : 起始抽樣數
% p : quantile
% k : 系統數
% c,v : CFD bandwidth parameter b = c / r ^(v)
% nb : 之後階段要抽的樣本數
% delta : 此系統比其他系統多delta，IZ parameter
% va : 變異數
% clear;
%rng(10)
%------------------------------------------------------------setting
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
options = optimset('disp', 'off','TolX',1e-10);
delta = delta;
c = c; %fd參數
v = v;
ak = ak;
w = 1;
L=nb;
Xp=norminv(p,0,va^(1/2));
s=1;
%------------------------------------------------------------ set k system
for i =1:k
    MU_Y(i) = mu0;%%平均數
    VAR_Y(i) = va;%%變異數
end

MU_Y(1) = MU_Y(1)+delta;%%%比其他系統多一個delta



%有個a參數,在function裡面,信心水準=0.05
correct=zeros(1,k);%正確選擇次數
numwsample=[];%#wholesample
trial=trial;%循環測試次數
FileName =  ['IS_SA_N(0,',num2str(VAR_Y(1)^(1/2)),') , p = ',num2str(p),' , n0 = ',num2str(n0),' , nb = ',num2str(nb),' , k = ',num2str(k),' , trial = ',num2str(trial),' , ak = ',num2str(ak),'.mat'];

%-----trial start------
for t = 1:trial
    disp(t)
    clear obser count_num stage MU_C k_count
    MU_C = MU_Y ;
    count_num =0;
    stage = 1;
    psi_all_previous = zeros(k,1);
    index(1:k) = 0;%第幾個系統是否已剔除
    k_count = 0;%已刪除系統個數;
    r = n0; % sample counter
    MU_New=MU_C;
    allobser=[];
%-------------------------------產生初始階段樣本與計算起始解-------------------------obser(system_code,sample,LR)
    b = c / r^(v);
s=1;

    for i =1:k % 起始解抽樣 用順序統計算
        timer5=tic;
        obser(i,1:n0,1) = normrnd(MU_Y(i), VAR_Y(i)^(1/2), [1 n0]); % %每個系統抽取的觀察值個數n0，大小1*n0，obser(系統數，抽取樣本數，1)
        sample_TIME = sample_TIME + toc(timer5);
        
        timer6=tic;
        obser(i,1:n0,2) = normpdf(obser(i,1:n0,1),MU_Y(i),VAR_Y(i)^(1/2)) ./ normpdf(obser(i,1:n0,1),MU_C(i),VAR_Y(i)^(1/2));%%%%LR
        count_num = count_num + n0;
        LR_TIME = LR_TIME + toc(timer6);

        timer3=tic;
        quantile_first(i) = quantile(obser(i,:,1),p);   % 計算起始解
        quantile_initial(i) = quantile_first(i);  

        if(p+b > 1)%  計算起始解 bandwidth
            quantile_first_bw(i,1) = quantile(obser(i,:,1),1-(1-p)/30);% 初始階段如果很罕見 這邊會 = 0
            quantile_first_bw(i,2) = quantile(obser(i,:,1),2*p-1+(1-p)/30);

        else
            quantile_first_bw(i,1) = quantile(obser(i,:,1),p+b);   % 計算起始解 bandwidth_1 (p+b)
            quantile_first_bw(i,2) = quantile(obser(i,:,1),p-b);   % 計算起始解 bandwidth_2 (p-b)

        end
    quantile_TIME = quantile_TIME + toc(timer3);

    allobser(i,1:n0,:)=obser(i,1:n0,:);
    end
    
%-------------------------------產生初始階段樣本與計算起始解-------------------------obser(system_code,sample,LR)
%-------刪系統-------
while true
    
    
    
%------------------------------產生下一階段的抽樣分配-------------------------------- Algorithm 3.2
 timer2=tic;% 計算最佳化時間   
 freq = (r-n0)/L;%%%% 更新頻率
         
    if fix (freq) == freq && freq >= 0 % update parameter，fix 取至整數
        for i = 1:k
            if index(i) == 0
                band = c / r^(v);
                temp1 = squeeze(allobser(i,:,:));%%%dimension length是1的刪掉，也就是第一項，系統數，只留下sample LR
                fun = @(MU_New) varfun(MU_New,temp1,s,L,n0,r,p,MU_Y(i),VAR_Y(i),quantile_initial(i),band,w);%%%MU_New為變數，要求解          
                [x,var1] = fmincon(fun,[MU_C(i)],[],[],[],[],[],[],[],options);
                MU_C(i)=x;
            end
        end        
        s = s + 1;
    end 

 optmize_TIME = optmize_TIME + toc(timer2);   
%------------------------------產生下一階段的抽樣分配-------------------------------- Algorithm 3.2 
    
%------------------------------抽下一階段樣本---------------------------------------
    obser = [];
    
    for i = 1 : k
        if (index(i))==0
            timer5=tic;
            obser(i,1:nb,1) = normrnd(MU_C(i), VAR_Y(i)^(1/2), [1 nb]);% 抽樣
            sample_TIME = sample_TIME + toc(timer5);
            
            timer6=tic;
            obser(i,1:nb,2) = normpdf(obser(i,1:nb,1),MU_Y(i),VAR_Y(i)^(1/2)) ./ normpdf(obser(i,1:nb,1),MU_C(i),VAR_Y(i)^(1/2));
            count_num = count_num + nb;% 計算ANS
            LR_TIME = LR_TIME + toc(timer6);

            allobser(i,r+1:r+nb,1)=obser(i,1:nb,1);
            allobser(i,r+1:r+nb,2)=obser(i,1:nb,2);
        
        else
            timer5=tic;
            obser(i,1:nb,:) = 0;
            sample_TIME = sample_TIME + toc(timer5);
        end  
    end
    
    r = r + nb;
    b = c / r^(v);
%------------------------------抽下一階段樣本---------------------------------------    

%-------------------------------------刪系統---------------------------------------
    [delesys,temp_quan,temp_quanbw,variance,temp_all_stage,temp_quantile_TIME,temp_dele_TIME] = dele_system_IS(index,obser,allobser,k,r,delta,p,b,w,quantile_initial,ak,stage,psi_all_previous,n0,C,alpha,quantile_TIME,dele_TIME);
 %%%%%%%%%%%%%%%%%%Projection%%%%%%%%%%%%
    if (index(1))==0
        if  temp_quan(1)>(Xp+delta)*1.3
            temp_quan(1)=(Xp+delta)*1.3 ;
        end
        
        if  temp_quan(1)<(Xp+delta)*0.7
            temp_quan(1)=(Xp+delta)*0.7 ;
        end
   end
    
   for i=2:k
       if (index(i))==0
          if  temp_quan(i)>Xp*1.3
              temp_quan(i)=Xp*1.3 ;
          end
          
          if  temp_quan(i)<Xp*0.7
              temp_quan(i)=Xp*0.7 ;
          end
      end
   end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
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
        a = unique(a); % 重複出現的刪掉
        
            for i =1:length(a)
                index(a(i)) = index(a(i)) + 1; % 紀錄被刪除的系統
            end
        k_count = k_count+length(a); % 紀錄被刪系統數
        
        if(k_count >= k-1) % 剩一個系統則停止
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

function [variance] = varfun(MU_New,obs,s,L,n0,r,p,MU_Y,VAR_Y,quan,b,w)
    obs = squeeze(obs);
    phi = zeros(1,length(b));
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
        temp2(i+1) = 1/length(temp1)*((temp1(:,1)>quan)'*temp1(:,2)) - (1-p)^2; %%%psi式子
        %%%%%%%%%%%%%%%%%%%%%%%%%把上面temp1平方拿掉。
    end
    Psi = r * sum(temp2) / (s^2);%%% 計算psi
    variance = (Psi*Phi^2)/r;
    clear temp1 temp2 obs_new
% end

end
