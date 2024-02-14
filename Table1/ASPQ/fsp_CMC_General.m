function [PCS, ANS,CPU_TIME,quantile_TIME,dele_TIME,sample_TIME]=fsp_CMC_General(n0,p,k,c,v,nb,delta,mu,va,trial)%%块n0 p k c v nb delta va
% n0 : 癬﹍┾妓计
% p : quantile
% k : ╰参计
% c,v : CFD bandwidth parameter b = c / r ^(v)
% nb : ぇ顶琿璶┾妓セ计
% delta : ╰参ゑㄤ╰参deltaIZ parameter
% va : 跑钵计

%------璹-----

rng(13) %%%北睹计贺
timer1=tic;
quantile_TIME = 0;
dele_TIME = 0;
sample_TIME = 0;
store_TIME = 0;


n0 = n0;%ノ材碭妓セ秨﹍setp2
VAR_Y = []; % Original state parameter
MU_Y = []; % Original state parameter
MU_C = []; % Current sampling parameter
final_quantile = [];
p = p;

delta = delta ;
k = k; %╰参计
c = c; %fd把计
v = v;

nb = nb;%┾计秖



for i =1:k
    MU_Y(i) = mu;%%キА计
    VAR_Y(i) = va;%%跑钵计
end

MU_Y(1) = MU_Y(1)+delta;%%%ゑㄤ╰参delta



%Τa把计,function柑,獺み非=0.05
correct=0;%タ絋匡拒Ω计
numwsample=[];%#wholesample
trial=trial;%碻吏代刚Ω计
FileName = ['CMC_N(0,',num2str(VAR_Y(1)^(1/2)),') , p = ',num2str(p),' , n0 = ',num2str(n0),' , nb = ',num2str(nb),' , k = ',num2str(k),' , trial = ',num2str(trial),'.mat'];

%-----碻吏代刚------
for t = 1:trial
clear obser MU_C

stage = 1;
index(1:k) = 0;%材碭╰参琌埃
k_count = 0;%埃╰参计;
r = n0; % sample counter

%---Generate sample 
timer2=tic;
    for i =1:k
        obser(i,1:n0,1) = normrnd(MU_Y(i), VAR_Y(i)^(1/2), [1 n0]);%%%%%%%%–╰参┾芠诡计n01*n0obser(╰参计┾妓セ计1)
    end
sample_TIME = sample_TIME + toc(timer2);
%-------╰参-------
while true
    b = c / r^(v); %CFD bandwidth
    
    [delesys,quantile,variance,temp_quantile_TIME,temp_dele_TIME] = dele_system(index,obser,k,r,delta,p,b,quantile_TIME,dele_TIME);
    quantile_TIME = temp_quantile_TIME;
    dele_TIME = temp_dele_TIME;
    clear temp_quantile_TIME temp_dele_TIME
    
    
    timer4=tic;
    if(any(any(delesys)) == true) %determine whether there's a system to be deleted1
        [a,b] = find(delesys == true);%%a╰参
        a = unique(a); %%%狡瞷奔
        
            for i =1:length(a)
                index(a(i)) = index(a(i)) + 1;%魁砆埃╰参
            end
        k_count = k_count+length(a);%魁砆╰参计
        
        if(k_count >= k-1)%%%%逞╰参玥氨ゎ
            break;
        end
    end
    clear a b delesys;
    dele_TIME = dele_TIME + toc(timer4);
% var_first(stage,1:k) = variance; 
stage = stage + 1;

%-------┾妓セ材˙--------
timer2=tic;
    for i =1 : k
       if index(i) == 0 %%%临⊿砆奔
            
           
           store_temp = normrnd(MU_Y(i), VAR_Y(i)^(1/2),[1 nb]);
           
           timer3=tic;
           obser(i,r+1:r+nb,1) = store_temp;
           store_TIME = store_TIME + + toc(timer3);
           
           clear store_temp
       end
    end
    r = r + nb;%┾妓セ
sample_TIME = sample_TIME + toc(timer2);

end
if (index(1))==0%SC
    correct = correct+1;
end
numwsample(t)=length(find(obser));%%%тぃ0
end
PCS=correct/trial;
ANS=sum(numwsample)/trial/k;

CPU_TIME=toc(timer1);
save(FileName)
end
