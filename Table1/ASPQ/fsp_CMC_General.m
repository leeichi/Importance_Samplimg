function [PCS, ANS,CPU_TIME,quantile_TIME,dele_TIME,sample_TIME]=fsp_CMC_General(n0,p,k,c,v,nb,delta,mu,va,trial)%%��Jn0 p k c v nb delta va
% n0 : �_�l��˼�
% p : quantile
% k : �t�μ�
% c,v : CFD bandwidth parameter b = c / r ^(v)
% nb : ���ᶥ�q�n�⪺�˥���
% delta : ���t�Τ��L�t�Φhdelta�AIZ parameter
% va : �ܲ���

%------�ۭq-----

rng(13) %%%����üƺؤl
timer1=tic;
quantile_TIME = 0;
dele_TIME = 0;
sample_TIME = 0;
store_TIME = 0;


n0 = n0;%�Ψ�ĴX�Ӽ˥��}�lsetp2
VAR_Y = []; % Original state parameter
MU_Y = []; % Original state parameter
MU_C = []; % Current sampling parameter
final_quantile = [];
p = p;

delta = delta ;
k = k; %�t�έӼ�
c = c; %fd�Ѽ�
v = v;

nb = nb;%�h�⪺�ƶq



for i =1:k
    MU_Y(i) = mu;%%������
    VAR_Y(i) = va;%%�ܲ���
end

MU_Y(1) = MU_Y(1)+delta;%%%���L�t�Φh�@��delta



%����a�Ѽ�,�bfunction�̭�,�H�ߤ���=0.05
correct=0;%���T��ܦ���
numwsample=[];%#wholesample
trial=trial;%�`�����զ���
FileName = ['CMC_N(0,',num2str(VAR_Y(1)^(1/2)),') , p = ',num2str(p),' , n0 = ',num2str(n0),' , nb = ',num2str(nb),' , k = ',num2str(k),' , trial = ',num2str(trial),'.mat'];

%-----�`������------
for t = 1:trial
clear obser MU_C

stage = 1;
index(1:k) = 0;%�ĴX�Өt�άO�_�w�簣
k_count = 0;%�w�R���t�έӼ�;
r = n0; % sample counter

%---Generate sample 
timer2=tic;
    for i =1:k
        obser(i,1:n0,1) = normrnd(MU_Y(i), VAR_Y(i)^(1/2), [1 n0]);%%%%%%%%�C�Өt�Ω�����[��ȭӼ�n0�A�j�p1*n0�Aobser(�t�μơA����˥��ơA1)
    end
sample_TIME = sample_TIME + toc(timer2);
%-------�R�t��-------
while true
    b = c / r^(v); %CFD bandwidth
    
    [delesys,quantile,variance,temp_quantile_TIME,temp_dele_TIME] = dele_system(index,obser,k,r,delta,p,b,quantile_TIME,dele_TIME);
    quantile_TIME = temp_quantile_TIME;
    dele_TIME = temp_dele_TIME;
    clear temp_quantile_TIME temp_dele_TIME
    
    
    timer4=tic;
    if(any(any(delesys)) == true) %determine whether there's a system to be deleted�A1
        [a,b] = find(delesys == true);%%a�����Өt��
        a = unique(a); %%%���ƥX�{���R��
        
            for i =1:length(a)
                index(a(i)) = index(a(i)) + 1;%�����Q�R�����t��
            end
        k_count = k_count+length(a);%�����Q�R�t�μ�
        
        if(k_count >= k-1)%%%%�Ѥ@�Өt�Ϋh����
            break;
        end
    end
    clear a b delesys;
    dele_TIME = dele_TIME + toc(timer4);
% var_first(stage,1:k) = variance; 
stage = stage + 1;

%-------�h��Ӽ˥��^�ĤG�B--------
timer2=tic;
    for i =1 : k
       if index(i) == 0 %%%�٨S�Q�R��
            
           
           store_temp = normrnd(MU_Y(i), VAR_Y(i)^(1/2),[1 nb]);
           
           timer3=tic;
           obser(i,r+1:r+nb,1) = store_temp;
           store_TIME = store_TIME + + toc(timer3);
           
           clear store_temp
       end
    end
    r = r + nb;%�h��Ӽ˥�
sample_TIME = sample_TIME + toc(timer2);

end
if (index(1))==0%SC
    correct = correct+1;
end
numwsample(t)=length(find(obser));%%%�䤣��0������
end
PCS=correct/trial;
ANS=sum(numwsample)/trial/k;

CPU_TIME=toc(timer1);
save(FileName)
end
