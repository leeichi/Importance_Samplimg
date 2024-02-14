function[empcdf] = empCDF(x,obs,s,r,L,n0,s_previously,amount)
    tt = [];
    if s_previously > amount
        for z = 0: s-1
            if(z == 0)
                t = obs((1:L),:);
            else
                t = obs((L+L*(z-1)+1:min(L +L*(z),r)),:);
            end %%%Z+1?? –顶琿だ秨
            tt(z+1) = (t(:,1) >= x)' * t(:,2) /length(t) ;%%%%P33,程F场だぃ1-
        end
        empcdf = 1-mean(tt);%%%%1-
   
    else
        for z = 0: s-1
            if(z == 0)
                t = obs((1:n0),:);
            else
                t = obs((n0+L*(z-1)+1:min(n0 +L*(z),r)),:);
            end%%%Z+1?? –顶琿だ秨
            tt(z+1) = (t(:,1) >= x)' * t(:,2) /length(t) ;%%%%P33,程F场だぃ1-
        end
        empcdf = 1-mean(tt);%%%%1-
   end

end