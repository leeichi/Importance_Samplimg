function[empcdf] = empCDF(x,obs,s,r,L,n0)
    tt = [];
    for z = 0: s-1
        if(z == 0)
            t = obs((1:n0),:);
        else
            t = obs((n0+L*(z-1)+1:min(n0 +L*(z),r)),:);
        end%%%Z+1?? 每一階段分開
        tt(z+1) = (t(:,1) >= x)' * t(:,2) /length(t) ;%%%%P33,最上面F後面的部分，不包刮1-
    end
    empcdf = 1-mean(tt);%%%%1-加上去
end