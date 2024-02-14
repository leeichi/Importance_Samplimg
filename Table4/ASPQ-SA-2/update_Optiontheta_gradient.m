function [theta_next] = update_Optiontheta_gradient(obser,p,stage,theta_initial,quantile_initial,tk,MU_Y,VAR_Y,MU_C,C,alpha,sample_path)
    
    step_size = tk / (C+ (stage)^(alpha));
    n = size(obser,2);
    theta_next = 0;
    theta_next = theta_initial - step_size * ...
        (-2*(1/n)* ( (obser(1,1:n,1) > quantile_initial) * ...
        ((normpdf(sample_path(1,:),MU_Y(1),VAR_Y(1)^(1/2))).^(2) ./ (normpdf(sample_path(1,:),MU_C(1),VAR_Y(1)^(1/2))).^(3) .*...
        (exp(-((MU_C(1)-sample_path(1,:)).^2) /(2*(VAR_Y^(1/2))^2)).*(sample_path(1,:)- MU_C(1)))./...
        ((2*pi)^(1/2)*(VAR_Y^(1/2))^3)  )'  )  );
          
%%%%thetaªº­pºâ
end

