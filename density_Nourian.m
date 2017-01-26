%Density from Nourian_Perturbation paper, beta=0, c=3, lambda=1/3, Nourian
%cost, steady state

normalize=true
sigma=0.1
s_2=sigma^2/2

num_y=41 %needs to be odd

delta_y=0.05

y_min=-(num_y-1)/2*delta_y;
y_max=(num_y-1)/2*delta_y;
y_grid=linspace(y_min,y_max,num_y);

mu=1/(2*pi*s_2)^(0.5)*exp(-(y_grid).^2/(2*s_2))*delta_y;

if normalize
    mu=mu/sum(mu);
end

plot(y_grid,mu)
title('mu(y1) stationary, beta=0')
xlabel('y1')
ylabel('mu(y1)')