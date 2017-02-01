%Calculate the variance

num_y=121 %needs to be odd

delta_y=0.05

y_min=-(num_y-1)/2*delta_y;
y_max=(num_y-1)/2*delta_y;
y_grid=linspace(y_min,y_max,num_y);

f(:,:)=final_mu;
integral_value=sum(sum(f));
f=f/integral_value;
f2=rot90(f);
mu_v=sum(f2,2)';

variance=0;
num_moments=100;
moments=zeros(num_moments,1);
for k=1:num_moments
for i=1:num_y
    y=y_grid(i);
    moments(k)=moments(k)+mu_v(i)*y^k;
end
end

moments;


plot(moments_jan_31(2,:))
title('Variance of marginal in v')
xlabel('beta')
ylabel('Variance')