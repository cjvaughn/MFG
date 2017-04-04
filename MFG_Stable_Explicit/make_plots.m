num_time_points=15001;
num_x=87; %127; %127; %495; %495; %573; %399; %1567; %787;
num_y=61; %21; %41;

delta_x=0.05;
delta_y=0.05;

x_max=(num_x-1)/2*delta_x;
y_max=(num_y-1)/2*delta_y;


%For 1d plots in x and v
%f(:,:)=squeeze(mu(n,left:right,left:right))*delta_x^2;
%For 2d plots in x1 and x2
% mu_curr=squeeze(mu(num_time_points,:,:,:,:));
% mu_12=sum(sum(mu_curr,3),4).*delta_x^4;
% mu_13=squeeze(sum(sum(mu_curr,2),4)).*delta_x^4;
% mu_24=squeeze(sum(sum(mu_curr,1),3)).*delta_x^4;
% mu_34=squeeze(sum(sum(mu_curr,1),2)).*delta_x^4;
%f(:,:)=final_mu;

mu_plot=true;
mu_short_plot=false; %use make_plots_movie instead
V_plot=false;

if mu_plot
    f(:,:)=final_mu;
    integral_value=sum(sum(f));
    f=f/integral_value;
elseif mu_short_plot
    f(:,:)=mu_short(:,:,45);
    integral_value=sum(sum(f));
    f=f/integral_value; 
elseif V_plot
    f(:,:)=initial_V;
else
    n=10000
    V_curr=squeeze(V(n,:,:));
    left(:,:)=V_curr-shift(V_curr,-1,2);
    right(:,:)=shift(V_curr,1,2)-V_curr;
    V_y=zeros(num_x,num_y);
    V_y(left<0 & right<0)=left(left<0 & right<0);
    V_y(left>0 & right>0)=right(left>0 & right>0);
    for i=1:num_x
        for j=1:num_y
            if i>(187-n+1) && i<(187+n-1) && j>(26-n+1) && j<(26+n-1)
                dummy=0;
            elseif i>(387-n+1) && i<(387+n-1) && j>(26-n+1) && j<(26+n-1)
                dummy=0;
            else
                V_y(i,j)=0;
            end
        end
    end
    f(:,:)=V_y;
end

discard=0; %390; %300+390;
x_left=discard+1;
x_right=num_x-discard;
f=f(x_left:x_right,:);
a=-x_max+(x_left-1)*delta_x;
b=-x_max+(x_right-1)*delta_x;

f2=rot90(f);
colormap(flipud(gray))
x=[a b];
y=[-y_max y_max];
imagesc(x,y,f2);
colorbar;
xlabel('x1')
ylabel('v1')
if mu_plot
    title('mu(T,x,v), params=[1,1.1,1.1,1], c=2, lambda=0.5')
elseif mu_short_plot
    title('mu(t1,x,v), beta=1.5, c=2, lambda=0.5')
elseif V_plot
    title('V(0,x,v), beta=1.5')
else
    title('grad_v V(10000,x,v), beta=1.5')
end

% % Marginals in Velocity
% v1_values=sum(f2,2)';
% v2_values=sum(f2,1);
% y=linspace(-y_max,y_max,num_y);
% plot(y,v1_values)
% % bar(y,v1_values)
% title('mu(T,y1), beta=0')
% xlabel('y1')
% ylabel('mu(T,y1)')


%% Attic

% %Marginal in Velocity
% for i=1:num_x_points
% for j=1:num_y_points
%     mu_marginal(j)=mu_marginal(j)+mu_prime(n,i,j);
% end
% end

%mu_t(:,:)=mu_t(:,:)/sum(sum(mu_t(:,:)));

% x_min=-5;
% x_max=5;
% y_min=-5;
% y_max=5;
% 
% delta_x=(x_max-x_min)/(num_x_points-1);
% delta_y=(y_max-y_min)/(num_y_points-1);
% 
% x=linspace(-5,5,num_x_points);
% mu_marginal=mu_marginal*delta_x*delta_y;
% bar(x',mu_marginal)
% title('mu_0(T,x)')
% xlabel('x')
% ylabel('mu(T,x)')