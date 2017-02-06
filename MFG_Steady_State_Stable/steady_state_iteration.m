clearvars
tic
jobstring='test' %'october_25_steady_Ex2'

%WARNING: we have a lambda here which is not the weighting the cost lambda!

%August 28th: NewBC's!!! (From notes Carmona sent me)

%August 26th: Poisson_perturbed to follow Achdou_MFG_Numerical existence
%proof

%August 22nd: Solving for lambda and V all at once in HJB.

%July 15th: Solving infinite horizon/steady state problem using an
%iteration method

c=2
lambda_cost=1/2

initial_uniform=false
initial_center=true

num_iterations=10;

epsilon=0 %for making it a viscosity solution (noise in x too)
sigma=0.1
beta=0
ro=10

num_x=41 %101;
num_y=41
delta_x=0.5
delta_y=0.05

x_min=-(num_x-1)/2*delta_x;
x_max=(num_x-1)/2*delta_x;
y_min=-(num_y-1)/2*delta_y;
y_max=(num_y-1)/2*delta_y;

x_grid=linspace(x_min,x_max,num_x);
y_grid=linspace(y_min,y_max,num_y);

v=get_v_matrix(num_x,num_y,delta_x,delta_y,beta);
vpad=zeros(3*num_x-2,3*num_y-2);
vpad(1:2*num_x-1,1:2*num_y-1)=v;

%initial guesses
V=zeros(num_x,num_y);
old_V_y=zeros(num_x,num_y);
if initial_uniform
    mu=zeros(num_x,num_y);
    mu(:,:)=1/(num_x*num_y*delta_x*delta_y);
elseif initial_center
    mu=zeros(num_x,num_y);
    mu(ceil(num_x/2),ceil(num_y/2))=1/(delta_x*delta_y);
else
    initial_mu_file='october_24_Ex4_final_mu.mat'
    mu_guess_struct=load(initial_mu_file);
    mu=mu_guess_struct.final_mu;
end
lambda=1; %not needed (only needed to calculate the first difference between lambda and old_lambda

lambda_history=zeros(num_iterations,1);
lambda2_history=zeros(num_iterations,1);
lambda_max_history=zeros(num_iterations,1);
lambda_min_history=zeros(num_iterations,1);
r_history=zeros(num_iterations,1);
r2_history=zeros(num_iterations,1);
diff_history=zeros(num_x*num_y+2,num_iterations);
diff2_history=zeros(num_x*num_y,num_iterations);
diff3_history=zeros(num_x*num_y,num_iterations);
integral_history=zeros(num_iterations,1);
for k=1:num_iterations
old_V=V;
old_mu=mu;
old_lambda=lambda;

% 'Finding Lambda'
% %Use V and mu to estimate lambda
% [lambda,lambda_max,lambda_min]=find_lambda_avg(V,mu,vpad,num_x,num_y,delta_x,delta_y,sigma,y_grid);


'Solving HJB'
%Use mu to solve for V and lambda
[V,r,lambda,diff,M]=solve_HJB_extra_eq(mu,old_V_y,vpad,num_x,num_y,delta_x,delta_y,epsilon,sigma,y_grid,c,lambda_cost);

lambda2=eval_lambda(V,mu,vpad,num_x,num_y,delta_x,delta_y,sigma,y_grid) %TODO: eval_lambda is wrapping around!

%Calculate old_V_y from V
right(:,:)=shift(V,1,2)-V;
left(:,:)=V-shift(V,-1,2);
old_V_y=zeros(num_x,num_y);
old_V_y(left<0 & right<0)=right(left<0 & right<0);
old_V_y(left>0 & right>0)=left(left>0 & right>0);
%New Scheme BC:
old_V_y(:,1)=0;
old_V_y(:,num_y)=0;

old_V_y=old_V_y/delta_y;
  
%Calculate alpha from old_V_y
alpha=-old_V_y/(c*(1-lambda_cost));

'Solving Poisson'
%Uses perturbed equation from Achdou_MFG_Numerical paper
[mu,r2,diff2,M2,b,diff3]=solve_Poisson_perturbed(ro,old_mu,V,num_x,num_y,delta_x,delta_y,epsilon,sigma,y_grid,alpha);

k
abs(lambda-old_lambda)
max(max(abs(mu-old_mu)*delta_x*delta_y))
max(max(abs(V-old_V)))


lambda_history(k,1)=lambda;
lambda2_history(k,1)=lambda2;
% lambda_max_history(k,1)=lambda_max;
% lambda_min_history(k,1)=lambda_min;

r_history(k)=r;
r2_history(k)=r2;
diff_history(:,k)=diff;
diff2_history(:,k)=diff2;
diff3_history(:,k)=diff3;
'max diff histories:'
max(diff)
max(diff2)
max(diff3)
integral=sum(sum(mu))*delta_x*delta_y;
'sum of mu'
integral
integral_history(k,1)=integral;
'max diff mu'
max(max(abs((mu-old_mu)*delta_x*delta_y)))
end
toc
save(strcat(jobstring,'_integral_history.mat'),'integral_history')
save(strcat(jobstring,'_diff3_history.mat'),'diff3_history') %TODO
save(strcat(jobstring,'_V.mat'),'V','-v7.3')
save(strcat(jobstring,'_mu.mat'),'mu','-v7.3')
save(strcat(jobstring,'_lambda.mat'),'lambda')











