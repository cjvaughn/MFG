clearvars
tic
jobstring='january_31_Ex29'

%January 27th: added rho_0 to HJB (can set to 0 to not have it)

%January 26th: putting all versions of cost in this folder, with booleans
%cost_unsq_norm, cost_unsq_norm_2, cost_unsq_norm_3, and cost_Nourian

%January 25th: made different version of the cost in other folders:
%MFG_Stable_Explicit_Nourian, ...

%December 16th: changing the scheme for the boundaries.

%December 7th: continuing to correct the use of alpha

%December 6th: rewriting Kolmogorov in terms of alpha only, as in
%Flocking_Derivation

%November 28th: bounding alpha

%November 16th: calculate fftn(v_pad) once before the while loop, replaced
%loop at 383, replaced inside while loop with mu_diff_max

%October 28th: added extra cost term for flight path, and lambda2 to weight
%it

%October 19th: changed stability condition to allow for delta_x!=delta_y
% changed box_r to depend on delta_x

%October 18th: checked for typos and that delta_x and delta_y can be
%different

%October 17th: added new initial density: initial_skew, initial_skew2

%October 14th: added new initial densities: initial_5_points and
%initial_5_points_xandv

%October 4th: adding parameters c and lambda to weight the two costs:
% cost=c*[(1-lambda) 1/2 alpha^2+lambda*integral dmu]
% c=2 and lamda=1/2 corresponds to what we had before
% keeping how we set the stability condition the same, but checking        
% V_y/(c*(1-lambda)) against it instead of just V_y

%September 20th: adding initial_2_points boolean to initialize with 2 point
%masses located box_r away from the origin in x (and v=0)

%September 19th: this is a copy of main_vectorized_NewBCs_Extend_x.m
% removing warm_start and Part A
% removing text file for writing (fprint(fileID...)
% removing throw_away
% removing more_room boolean
% added calculation to extend_x further for initial_box or initial_2_boxes
% initial_2_boxes puts 2 boxes in quadrants 2 and 4 (with overlap at the
% origin!)

cost_unsq_norm=false
cost_unsq_norm_2=false
cost_unsq_norm_3=false
cost_Nourian=true

threshold=10^(-5) %for checking if sum is 1, and alpha<alpha_max, V>0
normalize=false
bound_alpha=true
c=3
lambda=1/3 %ToDo
lambda2=0

initial_mu_guess=false

if initial_mu_guess
    initial_mu_file='august_29_Ex1_final_mu.mat'
    mu_guess_struct=load(initial_mu_file);
    mu_guess=mu_guess_struct.mu;
end

initial_box=false %%birds are initially distributed with x=0 and v ranging from -box_r*delta_y:box_r*delta_y
initial_box_x=false %birds are initially distributed with v=0 and x ranging from -box_r*delta_x:box_r*delta_x
initial_2_boxes=false %birds are in 1 of two boxes in quadrants 2 and 4
initial_2_points=false %birds are either at (-1,0) or (1,0) (when box_r=100)
initial_5_points=false %birds are either at (-1,0) (-0.5,0) (0,0) (.5,0) or (1,0) (when box_r=100)
initial_5_points_xandv=false %birds are either at (-1,.1) (-1,-.1) (0,0) (1,.1) or (1,-.1) (when box_r=100)
initial_skew=false %birds are either at (-1,-.1) (0,0) or (1,.1) (when box_r=100)
initial_skew2=false %birds are either at (-1,.1) (0,0) or (1,-.1) (when box_r=100)
initial_skew3=false %birds are either at (0,-y) or (0,y) (when box_r_y=y/delta_y)
initial_skew4=true

num_iterations=20 %TODO

alpha_max=0.1         %previously more_room_factor*sqrt(2)*y_max
alpha_min=-alpha_max

sigma=0.1
rho_0=0; %ToDo, make 0
beta=100

num_time_points=1201 %todo
num_y=121 %needs to be odd

delta_x=0.5
delta_y=0.05

box_r=round(1.0/delta_x)
box_r_y=round(0.6/delta_y)

y_min=-(num_y-1)/2*delta_y;
y_max=(num_y-1)/2*delta_y;


delta_t=0.5*1/(sigma^2/(delta_y)^2+(alpha_max/(delta_y)+y_max/(delta_x))); %removed extra 1/2
delta_t
T=(num_time_points-1)*delta_t;
T

x_max_desired=T*y_max;
num_x_one_side=ceil(x_max_desired/delta_x)+3;
%the plus three gives extra room to make sure the birds can't reach the boundary
if initial_box
    num_x_one_side=num_x_one_side;
elseif initial_2_boxes
    num_x_one_side=num_x_one_side+box_r;
elseif initial_box_x
    num_x_one_side=num_x_one_side+box_r;
elseif initial_2_points
    num_x_one_side=num_x_one_side+box_r;
elseif initial_5_points
    num_x_one_side=num_x_one_side+box_r;
elseif initial_5_points_xandv
    num_x_one_side=num_x_one_side+box_r;
elseif initial_skew
    num_x_one_side=num_x_one_side+box_r;
elseif initial_skew2
    num_x_one_side=num_x_one_side+box_r;
elseif initial_skew3
    num_x_one_side=num_x_one_side;
elseif initial_skew4
    num_x_one_side=num_x_one_side+box_r;
end
num_x=num_x_one_side*2+1;
num_x
x_max=(num_x-1)/2*delta_x;
x_min=-x_max;
x_max
t_grid=linspace(0,T,num_time_points);
x_grid=linspace(x_min,x_max,num_x);
y_grid=linspace(y_min,y_max,num_y);

v_init_mean=0; %since v_init_mean=0, could just ignore chi term in V_y since it will be 0
if v_init_mean ~=0
    'v_init_mean is not 0 so uncomment chi!!!!'
end
y_j=repmat(y_grid,num_x,1);
x_i=repmat(x_grid',1,num_y);



%Initializing iterating:
K=1;
mu=zeros(num_time_points,num_x,num_y);
if initial_mu_guess
    mu_3D=reshape(mu_guess,[1,size(mu_guess)]);
    mu(:,:,:)=repmat(mu_3D,num_time_points,1,1);
elseif initial_box
    mu(:,ceil(num_x/2),ceil(num_y/2)-box_r:ceil(num_y/2)+box_r)=1/((2*box_r+1)*delta_x*delta_y);
elseif initial_box_x
    mu(:,ceil(num_x/2)-box_r:ceil(num_x/2)+box_r,ceil(num_y/2))=1/((2*box_r+1)*delta_x*delta_y);
elseif initial_2_boxes
    mu(:,ceil(num_x/2)-box_r:ceil(num_x/2),ceil(num_y/2):ceil(num_y/2)+box_r)=1/((2*(box_r+1)^2-1)*delta_x*delta_y);
    mu(:,ceil(num_x/2):ceil(num_x/2)+box_r,ceil(num_y/2)-box_r:ceil(num_y/2))=1/((2*(box_r+1)^2-1)*delta_x*delta_y);
elseif initial_2_points
    mu(:,ceil(num_x/2)-box_r,ceil(num_y/2))=1/(2*delta_x*delta_y);
    mu(:,ceil(num_x/2)+box_r,ceil(num_y/2))=1/(2*delta_x*delta_y);
elseif initial_5_points
    mu(:,ceil(num_x/2)-box_r,ceil(num_y/2))=1/(5*delta_x*delta_y);
    mu(:,ceil(num_x/2)+box_r,ceil(num_y/2))=1/(5*delta_x*delta_y);
    mu(:,ceil(num_x/2)-floor(box_r/2),ceil(num_y/2))=1/(5*delta_x*delta_y);
    mu(:,ceil(num_x/2)+floor(box_r/2),ceil(num_y/2))=1/(5*delta_x*delta_y);
    mu(:,ceil(num_x/2),ceil(num_y/2))=1/(5*delta_x*delta_y);
elseif initial_5_points_xandv
    mu(:,ceil(num_x/2)-box_r,ceil(num_y/2)-box_r_y)=1/(5*delta_x*delta_y);
    mu(:,ceil(num_x/2)+box_r,ceil(num_y/2)-box_r_y)=1/(5*delta_x*delta_y);
    mu(:,ceil(num_x/2)-box_r,ceil(num_y/2)+box_r_y)=1/(5*delta_x*delta_y);
    mu(:,ceil(num_x/2)+box_r,ceil(num_y/2)+box_r_y)=1/(5*delta_x*delta_y);
    mu(:,ceil(num_x/2),ceil(num_y/2))=1/(5*delta_x*delta_y);
elseif initial_skew
    mu(:,ceil(num_x/2)-box_r,ceil(num_y/2)-box_r_y)=1/(3*delta_x*delta_y);
    mu(:,ceil(num_x/2)+box_r,ceil(num_y/2)+box_r_y)=1/(3*delta_x*delta_y);
    mu(:,ceil(num_x/2),ceil(num_y/2))=1/(3*delta_x*delta_y);
elseif initial_skew2
    mu(:,ceil(num_x/2)-box_r,ceil(num_y/2)+box_r_y)=1/(3*delta_x*delta_y);
    mu(:,ceil(num_x/2)+box_r,ceil(num_y/2)-box_r_y)=1/(3*delta_x*delta_y);
    mu(:,ceil(num_x/2),ceil(num_y/2))=1/(3*delta_x*delta_y);
elseif initial_skew3
    mu(:,ceil(num_x/2),ceil(num_y/2)-box_r_y)=1/(2*delta_x*delta_y);
    mu(:,ceil(num_x/2),ceil(num_y/2)+box_r_y)=1/(2*delta_x*delta_y);
elseif initial_skew4
    mu(:,ceil(num_x/2)-box_r,ceil(num_y/2)+box_r_y)=1/(2*delta_x*delta_y);
    mu(:,ceil(num_x/2)+box_r,ceil(num_y/2)-box_r_y)=1/(2*delta_x*delta_y);
else
    mu(:,ceil(num_x/2),ceil(num_y/2))=1/(delta_x*delta_y); %puts everything at the origin.
end
value=sum(sum(mu(1,:,:,:,:)))*delta_x*delta_y;
if value~=1
    'Oh no!!!!! Initial mu does not sum to 1'
    value
end

if cost_unsq_norm
    v=get_v_matrix_unsq_norm(num_x,num_y,delta_x,delta_y,beta);
elseif cost_unsq_norm_2
    v=get_v_matrix_unsq_norm_2(num_x,num_y,delta_x,delta_y,beta);
elseif cost_unsq_norm_3
    v=get_v_matrix_unsq_norm_3(num_x,num_y,delta_x,delta_y,beta);
elseif cost_Nourian
    v=get_v_matrix_Nourian(num_x,num_y,delta_x,delta_y,beta);
else
    v=get_v_matrix(num_x,num_y,delta_x,delta_y,beta);
end
vpad=zeros(3*num_x-2,3*num_y-2);
vpad(1:2*num_x-1,1:2*num_y-1)=v;
fftn_vpad=fftn(vpad);

mu_diff_max=1;
old_alpha=zeros(num_time_points,num_x,num_y);
cost_alpha=zeros(num_time_points-1,1);
cost_integral=zeros(num_time_points-1,1);
%Iteration:
while(mu_diff_max>0 && K<num_iterations+1) %K<2 means 1 iteration
old_mu=mu;
'Part B: HJB'
%% Given mu, solve for V

max_alpha=0;
V=zeros(num_time_points,num_x,num_y);
V(num_time_points,:,:)=0;

max_diff_alpha=0;
for counter=1:num_time_points-1
    n=num_time_points-counter;
    t_n=t_grid(n);
    u(:,:)=mu(n,:,:)*delta_x*delta_y; %u is (num_x,num_y), v is (2*num_x-1,2*num_y-1), output=(3*num_x-2,3*num_y-2)
    
    upad=zeros(3*num_x-2,3*num_y-2);
    upad(1:num_x,1:num_y)=u;

    F2=ifftn(fftn(upad).*fftn_vpad);
    F=F2(num_x:2*num_x-1,num_y:2*num_y-1);
    if cost_Nourian
        F=F.^2;
    end

    mu_curr=squeeze(mu(n,:,:));
    V_curr=squeeze(V(n+1,:,:));
    
    V_yy(:,:)=shift(V_curr,1,2)-2*V_curr+shift(V_curr,-1,2); %centered
    %New Scheme BC:
    V_yy(:,1)=2*V_curr(:,2)-2*V_curr(:,1);
    V_yy(:,num_y)=2*V_curr(:,num_y-1)-2*V_curr(:,num_y);
    V_x=zeros(num_x,num_y);
    V_x(:,1:ceil(num_y/2)-1)=V_curr(:,1:ceil(num_y/2)-1)-shift(V_curr(:,1:ceil(num_y/2)-1),-1,1); %if v_j<0, one sided backwards
    V_x(:,ceil(num_y/2)+1:num_y)=shift(V_curr(:,ceil(num_y/2)+1:num_y),1,1)-V_curr(:,ceil(num_y/2)+1:num_y); %if v_j>0, one sided forward
    %New Scheme BC:
    V_x(1,:)=0;
    V_x(num_x,:)=0;
    
    %alpha_approx(:,:)=-(shift(V_curr,1,2)-shift(V_curr,-1,2)); %centered, Could have done other things here. Such as using last iterate for V, or Beta=0 solution.
    right(:,:)=shift(V_curr,1,2)-V_curr;
    left(:,:)=V_curr-shift(V_curr,-1,2);
    V_y=zeros(num_x,num_y);
    V_y(left<0 & right<0)=right(left<0 & right<0);
    V_y(left>0 & right>0)=left(left>0 & right>0);
    %New Scheme BC:
    V_y(:,1)=0;
    V_y(:,num_y)=0;
    
    alpha=-V_y/(c*(1-lambda-lambda2)*delta_y);
    if bound_alpha
        alpha=min(alpha,alpha_max);
        alpha=max(alpha,alpha_min);
    end
    
    if lambda==1
        if V_y>0
            alpha=alpha_min;
        elseif V_y<0
            alpha=alpha_max;
        else
            alpha=0;
        end
    end
    
    %theta_t=0.1*delta_t*n;
    theta_t=0.1*sin(delta_t*n);
    
%     cost=(x_i-theta_t).^2;
%     'here'
%     cost(ceil(num_x/2),ceil(num_y/2))
%     F(ceil(num_x/2),ceil(num_y/2))
    
    %%%% Using convolution
    %Linear, using the previous estimate for V to approximate gradient V
    if K>1
        V(n,:,:)=V_curr+delta_t*(sigma^2/2*V_yy/(delta_y)^2+y_j.*V_x/delta_x+alpha.*V_y/delta_y+1/2*c*(1-lambda-lambda2)*(squeeze(old_alpha(n,:,:)).*alpha)+c*lambda*F(:,:)+c*lambda2*100*(y_j-theta_t).^2-rho_0);
    else
        V(n,:,:)=V_curr+delta_t*(sigma^2/2*V_yy/(delta_y)^2+y_j.*V_x/delta_x+alpha.*V_y/delta_y+1/2*c*(1-lambda-lambda2)*(alpha.*alpha)+c*lambda*F(:,:)+c*lambda2*100*(y_j-theta_t).^2-rho_0);
    end
    
    if K>1
        alpha_diff=abs(squeeze(old_alpha(n,:,:))-alpha);
        alpha_diff_real=zeros(num_x,num_y);
        alpha_diff_real(mu_curr>0)=alpha_diff(mu_curr>0);
        max_curr=max(max(alpha_diff_real));
        max_diff_alpha=max(max_curr,max_diff_alpha);
        cost_alpha(n)=sum(sum(1/2*c*(1-lambda-lambda2)*(squeeze(old_alpha(n,:,:)).*alpha)));
        cost_integral(n)=sum(sum(c*lambda*F(:,:)));
    end

    max_alpha=max(max_alpha,max(max(abs(alpha))));
    old_alpha(n,:,:)=alpha;
end

%Checking if the solution is valid:
if K>1
    'Difference in alpha'
    max_diff_alpha
end
'Part B Max Alpha'
max_alpha
if max_alpha>alpha_max+threshold
    'Oh no!!!!! Part B max_alpha>alpha_max'
    max_alpha
end
'Minimum of V'
value=min(min(min(V(:,:,:))))
if value<-threshold
    'Oh no!!!!! Negatives in V'
    value
end
'Maximum of V'
value=max(max(max(V(:,:,:))))


'Part C: Kolmogorov'
%% Given V, solve for mu
max_alpha=0;
mu=zeros(num_time_points,num_x,num_y);

if initial_mu_guess
    mu(1,:,:)=mu_guess;
elseif initial_box
    mu(1,ceil(num_x/2),ceil(num_y/2)-box_r:ceil(num_y/2)+box_r)=1/((2*box_r+1)*delta_x*delta_y);
elseif initial_box_x
    mu(1,ceil(num_x/2)-box_r:ceil(num_x/2)+box_r,ceil(num_y/2))=1/((2*box_r+1)*delta_x*delta_y);
elseif initial_2_boxes
    mu(1,ceil(num_x/2)-box_r:ceil(num_x/2),ceil(num_y/2):ceil(num_y/2)+box_r)=1/((2*(box_r+1)^2-1)*delta_x*delta_y);
    mu(1,ceil(num_x/2):ceil(num_x/2)+box_r,ceil(num_y/2)-box_r:ceil(num_y/2))=1/((2*(box_r+1)^2-1)*delta_x*delta_y);
elseif initial_2_points
    mu(1,ceil(num_x/2)-box_r,ceil(num_y/2))=1/(2*delta_x*delta_y);
    mu(1,ceil(num_x/2)+box_r,ceil(num_y/2))=1/(2*delta_x*delta_y);
elseif initial_5_points
    mu(1,ceil(num_x/2)-box_r,ceil(num_y/2))=1/(5*delta_x*delta_y);
    mu(1,ceil(num_x/2)+box_r,ceil(num_y/2))=1/(5*delta_x*delta_y);
    mu(1,ceil(num_x/2)-floor(box_r/2),ceil(num_y/2))=1/(5*delta_x*delta_y);
    mu(1,ceil(num_x/2)+floor(box_r/2),ceil(num_y/2))=1/(5*delta_x*delta_y);
    mu(1,ceil(num_x/2),ceil(num_y/2))=1/(5*delta_x*delta_y);
elseif initial_5_points_xandv
    mu(1,ceil(num_x/2)-box_r,ceil(num_y/2)-box_r_y)=1/(5*delta_x*delta_y);
    mu(1,ceil(num_x/2)+box_r,ceil(num_y/2)-box_r_y)=1/(5*delta_x*delta_y);
    mu(1,ceil(num_x/2)-box_r,ceil(num_y/2)+box_r_y)=1/(5*delta_x*delta_y);
    mu(1,ceil(num_x/2)+box_r,ceil(num_y/2)+box_r_y)=1/(5*delta_x*delta_y);
    mu(1,ceil(num_x/2),ceil(num_y/2))=1/(5*delta_x*delta_y);
elseif initial_skew
    mu(1,ceil(num_x/2)-box_r,ceil(num_y/2)-box_r_y)=1/(3*delta_x*delta_y);
    mu(1,ceil(num_x/2)+box_r,ceil(num_y/2)+box_r_y)=1/(3*delta_x*delta_y);
    mu(1,ceil(num_x/2),ceil(num_y/2))=1/(3*delta_x*delta_y);
elseif initial_skew2
    mu(1,ceil(num_x/2)-box_r,ceil(num_y/2)+box_r_y)=1/(3*delta_x*delta_y);
    mu(1,ceil(num_x/2)+box_r,ceil(num_y/2)-box_r_y)=1/(3*delta_x*delta_y);
    mu(1,ceil(num_x/2),ceil(num_y/2))=1/(3*delta_x*delta_y);
elseif initial_skew3
    mu(1,ceil(num_x/2),ceil(num_y/2)-box_r_y)=1/(2*delta_x*delta_y);
    mu(1,ceil(num_x/2),ceil(num_y/2)+box_r_y)=1/(2*delta_x*delta_y);
elseif initial_skew4
    mu(1,ceil(num_x/2)-box_r,ceil(num_y/2)+box_r_y)=1/(2*delta_x*delta_y);
    mu(1,ceil(num_x/2)+box_r,ceil(num_y/2)-box_r_y)=1/(2*delta_x*delta_y);
else
    mu(1,ceil(num_x/2),ceil(num_y/2))=1/(delta_x*delta_y); %puts everything at the origin.
end

for n=1:num_time_points-1
    t_n=t_grid(n);
        
    mu_curr=squeeze(mu(n,:,:));
    V_curr=squeeze(V(n,:,:));
    
    mu_yy(:,:)=shift(mu_curr,1,2)-2*mu_curr+shift(mu_curr,-1,2); %centered
    %New Scheme BC:
    mu_yy(:,2)=mu_yy(:,2)+mu_curr(:,1);
    mu_yy(:,num_y-1)=mu_yy(:,num_y-1)+mu_curr(:,num_y);
    mu_yy(:,1)=mu_yy(:,1)-mu_curr(:,num_y);
    mu_yy(:,num_y)=mu_yy(:,num_y)-mu_curr(:,1);
    V_yy(:,:)=shift(V_curr,1,2)-2*V_curr+shift(V_curr,-1,2); %centered
    mu_x=zeros(num_x,num_y);
    mu_x(:,1:ceil(num_y/2)-1)=shift(mu_curr(:,1:ceil(num_y/2)-1),1,1)-mu_curr(:,1:ceil(num_y/2)-1); %if v_j<0, one sided forward
    mu_x(:,ceil(num_y/2)+1:num_y)=mu_curr(:,ceil(num_y/2)+1:num_y)-shift(mu_curr(:,ceil(num_y/2)+1:num_y),-1,1); %if v_j>0, one sided backwards (if v_j=0, 0)
    %New Scheme BC:
    mu_x(2,ceil(num_y/2)+1:num_y)=mu_curr(2,ceil(num_y/2)+1:num_y);
    mu_x(num_x-1,1:ceil(num_y/2)-1)=-mu_curr(num_x-1,1:ceil(num_y/2)-1);
    mu_x(1,:)=0;
    mu_x(1,1:ceil(num_y/2)-1)=mu_curr(2,1:ceil(num_y/2)-1);
    mu_x(num_x,:)=0;
    mu_x(num_x,ceil(num_y/2)+1:num_y)=-mu_curr(num_x-1,ceil(num_y/2)+1:num_y);
    
    left(:,:)=V_curr-shift(V_curr,-1,2);
    right(:,:)=shift(V_curr,1,2)-V_curr;
    
    V_y=zeros(num_x,num_y);
    V_y(left<0 & right<0)=right(left<0 & right<0);
    V_y(left>0 & right>0)=left(left>0 & right>0);
    V_y(:,1)=0;
    V_y(:,num_y)=0;
    
    alpha=-V_y/(c*(1-lambda-lambda2)*delta_y);
    if bound_alpha
        alpha=min(alpha,alpha_max);
        alpha=max(alpha,alpha_min);
    end
    
    if lambda==1
        if V_y>0
            alpha=alpha_min;
        elseif V_y<0
            alpha=alpha_max;
        else
            alpha=0;
        end
    end
    
    alpha_minus=zeros(num_x,num_y);
    alpha_plus=zeros(num_x,num_y);
    alpha_minus(alpha<0)=-alpha(alpha<0);
    alpha_plus(alpha>0)=alpha(alpha>0);
    
    alpha_term_a=shift(alpha_minus,1,2);
    alpha_term_b=shift(alpha_plus,-1,2);
    alpha_term_c=-(alpha_plus+alpha_minus);
    %New Scheme BC:
    alpha_term_b(:,2)=0;
    alpha_term_a(:,num_y-1)=0;
    alpha_term_b(:,1)=0;
    alpha_term_a(:,num_y)=0;
    alpha_term_c(:,1)=0;
    alpha_term_c(:,num_y)=0;
    
    mu(n+1,:,:)=mu_curr+delta_t*(sigma^2/2*mu_yy/(delta_y)^2-y_j.*mu_x/delta_x+alpha_term_a.*shift(mu_curr,1,2)/delta_y+alpha_term_b.*shift(mu_curr,-1,2)/delta_y+alpha_term_c.*mu_curr/delta_y);

    max_alpha=max(max_alpha,max(max(abs(alpha))));
end

%Checking if the solution is valid:
integral_values_2=zeros(num_time_points,1);
for n=1:num_time_points
    integral=sum(sum(mu(n,:,:)))*delta_x*delta_y;
    integral_values_2(n)=integral;
    if integral>1+threshold || integral<1-threshold
        'Oh no!!!!! Sum is not 1!!!!'
        integral
    end
end
if normalize
    for n=1:num_time_points
        mu(n,:,:)=mu(n,:,:)/integral_values_2(n);
    end
end
'Part C Max Alpha'
max_alpha
if max_alpha>alpha_max+threshold
    'Oh no!!!!! Part C max_alpha>alpha_max'
    max_alpha
end
%'Minimum Density in Part C'
value=min(min(min(mu(:,:,:))));
if value<-threshold
    'Oh no!!!!! Negatives in mu'
    value
end
mu_diff=abs(mu-old_mu);
mu_diff_max=max(max(max(mu_diff(:,:,:))))*delta_x*delta_y;
value=mu_diff_max;
K
'Difference in mu'
value
mu_diff_frac=mu_diff./abs(mu);
mu_diff_frac(abs(mu)<10^(-10) & abs(old_mu)<10^(-10))=-1;
value=max(max(max(mu_diff_frac(:,:,:))));
'Largest Fractional Difference'
value
K=K+1;
end %This ends the while loop

%Final Calculations
timer=toc
final_mu=squeeze(mu(num_time_points,:,:)).*delta_x*delta_y;
save(strcat(jobstring,'_final_mu.mat'),'final_mu')
initial_V=squeeze(V(1,:,:));
save(strcat(jobstring,'_initial_V.mat'),'initial_V')
%save(strcat(jobstring,'_V.mat'),'V','-v7.3')
save(strcat(jobstring,'_cost_alpha.mat'),'cost_alpha')
save(strcat(jobstring,'_cost_integral.mat'),'cost_integral')

output_freq=1000;
num_times=floor(num_time_points/output_freq)+1;
mu_short=zeros(num_x,num_y,num_times);
for i=1:num_times
    mu_short(:,:,i)=mu((i-1)*output_freq+1,:,:);
end
%save(strcat(jobstring,'_mu.mat'),'mu','-v7.3')
save(strcat(jobstring,'_mu_short.mat'),'mu_short','-v7.3')
save(strcat(jobstring,'_integral_values_2.mat'),'integral_values_2')
'Done'
