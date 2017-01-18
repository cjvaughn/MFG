clearvars
tic
jobstring='january_18_Ex1_2d'

%December 16th: ToDo: make this code new stable scheme like 1d code

%December 8th: making corresponding changes from 1d code:
    %December 7th: continuing to correct the use of alpha
    %December 6th: rewriting Kolmogorov in terms of alpha only, as in
    %Flocking_Derivation

%November 28th: bounding alpha

%November 16th: making some changes (see main.m)

%October 19th: copied this from main.m. going to make 2d
% this time: (time,x1,y1,x2,y2)

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
% removing more_room
% added calculation to extend_x further for initial_box or initial_2_boxes
% initial_2_boxes puts 2 boxes in quadrants 2 and 4 (with overlap at the
% origin!)

threshold=10^(-5) %for checking if sum is 1, and alpha<alpha_max, V>0
normalize=false
bound_alpha=true
c=2
lambda=0.5

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

num_iterations=10 %ToDo

%more_room_factor=15
alpha_max=0.1 %todo, 1
alpha_min=-alpha_max

sigma=0.1
beta=1.5

num_time_points=1501
num_y=21 %needs to be odd

delta_x=0.5
delta_y=0.025

box_r=round(0.25/delta_x)
box_r_y=round(0.2/delta_y)

y_min=-(num_y-1)/2*delta_y;
y_max=(num_y-1)/2*delta_y;



delta_t=0.5*1/(2*sigma^2/(delta_y)^2+2*(alpha_max/(delta_y)+y_max/(delta_x))); %removed extra 1/2
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

y_j1=repmat(reshape(y_grid,1,num_y,1,1),num_x,1,num_x,num_y);
y_j2=repmat(reshape(y_grid,1,1,1,num_y),num_x,num_y,num_x,1);

%Initializing iterating:
K=1;
mu=zeros(num_time_points,num_x,num_y,num_x,num_y);
if initial_mu_guess
    mu_5D=reshape(mu_guess,[1,size(mu_guess)]);
    mu(:,:,:,:,:)=repmat(mu_3D,num_time_points,1,1,1,1);
elseif initial_box
    mu(:,ceil(num_x/2),ceil(num_y/2)-box_r:ceil(num_y/2)+box_r,ceil(num_x/2),ceil(num_y/2)-box_r:ceil(num_y/2)+box_r)=1/((2*box_r+1)^2*(delta_x*delta_y)^2);
elseif initial_box_x
    mu(:,ceil(num_x/2)-box_r:ceil(num_x/2)+box_r,ceil(num_y/2),ceil(num_x/2)-box_r:ceil(num_x/2)+box_r,ceil(num_y/2))=1/((2*box_r+1)^2*(delta_x*delta_y)^2);
elseif initial_2_boxes
    mu(:,ceil(num_x/2)-box_r:ceil(num_x/2),ceil(num_y/2):ceil(num_y/2)+box_r,ceil(num_x/2)-box_r:ceil(num_x/2),ceil(num_y/2):ceil(num_y/2)+box_r)=1/((2*(box_r+1)^2-1)^2*(delta_x*delta_y)^2);
    mu(:,ceil(num_x/2):ceil(num_x/2)+box_r,ceil(num_y/2)-box_r:ceil(num_y/2),ceil(num_x/2):ceil(num_x/2)+box_r,ceil(num_y/2)-box_r:ceil(num_y/2))=1/((2*(box_r+1)^2-1)^2*(delta_x*delta_y)^2);
elseif initial_2_points
    mu(:,ceil(num_x/2)-box_r,ceil(num_y/2),ceil(num_x/2)-box_r,ceil(num_y/2))=1/(2*(delta_x*delta_y)^2);
    mu(:,ceil(num_x/2)+box_r,ceil(num_y/2),ceil(num_x/2)+box_r,ceil(num_y/2))=1/(2*(delta_x*delta_y)^2);
elseif initial_5_points
    mu(:,ceil(num_x/2)-box_r,ceil(num_y/2),ceil(num_x/2)-box_r,ceil(num_y/2))=1/(5*(delta_x*delta_y)^2);
    mu(:,ceil(num_x/2)+box_r,ceil(num_y/2),ceil(num_x/2)+box_r,ceil(num_y/2))=1/(5*(delta_x*delta_y)^2);
    mu(:,ceil(num_x/2)-floor(box_r/2),ceil(num_y/2),ceil(num_x/2)-floor(box_r/2),ceil(num_y/2))=1/(5*(delta_x*delta_y)^2);
    mu(:,ceil(num_x/2)+floor(box_r/2),ceil(num_y/2),ceil(num_x/2)+floor(box_r/2),ceil(num_y/2))=1/(5*(delta_x*delta_y)^2);
    mu(:,ceil(num_x/2),ceil(num_y/2),ceil(num_x/2),ceil(num_y/2))=1/(5*(delta_x*delta_y)^2);
elseif initial_5_points_xandv
    mu(:,ceil(num_x/2)-box_r,ceil(num_y/2)-box_r_y,ceil(num_x/2)-box_r,ceil(num_y/2)-box_r_y)=1/(5*(delta_x*delta_y)^2);
    mu(:,ceil(num_x/2)+box_r,ceil(num_y/2)-box_r_y,ceil(num_x/2)+box_r,ceil(num_y/2)-box_r_y)=1/(5*(delta_x*delta_y)^2);
    mu(:,ceil(num_x/2)-box_r,ceil(num_y/2)+box_r_y,ceil(num_x/2)-box_r,ceil(num_y/2)+box_r_y)=1/(5*(delta_x*delta_y)^2);
    mu(:,ceil(num_x/2)+box_r,ceil(num_y/2)+box_r_y,ceil(num_x/2)+box_r,ceil(num_y/2)+box_r_y)=1/(5*(delta_x*delta_y)^2);
    mu(:,ceil(num_x/2),ceil(num_y/2),ceil(num_x/2),ceil(num_y/2))=1/(5*(delta_x*delta_y)^2);
elseif initial_skew
    mu(:,ceil(num_x/2)-box_r,ceil(num_y/2)-box_r_y,ceil(num_x/2)-box_r,ceil(num_y/2)-box_r_y)=1/(3*(delta_x*delta_y)^2);
    mu(:,ceil(num_x/2)+box_r,ceil(num_y/2)+box_r_y,ceil(num_x/2)+box_r,ceil(num_y/2)+box_r_y)=1/(3*(delta_x*delta_y)^2);
    mu(:,ceil(num_x/2),ceil(num_y/2),ceil(num_x/2),ceil(num_y/2))=1/(3*(delta_x*delta_y)^2);
elseif initial_skew2
    mu(:,ceil(num_x/2)-box_r,ceil(num_y/2)+box_r_y,ceil(num_x/2)-box_r,ceil(num_y/2)+box_r_y)=1/(3*(delta_x*delta_y)^2);
    mu(:,ceil(num_x/2)+box_r,ceil(num_y/2)-box_r_y,ceil(num_x/2)+box_r,ceil(num_y/2)-box_r_y)=1/(3*(delta_x*delta_y)^2);
    mu(:,ceil(num_x/2),ceil(num_y/2),ceil(num_x/2),ceil(num_y/2))=1/(3*(delta_x*delta_y)^2);
else
    mu(:,ceil(num_x/2),ceil(num_y/2),ceil(num_x/2),ceil(num_y/2))=1/(delta_x*delta_y)^2; %puts everything at the origin.
end
value=sum(sum(sum(sum(mu(1,:,:,:,:)))))*(delta_x*delta_y)^2;
if value~=1
    'Oh no!!!!! Initial mu does not sum to 1'
    value
end
    
% old_mu=mu-ones(num_time_points,num_x,num_y,num_x,num_y);
v=get_v_matrix_2d(num_x,num_y,delta_x,delta_y,beta);
vpad=zeros(3*num_x-2,3*num_y-2,3*num_x-2,3*num_y-2);
vpad(1:2*num_x-1,1:2*num_y-1,1:2*num_x-1,1:2*num_y-1)=v;
fftn_vpad=fftn(vpad);

mu_diff_max=1;
old_alpha=zeros(num_time_points,num_x,num_y);
%Iteration:
while(mu_diff_max>0 && K<num_iterations+1) %K<2 means 1 iteration
old_mu=mu;
'Part B: HJB'
%% Given mu, solve for V

max_alpha=0;
V=zeros(num_time_points,num_x,num_y,num_x,num_y);
V(num_time_points,:,:,:,:)=0;

max_diff_alpha_1=0;
max_diff_alpha_2=0;
for counter=1:num_time_points-1
    n=num_time_points-counter;
    t_n=t_grid(n);
    u(:,:,:,:)=mu(n,:,:,:,:)*(delta_x*delta_y)^2; %u is (num_x,num_y), v is (2*num_x-1,2*num_y-1), output=(3*num_x-2,3*num_y-2)
    
    upad=zeros(3*num_x-2,3*num_y-2,3*num_x-2,3*num_y-2);
    upad(1:num_x,1:num_y,1:num_x,1:num_y)=u;

    F2=ifftn(fftn(upad).*fftn_vpad);
    F=F2(num_x:2*num_x-1,num_y:2*num_y-1,num_x:2*num_x-1,num_y:2*num_y-1);

    mu_curr=squeeze(mu(n,:,:,:,:));
    V_curr=squeeze(V(n+1,:,:,:,:));
    
    V_yy1(:,:,:,:)=shift(V_curr,1,2)-2*V_curr+shift(V_curr,-1,2); %centered
    %New Scheme BC:
    V_yy1(:,1,:,:)=2*V_curr(:,2,:,:)-2*V_curr(:,1,:,:);
    V_yy1(:,num_y,:,:)=2*V_curr(:,num_y-1,:,:)-2*V_curr(:,num_y,:,:);
    V_x1=zeros(num_x,num_y,num_x,num_y);
    V_x1(:,1:ceil(num_y/2)-1,:,:)=V_curr(:,1:ceil(num_y/2)-1,:,:)-shift(V_curr(:,1:ceil(num_y/2)-1,:,:),-1,1); %if v_j<0, one sided backwards
    V_x1(:,ceil(num_y/2)+1:num_y,:,:)=shift(V_curr(:,ceil(num_y/2)+1:num_y,:,:),1,1)-V_curr(:,ceil(num_y/2)+1:num_y,:,:); %if v_j>0, one sided forward
    %New Scheme BC:
    V_x1(1,:,:,:)=0;
    V_x1(num_x,:,:,:)=0;
    
    right1=shift(V_curr,1,2)-V_curr;
    left1(:,:,:,:)=V_curr-shift(V_curr,-1,2);
    V_y1=zeros(num_x,num_y,num_x,num_y);
    V_y1(left1<0 & right1<0)=right1(left1<0 & right1<0);
    V_y1(left1>0 & right1>0)=left1(left1>0 & right1>0);
    %New Scheme BC:
    V_y1(:,1,:,:)=0;
    V_y1(:,num_y,:,:)=0;
    
    V_yy2(:,:,:,:)=shift(V_curr,1,4)-2*V_curr+shift(V_curr,-1,4); %centered
    %New Scheme BC:
    V_yy2(:,:,:,1)=2*V_curr(:,:,:,2)-2*V_curr(:,:,:,1);
    V_yy2(:,:,:,num_y)=2*V_curr(:,:,:,num_y-1)-2*V_curr(:,:,:,num_y);
    V_x2=zeros(num_x,num_y,num_x,num_y);
    V_x2(:,:,:,1:ceil(num_y/2)-1)=V_curr(:,:,:,1:ceil(num_y/2)-1)-shift(V_curr(:,:,:,1:ceil(num_y/2)-1),-1,3); %if v_j<0, one sided backwards
    V_x2(:,:,:,ceil(num_y/2)+1:num_y)=shift(V_curr(:,:,:,ceil(num_y/2)+1:num_y),1,3)-V_curr(:,:,:,ceil(num_y/2)+1:num_y); %if v_j>0, one sided forward
    %New Scheme BC:
    V_x2(:,:,1,:)=0;
    V_x2(:,:,num_x,:)=0;
    
    right2(:,:,:,:)=shift(V_curr,1,4)-V_curr;
    left2(:,:,:,:)=V_curr-shift(V_curr,-1,4);
    V_y2=zeros(num_x,num_y,num_x,num_y);
    V_y2(left2<0 & right2<0)=right2(left2<0 & right2<0);
    V_y2(left2>0 & right2>0)=left2(left2>0 & right2>0);
    %New Scheme BC:
    V_y2(:,:,:,1)=0;
    V_y2(:,:,:,num_y)=0;
    
    alpha_1=-V_y1/(c*(1-lambda)*delta_y);
    alpha_2=-V_y2/(c*(1-lambda)*delta_y);
    alpha_norm=sqrt(alpha_1.*alpha_1+alpha_2.*alpha_2);
    if bound_alpha
        alpha_1(alpha_norm>alpha_max)=alpha_1(alpha_norm>alpha_max)./alpha_norm(alpha_norm>alpha_max)*alpha_max;
        alpha_2(alpha_norm>alpha_max)=alpha_2(alpha_norm>alpha_max)./alpha_norm(alpha_norm>alpha_max)*alpha_max;
    end
    
    
    %%%% Using convolution 
    %Linear, using the previous estimate for V to approximate gradient V
    if K>1
        V(n,:,:,:,:)=V_curr+delta_t*(sigma^2/2*(V_yy1+V_yy2)/(delta_y)^2+(y_j1.*V_x1+y_j2.*V_x2)/delta_x+alpha_1.*V_y1/delta_y+alpha_2.*V_y2/delta_y+1/2*c*(1-lambda)*(squeeze(old_alpha_1(n,:,:,:,:)).*alpha_1+squeeze(old_alpha_2(n,:,:,:,:)).*alpha_2)+c*lambda*F(:,:,:,:));
    else
        V(n,:,:,:,:)=V_curr+delta_t*(sigma^2/2*(V_yy1+V_yy2)/(delta_y)^2+(y_j1.*V_x1+y_j2.*V_x2)/delta_x+alpha_1.*V_y1/delta_y+alpha_2.*V_y2/delta_y+1/2*c*(1-lambda)*(alpha_1.*alpha_1+alpha_2.*alpha_2)+c*lambda*F(:,:,:,:));
    end
    
    if K>1
        alpha_1_diff=abs(squeeze(old_alpha_1(n,:,:,:,:))-alpha_1);
        alpha_1_diff_real=zeros(num_x,num_y,num_x,num_y);
        alpha_1_diff_real(mu_curr>0)=alpha_1_diff(mu_curr>0);
        max_curr=max(max(max(max(alpha_1_diff_real))));
        max_diff_alpha_1=max(max_curr,max_diff_alpha_1);
        alpha_2_diff=abs(squeeze(old_alpha_2(n,:,:,:,:))-alpha_2);
        alpha_2_diff_real=zeros(num_x,num_y,num_x,num_y);
        alpha_2_diff_real(mu_curr>0)=alpha_2_diff(mu_curr>0);
        max_curr=max(max(max(max(alpha_2_diff_real))));
        max_diff_alpha_2=max(max_curr,max_diff_alpha_2);
    end

    alpha_norm=sqrt(alpha_1.*alpha_1+alpha_2.*alpha_2);
    max_alpha=max(max_alpha,max(max(max(max(alpha_norm)))));
    old_alpha_1(n,:,:,:,:)=alpha_1;
    old_alpha_2(n,:,:,:,:)=alpha_2;
end

%Checking if the solution is valid:
if K>1
    'Difference in alpha_1'
    max_diff_alpha_1
    'Difference in alpha_2'
    max_diff_alpha_2
end
'Part B Max Alpha'
max_alpha
if max_alpha>alpha_max+threshold
    'Oh no!!!!! Part B max_alpha>alpha_max'
    max_alpha
end
'Minimum of V'
value=min(min(min(min(min(V(:,:,:,:,:))))))
if value<-threshold
    'Oh no!!!!! Negatives in V'
    value
end
'Maximum of V'
value=max(max(max(max(max(V(:,:,:,:,:))))))


'Part C: Kolmogorov'
%% Given V, solve for mu
% max_V_y=0;
max_alpha=0;
mu=zeros(num_time_points,num_x,num_y,num_x,num_y);

if initial_mu_guess
    mu_5D=reshape(mu_guess,[1,size(mu_guess)]);
    mu(1,:,:,:,:)=repmat(mu_3D,num_time_points,1,1,1,1);
elseif initial_box
    mu(1,ceil(num_x/2),ceil(num_y/2)-box_r:ceil(num_y/2)+box_r,ceil(num_x/2),ceil(num_y/2)-box_r:ceil(num_y/2)+box_r)=1/((2*box_r+1)^2*(delta_x*delta_y)^2);
elseif initial_box_x
    mu(1,ceil(num_x/2)-box_r:ceil(num_x/2)+box_r,ceil(num_y/2),ceil(num_x/2)-box_r:ceil(num_x/2)+box_r,ceil(num_y/2))=1/((2*box_r+1)^2*(delta_x*delta_y)^2);
elseif initial_2_boxes
    mu(1,ceil(num_x/2)-box_r:ceil(num_x/2),ceil(num_y/2):ceil(num_y/2)+box_r,ceil(num_x/2)-box_r:ceil(num_x/2),ceil(num_y/2):ceil(num_y/2)+box_r)=1/((2*(box_r+1)^2-1)^2*(delta_x*delta_y)^2);
    mu(1,ceil(num_x/2):ceil(num_x/2)+box_r,ceil(num_y/2)-box_r:ceil(num_y/2),ceil(num_x/2):ceil(num_x/2)+box_r,ceil(num_y/2)-box_r:ceil(num_y/2))=1/((2*(box_r+1)^2-1)^2*(delta_x*delta_y)^2);
elseif initial_2_points
    mu(1,ceil(num_x/2)-box_r,ceil(num_y/2),ceil(num_x/2)-box_r,ceil(num_y/2))=1/(2*(delta_x*delta_y)^2);
    mu(1,ceil(num_x/2)+box_r,ceil(num_y/2),ceil(num_x/2)+box_r,ceil(num_y/2))=1/(2*(delta_x*delta_y)^2);
elseif initial_5_points
    mu(1,ceil(num_x/2)-box_r,ceil(num_y/2),ceil(num_x/2)-box_r,ceil(num_y/2))=1/(5*(delta_x*delta_y)^2);
    mu(1,ceil(num_x/2)+box_r,ceil(num_y/2),ceil(num_x/2)+box_r,ceil(num_y/2))=1/(5*(delta_x*delta_y)^2);
    mu(1,ceil(num_x/2)-floor(box_r/2),ceil(num_y/2),ceil(num_x/2)-floor(box_r/2),ceil(num_y/2))=1/(5*(delta_x*delta_y)^2);
    mu(1,ceil(num_x/2)+floor(box_r/2),ceil(num_y/2),ceil(num_x/2)+floor(box_r/2),ceil(num_y/2))=1/(5*(delta_x*delta_y)^2);
    mu(1,ceil(num_x/2),ceil(num_y/2),ceil(num_x/2),ceil(num_y/2))=1/(5*(delta_x*delta_y)^2);
elseif initial_5_points_xandv
    mu(1,ceil(num_x/2)-box_r,ceil(num_y/2)-box_r_y,ceil(num_x/2)-box_r,ceil(num_y/2)-box_r_y)=1/(5*(delta_x*delta_y)^2);
    mu(1,ceil(num_x/2)+box_r,ceil(num_y/2)-box_r_y,ceil(num_x/2)+box_r,ceil(num_y/2)-box_r_y)=1/(5*(delta_x*delta_y)^2);
    mu(1,ceil(num_x/2)-box_r,ceil(num_y/2)+box_r_y,ceil(num_x/2)-box_r,ceil(num_y/2)+box_r_y)=1/(5*(delta_x*delta_y)^2);
    mu(1,ceil(num_x/2)+box_r,ceil(num_y/2)+box_r_y,ceil(num_x/2)+box_r,ceil(num_y/2)+box_r_y)=1/(5*(delta_x*delta_y)^2);
    mu(1,ceil(num_x/2),ceil(num_y/2),ceil(num_x/2),ceil(num_y/2))=1/(5*(delta_x*delta_y)^2);
elseif initial_skew
    mu(1,ceil(num_x/2)-box_r,ceil(num_y/2)-box_r_y,ceil(num_x/2)-box_r,ceil(num_y/2)-box_r_y)=1/(3*(delta_x*delta_y)^2);
    mu(1,ceil(num_x/2)+box_r,ceil(num_y/2)+box_r_y,ceil(num_x/2)+box_r,ceil(num_y/2)+box_r_y)=1/(3*(delta_x*delta_y)^2);
    mu(1,ceil(num_x/2),ceil(num_y/2),ceil(num_x/2),ceil(num_y/2))=1/(3*(delta_x*delta_y)^2);
elseif initial_skew2
    mu(1,ceil(num_x/2)-box_r,ceil(num_y/2)+box_r_y,ceil(num_x/2)-box_r,ceil(num_y/2)+box_r_y)=1/(3*(delta_x*delta_y)^2);
    mu(1,ceil(num_x/2)+box_r,ceil(num_y/2)-box_r_y,ceil(num_x/2)+box_r,ceil(num_y/2)-box_r_y)=1/(3*(delta_x*delta_y)^2);
    mu(1,ceil(num_x/2),ceil(num_y/2),ceil(num_x/2),ceil(num_y/2))=1/(3*(delta_x*delta_y)^2);
else
    mu(1,ceil(num_x/2),ceil(num_y/2),ceil(num_x/2),ceil(num_y/2))=1/(delta_x*delta_y)^2; %puts everything at the origin.
end

for n=1:num_time_points-1
    t_n=t_grid(n);
        
    mu_curr=squeeze(mu(n,:,:,:,:));
    V_curr=squeeze(V(n,:,:,:,:));
    
    mu_yy1(:,:,:,:)=shift(mu_curr,1,2)-2*mu_curr+shift(mu_curr,-1,2); %centered
    %New Scheme BC:
    mu_yy1(:,2,:,:)=mu_yy1(:,2,:,:)+mu_curr(:,1,:,:);
    mu_yy1(:,num_y-1,:,:)=mu_yy1(:,num_y-1,:,:)+mu_curr(:,num_y,:,:);
    mu_yy1(:,1,:,:)=mu_yy1(:,1,:,:)-mu_curr(:,num_y,:,:);
    mu_yy1(:,num_y,:,:)=mu_yy1(:,num_y,:,:)-mu_curr(:,1,:,:);
    V_yy1(:,:,:,:)=shift(V_curr,1,2)-2*V_curr+shift(V_curr,-1,2); %centered
    mu_x1=zeros(num_x,num_y,num_x,num_y);
    mu_x1(:,1:ceil(num_y/2)-1,:,:)=shift(mu_curr(:,1:ceil(num_y/2)-1,:,:),1,1)-mu_curr(:,1:ceil(num_y/2)-1,:,:); %if v_j<0, one sided forward
    mu_x1(:,ceil(num_y/2)+1:num_y,:,:)=mu_curr(:,ceil(num_y/2)+1:num_y,:,:)-shift(mu_curr(:,ceil(num_y/2)+1:num_y,:,:),-1,1); %if v_j>0, one sided backwards (if v_j=0, 0)
    %New Scheme BC:
    mu_x1(2,ceil(num_y/2)+1:num_y,:,:)=mu_curr(2,ceil(num_y/2)+1:num_y,:,:);
    mu_x1(num_x-1,1:ceil(num_y/2)-1,:,:)=-mu_curr(num_x-1,1:ceil(num_y/2)-1,:,:);
    mu_x1(1,:,:,:)=0;
    mu_x1(1,1:ceil(num_y/2)-1,:,:)=mu_curr(2,1:ceil(num_y/2)-1,:,:);
    mu_x1(num_x,:,:,:)=0;
    mu_x1(num_x,ceil(num_y/2)+1:num_y,:,:)=-mu_curr(num_x-1,ceil(num_y/2)+1:num_y,:,:);
    
    left1(:,:,:,:)=V_curr-shift(V_curr,-1,2);
    right1(:,:,:,:)=shift(V_curr,1,2)-V_curr;

    V_y1=zeros(num_x,num_y,num_x,num_y);
    V_y1(left1<0 & right1<0)=right1(left1<0 & right1<0);
    V_y1(left1>0 & right1>0)=left1(left1>0 & right1>0);
    V_y1(:,1,:,:)=0;
    V_y1(:,num_y,:,:)=0;
    
    mu_yy2(:,:,:,:)=shift(mu_curr,1,4)-2*mu_curr+shift(mu_curr,-1,4); %centered
    %New Scheme BC:
    mu_yy2(:,:,:,2)=mu_yy2(:,:,:,2)+mu_curr(:,:,:,1);
    mu_yy2(:,:,:,num_y-1)=mu_yy2(:,:,:,num_y-1)+mu_curr(:,:,:,num_y);
    mu_yy2(:,:,:,1)=mu_yy2(:,:,:,1)-mu_curr(:,:,:,num_y);
    mu_yy2(:,:,:,num_y)=mu_yy2(:,:,:,num_y)-mu_curr(:,:,:,1);
    V_yy2(:,:,:,:)=shift(V_curr,1,4)-2*V_curr+shift(V_curr,-1,4); %centered
    mu_x2=zeros(num_x,num_y,num_x,num_y);
    mu_x2(:,:,:,1:ceil(num_y/2)-1)=shift(mu_curr(:,:,:,1:ceil(num_y/2)-1),1,3)-mu_curr(:,:,:,1:ceil(num_y/2)-1); %if v_j<0, one sided forward
    mu_x2(:,:,:,ceil(num_y/2)+1:num_y)=mu_curr(:,:,:,ceil(num_y/2)+1:num_y)-shift(mu_curr(:,:,:,ceil(num_y/2)+1:num_y),-1,3); %if v_j>0, one sided backwards (if v_j=0, 0)
    %New Scheme BC:
    mu_x2(:,:,2,ceil(num_y/2)+1:num_y)=mu_curr(:,:,2,ceil(num_y/2)+1:num_y);
    mu_x2(:,:,num_x-1,1:ceil(num_y/2)-1)=-mu_curr(:,:,num_x-1,1:ceil(num_y/2)-1);
    mu_x2(:,:,1,:)=0;
    mu_x2(:,:,1,1:ceil(num_y/2)-1)=mu_curr(:,:,2,1:ceil(num_y/2)-1);
    mu_x2(:,:,num_x,:)=0;
    mu_x2(:,:,num_x,ceil(num_y/2)+1:num_y)=-mu_curr(:,:,num_x-1,ceil(num_y/2)+1:num_y);
    
    left2(:,:,:,:)=V_curr-shift(V_curr,-1,4);
    right2(:,:,:,:)=shift(V_curr,1,4)-V_curr;

    V_y2=zeros(num_x,num_y,num_x,num_y);
    V_y2(left2<0 & right2<0)=right2(left2<0 & right2<0);
    V_y2(left2>0 & right2>0)=left2(left2>0 & right2>0);
    V_y2(:,:,:,1)=0;
    V_y2(:,:,:,num_y)=0;
    
    alpha_1=-V_y1/(c*(1-lambda)*delta_y);
    alpha_2=-V_y2/(c*(1-lambda)*delta_y);
    alpha_norm=sqrt(alpha_1.*alpha_1+alpha_2.*alpha_2);
    if bound_alpha
        alpha_1(alpha_norm>alpha_max)=alpha_1(alpha_norm>alpha_max)./alpha_norm(alpha_norm>alpha_max)*alpha_max;
        alpha_2(alpha_norm>alpha_max)=alpha_2(alpha_norm>alpha_max)./alpha_norm(alpha_norm>alpha_max)*alpha_max;
    end

    alpha_1_minus=zeros(num_x,num_y,num_x,num_y);
    alpha_1_plus=zeros(num_x,num_y,num_x,num_y);
    alpha_1_minus(alpha_1<0)=-alpha_1(alpha_1<0);
    alpha_1_plus(alpha_1>0)=alpha_1(alpha_1>0);
    alpha_2_minus=zeros(num_x,num_y,num_x,num_y);
    alpha_2_plus=zeros(num_x,num_y,num_x,num_y);
    alpha_2_minus(alpha_2<0)=-alpha_2(alpha_2<0);
    alpha_2_plus(alpha_2>0)=alpha_2(alpha_2>0);
    
    alpha_term_a_1=shift(alpha_1_minus,1,2);
    alpha_term_b_1=shift(alpha_1_plus,-1,2);
    alpha_term_c_1=-(alpha_1_plus+alpha_1_minus);
    %New Scheme BC:
    alpha_term_b_1(:,2,:,:)=0;
    alpha_term_a_1(:,num_y-1,:,:)=0;
    alpha_term_b_1(:,1,:,:)=0;
    alpha_term_a_1(:,num_y,:,:)=0;
    alpha_term_c_1(:,1,:,:)=0;
    alpha_term_c_1(:,num_y,:,:)=0;
    
    alpha_term_a_2=shift(alpha_2_minus,1,4);
    alpha_term_b_2=shift(alpha_2_plus,-1,4);
    alpha_term_c_2=-(alpha_2_plus+alpha_2_minus);
    %New Scheme BC:
    alpha_term_b_2(:,:,:,2)=0;
    alpha_term_a_2(:,:,:,num_y-1)=0;
    alpha_term_b_2(:,:,:,1)=0;
    alpha_term_a_2(:,:,:,num_y)=0;
    alpha_term_c_2(:,:,:,1)=0;
    alpha_term_c_2(:,:,:,num_y)=0;
    
    mu(n+1,:,:,:,:)=mu_curr+delta_t*(sigma^2/2*(mu_yy1+mu_yy2)/(delta_y)^2-(y_j1.*mu_x1+y_j2.*mu_x2)/delta_x+alpha_term_a_1.*shift(mu_curr,1,2)/delta_y+alpha_term_b_1.*shift(mu_curr,-1,2)/delta_y+alpha_term_a_2.*shift(mu_curr,1,4)/delta_y+alpha_term_b_2.*shift(mu_curr,-1,4)/delta_y+(alpha_term_c_1+alpha_term_c_2).*mu_curr/delta_y);

    alpha_norm=sqrt(alpha_1.*alpha_1+alpha_2.*alpha_2);
    max_alpha=max(max_alpha,max(max(max(max(alpha_norm)))));
end

%Checking if the solution is valid:
integral_values_2=zeros(num_time_points,1);
for n=1:num_time_points
    integral=sum(sum(sum(sum(mu(n,:,:,:,:)))))*(delta_x*delta_y)^2;
    integral_values_2(n)=integral;
    if integral>1+threshold || integral<1-threshold
        'Oh no!!!!! Sum is not 1!!!!'
        integral
    end
end
if normalize
    for n=1:num_time_points
        mu(n,:,:,:,:)=mu(n,:,:,:,:)/integral_values_2(n);
    end
end
'Part C Max Alpha'
max_alpha
if max_alpha>alpha_max+threshold
    'Oh no!!!!! Part C max_alpha>alpha_max'
    max_alpha
end
%'Minimum Density in Part C'
value=min(min(min(min(min(mu(:,:,:,:,:))))));
if value<-threshold
    'Oh no!!!!! Negatives in mu'
    value
end
mu_diff=abs(mu-old_mu);
mu_diff_max=max(max(max(max(max(mu_diff(:,:,:,:,:))))))*(delta_x*delta_y)^2;
value=mu_diff_max;
K
'Difference in mu'
value
mu_diff_frac=mu_diff./abs(mu);
mu_diff_frac(abs(mu)<10^(-10) & abs(old_mu)<10^(-10))=-1;
value=max(max(max(max(max(mu_diff_frac(:,:,:,:,:))))));
'Largest Fractional Difference'
value
K=K+1;
end %This ends the while loop

%Final Calculations
timer=toc
final_mu=squeeze(mu(num_time_points,:,:,:,:)).*(delta_x*delta_y)^2;

mu_12=sum(sum(final_mu,3),4).*(delta_x*delta_y)^2;
mu_13=squeeze(sum(sum(final_mu,2),4)).*(delta_x*delta_y)^2;
mu_24=squeeze(sum(sum(final_mu,1),3)).*(delta_x*delta_y)^2;
mu_34=squeeze(sum(sum(final_mu,1),2)).*(delta_x*delta_y)^2;

save(strcat(jobstring,'_mu_12.mat'),'mu_12')
save(strcat(jobstring,'_mu_13.mat'),'mu_13')
save(strcat(jobstring,'_mu_24.mat'),'mu_24')
save(strcat(jobstring,'_mu_34.mat'),'mu_34')

output_freq=1000;
num_times=floor(num_time_points/output_freq)+1;
mu_short=zeros(num_x,num_y,num_times);
for i=1:num_times
    mu_curr=squeeze(mu((i-1)*output_freq+1,:,:,:,:)).*(delta_x*delta_y)^2;
    mu12_short(:,:,i)=sum(sum(mu_curr,3),4).*(delta_x*delta_y)^2;
    mu13_short(:,:,i)=squeeze(sum(sum(mu_curr,2),4)).*(delta_x*delta_y)^2;
    mu24_short(:,:,i)=squeeze(sum(sum(mu_curr,1),3)).*(delta_x*delta_y)^2;
    mu34_short(:,:,i)=squeeze(sum(sum(mu_curr,1),2)).*(delta_x*delta_y)^2;
end
save(strcat(jobstring,'_mu12_short.mat'),'mu12_short','-v7.3')
save(strcat(jobstring,'_mu13_short.mat'),'mu13_short','-v7.3')
save(strcat(jobstring,'_mu24_short.mat'),'mu24_short','-v7.3')
save(strcat(jobstring,'_mu34_short.mat'),'mu34_short','-v7.3')


% save(strcat(jobstring,'_final_mu.mat'),'final_mu')
% initial_V=squeeze(V(1,:,:,:,:));
% save(strcat(jobstring,'_initial_V.mat'),'initial_V')
% %save(strcat(jobstring,'_V.mat'),'V','-v7.3')
% 
% num_times=floor(num_time_points/1000);
% mu_short=zeros(num_x,num_y,num_times);
% for i=1:num_times
%     mu_short(:,:,i)=mu(i*1000,:,:);
% end
% %save(strcat(jobstring,'_mu.mat'),'mu','-v7.3')
% save(strcat(jobstring,'_mu_short.mat'),'mu_short','-v7.3')
save(strcat(jobstring,'_integral_values_2.mat'),'integral_values_2')
'Done'
