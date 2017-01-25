num_time_points=2501;
num_x=33; %495; %495; %573; %399; %1567; %787;
num_y=41;

delta_x=0.5;
delta_y=0.05;

x_max=(num_x-1)/2*delta_x;
y_max=(num_y-1)/2*delta_y;


time=1;
f(:,:)=mu_short(:,:,time);
integral_value=sum(sum(f));
f=f/integral_value; 

discard=0;
x_left=discard+1;
x_right=num_x-discard;
f=f(x_left:x_right,:);
a=-x_max+(x_left-1)*delta_x;
b=-x_max+(x_right-1)*delta_x;

f2=rot90(f);

colormap(flipud(gray))
x=[a b];
y=[-y_max y_max];
colorbar;
imagesc(x,y,f2);
xlabel('x1')
ylabel('v1')
title('mu(t,x,v) c=2, lambda=0.5') %beta=1.5,

size_mu_short=size(mu_short);
max_time=size_mu_short(3);
while true
    [~,~,b] = ginput(1);
    if b==29
        time=time+1;
    elseif b==28
        time=time-1;
    elseif b==27
        break
    end
    if time>0 && time<=max_time
        f(:,:)=mu_short(:,:,time);
        integral_value=sum(sum(f));
        f=f/integral_value; 
        discard=0;
        x_left=discard+1;
        x_right=num_x-discard;
        f=f(x_left:x_right,:);
        a=-x_max+(x_left-1)*delta_x;
        b=-x_max+(x_right-1)*delta_x;
        f2=rot90(f);
        imagesc(x,y,f2);
        xlabel('x1')
        ylabel('v1')
        title('mu(t,x,v), beta=1.5, c=2, lambda=0.5')
    end
end