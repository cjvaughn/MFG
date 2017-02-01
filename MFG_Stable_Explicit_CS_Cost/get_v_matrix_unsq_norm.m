function[v]=get_v_matrix_unsq_norm(num_x,num_y,delta_x,delta_y,beta)
i=1:num_x;
j=1:num_y;
a=1./(1+((i-1)*delta_x).^2).^beta;
b=(j-1)*delta_y;
v2=a'*b;

% v=zeros(2*num_x-1,2*num_y-1);
% v(num_x:2*num_x-1,num_y:2*num_y-1)=v2;
% v(num_x:2*num_x-1,1:num_y)=fliplr(v2);
% v(1:num_x,num_y:2*num_y-1)=flipud(v2);
% v(1:num_x,1:num_y)=fliplr(flipud(v2));

v3=cat(1,flip(v2,1),v2);
v4=cat(2,flip(v3,2),v3);
v=v4([1:num_x-1,num_x+1:end],[1:num_y-1,num_y+1:end]);
end