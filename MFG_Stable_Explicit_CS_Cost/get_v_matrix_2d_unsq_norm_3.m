function[v]=get_v_matrix_2d_unsq_norm_3(num_x,num_y,delta_x,delta_y,beta)

i1=1:num_x;
i2=1:num_x;
j1=1:num_y;
j2=1:num_y;
a1=((i1-1)*delta_x).^2;
a2=((i2-1)*delta_x).^2;
a=1./(1+(bsxfun(@plus,a1,a2')).^0.5).^beta; %+ here needs to iterate over all combinations
b1=((j1-1)*delta_y).^2;
b2=((j2-1)*delta_y).^2;
b=(bsxfun(@plus,b1,b2')).^0.5; %+ here needs to iterate over all combinations
a_vector=reshape(a,num_x^2,1);
b_vector=reshape(b,1,num_y^2);
c_matrix=a_vector*b_vector;
v2=reshape(c_matrix,num_x,num_x,num_y,num_y); %* here needs to multiply over all combinations

v3=cat(1,flip(v2,1),v2);
v4=cat(2,flip(v3,2),v3);
v5=cat(3,flip(v4,3),v4);
v6=cat(4,flip(v5,4),v5);
v=v6([1:num_x-1,num_x+1:end],[1:num_x-1,num_x+1:end],[1:num_y-1,num_y+1:end],[1:num_y-1,num_y+1:end]);

v=permute(v,[1 3 2 4]);

end