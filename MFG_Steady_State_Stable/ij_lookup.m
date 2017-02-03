function[n]=ij_lookup(i,j,num_x,num_y)
i=mod(i-1,num_x)+1;
j=mod(j-1,num_y)+1;
n=i+(j-1)*num_x;
% if i>num_x
%     n=false
% end
% if j>num_y
%     n=false
% end
end