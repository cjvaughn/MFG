function[lambda]=eval_lambda(V,mu,vpad,num_x,num_y,delta_x,delta_y,sigma,y_grid)

u(:,:)=mu(:,:)*delta_x*delta_y; %u is (num_x,num_y), v is (2*num_1,2*num_y-1), output=(3*num_x-2,3*num_y-2)
    
upad=zeros(3*num_x-2,3*num_y-2);
upad(1:num_x,1:num_y)=u;

F2=ifftn(fftn(upad).*fftn(vpad));
F=F2(num_x:2*num_x-1,num_y:2*num_y-1);

lambda=0;
for j=1:num_y-1
    if j==1
        jminus=num_y;
        jplus=2;
    elseif j==num_y
        jminus=num_y-1;
        jplus=1;
    else
        jminus=j-1;
        jplus=j+1;
    end
    for i=1:num_x-1
        if j<ceil(num_y/2) %forwards
            grad_v_V=(V(i,jplus)-V(i,j))/delta_y;
        elseif j>ceil(num_y/2) %backwards
            grad_v_V=(V(i,j)-V(i,jminus))/delta_y;
        end
        lambda=lambda+(1/2*grad_v_V^2+F(i,j))*mu(i,j)*delta_x*delta_y;
    end
end

end