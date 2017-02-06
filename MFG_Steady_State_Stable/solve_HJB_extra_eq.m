function[V,r,lambda,diff,M]=solve_HJB_extra_eq(mu,old_V_y,vpad,num_x,num_y,delta_x,delta_y,epsilon,sigma,y_grid,c,lambda_cost)

M=zeros(num_x*num_y+2,num_x*num_y+1);
b=zeros(num_x*num_y+2,1);
b(num_x*num_y+1,1)=0; %redundant, but emphasis on sum of the V_ij is 0

u(:,:)=mu(:,:)*delta_x*delta_y; %u is (num_x,num_y), v is (2*num_1,2*num_y-1), output=(3*num_x-2,3*num_y-2)
    
upad=zeros(3*num_x-2,3*num_y-2);
upad(1:num_x,1:num_y)=u;

F2=ifftn(fftn(upad).*fftn(vpad));
F=F2(num_x:2*num_x-1,num_y:2*num_y-1);

for j=1:num_y
    for i=1:num_x
        y_j=y_grid(j);
        
        ij=ij_lookup(i,j,num_x,num_y);
        ilj=ij_lookup(i-1,j,num_x,num_y);
        irj=ij_lookup(i+1,j,num_x,num_y);
        ijl=ij_lookup(i,j-1,num_x,num_y);
        ijr=ij_lookup(i,j+1,num_x,num_y);
        
        
        %This part is the HJB equation
        eq=ij;
        
%         %This makes it a viscosity solution
%         %from laplacian_x V
%         M(eq,irj)=M(eq,irj)+epsilon^2/(2*delta_x^2);
%         M(eq,ij)=M(eq,ij)-epsilon^2/(delta_x^2);
%         M(eq,ilj)=M(eq,ilj)+epsilon^2/(2*delta_x^2);
        
        %from laplacian_v V
        if j>1 && j<num_y
            M(eq,ijr)=M(eq,ijr)+sigma^2/(2*delta_y^2);
            M(eq,ij)=M(eq,ij)-sigma^2/(delta_y^2);
            M(eq,ijl)=M(eq,ijl)+sigma^2/(2*delta_y^2);
        elseif j==1
            M(eq,ijr)=M(eq,ijr)+sigma^2/(delta_y^2);
            M(eq,ij)=M(eq,ij)-sigma^2/(delta_y^2);
        elseif j==num_y            
            M(eq,ijl)=M(eq,ijl)+sigma^2/(delta_y^2);
            M(eq,ij)=M(eq,ij)-sigma^2/(delta_y^2); 
        end
        
        %next part is v* grad_x V
        if i>1 && i<num_x
        if j<ceil(num_y/2) %if v_j<0, backwards
            M(eq,ij)=M(eq,ij)+y_j/delta_x;
            M(eq,ilj)=M(eq,ilj)-y_j/delta_x;
        elseif j>ceil(num_y/2) %if v_j>0, forwards
            M(eq,irj)=M(eq,irj)+y_j/delta_x;
            M(eq,ij)=M(eq,ij)-y_j/delta_x;
        end
        end
        
        %next part is -1/2*(grad_v V)^2
        if j>1 && j<num_y
        if old_V_y(i,j)<0 %forwards
            M(eq,ijr)=M(eq,ijr)+1/delta_y*(-1/(2*c*(1-lambda_cost))*old_V_y(i,j));
            M(eq,ij)=M(eq,ij)-1/delta_y*(-1/(2*c*(1-lambda_cost))*old_V_y(i,j));
        elseif old_V_y(i,j)>0 %backwards
            M(eq,ij)=M(eq,ij)+1/delta_y*(-1/(2*c*(1-lambda_cost))*old_V_y(i,j));
            M(eq,ijl)=M(eq,ijl)-1/delta_y*(-1/(2*c*(1-lambda_cost))*old_V_y(i,j));
        end
        end
        
        %this is the lambda term
        M(eq,num_x*num_y+1)=M(eq,num_x*num_y+1)-1;
        
        %this is the F term
        b(eq)=b(eq)-c*lambda_cost*F(i,j);
        
        M(num_x*num_y+1,ij)=M(num_x*num_y+1,ij)+1; %sum of the V_ij is 0
    end
end
%add symmetry constraint to buff up rank of M
center=ij_lookup(ceil(num_x/2),ceil(num_y/2),num_x,num_y);
first=ij_lookup(1,1,num_x,num_y);
last=ij_lookup(num_x,num_y,num_x,num_y);
M(num_x*num_y+2,first)=1;
M(num_x*num_y+2,last)=-1;

r=rank(M);
V=M\b;
diff=M*V-b;
%V=linsolve(M,b);
lambda=V(num_x*num_y+1);
V=reshape(V(1:num_x*num_y),[num_x,num_y]);
end