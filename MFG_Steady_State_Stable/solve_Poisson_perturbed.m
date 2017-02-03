function[mu,r,diff,M,b,diff_orig]=solve_Poisson_perturbed(ro,old_mu,V,num_x,num_y,delta_x,delta_y,epsilon,sigma,y_grid,alpha)

M=zeros(num_x*num_y,num_x*num_y);
b=zeros(num_x*num_y,1);

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

for j=1:num_y
    for i=1:num_x
        jminus=j-1;
        jplus=j+1;
        y_j=y_grid(j);
        
        ij=ij_lookup(i,j,num_x,num_y);
        ilj=ij_lookup(i-1,j,num_x,num_y);
        irj=ij_lookup(i+1,j,num_x,num_y);
        ijl=ij_lookup(i,j-1,num_x,num_y);
        ijr=ij_lookup(i,j+1,num_x,num_y);
        
        eq=ij;
        
%         %This makes it a viscosity solution
%         %from - laplacian_x mu
%         M(eq,irj)=M(eq,irj)-epsilon^2/(2*delta_x^2);
%         M(eq,ij)=M(eq,ij)+epsilon^2/(delta_x^2);
%         M(eq,ilj)=M(eq,ilj)-epsilon^2/(2*delta_x^2);
        
        %from -laplacian_v mu
        M(eq,ijr)=M(eq,ijr)-sigma^2/(2*delta_y^2);
        M(eq,ij)=M(eq,ij)+sigma^2/(delta_y^2);
        M(eq,ijl)=M(eq,ijl)-sigma^2/(2*delta_y^2);
        if j==2
            M(eq,ijl)=M(eq,ijl)-sigma^2/(2*delta_y^2);
        elseif j==num_y-1
            M(eq,ijr)=M(eq,ijr)-sigma^2/(2*delta_y^2);
        elseif j==1
            M(eq,ijl)=M(eq,ijl)+sigma^2/(2*delta_y^2);
        elseif j==num_y
            M(eq,ijr)=M(eq,ijr)+sigma^2/(2*delta_y^2);
        end
        
        %from v* grad_x *mu
        if i>2 && i<num_x-1
            if j<ceil(num_y/2) %if v_j<0, forwards
                M(eq,irj)=M(eq,irj)+y_j/delta_x;
                M(eq,ij)=M(eq,ij)-y_j/delta_x;
            elseif j>ceil(num_y/2) %if v_j>0, backwards
                M(eq,ij)=M(eq,ij)+y_j/delta_x;
                M(eq,ilj)=M(eq,ilj)-y_j/delta_x;
            end
        elseif i==2
            if j<ceil(num_y/2)
                M(eq,irj)=M(eq,irj)+y_j/delta_x;
                M(eq,ij)=M(eq,ij)-y_j/delta_x;
            elseif j>ceil(num_y/2)
                M(eq,ij)=M(eq,ij)+y_j/delta_x;
            end
        elseif i==num_x-1
            if j<ceil(num_y/2)
                M(eq,ij)=M(eq,ij)-y_j/delta_x;
            elseif j>ceil(num_y/2)
                M(eq,ij)=M(eq,ij)+y_j/delta_x;
                M(eq,ilj)=M(eq,ilj)-y_j/delta_x;
            end
        elseif i==1
            if j<ceil(num_y/2)
                M(eq,irj)=M(eq,irj)+y_j/delta_x;
            end
        elseif i==num_x
            if j>ceil(num_y/2)
                M(eq,ilj)=M(eq,ilj)-y_j/delta_x;
            end
        end
        
        M(eq,ijr)=M(eq,ijr)-alpha_term_a(i,j)/delta_y;
        M(eq,ijl)=M(eq,ijl)-alpha_term_b(i,j)/delta_y;
        M(eq,ij)=M(eq,ij)-alpha_term_c(i,j)/delta_y;
        
        
        %This is the perturbed part
        M(eq,ij)=M(eq,ij)+ro;
        b(eq,1)=b(ij,1)+old_mu(i,j)*ro;
    end
end

r=rank(M);
mu=M\b;
diff=M*mu-b;

old_M=M-ro*eye(num_x*num_y);

diff_orig=old_M*mu;

%mu=linsolve(M,b);

% Aeq=ones(1,num_x*num_y)*delta_x*delta_y;
% lb=zeros(num_x*num_y,1);
% ub=Inf(num_x*num_y,1);
% mu=lsqlin(M,b,[],[],Aeq,1,lb,ub);

mu=reshape(mu,[num_x,num_y]);
end