function[f_shifted]=shift(f,difference_type,axis)
if difference_type==1;
    f_shifted=circshift(f,-1,axis);
else difference_type==-1;
    f_shifted=circshift(f,1,axis);
end
end