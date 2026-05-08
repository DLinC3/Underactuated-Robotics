function param=init_param()

    param.l1 = 1;
    param.l2 = 2;  
    param.m1 = 1; 
    param.m2 = 1;  
    param.g = 9.81;
    param.b1=.2;  
    param.b2=.2;
    param.w1=0;
    param.w2=1000;
%    b1=0; b2=0;
    param.lc1 = .5; 
    param.lc2 = 1; 
    param.I1 = 0.083 + param.m1*param.lc1^2;  
    param.I2 = 0.33 + param.m2*param.lc2^2;  
end