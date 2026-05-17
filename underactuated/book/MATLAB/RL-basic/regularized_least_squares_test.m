function regularized_least_squares_test    
    x = [-1:0.01:1];%specify a set of sample points
    N = length(x);
    for i=1:N
        ytrue(i) = f(x(i));%evaluate the true function
        y(i) = f(x(i))+0.1*randn(1,1);%evaluate the true function with noise
    end
    param=[];
    theta=regularized_least_squares(@cubic,x,y',param) %apply regularized least-squares
    yest=eval_cubic(x,theta);%eval the estimated output

    %plot the results
    figure(1)
    plot(x,y,'r',x,ytrue,'--k',x,yest,'-b')
    xlabel('x')
    ylabel('y')
end

%actual cubic function
function y=f(x)
    y = 0.1*x^3-0.2*x^2+0.3*x+0.4;
end

%cubic basis functions
function phi=cubic(x,param)
    phi=[x(1)^3 x(1)^2 x(1) 1];
    phi=phi';
end

%evaluate the cubic function
function y=eval_cubic(x,theta)
    param=[];
    for i=1:length(x)
        phi=cubic(x(i),param);
        y(i)= theta'*phi;
    end
end

%compute the parameters using regularized least squares
function theta=regularized_least_squares(BASIS,x,y,param)
    N=length(x);
    Phi=[];
    for k=1:N
        phi=BASIS(x(:,k),param);
        Phi=[Phi phi];
    end
    gamma=0.00;
    Phi=Phi';
    P=(Phi'*Phi);
    I=eye(length(P));
    B=Phi'*y;
    theta=(P+gamma*I)\B;
end