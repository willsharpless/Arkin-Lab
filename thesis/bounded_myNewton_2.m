function [lambda, its, hist] = bounded_myNewton_2(f, df, x0, TOL, maxit)
%Bounded Newton's Method - Will Sharpless
    %willsharpless@berkeley.edu
    %Newton's method with a predetermined gradient f and hessian df. This
    %version uses function arrays in order to rapidly calculate.

its = 0; hist(:,1) = cell2mat(x0); lambda = x0;
while norm(f(lambda{:})) > TOL && its < maxit
    lambda = cell2mat(lambda) - df(lambda{:})\f(lambda{:});
    for i=1:3
        if lambda(i) < 0
            lambda(i) = 0;
        end
    end
    hist = [hist lambda]; lambda = num2cell(lambda);
    its = its + 1;
end
end