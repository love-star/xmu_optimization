function [alpha] = Armjio(fun, grid, x0, dk)
    beta = 0.333; 		% 步长 alpha 的迭代系数，小于 1
    rho = 1e-3; 		% 泰勒展开式补足系数，0 < rho < 1/2
    alpha = 1; 			% 初始步长为 1
    k = 0; 				% 统计迭代次数
    gk = feval(grid, x0);	% x0处的梯度值
    fd = feval(fun, x0 + alpha * dk); 	% 函数在下一个迭代点处的目标函数值
    fk = feval(fun, x0) + alpha * rho * gk' * dk; 	% 函数在下一个迭代点处的泰勒展开值
    while fd > fk
        alpha = beta * alpha;
        fd = feval(fun, x0 + alpha * dk);
        fk = feval(fun, x0) + alpha * rho * gk' * dk;
        k = k + 1;
    end
end

