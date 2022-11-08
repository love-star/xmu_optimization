syms x1 x2;
%f = xs^2+2*ys^2-2*xs*ys + 2*ys + 2;
f = 100 * (x2-x1.^2).^2 + (1 - x1).^2;
%f = (1-x1)^2 + 2*(x2 - x1^2)^2;
x = {x1, x2};

% initial value
x0 = [-2 2];
% tolerance
epsilon = 1e-5;
%% call conjungate gradient method
show_detail = true;
[bestf, bestx, count] = conjungate_gradient(f, x, x0, epsilon, show_detail);
% print result
fprintf('bestx = %s, bestf = %f, count = %d\n', num2str(bestx), bestf, count);