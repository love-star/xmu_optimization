function [fv, bestx, iter_num] = conjungate_gradient(f, x, x0, epsilon, show_detail)

syms lambdas 

n = length(x);

nf = cell(1, n); 
for i = 1 : n
    nf{i} = diff(f, x{i});
end

nfv = subs(nf, x, x0);

nfv_pre = nfv;
count = 0;
k = 0;
xv = x0;
d = - nfv; 

if show_detail
    fprintf('Initial:\n');
    fprintf('f = %s, x0 = %s, epsilon = %f\n\n', char(f), num2str(x0), epsilon);
end

while (norm(nfv) > epsilon)
    xv = xv+lambdas*d;
    phi = subs(f, x, xv);
    nphi = diff(phi);
    lambda = solve(nphi);
    lambda = double(lambda);  
    if length(lambda) > 1
        lambda = lambda(abs(imag(lambda)) < 1e-5);
        lambda = lambda(lambda > 0);
        lambda = lambda(1);
    end
    if lambda < 1e-5
        break;
    end

    xv = subs(xv, lambdas, lambda); 
    xv = double(xv);
    nfv = subs(nf, x, xv);   
    count = count + 1;
    k = k + 1; 
    alpha = sum(nfv(:).*nfv(:)) / sum(nfv_pre(:).*nfv_pre(:));

    if show_detail
        fprintf('epoch: %d\n', count);
        fprintf('x(%d) = %s, lambda = %f\n', count, num2str(xv), lambda);
        fprintf('nf(x) = %s, norm(nf) = %f\n', num2str(double(nfv)), norm(double(nfv)));
        fprintf('d = %s, alpha = %f\n', num2str(double(d)), double(alpha));
        fprintf('\n');
    end

    d = -nfv + alpha .* d;
    nfv_pre = nfv;
    if k >= 6
        k = 0;
        d = - nfv;
    end
end % while

fv = double(subs(f, x, xv));
bestx = double(xv);
iter_num = count;

end