function [f, xk, k] = BFGS(x0, fun, grid, eps, kmax)
  k = 0;
  n = length(x0);
  H0 = eye(n); % 初始选取单位阵作为Hessen矩阵的逆的近似阵
  Hk = H0;
  xk = x0;
  gk = feval(grid, xk);
  while k <= kmax
      if norm(gk) < eps
          break;
      end
      dk = -Hk * gk; % 拟牛顿下降方向
      alpha = Armjio(fun, grid, xk, dk);
      x_ = xk; % x_ 保存上一个点坐标
      xk = x_ + alpha * dk; % 更新 xk
      gk_ = gk; % gk_ 保存上一个点的梯度值
      gk = feval(grid, xk); % 更新 gk
      sk = xk - x_; % 记 xk - x_ 为 sk
      yk = gk - gk_; % 记 gk - gk_ 为 yk
      if sk' * yk > 0
          v = yk' * sk;
          % BFGS公式
          Hk = Hk + (1 + (yk' * Hk * yk) / v) * (sk * sk') / v - (sk * yk' * Hk + Hk * yk * sk') / v;
      end
      k = k + 1;
  end
  f = feval(fun, xk);
end
