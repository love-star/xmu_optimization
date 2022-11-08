function [f, xk, k] = BFGS(x0, fun, grid, eps, kmax)
  k = 0;
  n = length(x0);
  H0 = eye(n); % ��ʼѡȡ��λ����ΪHessen�������Ľ�����
  Hk = H0;
  xk = x0;
  gk = feval(grid, xk);
  while k <= kmax
      if norm(gk) < eps
          break;
      end
      dk = -Hk * gk; % ��ţ���½�����
      alpha = Armjio(fun, grid, xk, dk);
      x_ = xk; % x_ ������һ��������
      xk = x_ + alpha * dk; % ���� xk
      gk_ = gk; % gk_ ������һ������ݶ�ֵ
      gk = feval(grid, xk); % ���� gk
      sk = xk - x_; % �� xk - x_ Ϊ sk
      yk = gk - gk_; % �� gk - gk_ Ϊ yk
      if sk' * yk > 0
          v = yk' * sk;
          % BFGS��ʽ
          Hk = Hk + (1 + (yk' * Hk * yk) / v) * (sk * sk') / v - (sk * yk' * Hk + Hk * yk * sk') / v;
      end
      k = k + 1;
  end
  f = feval(fun, xk);
end
