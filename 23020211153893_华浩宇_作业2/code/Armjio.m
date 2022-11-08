function [alpha] = Armjio(fun, grid, x0, dk)
    beta = 0.333; 		% ���� alpha �ĵ���ϵ����С�� 1
    rho = 1e-3; 		% ̩��չ��ʽ����ϵ����0 < rho < 1/2
    alpha = 1; 			% ��ʼ����Ϊ 1
    k = 0; 				% ͳ�Ƶ�������
    gk = feval(grid, x0);	% x0�����ݶ�ֵ
    fd = feval(fun, x0 + alpha * dk); 	% ��������һ�������㴦��Ŀ�꺯��ֵ
    fk = feval(fun, x0) + alpha * rho * gk' * dk; 	% ��������һ�������㴦��̩��չ��ֵ
    while fd > fk
        alpha = beta * alpha;
        fd = feval(fun, x0 + alpha * dk);
        fk = feval(fun, x0) + alpha * rho * gk' * dk;
        k = k + 1;
    end
end

