%下一步用蒙特卡洛方法求统计平均，也就是将产生的1000个点逐一迭代做平均
%这里假定所有的功率以分贝表示并且都是整数形式的分贝
%下面的第一步就是将优化暴力搜索的函数找到
%The num of simulation is set to 1000.
%The produce of Raylay distribution(using Gaussian distribution)
%The expectation is 0, while varition is 1, so E[|alpha_xy|^2] = 1
alpha_AB = sqrt(1/2)*(randn(1,1000) +1i*randn(1,1000));
alpha_AE = sqrt(1/2)*(randn(1,1000) +1i*randn(1,1000));
alpha_BE = sqrt(1/2)*(randn(1,1000) +1i*randn(1,1000));
alpha_AH = sqrt(1/2)*(randn(10,1000) +1i*randn(10,1000));
alpha_BH = sqrt(1/2)*(randn(10,1000) +1i*randn(10,1000));
alpha_HB = sqrt(1/2)*(randn(10,1000) +1i*randn(10,1000));
alpha_HE = sqrt(1/2)*(randn(10,1000) +1i*randn(10,1000));
 
%definition of symbol variation(for optimization)
%new simulation is 4 and 10(nums of helper nodes)
syms kexi Px Py Pz Pn Pe;
%num of helper nodes is 10
beta = sqrt(kexi/2)*(randn(1,10) + 1i*randn(1,10));
%节点数为10的保密信道
Csec = vpa(zeros(1, 1000));
%num of helper nodes is 4
beta_4 = sqrt(kexi/2)*(randn(1,10) + 1i*randn(1,10));
Csec_4 = vpa(zeros(1, 1000));
%蒙特卡洛仿真1000次
for j = 1:1000
    %num of helper nodes is 10
    gama_s = vpa(zeros(1, 2));
    yita_iter = 0;
    k_BE = 0;
    sigma_nb2 = 0;
    for i = 1:10
        gama_s(1) = gama_s(1) + beta(i) * alpha_AH(i, j) * alpha_HE(i, j);
        gama_s(2) = gama_s(2) + beta(i) * alpha_AH(i, j) * alpha_HB(i, j);
        k_BE = k_BE + beta(i) * alpha_BH(i, j) * alpha_HE(i, j);
        yita_iter = yita_iter + abs(alpha_HE(i, j))^2;
        sigma_nb2 = sigma_nb2 + abs(alpha_HB(i, j))^2;
    end
    hz = [0; alpha_AE(j)];
    gama = alpha_AB(j) * gama_s(1) - alpha_AE(j) * gama_s(2);
    Hxy = [alpha_AB(j) * alpha_AE(j), alpha_BE(j); gama, k_BE];
    sigma_nb2 = sigma_nb2 * kexi * Pn + Pe;
    yita = abs(Hxy(2, 1))^2 * Px + abs(Hxy(2, 2))^2 * Py + yita_iter * kexi * Pn + Pe;
    K = [abs(Hxy(1, 1))^2 * Px + abs(Hxy(1, 2))^2 * Py + Pe, 0; 0, yita];
    %保密信道
    Csec(j) = log(1 + abs(alpha_AB(j))^2 * Pz/sigma_nb2) - log(det(hz * hz' * Pz + K) / det(K));
    
    %num of helper nodes is 4
    gama_s_4 = vpa(zeros(1, 2));
    yita_iter_4 = 0;
    k_BE_4 = 0;
    sigma_nb2_4 = 0;
    for k = 1:4
        gama_s_4(1) = gama_s_4(1) + beta_4(k) * alpha_AH(k, j) * alpha_HE(k, j);
        gama_s_4(2) = gama_s_4(2) + beta_4(k) * alpha_AH(k, j) * alpha_HB(k, j);
        k_BE_4 = k_BE_4 + beta_4(k) * alpha_BH(k, j) * alpha_HE(k, j);
        yita_iter_4 = yita_iter_4 + abs(alpha_HE(k, j))^2;
        sigma_nb2_4 = sigma_nb2_4 + abs(alpha_HB(k, j))^2;
    end
    hz_4 = [0; alpha_AE(j)];
    gama_4 = alpha_AB(j) * gama_s_4(1) - alpha_AE(j) * gama_s_4(2);
    Hxy_4 = [alpha_AB(j) * alpha_AE(j), alpha_BE(j); gama_4, k_BE_4];
    sigma_nb2_4 = sigma_nb2_4 * kexi * Pn + Pe;
    yita_4 = abs(Hxy_4(2, 1))^2 * Px + abs(Hxy_4(2, 2))^2 * Py + yita_iter_4 * kexi * Pn + Pe;
    K_4 = [abs(Hxy_4(1, 1))^2 * Px + abs(Hxy_4(1, 2))^2 * Py + Pe, 0; 0, yita_4];
    %保密信道
    Csec_4(j) = log(1 + abs(alpha_AB(j))^2 * Pz/sigma_nb2_4) - log(det(hz_4 * hz_4' * Pz + K_4) / det(K_4));
end
%计算保密信道均值
mean_Csec = mean(Csec);
mean_Csec
mean_Csec_4 = mean(Csec_4);
mean_Csec_4
P0 = 50;
%优化函数f2 and capacity
f2 = (1 + 2 * 10 * kexi) * Px + (1 + 10 * kexi) * Py + Pz + 10 * kexi * Pn;
f2_4 = (1 + 2 * 4 * kexi) * Px + (1 + 4 * kexi) * Py + Pz + 4 * kexi * Pn;
%search
SNR_P0_Pe = [20, 22, 24, 26, 28, 30];
capa_res = zeros(1, 6);
mc_res = zeros(1, 6);
mc_res2 = zeros(1, 6);
for i = 1:6
    SNR = 10^(SNR_P0_Pe(i) / 10);
    P0_test = 10^(P0 / 10);
    Pe_test = P0_test / SNR;
    Pn_test = 0;
    mc = subs(mean_Csec, [Pn, Pe], [Pn_test, Pe_test]);
    mc4 = subs(mean_Csec_4, [Pn, Pe], [Pn_test, Pe_test]);
    f2_test = subs(f2, Pn, Pn_test);
    f2_4_test = subs(f2_4, Pn, Pn_test);
    max_Csec = 0;
    max_Csec4 = 0;
    capa_res(i) = log(1 + SNR);
    for Px_test_db = 1 :3 : 15
        for Py_test_db = 1 : 3 : 15
            for Pz_test_db = 20 : 3 : P0 - Px_test_db - Py_test_db
                Px_test = 10^(Px_test_db / 10);
                Py_test = 10^(Py_test_db / 10);
                Pz_test = 10^(Pz_test_db / 10);
                f2_test_1 = subs(f2_test, [Px, Py, Pz], [Px_test, Py_test, Pz_test]);
                f2_test_1 = f2_test_1 - 10^(P0 / 10);
                kexi_test_1 = solve(f2_test_1, kexi); 
                f2_4_test_1 = subs(f2_4_test, [Px, Py, Pz], [Px_test, Py_test, Pz_test]);
                f2_4_test_1 = f2_4_test_1 - 10^(P0 / 10);
                kexi_4_test_1 = solve(f2_4_test_1, kexi); 
                if kexi_test_1 > 0
                    res = subs(mc, [Px, Py, Pz, kexi], [Px_test, Py_test, Pz_test, kexi_test_1]);
                    res = abs(vpa(res));
                    res
                    if res > max_Csec
                        max_Csec = res;
                        max_Csec
                    end
                end
                if kexi_4_test_1 > 0
                    res4 = subs(mc4, [Px, Py, Pz, kexi], [Px_test, Py_test, Pz_test, kexi_4_test_1]);
                    res4 = abs(vpa(res4));
                    if res4 > max_Csec4
                        max_Csec4 = res4;
                    end
                end
            end
        end
        max_Csec
    end
    mc_res(i) = max_Csec;
    mc_res2(i) = max_Csec4;
end
%画图
plot(SNR_P0_Pe, capa_res, '-rx', SNR_P0_Pe, mc_res, '-bo',  SNR_P0_Pe, mc_res2, '-g*');
legend('E(capacity)', 'E(c):10 nodes', 'E(c):4 nodes');
title('E(C):variation with P0'); 
xlabel('SNR:P0/Pe'); 
ylabel('E[c](nats/symbol)');


