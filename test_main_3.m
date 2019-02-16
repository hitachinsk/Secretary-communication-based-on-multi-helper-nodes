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
Csec_t2 = vpa(zeros(1, 1000));
yita = 0;
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
    Csec_t2(j) = log(det(hz * hz' * Pz + K) / det(K));
    Csec(j) = log(1 + abs(alpha_AB(j))^2 * Pz/sigma_nb2) - Csec_t2(j);
end
%计算保密信道均值
mean_Csec = mean(Csec);
mean_Csec
mean_Csec_t2 = mean(Csec_t2);

%优化函数f2
f2 = (1 + 2 * 10 * kexi) * Px + (1 + 10 * kexi) * Py + Pz + 10 * kexi * Pn;
SNRatE = [-10, 0, 10, 20, 30, 40];
P0 = 50;
mc_res = zeros(1, 6);
mc_res2 = zeros(1, 6);
C_B1 = log(1 + 10^(30/10)) - mean_Csec_t2;
C_B2 = log(1 + 10^(20/10)) - mean_Csec_t2;
Px_test_db = 2;
Py_test_db = 2;
Pz_test_db = 45;
Px_test = 10^(Px_test_db / 10);
Py_test = 10^(Py_test_db / 10);
Pz_test = 10^(Pz_test_db / 10);
for i = 1:6
    SNRE = 10^(SNRatE(i) / 10);
    Pe_test = 10^(P0 / 10) / SNRE;
    Pn_1 = (Pz_test / 1000 - Pe_test) / (10 * kexi);
    Pn_2 = (Pz_test / 100 - Pe_test) / (10 * kexi);
    f2_test_1 = subs(f2, [Px, Py, Pz, Pn], [Px_test, Py_test, Pz_test, Pn_1]);
    f2_test_1 = f2_test_1 - 10^(P0 / 10);
    kexi_1 = solve(f2_test_1, kexi);
    Pn_1_test = subs(Pn_1, kexi, kexi_1);
    f2_test_2 = subs(f2, [Px, Py, Pz, Pn], [Px_test, Py_test, Pz_test, Pn_2]);
    f2_test_2 = f2_test_2 - 10^(P0 / 10);
    kexi_2 = solve(f2_test_2, kexi);
    Pn_2_test = subs(Pn_2, kexi, kexi_2);
    C_B_1 = abs(subs(C_B1, [Px, Py, Pz, Pn, Pe, kexi], [Px_test, Py_test, Pz_test, Pn_1_test, Pe_test, kexi_1]));
    C_B_2 = abs(subs(C_B2, [Px, Py, Pz, Pn, Pe, kexi], [Px_test, Py_test, Pz_test, Pn_2_test, Pe_test, kexi_2]));
    mc_res(i) = C_B_1;
    mc_res2(i) = C_B_2;
end
%画图
kg = [log(1 + 10^(30/10)), log(1 + 10^(30/10)), log(1 + 10^(30/10)), log(1 + 10^(30/10)), log(1 + 10^(30/10)), log(1 + 10^(30/10))];
ed = [log(1 + 10^(20/10)), log(1 + 10^(20/10)), log(1 + 10^(20/10)), log(1 + 10^(20/10)), log(1 + 10^(20/10)), log(1 + 10^(20/10))];
plot(SNRatE, kg, 'b*--', SNRatE, mc_res, 'b*-', SNRatE, ed, 'r<--', SNRatE, mc_res2, 'r<-'); 
legend('E(capacity[30db])', 'C_sec_h[30db]', 'E(capacity[20db]', 'C_sec_h[20db]'); 
xlabel('SNR at eavesdropper'); 
ylabel('C_sec_h(nats/symbol)');
title('Fig5:C_sec_h:variation with distance');