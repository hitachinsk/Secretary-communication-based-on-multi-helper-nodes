NodesNum = [2, 4, 6];
alpha_AB = sqrt(1/2)*(randn(1,8000) +1i*randn(1,8000));
alpha_AE = sqrt(1/2)*(randn(1,8000) +1i*randn(1,8000));
alpha_BE = sqrt(1/2)*(randn(1,8000) +1i*randn(1,8000));
%记录停机出现的次数
outrage_num1 = zeros(1, 3);
outrage_num2 = zeros(1, 3);
outrage_num3 = zeros(1, 3);
outrage_num4 = zeros(1, 3);
%记录停机限制信道
cLimit1 = 2;
cLimit2 = 1;
cLimit3 = 0.5;
cLimit4 = 0.1;

for i = 1:3
    k = NodesNum(i);
    [PxOpt, PyOpt, PzOpt, kxopt, PnOpt, PeOpt] = optCesc(k);
    PxOpt = 10^(PxOpt / 10);
    PyOpt = 10^(PyOpt / 10);
    PzOpt = 10^(PzOpt / 10);
    beta = sqrt(kxopt/2)*(randn(1,k) + 1i*randn(1,k));
    gama_s = vpa(zeros(1, 2));
    yita_iter = 0;
    k_BE = 0;
    sigma_nb2 = 0;
    alpha_AH = sqrt(1/2)*(randn(k,8000) +1i*randn(k,8000));
    alpha_BH = sqrt(1/2)*(randn(k,8000) +1i*randn(k,8000));
    alpha_HB = sqrt(1/2)*(randn(k,8000) +1i*randn(k,8000));
    alpha_HE = sqrt(1/2)*(randn(k,8000) +1i*randn(k,8000));
    for j = 1:8000
        for nfn = 1:k
            gama_s(1) = gama_s(1) + beta(nfn) * alpha_AH(nfn, j) * alpha_HE(nfn, j);
            gama_s(2) = gama_s(2) + beta(nfn) * alpha_AH(nfn, j) * alpha_HB(nfn, j);
            k_BE = k_BE + beta(k) * alpha_BH(nfn, j) * alpha_HE(nfn, j);
            yita_iter = yita_iter + abs(alpha_HE(nfn, j))^2;
            sigma_nb2 = sigma_nb2 + abs(alpha_HB(nfn, j))^2;
        end
        hz = [0; alpha_AE(j)];
        gama = alpha_AB(j) * gama_s(1) - alpha_AE(j) * gama_s(2);
        Hxy = [alpha_AB(j) * alpha_AE(j), alpha_BE(j); gama, k_BE];
        sigma_nb2 = sigma_nb2 * kxopt * PnOpt + PeOpt;
        yita = abs(Hxy(2, 1))^2 * PxOpt + abs(Hxy(2, 2))^2 * PyOpt + yita_iter * kxopt * PnOpt + PeOpt;
        K = [abs(Hxy(1, 1))^2 * PxOpt + abs(Hxy(1, 2))^2 * PyOpt + PeOpt, 0; 0, yita];
        %保密信道
        Csec = vpa(log(1 + abs(alpha_AB(j))^2 * PzOpt/sigma_nb2) - log(det(hz * hz' * PzOpt + K) / det(K)));
        Csec
        if Csec < cLimit1
            outrage_num1(i) = outrage_num1(i) + 1;
        end
        if Csec < cLimit2
            outrage_num2(i) = outrage_num2(i) + 1;
        end
        if Csec < cLimit3
            outrage_num3(i) = outrage_num3(i) + 1;
        end
        if Csec < cLimit4
            outrage_num4(i) = outrage_num4(i) + 1;
        end
        i
        j
    end
end
%计算停机概率
outrage_percent1 = outrage_num1 / 8000;
outrage_percent2 = outrage_num2 / 8000;
outrage_percent3 = outrage_num3 / 8000;
outrage_percent4 = outrage_num4 / 8000;

plot(NodesNum, outrage_percent1, '-rx', NodesNum, outrage_percent2, '-bo',  NodesNum, outrage_percent3, '-g*', NodesNum, outrage_percent4, '-r*');
legend('C outrage = 2', 'C outrage = 1', 'C outrage = 0.5', 'C outrage = 0.1');
title('Outrage Probability in Multiple Relays Scenario'); 
xlabel('Number of helper nodes'); 
ylabel('Outrage Probability');
        