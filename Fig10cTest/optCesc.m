function [PxRes, PyRes, PzRes, kxres, Pnres, Peres] = optCesc(NumofNodes)
alpha_AB = sqrt(1/2)*(randn(1,1000) +1i*randn(1,1000));
alpha_AE = sqrt(1/2)*(randn(1,1000) +1i*randn(1,1000));
alpha_BE = sqrt(1/2)*(randn(1,1000) +1i*randn(1,1000));
alpha_AH = sqrt(1/2)*(randn(NumofNodes,1000) +1i*randn(NumofNodes,1000));
alpha_BH = sqrt(1/2)*(randn(NumofNodes,1000) +1i*randn(NumofNodes,1000));
alpha_HB = sqrt(1/2)*(randn(NumofNodes,1000) +1i*randn(NumofNodes,1000));
alpha_HE = sqrt(1/2)*(randn(NumofNodes,1000) +1i*randn(NumofNodes,1000));

syms kexi Px1 Py1 Pz1 Pn1 Pe1;
beta = sqrt(kexi/2)*(randn(1,NumofNodes) + 1i*randn(1,NumofNodes));

Csec = vpa(zeros(1, 1000));
for j = 1:1000
    gama_s = vpa(zeros(1, 2));
    yita_iter = 0;
    k_BE = 0;
    sigma_nb2 = 0;
    for i = 1:NumofNodes
        gama_s(1) = gama_s(1) + beta(i) * alpha_AH(i, j) * alpha_HE(i, j);
        gama_s(2) = gama_s(2) + beta(i) * alpha_AH(i, j) * alpha_HB(i, j);
        k_BE = k_BE + beta(i) * alpha_BH(i, j) * alpha_HE(i, j);
        yita_iter = yita_iter + abs(alpha_HE(i, j))^2;
        sigma_nb2 = sigma_nb2 + abs(alpha_HB(i, j))^2;
    end
    hz = [0; alpha_AE(j)];
    gama = alpha_AB(j) * gama_s(1) - alpha_AE(j) * gama_s(2);
    Hxy = [alpha_AB(j) * alpha_AE(j), alpha_BE(j); gama, k_BE];
    sigma_nb2 = sigma_nb2 * kexi * Pn1 + Pe1;
    yita = abs(Hxy(2, 1))^2 * Px1 + abs(Hxy(2, 2))^2 * Py1 + yita_iter * kexi * Pn1 + Pe1;
    K = [abs(Hxy(1, 1))^2 * Px1 + abs(Hxy(1, 2))^2 * Py1 + Pe1, 0; 0, yita];
    %±£ÃÜÐÅµÀ
    Csec(j) = log(1 + abs(alpha_AB(j))^2 * Pz1/sigma_nb2) - log(det(hz * hz' * Pz1 + K) / det(K));
end
meanCsec = mean(Csec);
P0 = 50;
%ÓÅ»¯º¯Êýf2 and capacity
f2 = (1 + 2 * NumofNodes * kexi) * Px1 + (1 + NumofNodes * kexi) * Py1 + Pz1 + NumofNodes * kexi * Pn1;
P0Test = 10^(P0 / 10);
PeTest = P0Test / 100;
meanCsecTest = subs(meanCsec, [Pe1, Pn1], [PeTest, 0]);
f2Test = subs(f2, Pn1, 0);
maxCsec = 0;
Pnres = 0;
Peres = PeTest;
for PxTestDb = 1 :3 : 15
    for PyTestDb = 1 : 3 : 15
        for PzTestDb = 20 : 3 : P0 - PxTestDb - PyTestDb
            PxTest = 10^(PxTestDb / 10);
            PyTest = 10^(PyTestDb / 10);
            PzTest = 10^(PzTestDb / 10);
            f2Test1 = subs(f2Test, [Px1, Py1, Pz1], [PxTest, PyTest, PzTest]);
            f2Test1 = f2Test1 - 10^(P0 / 10);
            kexiTest1 = solve(f2Test1, kexi);  
            if kexiTest1 > 0
                res = subs(meanCsecTest, [Px1, Py1, Pz1, kexi], [PxTest, PyTest, PzTest, kexiTest1]);
                res = abs(vpa(res));
                res
                if res > maxCsec
                    PxRes = PxTestDb;
                    PyRes = PyTestDb;
                    PzRes = PzTestDb;
                    kxres = kexiTest1;
                    maxCsec = res;
                    maxCsec
                end
            end
        end
    end
end
maxCsec
end
