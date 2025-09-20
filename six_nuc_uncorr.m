clear;
% Всё считаем в герцах
% Константы
beta = 5.05; % 1e-24
g = 5.585;
gn=-0.567;
h = 1.054; % 1e-27

% Количество спинов
n_spins = 6;
dim = 2^n_spins;

% Параметры основной пары
sigma1 = 509.94e-6;
sigma2 = 509.94e-6;
r12 = 1.248; % Å
J12 = 16 * 2 * pi; % Гц переводим в цикл. частоту
dJ12 = 0; % Гц

% Параметры дополнительных ядер
% Расстояния до основной пары (в Å)
all_pairs = nchoosek(1:n_spins, 2);
n_pairs = size(all_pairs, 1);

pairs = {[1, 2],[1,3], [1,4], [1,5], [1,6], [2,3], [2,4], [2,5], [2,6], ...
                 [3,4], [3,5], [3,6], [4,5], [4,6], [5,6]};

r13 = 2.539; r14 = 2.470; r15 = 3.779; r16 = 2.734;
r23 = 3.779; r24 = 2.734; r25 = 2.539; r26 = 2.470;
r34 = 3.813; r35 = 6.317; r36 = 4.267;
r45 = 4.267; r46 = 5.060;
r56 = 3.813;
r_mas=[r12 r13 r14 r15 r16 r23 r24 r25 r26 r34 r35 r36 r45 r46 r56];
% Химические сдвиги дополнительных ядер (ppm)
sigma3 = 8.494e-6; sigma4 = 8.373e-6; sigma5 = 8.494e-6; sigma6 = 8.373e-6;

% Скалярные константы связи (Гц)
J13 = 1.67 * 2 * pi; J14 = -0.42 * 2 * pi; J15 = -0.42 * 2 * pi; J16 = 1.67 * 2 * pi;
J23 = -0.42 * 2 * pi; J24 = 1.67 * 2 * pi; J25 = 1.67 * 2 * pi; J26 = -0.42 * 2 * pi;
J34 = 0 * 2 * pi; J35 = 0 * 2 * pi; J36 = 2.11 * 2 * pi;
J45 = 2.11 * 2 * pi; J46 = 0 * 2 * pi;
J56 = 0 * 2 * pi;
J_mas = [J12 J13 J14 J15 J16 J23 J24 J25 J26 J34 J35 J36 J45 J46 J56];
% Параметры расчёта
NB = 1;
B = linspace(5, 5, NB);
B = 1 * 10.^(B);
tau = [9e-11]; 
% Константы для диполь-дипольного взаимодействия
const_HH = 1e6 * g^4 * beta^4 / (h^2);
const_HN = 1e6 * g^2*gn^2 * beta^4 / (h^2);
const_NN = 1e6 * gn^4 * beta^4 / (h^2);

% Создание операторов для каждого спина
up=[0 1; 0 0]; dn=[0 0; 1 0]; z=[0.5 0; 0 -0.5];
for i=1:n_spins  
    Iup{i}=kron(eye(2^(i-1)),kron(up,eye(2^(n_spins-i))));
    Idn{i}=kron(eye(2^(i-1)),kron(dn,eye(2^(n_spins-i))));
    Iz{i}=kron(eye(2^(i-1)),kron(z,eye(2^(n_spins-i))));
%     Ii{i}=Ix{i}+1i*Iy{i};
%     Id{i}=Ix{i}-1i*Iy{i};
end
% Начальная плотность матрицы (только первая пара в синглете)
alpha=[1 0; 0 0];
bita=[0 0; 0 1];
equil = [1 0; 0 1]/2;
singlet=[0 0 0 0;
    0 0.5 -0.5 0;
    0 -0.5 0.5 0;
    0 0 0 0];

ro0 = kron(singlet, kron(equil, kron(equil, kron(equil, equil))));
% Проектор на синглетное состояние первых двух спинов
PS = (eye(dim, dim)-4*Iz{1}*Iz{2}-2*(Iup{1}*Idn{2}+Idn{1}*Iup{2}))/4;
% Основной цикл по времени корреляции
for p = 1:length(tau)
    tau_S = zeros(NB, 1);
    ttau_S = zeros(NB, 1);
    
    for l = 1:NB
        % Гамильтониан Зеемана
        H_zeeman = zeros(dim, dim);
        H_zeeman = H_zeeman - 1e3 * gn * beta * B(l) * (1 - sigma1) / h .* (Iz{1}+Iz{2});
        for k = 3:n_spins
            sigma_k = eval(['sigma' num2str(k)]);            
            H_zeeman = H_zeeman - 1e3 * g * beta * B(l) * (1 - sigma_k) / h .* Iz{k};
        end        
        % Гамильтониан скалярного взаимодействия
        H_J = zeros(dim, dim);                
        
        % Взаимодействия с дополнительными ядрами
        pairs = {[1, 2], [1,3], [1,4], [1,5], [1,6], [2,3], [2,4], [2,5], [2,6], ...
                 [3,4], [3,5], [3,6], [4,5], [4,6], [5,6]};
        
        for pair = pairs
            i = pair{1}(1); j = pair{1}(2);
            J_val = eval(['J' num2str(i) num2str(j)]);
            H_J = H_J + J_val*(Iz{i}*Iz{j} + 0.5*(Iup{i}*Idn{j} + Idn{i}*Iup{j}));
        end
        
        % Полный гамильтониан
        H_total = H_zeeman + H_J;        
        % Диагонализация
        [V, D] = eig(H_total);
        lam = kron(D, eye(dim)) - kron(eye(dim), conj(D));
        U = kron(V, conj(V));
        i_U = kron(conj(inv(V)), inv(V));
        % Релаксационный оператор Редфилда
        Rrf = zeros(dim^2, dim^2);
        % Спектральная плотность            
        Jlam=zeros(dim^2, dim^2);
        for k1 = 1:dim^2
            for k2 = 1:dim^2
                %Jlam(i,k)=2*tau/(1+tau^2*(lam(i,i)-lam(k,k))^2);   %-inf->inf
                Jlam(k1,k2)=1/(1/tau(p)+1i*(lam(k1,k1)-lam(k2,k2))); %0->inf 2xtimes slower
            end    
        end 
        % Диполь-дипольные взаимодействия между всеми парами
        all_pairs = nchoosek(1:n_spins, 2);    
        i_idx = all_pairs(:, 1);
        j_idx = all_pairs(:, 2);        
        parfor idx = 1:size(all_pairs, 1)
            i = i_idx(idx);
            j = j_idx(idx);                                     
            % Операторы диполь-дипольного взаимодействия 
            A_cs2 = Iup{i}*Idn{j}+Idn{i}*Iup{j}-4*Iz{i}*Iz{j};
            A_up = Iz{i}*Iup{j}+Iup{i}*Iz{j};
            A_dn = Iz{i}*Idn{j}+Idn{i}*Iz{j};
            A4 = Iup{i}*Iup{j};
            A5 = Idn{i}*Idn{j};
            % размерность
            A_cs2_m=kron(A_cs2, eye(dim)) - kron(eye(dim), A_cs2');
            A_up_m=kron(A_up, eye(dim)) - kron(eye(dim), A_up');
            A_dn_m=kron(A_dn, eye(dim)) - kron(eye(dim), A_dn');
            A_4_m=kron(A4, eye(dim)) - kron(eye(dim), A4');
            A_5_m=kron(A5, eye(dim)) - kron(eye(dim), A5');                        
            % Расстояние между ядрами
            r_ij = r_mas(idx);
            if (i==1)&&(j==2)
                const_rel=const_NN/r_ij^6;                
            else
                if (i==1)||(i==2)
                    const_rel=const_HN/r_ij^6;                    
                else
                    const_rel=const_HH/r_ij^6;                    
                end
            end            
            % Вклад в релаксационный оператор            
            Rrf = Rrf -const_rel*0.05*A_cs2_m'*U*((U\A_cs2_m*U).*Jlam)*i_U;
            Rrf = Rrf -const_rel*0.3*A_up_m'*U*((U\A_up_m*U).*Jlam)*i_U;
            Rrf = Rrf -const_rel*0.3*A_dn_m'*U*((U\A_dn_m*U).*Jlam)*i_U;
            Rrf = Rrf -const_rel*0.3*A_4_m'*U*((U\A_4_m*U).*Jlam)*i_U;
            Rrf = Rrf -const_rel*0.3*A_5_m'*U*((U\A_5_m*U).*Jlam)*i_U;
        end            
        % Оператор эволюции
        diff_M = -1i*U*lam*i_U+Rrf;        
        % Преобразование начальной матрицы плотности
        rv0 = reshape(ro0, [dim^2, 1]);        
        % Вычисление времени жизни через преобразование Лапласа
        s = 1e-5;
        I0 = eye(dim^2);
        
        rv_s = (I0 * s - diff_M) \ rv0;
        rho_s = reshape(rv_s, [dim, dim]);
        prob_S_s = trace(PS * rho_s);
        tau_S(l) = prob_S_s - 1/(4 * s);
        
        rv_s = (I0 * s - diff_M) \ ((I0 * s - diff_M) \ rv0);
        rho_s = reshape(rv_s, [dim, dim]);
        prob_S_s = trace(PS * rho_s);
        ttau_S(l) = prob_S_s - 1/(4 * s^2);
        disp(l);
    end
    disp(ttau_S ./ tau_S);
    % Визуализация результатов
    hold on;
    plot(B, real(ttau_S ./ tau_S), 'DisplayName', ['\tau_c = ' num2str(tau(p)*1e9) ' ns'], 'LineWidth', 2);
    %сохранение в файл
    to_print=[B; real(ttau_S ./ tau_S)'];
    timestamp = datestr(now, 'yyyy_mm_dd__HH_MM_SS');
    fileID = fopen(['data_' num2str(tau(p)*1e9) '_' timestamp '.txt'],'w');
    fprintf(fileID,'%6.6f %4.4f\r\n',to_print);
    fclose(fileID);
end

xlabel('Магнитное поле, Гс');
ylabel('\tau_S, с');
title('Время жизни синглетного состояния с дополнительными ядрами');
legend;
grid on;
set(gca, 'XScale', 'log');
%set(gca, 'YScale', 'log');
