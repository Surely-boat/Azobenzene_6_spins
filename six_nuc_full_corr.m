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
r21=r12;
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
r31 = 2.539; r41 = 2.470; r51 = 3.779; r161 = 2.734;
r32 = 3.779; r42 = 2.734; r52 = 2.539; r62 = 2.470;
r43 = 3.813; r53 = 6.317; r63 = 4.267;
r54 = 4.267; r64 = 5.060;
r65 = 3.813;
r_mas=[r12 r13 r14 r15 r16 r23 r24 r25 r26 r34 r35 r36 r45 r46 r56];
% Углы между связями HN-HN прилегающие
phi_3114=99.2*pi/180;    phi_3116=108*pi/180;     phi_4116=152.9*pi/180;   phi_1442=27.1*pi/180;   
phi_5226=phi_3114;        phi_4225=phi_3116;        phi_4226=phi_4116;        phi_1662=phi_1442;
%NN-HN
phi_2113=172.5*pi/180;   phi_2114=88.3*pi/180;    phi_2116=64.6*pi/180;
phi_1225=172.5*pi/180;   phi_1226=88.3*pi/180;    phi_1224=64.6*pi/180;
%HN-HN нет взаимного атома
phi_1324=72.02*pi/180;  phi_1326=80.85*pi/180;  phi_1325=pi;
phi_1625=72.02*pi/180;  phi_1624=pi;
phi_2514=80.85*pi/180; phi_2614=pi;
% Химические сдвиги дополнительных ядер (ppm)
sigma3 = 8.494e-6; sigma4 = 8.373e-6; sigma5 = 8.494e-6; sigma6 = 8.373e-6;

% Скалярные константы связи (Гц)
J13 = 1.7 * 2 * pi; J14 = -0.4 * 2 * pi; J15 = -0.42 * 2 * pi; J16 = 1.67 * 2 * pi;
J23 = -0.42 * 2 * pi; J24 = 1.67 * 2 * pi; J25 = 1.67 * 2 * pi; J26 = -0.42 * 2 * pi;
J34 = 0 * 2 * pi; J35 = 0 * 2 * pi; J36 = 2.11 * 2 * pi;
J45 = 2.11 * 2 * pi; J46 = 0 * 2 * pi;
J56 = 0 * 2 * pi;
J_mas = [J12 J13 J14 J15 J16 J23 J24 J25 J26 J34 J35 J36 J45 J46 J56];
% Параметры расчёта
NB = 15;
B = linspace(-3, -2, NB);
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
        % Учёт корреляции
        angles = {[3, 1, 1, 4], [3, 1, 1, 6], [4, 1, 1, 6], [1, 4, 4, 2], [5, 2, 2, 6], [4, 2, 2, 5], [4, 2, 2, 6], [1, 6, 6, 2], ...
            [1, 3, 2, 4], [1, 3, 2, 6], [1, 3, 2, 5], [1, 6, 2, 5], [1, 6, 2, 4], [2, 5, 1, 4], [2, 6, 1, 4], ...
        [2, 1, 1, 3], [2, 1, 1, 4], [2, 1, 1, 6], [1, 2, 2, 5], [1, 2, 2, 6], [1, 2, 2, 4]};%NH-NN
        r_1_mas=[r13 r13 r14 r14 r25 r24 r24 r16 r13 r13 r13 r16 r16 r25 r26 r12 r12 r12 r12 r12 r12];
        r_2_mas=[r14 r16 r16 r24 r26 r25 r26 r26 r24 r26 r25 r25 r24 r14 r14 r13 r14 r16 r25 r26 r24];
        phi_mas=[phi_3114 phi_3116 phi_4116 phi_1442 phi_5226 phi_4225 phi_4226 phi_1662 phi_1324 phi_1326 phi_1325 phi_1625 ...
            phi_1624 phi_2514 phi_2614 phi_2113 phi_2114 phi_2116 phi_1225 phi_1226 phi_1224];
        i_idx=zeros(1, length(angles));
        j_idx=zeros(1, length(angles));
        k_idx=zeros(1, length(angles));
        m_idx=zeros(1, length(angles));
        for o=1:length(angles)
            i_idx(o)=angles{o}(1);
            j_idx(o)=angles{o}(2);
            k_idx(o)=angles{o}(3);
            m_idx(o)=angles{o}(4);
        end
        parfor idx = 1:length(angles)
            i = i_idx(idx);  j = j_idx(idx);   k = k_idx(idx);  m = m_idx(idx);           
            NN_switch=0;
            %first pair
            A_cs2_1 = Iup{j}*Idn{i}+Idn{j}*Iup{i}-4*Iz{j}*Iz{i};
            A_up_1 = Iz{j}*Iup{i}+Iup{j}*Iz{i};
            A_dn_1 = Iz{j}*Idn{i}+Idn{j}*Iz{i};
            A4_1 = Iup{j}*Iup{i};
            A5_1 = Idn{j}*Idn{i};
            % dimentions
            A_cs2_m_1=kron(A_cs2_1, eye(dim)) - kron(eye(dim), A_cs2_1');
            A_up_m_1=kron(A_up_1, eye(dim)) - kron(eye(dim), A_up_1');
            A_dn_m_1=kron(A_dn_1, eye(dim)) - kron(eye(dim), A_dn_1');
            A_4_m_1=kron(A4_1, eye(dim)) - kron(eye(dim), A4_1');
            A_5_m_1=kron(A5_1, eye(dim)) - kron(eye(dim), A5_1');
            %second pair
            A_cs2_2 = Iup{m}*Idn{k}+Idn{m}*Iup{k}-4*Iz{m}*Iz{k};
            A_up_2 = Iz{m}*Iup{k}+Iup{m}*Iz{k};
            A_dn_2 = Iz{m}*Idn{k}+Idn{m}*Iz{k};
            A4_2 = Iup{m}*Iup{k};
            A5_2 = Idn{m}*Idn{k};
            % dimentions
            A_cs2_m_2=kron(A_cs2_2, eye(dim)) - kron(eye(dim), A_cs2_2');
            A_up_m_2=kron(A_up_2, eye(dim)) - kron(eye(dim), A_up_2');
            A_dn_m_2=kron(A_dn_2, eye(dim)) - kron(eye(dim), A_dn_2');
            A_4_m_2=kron(A4_2, eye(dim)) - kron(eye(dim), A4_2');
            A_5_m_2=kron(A5_2, eye(dim)) - kron(eye(dim), A5_2');
            % correlation constant
            r_1 = r_1_mas(idx);
            r_2 = r_2_mas(idx);
            phi = phi_mas(idx);  
            if (i==2)&&(j==1)&&(k==1)&&(m==3)
                NN_switch=1;                
            end            
            if (NN_switch==0)
                const_rel=(1+3*cos(2*phi))*const_HN/(r_1^3*r_2^3);
            else
                const_rel=(1+3*cos(2*phi))* const_HN*(gn/g)/(r_1^3*r_2^3);
            end
            % Вклад в релаксационный оператор
            Rrf = Rrf -const_rel*(1/80)*A_cs2_m_1'*U*((U\A_cs2_m_2*U).*Jlam)*i_U;
            Rrf = Rrf -const_rel*(3/40)*A_up_m_1'*U*((U\A_up_m_2*U).*Jlam)*i_U;
            Rrf = Rrf -const_rel*(3/40)*A_dn_m_1'*U*((U\A_dn_m_2*U).*Jlam)*i_U;
            Rrf = Rrf -const_rel*(3/40)*A_4_m_1'*U*((U\A_4_m_2*U).*Jlam)*i_U;
            Rrf = Rrf -const_rel*(3/40)*A_5_m_1'*U*((U\A_5_m_2*U).*Jlam)*i_U;
            
            Rrf = Rrf -const_rel*(1/80)*A_cs2_m_2'*U*((U\A_cs2_m_1*U).*Jlam)*i_U;
            Rrf = Rrf -const_rel*(3/40)*A_up_m_2'*U*((U\A_up_m_1*U).*Jlam)*i_U;
            Rrf = Rrf -const_rel*(3/40)*A_dn_m_2'*U*((U\A_dn_m_1*U).*Jlam)*i_U;
            Rrf = Rrf -const_rel*(3/40)*A_4_m_2'*U*((U\A_4_m_1*U).*Jlam)*i_U;
            Rrf = Rrf -const_rel*(3/40)*A_5_m_2'*U*((U\A_5_m_1*U).*Jlam)*i_U;            
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
    fileID = fopen(['data/data_' num2str(tau(p)*1e9) '_' timestamp '.txt'],'w');
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
