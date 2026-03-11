clear;
% Всё считаем в герцах
%% Константы
beta = 5.05; % 1e-24
g = 5.585;
gn=-0.567;
h = 1.054; % 1e-27
% Количество спинов
n_spins = 6;
dim = 2^n_spins;
%% Параметры молекулы азобензола
%{
    3   5
1     
  2
    4   6
%} 

% Химические (ppm)
sigmaXX=-789e-6;
sigmaYY=-146e-6;
sigmaZZ=136e-6;
sigma1 = 509.94e-6;
sigma2 = 509.94e-6;
sigma3 = 7.925e-6; sigma4 = 7.925e-6; sigma5 = 7.591e-6; sigma6 = 7.591e-6;
d_sig_34=0.2;
d_sig_56=0.1;
sigma_mas=[sigma1 sigma2 sigma3 sigma4 sigma5 sigma6];

r12 = 1.248; % Å
J12 = 16 * 2 * pi; % Гц переводим в цикл. частоту
dJ12 = 0; % Гц

% Параметры дополнительных ядер
% Расстояния до основной пары (в Å)
r13 = 2.495; r14 = 3.787; r15 = 4.818; r16 = 5.611;
r23 = 2.748; r24 = 2.544; r25 = 4.624; r26 = 4.520;
r34 = 4.279; r35 = 2.475; r36 = 4.956;
r45 = 4.945; r46 = 2.488;
r56 = 4.289;
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

dJ1_34=2;
dJ1_56=0.5;
dJ2_34=0;
dJ2_56=0;
% Скалярные константы связи (Гц)^M
J13 = (-0.42) * 2 * pi; J14 = (-0.42) * 2 * pi; J15 = (0.16) * 2 * pi; J16 = (0.16) * 2 * pi;
J23 = (1.67) * 2 * pi; J24 = (1.67) * 2 * pi; J25 = (0.2) * 2 * pi; J26 = (0.2) * 2 * pi;

J34 = 2.11 * 2 * pi; J35 = 7.96 * 2 * pi; J36 = 0.58 * 2 * pi;
J45 = 0.58 * 2 * pi; J46 = 7.96 * 2 * pi;
J56 = 1.53 * 2 * pi;
J_mas = [J12 J13 J14 J15 J16 J23 J24 J25 J26 J34 J35 J36 J45 J46 J56];
%параметры анизотроопии
psi=-37*pi/180; % угол между XX и 12
phi_12=0;
% отразил рисунок в ворде, чтобы была картинка, как в авогадро
phi_13=-7*pi/180;
phi_16=-114.9*pi/180;
phi_14=92.2*pi/180;
%phi_15=175.3*pi/180;
%phi_23=-4.7*pi/180;
phi_26=-87.8*pi/180;
phi_24=65.1*pi/180;
phi_25=173*pi/180;
%% Параметры расчёта
NB = 30;
B = linspace(3.3, 5.3, NB);
B = 1 * 10.^(B);
tau = [6e-11]; 
tau_flip=10e-2;% 0.1+H no picks
flip_flag=1;
D_D_flag = 1;
D_D_corr_flag=0;
CSA_D_corr_flag=0;
CSA_flag=1;
read_from_file_flag=0;
H_relax_flag=1;
T1_azo_flag=0;
kin_flag=0;
kin_1_flag=0;
if (kin_1_flag==1)
    NB=1;
    B = 2*10.^[3];
end
%% Создание операторов для каждого спина
up=[0 1; 0 0]; dn=[0 0; 1 0]; z=[0.5 0; 0 -0.5];
for i=1:n_spins  
    Iup{i}=kron(eye(2^(i-1)),kron(up,eye(2^(n_spins-i))));
    Idn{i}=kron(eye(2^(i-1)),kron(dn,eye(2^(n_spins-i))));
    Iz{i}=kron(eye(2^(i-1)),kron(z,eye(2^(n_spins-i))));
%     Ii{i}=Ix{i}+1i*Iy{i};
%     Id{i}=Ix{i}-1i*Iy{i};
end
% Начальная матрица плотности (только первая пара в синглете)
alpha=[1 0; 0 0];
bita=[0 0; 0 1];
equil = [1 0; 0 1]/2;
singlet=[0 0 0 0;
    0 0.5 -0.5 0;
    0 -0.5 0.5 0;
    0 0 0 0];

ro0 = kron(singlet, kron(equil, kron(equil, kron(equil, equil))));
%ro0 = kron(equil, kron(equil, kron(equil, kron(alpha, kron(equil, equil)))));
if (T1_azo_flag==1)
    ro0 = kron(kron(alpha, alpha), kron(equil, kron(equil, kron(equil, equil))));
end
% Проектор на синглетное состояние первых двух спинов
PS = (eye(dim, dim)-4*Iz{1}*Iz{2}-2*(Iup{1}*Idn{2}+Idn{1}*Iup{2}))/4; %Нужно поделить на 16, чтобы получить синглет+равновесие водородов
PT_p = (eye(dim, dim)+2*Iz{1}+2*Iz{2}+4*Iz{1}*Iz{2})/4;
PT_0 = (eye(dim, dim)-4*Iz{1}*Iz{2}+2*(Iup{1}*Idn{2}+Idn{1}*Iup{2}))/4;
PT_m = (eye(dim, dim)-2*Iz{1}-2*Iz{2}+4*Iz{1}*Iz{2})/4;
if (T1_azo_flag==1)
    PS = Iz{1}+Iz{2};
end
%% Собственная релаксация протонов
if (H_relax_flag==1)
    T1 = [0 0 10 10 10 10];
    T2 = T1;
    %T2 = [0 0 10 5 5 5];
    E=eye(dim, dim);
    R_fast_1 = zeros(dim^2, dim^2);
    R_fast_2 = zeros(dim^2, dim^2);
    for i=3:6
        R_fast_1=R_fast_1+(kron(Iup{i},Iup{i})+kron(Idn{i},Idn{i})-2*kron(Iz{i},Iz{i})-0.5*kron(E, E))/(2*T1(i));
        R_fast_2=R_fast_2+(2*kron(Iz{i},Iz{i})-0.5*kron(E, E))/(T2(i));    
    end
    R_fast=R_fast_1+R_fast_2;
end
%% Основные циклы по времени корреляции и магнитному полю
for p = 1:length(tau)
    tau_S = zeros(NB, 1);
    ttau_S = zeros(NB, 1);
    T1_kin_1 = zeros(NB, 1);
    T1_kin_2 = zeros(NB, 1);
    p1_kin = zeros(NB, 1);
    p2_kin = zeros(NB, 1);
    
    for l = 1:NB
        %% Гамильтониан Зеемана        
        H_zeeman = zeros(dim, dim);
        H_zeeman = H_zeeman - 1e3 * gn * beta * B(l) * (1 - sigma_mas(1)) / h .* (Iz{1}+Iz{2});
        for k = 3:n_spins                       
            H_zeeman = H_zeeman - 1e3 * g * beta * B(l) * (1 - sigma_mas(k)) / h .* Iz{k};
        end        
        %% Гамильтониан скалярного взаимодействия
        H_J = zeros(dim, dim);                
        
        % Все пары J-взаимодействия
        all_pairs = nchoosek(1:n_spins, 2);    
        i_idx = all_pairs(:, 1);
        j_idx = all_pairs(:, 2);         
        for idx = 1:size(all_pairs, 1)
            i = i_idx(idx);
            j = j_idx(idx);            
            H_J = H_J + J_mas(idx)*(Iz{i}*Iz{j} + 0.5*(Iup{i}*Idn{j} + Idn{i}*Iup{j}));
        end
        
        %% Полный гамильтониан и диагонализация
        H_total = H_zeeman + H_J;        
        % Диагонализация
        [V, D] = eig(H_total);        
        lam = kron(D, eye(dim)) - kron(eye(dim), conj(D));
        U = kron(V, conj(V));
        i_U = inv(U);                
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
        %% Диполь-дипольные взаимодействия между всеми парами
        if (D_D_flag == 1)&&(read_from_file_flag==0)
        % Константы для диполь-дипольного взаимодействия
            const_HH = 1e6 * g^4 * beta^4 / (h^2);
            const_HN = 1e6 * g^2*gn^2 * beta^4 / (h^2);
            const_NN = 1e6 * gn^4 * beta^4 / (h^2);

            all_pairs = nchoosek(1:n_spins, 2);  
            all_pairs = nchoosek(1:3, 2);
            i_idx = all_pairs(:, 1);
            j_idx = all_pairs(:, 2);  
            R_DD = zeros(dim^2, dim^2);
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
        end 
        %% CSA
        if (CSA_flag==1)&&(read_from_file_flag==0)
            % операторы
            Ap = Iup{1}+Iup{2};
            Am = Idn{1}+Idn{2};
            Az = Iz{1}+Iz{2};
            Ap_m=kron(Ap, eye(dim)) - kron(eye(dim), Ap');
            Am_m=kron(Am, eye(dim)) - kron(eye(dim), Am');
            Az_m=kron(Az, eye(dim)) - kron(eye(dim), Az');
            %константы
            sigma_const=(sigmaXX^2+sigmaYY^2+sigmaZZ^2-sigmaXX*sigmaYY-sigmaXX*sigmaZZ-sigmaZZ*sigmaYY);
            const_CSA=1e3*gn*beta/h;
            %вклад в релаксационный оператор
            Rrf = Rrf -(const_CSA*B(l))^2*(1/30)*(sigma_const)*Ap_m'*U*((U\Ap_m*U).*Jlam)*i_U;
            Rrf = Rrf -(const_CSA*B(l))^2*(1/30)*(sigma_const)*Am_m'*U*((U\Am_m*U).*Jlam)*i_U;
            Rrf = Rrf -(const_CSA*B(l))^2*(4/45)*(sigma_const)*Az_m'*U*((U\Az_m*U).*Jlam)*i_U;
        end            
        if (read_from_file_flag==1)
            load(['6_spin_side_Rel_mat_DD_CSA_3.3_5.3/' num2str(tau(p)*1e9) '_' num2str(log(B(l))/log(10), 3) '.mat']);
            Rrf = R_DD_CSA;
        end
        %% Учёт корреляции диполь-диполей
        if (D_D_corr_flag==1)
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
                if (idx<16)
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
        end
        %% CSA corr
        if (CSA_D_corr_flag==1)
            angles = {[1, 2], [1, 3], [1, 4], [1, 6], [2, 4], [2, 5], [2, 6]};%NH-NN
            r_CSA_mas=[r12 r13 r14 r16 r24 r25 r26];        
            phi_mas=[phi_12 phi_13 phi_14 phi_16 phi_24 phi_25 phi_26];
            i_idx=zeros(1, length(angles));
            j_idx=zeros(1, length(angles));        
            for o=1:length(angles)
                i_idx(o)=angles{o}(1);
                j_idx(o)=angles{o}(2);            
            end
            parfor idx = 1:length(angles)
                i = i_idx(idx);  j = j_idx(idx);                    
                %dipole
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
                %CSA
                Ap = Iup{1}+Iup{2};
                Am = Idn{1}+Idn{2};
                Az = Iz{1}+Iz{2};
                Ap_m=kron(Ap, eye(dim)) - kron(eye(dim), Ap');
                Am_m=kron(Am, eye(dim)) - kron(eye(dim), Am');
                Az_m=kron(Az, eye(dim)) - kron(eye(dim), Az');
                % correlation constant
                r = r_CSA_mas(idx);            
                phi = phi_mas(idx); 
                const_CSA=1e3*gn*beta/h;
                sigma_corr=2*sigmaZZ-sigmaXX-sigmaYY-3*(sigmaXX-sigmaYY)*cos(2*(phi-psi))

                if (idx<2)
                    const_rel=-sigma_corr*const_CSA*B(l)*1e3 *gn^2 * beta^2 / (h*r^3);                    
                else
                    const_rel=-sigma_corr*const_CSA*B(l)*1e3 *gn*g * beta^2 / (h*r^3);                    
                end                
                % Вклад в релаксационный оператор
                Rrf = Rrf -const_rel*(1/60)*A_cs2_m_1'*U*((U\Az_m*U).*Jlam)*i_U;
                Rrf = Rrf -const_rel*(1/40)*A_up_m_1'*U*((U\Ap_m*U).*Jlam)*i_U;
                Rrf = Rrf -const_rel*(1/40)*A_dn_m_1'*U*((U\Am_m*U).*Jlam)*i_U;

                Rrf = Rrf -const_rel*(1/60)*Az_m'*U*((U\A_cs2_m_1*U).*Jlam)*i_U;
                Rrf = Rrf -const_rel*(1/40)*Ap_m'*U*((U\A_up_m_1*U).*Jlam)*i_U;
                Rrf = Rrf -const_rel*(1/40)*Am_m'*U*((U\A_dn_m_1*U).*Jlam)*i_U;            
            end
        end
        %% flip C_N
        if (flip_flag==1)
            A_flip_1_34=dJ1_34*0.5*(Iup{1}*(Idn{3}-Idn{4})+Idn{1}*(Iup{3}-Iup{4})+Iz{1}*(Iz{3}-Iz{4}));
            A_flip_2_34=dJ2_34*0.5*(Iup{2}*(Idn{3}-Idn{4})+Idn{2}*(Iup{3}-Iup{4})+Iz{2}*(Iz{3}-Iz{4}));
            A_flip_H_34=1e3 * g * beta * B(l) * d_sig_34 / h .* (Iz{3}-Iz{4});     
            
            A_flip_1_56=dJ1_56*0.5*(Iup{1}*(Idn{5}-Idn{6})+Idn{1}*(Iup{5}-Iup{6})+Iz{1}*(Iz{5}-Iz{6}));
            A_flip_2_56=dJ2_56*0.5*(Iup{2}*(Idn{5}-Idn{6})+Idn{2}*(Iup{5}-Iup{6})+Iz{2}*(Iz{5}-Iz{6}));
            A_flip_H_56=1e3 * g * beta * B(l) * d_sig_56 / h .* (Iz{5}-Iz{6});
            
            A_flip = A_flip_1_34+A_flip_2_34+A_flip_H_34+A_flip_1_56+A_flip_2_56+A_flip_H_56;
            A_flip_m=kron(A_flip, eye(dim)) - kron(eye(dim), A_flip');            
            Jlam=zeros(dim^2, dim^2);
            for k1 = 1:dim^2
                for k2 = 1:dim^2            	    
            	    Jlam(k1,k2)=1/(1/tau_flip+1i*(lam(k1,k1)-lam(k2,k2))); %0->inf 2xtimes slower^M
                end
            end
            Rrf = Rrf -A_flip_m'*U*((U\A_flip_m*U).*Jlam)*i_U;                      
        end
        %% Оператор эволюции и среднее время жизни синглета
        if (H_relax_flag==1)
            Rrf=Rrf+R_fast;
        end
        diff_M = -1i*U*lam*i_U+Rrf;        
        % Преобразование начальной матрицы плотности
        rv0 = reshape(ro0, [dim^2, 1]);        
        % Вычисление времени жизни через преобразование Лапласа
        s = 1e-4;
        I0 = eye(dim^2);
        %PS = Iz{4};
        rv_s = (I0 * s - diff_M) \ rv0;
        rho_s = reshape(rv_s, [dim, dim]);
        prob_S_s = trace(PS * rho_s);
        tau_S(l) = prob_S_s - 1/(4 * s);
        %tau_S(l) = prob_S_s;
        if (T1_azo_flag==1)
            tau_S(l) = prob_S_s;
        end
        
        rv_s = (I0 * s - diff_M) \ ((I0 * s - diff_M) \ rv0);
        rho_s = reshape(rv_s, [dim, dim]);
        prob_S_s = trace(PS * rho_s);
        ttau_S(l) = prob_S_s - 1/(4 * s^2);
        %ttau_S(l) = prob_S_s;
        if (T1_azo_flag==1)
            ttau_S(l) = prob_S_s;
        end
        
        disp(l);
        %% расчёт кинетики
        if (kin_flag==1)||(kin_1_flag==1)
            dt = 0.1; %с
            Nt = 5000;            
            exp_M = expm(diff_M*dt);
            kin_S = zeros(Nt, 1);        
            t_mas=zeros(Nt, 1);        
            for w=1:Nt
                rho = reshape(rv0, [dim, dim]);               
                kin_S(w)=trace(PS * rho);            
                rv0=exp_M*rv0;
                t_mas(w)=(w-1)*dt;
            end
            %% аппроксимация кинетики
            exp_model = @(p, x) p(1)*exp(-x/p(2))+p(3)*exp(-x/p(4))+1/4; 
            if (T1_azo_flag==1)
                exp_model = @(p, x) p(1)*exp(-x/p(2))+p(3)*exp(-x/p(4));
            end
            initial_guess = [0.1, 50, 0.5, 300];
            lb = [0, 0, 0, 0]; % Нижние границы
            ub = [1, inf, 1, inf];   % Верхние границы
            params_fit = lsqcurvefit(exp_model, initial_guess, t_mas, real(kin_S), lb, ub);
            p1_kin(l)=params_fit(1);
            T1_kin_1(l)=params_fit(2) ;
            p2_kin(l)=params_fit(3);
            T1_kin_2(l) = params_fit(4);
            %% сохранение в файл
            if (T1_azo_flag==1)
                fileID = fopen(['data/kin_CSA_T1_' num2str(tau(p)*1e9) '_' num2str(log(B(l))/log(10), 3) '.txt'],'w');
            else
                fileID = fopen(['data/kin_CSA_S_' num2str(tau(p)*1e9) '_' num2str(log(B(l))/log(10), 3) '.txt'],'w');
            end
            fprintf(fileID,'%4.4f %4.4f %4.4f %4.4f %4.4f\n', ttau_S(l) ./ tau_S(l), p1_kin(l), T1_kin_1(l), p2_kin(l), T1_kin_2(l));
            for i = 1:Nt
                fprintf(fileID, '%4.4f %1.10f\n', t_mas(i), kin_S(i));
            end
            fclose(fileID);
        end
    end
    %% Визуализация результатов и сохранение
    disp(ttau_S ./ tau_S);
    if (kin_1_flag==1)
        plot(t_mas, real(real(kin_S)), 'DisplayName', 'calc, B = '+string(round(B(1)*1e-4, 3))+' Тл', 'LineWidth', 2);        
        plot(t_mas, exp_model(params_fit, t_mas), '--k', 'LineWidth', 1, 'DisplayName', 'fit');
        xlabel('Время, с');
        ylabel('Рs');
    end
    %hold on;
    plot(B, real(ttau_S ./ tau_S), 'DisplayName', ['\tau_c = ' num2str(tau(p)*1e9) ' ns'], 'LineWidth', 2);
    xlabel('Магнитное поле, Гс');
    ylabel('Ts, с');
    set(gca, 'XScale', 'log');
    % сохранение в файл
    to_print=[B; real(ttau_S ./ tau_S)'; p1_kin'; T1_kin_1'; p2_kin'; T1_kin_2'];
    timestamp = datestr(now, 'yyyy_mm_dd__HH_MM_SS');
    
    if (T1_azo_flag==1)
        fileID = fopen(['data/data_T1_' num2str(tau(p)*1e9) '_' timestamp '.txt'],'w');
    else
        fileID = fopen(['data/data_' num2str(tau(p)*1e9) '_' timestamp '.txt'],'w');
    end
    fprintf(fileID,'%6.6f %4.4f %4.4f %4.4f %4.4f %4.4f\r\n',to_print);
    fclose(fileID);
end

