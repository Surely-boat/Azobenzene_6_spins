clear;
tic;
% Всё считаем в герцах
%% Константы
beta = 5.05; % 1e-24
g = 5.585;
gn=-0.567;
h = 1.054; % 1e-27
% Количество спинов
n_spins = 7;
dim = 2^n_spins;
%% Параметры расчёта
OPTIMIZATION.use_sparse = true;      % Использовать разреженные матрицы
OPTIMIZATION.use_gpu = false;         % Использовать GPU (если доступно)

NB = 1;
B = linspace(3.3, 3.3, NB);
B = 1 * 10.^(B);
tau = [22e-12]; 
tau_flip=10e-2;
flip_flag=0;
D_D_flag = 1;
D_D_corr_flag=0;
CSA_D_corr_flag=0;
CSA_flag=1;
read_from_file_flag=0;
H_relax_flag=0;
T1_azo_flag=0;
kin_flag=0;
kin_1_flag=0;
if (kin_1_flag==1)
    NB=1;
    B = 1.485*10.^[4];
end
%% Проверка доступности GPU
if OPTIMIZATION.use_gpu
    try
        gpuDevice();
        OPTIMIZATION.gpu_available = true;
         fprintf('GPU доступен и будет использоваться\n');
    catch
        OPTIMIZATION.gpu_available = false;
        OPTIMIZATION.use_gpu = false;
        fprintf('GPU недоступен, используется CPU\n');
    end
else
    OPTIMIZATION.gpu_available = false;
end
%% Параметры молекулы азобензола
%{
    3   5
1     
  2         7
    4   6
%} 
% Химические (ppm)
sigmaXX=-789e-6;
sigmaYY=-146e-6;
sigmaZZ=136e-6;
sigma1 = 509.94e-6;
sigma2 = 509.94e-6;
sigma3 = 7.925e-6; sigma4 = 7.925e-6; sigma5 = 7.591e-6; sigma6 = 7.591e-6; sigma7 = 7.557e-6;
d_sig=0.2e-6;
sigma_mas=[sigma1 sigma2 sigma3 sigma4 sigma5 sigma6 sigma7];


r12 = 1.248; % Å
J12 = 16 * 2 * pi; % Гц переводим в цикл. частоту
dJ12 = 0; % Гц

% Параметры дополнительных ядер
% Расстояния до основной пары (в Å)
r13 = 2.495;
r23 = 2.748; r24 = 2.544;
r35 = 2.475; 
r46 = 2.488;
r57 = 2.474;
r67 = 2.477;

r_mas=[r12 r13 r23 r24 r35 r46 r57 r67];
DD_pairs = {[1, 2], [1, 3], [2, 3], [2, 4], [3, 5], [4, 6], [5, 7], [6, 7]};


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

dJ1=2;
dJ2=0;
% Скалярные константы связи (Гц)^M
J13 = (-0.42) * 2 * pi; J14 = (-0.42) * 2 * pi; J15 = (0.16) * 2 * pi; J16 = (0.16) * 2 * pi; J17 = (-0.32) * 2 * pi;
J23 = (1.67) * 2 * pi; J24 = (1.67) * 2 * pi; J25 = (0.2) * 2 * pi; J26 = (0.2) * 2 * pi; J27 = (0.31) * 2 * pi;

J34 = 2.11 * 2 * pi; J35 = 7.96 * 2 * pi; J36 = 0.58 * 2 * pi;J37 = 1.22 * 2 * pi;
J45 = 0.58 * 2 * pi; J46 = 7.96 * 2 * pi; J47 = 1.22 * 2 * pi;
J56 = 1.53 * 2 * pi; J57 = 7.37 * 2 * pi;
J67 = 7.37 * 2 * pi;

J_mas = [J12 J13 J14 J15 J16 J17 J23 J24 J25 J26 J27 J34 J35 J36 J37 J45 J46 J47 J56 J57 J67];

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
%% Создание операторов для каждого спина
up = [0 1; 0 0]; dn = [0 0; 1 0]; z = [0.5 0; 0 -0.5];
if OPTIMIZATION.use_sparse
    up = sparse(up);
    dn = sparse(dn);
    z = sparse(z);
end

for i = 1:n_spins
    if OPTIMIZATION.use_sparse
        Iup{i} = kron(kron(speye(2^(i-1)), up), speye(2^(n_spins-i)));
        Idn{i} = kron(kron(speye(2^(i-1)), dn), speye(2^(n_spins-i)));
        Iz{i} = kron(kron(speye(2^(i-1)), z), speye(2^(n_spins-i)));
    else
        Iup{i} = kron(kron(eye(2^(i-1)), up), eye(2^(n_spins-i)));
        Idn{i} = kron(kron(eye(2^(i-1)), dn), eye(2^(n_spins-i)));
        Iz{i} = kron(kron(eye(2^(i-1)), z), eye(2^(n_spins-i)));
    end
end
if OPTIMIZATION.use_sparse
    eye_dim = speye(dim);
    eye_dim2 = speye(dim^2);
else
    eye_dim = eye(dim);
    eye_dim2 = eye(dim^2);
end
% Начальная матрица плотности (только первая пара в синглете)
alpha=[1 0; 0 0];
bita=[0 0; 0 1];
equil = [1 0; 0 1]/2;
singlet=[0 0 0 0;
    0 0.5 -0.5 0;
    0 -0.5 0.5 0;
    0 0 0 0];
if OPTIMIZATION.use_sparse
    alpha = sparse(alpha);
    bita = sparse(bita);
    equil = sparse(equil);
    singlet = sparse(singlet);
end
ro0 = kron(singlet, kron(equil, kron(equil, kron(equil, kron(equil, equil)))));

%ro0 = kron(equil, kron(equil, kron(equil, kron(alpha, kron(equil, equil)))));
if (T1_azo_flag==1)
    ro0 = kron(kron(alpha, alpha), kron(equil, kron(equil, kron(equil, kron(equil, equil)))));
end
% Проектор на синглетное состояние первых двух спинов
PS = (eye_dim-4*Iz{1}*Iz{2}-2*(Iup{1}*Idn{2}+Idn{1}*Iup{2}))*0.25; %Нужно поделить на 16, чтобы получить синглет+равновесие водородов
PT_p = (eye_dim+2*Iz{1}+2*Iz{2}+4*Iz{1}*Iz{2})*0.25;
PT_0 = (eye_dim-4*Iz{1}*Iz{2}+2*(Iup{1}*Idn{2}+Idn{1}*Iup{2}))*0.25;
PT_m = (eye_dim-2*Iz{1}-2*Iz{2}+4*Iz{1}*Iz{2})*0.25;
if (T1_azo_flag==1)
    PS = Iz{1}+Iz{2};
end
% superoperators for CSA
Ap = Iup{1}+Iup{2};
Am = Idn{1}+Idn{2};
Az = Iz{1}+Iz{2};
Ap_m=kron(Ap, eye_dim) - kron(eye_dim, Ap');
Am_m=kron(Am, eye_dim) - kron(eye_dim, Am');
if OPTIMIZATION.use_gpu
    Ap_m=gpuArray(Ap_m);
    Am_m=gpuArray(Am_m);    
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
        H_zeeman = OPTIMIZATION.use_sparse * sparse(dim, dim) + ...
                   ~OPTIMIZATION.use_sparse * zeros(dim, dim);
        H_zeeman = H_zeeman - 1e3 * gn * beta * B(l) * (1 - sigma_mas(1)) / h .* (Iz{1}+Iz{2});
        for k = 3:n_spins                       
            H_zeeman = H_zeeman - 1e3 * g * beta * B(l) * (1 - sigma_mas(k)) / h .* Iz{k};
        end        
        %% Гамильтониан скалярного взаимодействия
        H_J = OPTIMIZATION.use_sparse * sparse(dim, dim) + ...
              ~OPTIMIZATION.use_sparse * zeros(dim, dim);                
        % Все пары диполь-дипольного взаимодействия
        all_pairs = nchoosek(1:n_spins, 2);    
        i_idx = all_pairs(:, 1);
        j_idx = all_pairs(:, 2);         
        for idx = 1:size(all_pairs, 1)
            spin_i = i_idx(idx);
            j = j_idx(idx);            
            H_J = H_J + J_mas(idx)*(Iz{spin_i}*Iz{j} + 0.5*(Iup{spin_i}*Idn{j} + Idn{spin_i}*Iup{j}));
        end
        %% Полный гамильтониан и диагонализация
        H_total = H_zeeman + H_J;        
        % Диагонализация
        if OPTIMIZATION.use_sparse
            [V, D] = eigs(H_total, dim);
        else
            [V, D] = eig(H_total);
        end

        lam = kron(D, eye_dim) - kron(eye_dim, conj(D));
        U = kron(V, conj(V));
        i_V=inv(V);
        i_U=kron(i_V, conj(i_V));              
        lam_diag = diag(D);
        lam_diff = lam_diag - lam_diag.';
        Jlam = 1 ./ (1/tau(p) + 1i * lam_diff);
        Jlam = Jlam(:);
        if OPTIMIZATION.use_sparse
            Jlam_sparse = spdiags(Jlam, 0, dim^2, dim^2);
            if OPTIMIZATION.use_gpu
                Jlam_sparse = gpuArray(Jlam_sparse);
            end
        end
        if OPTIMIZATION.use_gpu
            Jlam = gpuArray(Jlam);
            i_U = gpuArray(i_U);
            U = gpuArray(U);
        end
        
        
        % Релаксационный оператор Редфилда        
        Rrf = OPTIMIZATION.use_sparse * sparse(dim^2, dim^2) + ...
                   ~OPTIMIZATION.use_sparse * zeros(dim^2, dim^2);
        if OPTIMIZATION.use_gpu
            Rrf = gpuArray(Rrf);
        end
        
        %% Собственная релаксация протонов
        if (H_relax_flag==1)
            r_HH = 2.481;
            const_1 = 1e6*g^4*beta^4/(h^2*r_HH^6);
            sigma_meta = 7.591e-6;
            w1 = -1e3*g*beta*B(l)*(1-sigma3)/h;
            w2 = -1e3*g*beta*B(l)*(1-sigma_meta)/h;
            T1_const=1./(3/10*const_1*tau(p)*(1./(1+w2.^2*tau(p)^2)+4./(1+(w1+w2).^2*tau(p)^2)));
            T2_const=1./(3/20*const_1*tau(p)*(3+5./(1+w1.^2*tau(p)^2)+2./(1+(w1+w2).^2*tau(p)^2)));
            T1 = [0 0 T1_const T1_const T1_const T1_const];
            T2 = [0 0 T2_const T2_const T2_const T2_const];
            if OPTIMIZATION.use_gpu
                T1 = gpuArray(T1);
                T2 = gpuArray(T2);
            end
            %T2 = [0 0 10 5 5 5];            
            R_fast_1 = OPTIMIZATION.use_sparse * sparse(dim^2, dim^2) + ...
                   ~OPTIMIZATION.use_sparse * zeros(dim^2, dim^2);
            R_fast_2 = OPTIMIZATION.use_sparse * sparse(dim^2, dim^2) + ...
                   ~OPTIMIZATION.use_sparse * zeros(dim^2, dim^2);
            if OPTIMIZATION.use_gpu
                R_fast_1 = gpuArray(R_fast_1);
                R_fast_2 = gpuArray(R_fast_2);
            end
            for ii=3:6
                R_fast_1=R_fast_1+(kron(Iup{ii},Iup{ii})+kron(Idn{ii},Idn{ii})-2*kron(Iz{ii},Iz{ii})-0.5*kron(eye_dim, eye_dim))/(2*T1(ii));
                R_fast_2=R_fast_2+(2*kron(Iz{ii},Iz{ii})-0.5*kron(eye_dim, eye_dim))/(T2(ii));    
            end            
            Rrf=Rrf+R_fast_1+R_fast_2;
        end
        %% Диполь-дипольные взаимодействия между всеми парами
        if (D_D_flag == 1)&&(read_from_file_flag==0)
        % Константы для диполь-дипольного взаимодействия
            const_HH = 1e6 * g^4 * beta^4 / (h^2);
            const_HN = 1e6 * g^2*gn^2 * beta^4 / (h^2);
            const_NN = 1e6 * gn^4 * beta^4 / (h^2);
            for o=1:length(DD_pairs)
                i_idx(o)=DD_pairs{o}(1);
                j_idx(o)=DD_pairs{o}(2);
            end
            for idx = 1:length(DD_pairs)
                spin_i = i_idx(idx);
                j = j_idx(idx);                                     
                % Операторы диполь-дипольного взаимодействия 
                A_cs2 = Iup{spin_i}*Idn{j}+Idn{spin_i}*Iup{j}-4*Iz{spin_i}*Iz{j};
                A_up = Iz{spin_i}*Iup{j}+Iup{spin_i}*Iz{j};
                A_dn = Iz{spin_i}*Idn{j}+Idn{spin_i}*Iz{j};
                A4 = Iup{spin_i}*Iup{j};
                A5 = Idn{spin_i}*Idn{j};
                % размерность
                A_cs2_m=kron(A_cs2, eye_dim) - kron(eye_dim, A_cs2');
                A_up_m=kron(A_up, eye_dim) - kron(eye_dim, A_up');
                A_dn_m=kron(A_dn, eye_dim) - kron(eye_dim, A_dn');
                A_4_m=kron(A4, eye_dim) - kron(eye_dim, A4');
                A_5_m=kron(A5, eye_dim) - kron(eye_dim, A5');                        
                % Расстояние между ядрами
                r_ij = r_mas(idx);
                if (spin_i==1)&&(j==2)
                    const_rel=const_NN/r_ij^6;
                else
                    if (spin_i==1)||(spin_i==2)
                        const_rel=const_HN/r_ij^6;  
                    else
                        const_rel=const_HH/r_ij^6;                         
                    end
                end            
                %Вклад в релаксационный оператор
                if OPTIMIZATION.use_sparse
                    Rrf = Rrf -const_rel*0.05*A_cs2_m'*U*(Jlam_sparse *(i_U*A_cs2_m*U))*i_U;
                    Rrf = Rrf -const_rel*0.3*A_up_m'*U*(Jlam_sparse *(i_U*A_up_m*U))*i_U;
                    Rrf = Rrf -const_rel*0.3*A_dn_m'*U*(Jlam_sparse *(i_U*A_dn_m*U))*i_U;
                    Rrf = Rrf -const_rel*0.3*A_4_m'*U*(Jlam_sparse *(i_U*A_4_m*U))*i_U;
                    Rrf = Rrf -const_rel*0.3*A_5_m'*U*(Jlam_sparse *(i_U*A_5_m*U))*i_U;
                else
                    Rrf = Rrf -const_rel*0.05*A_cs2_m'*U*((i_U*A_cs2_m*U).*Jlam)*i_U;
                    Rrf = Rrf -const_rel*0.3*A_up_m'*U*((i_U*A_up_m*U).*Jlam)*i_U;
                    Rrf = Rrf -const_rel*0.3*A_dn_m'*U*((i_U*A_dn_m*U).*Jlam)*i_U;
                    Rrf = Rrf -const_rel*0.3*A_4_m'*U*((i_U*A_4_m*U).*Jlam)*i_U;
                    Rrf = Rrf -const_rel*0.3*A_5_m'*U*((i_U*A_5_m*U).*Jlam)*i_U;
                end
               disp(['DD numder ' num2str(idx) ' done']); 
            end
            R_DD = Rrf;
            save(['Seven_spin_DD_relmat_22ps/' num2str(log(B(l))/log(10),3) 'Gs.mat'], "R_DD", '-v7.3');
        end 
        %% CSA
        if (CSA_flag==1)&&(read_from_file_flag==0)
           %константы
            sigma_const=(sigmaXX^2+sigmaYY^2+sigmaZZ^2-sigmaXX*sigmaYY-sigmaXX*sigmaZZ-sigmaZZ*sigmaYY);
            const_CSA=1e3*gn*beta/h;
            %вклад в релаксационный оператор
            if OPTIMIZATION.use_sparse
                Rrf = Rrf -(const_CSA*B(l))^2*(1/30)*(sigma_const)*Ap_m'*U*(Jlam_sparse *(i_U*Ap_m*U))*i_U;
                Rrf = Rrf -(const_CSA*B(l))^2*(1/30)*(sigma_const)*Am_m'*U*(Jlam_sparse *(i_U*Am_m*U))*i_U;                
            else
                Rrf = Rrf -(const_CSA*B(l))^2*(1/30)*(sigma_const)*Ap_m'*U*((i_U*Ap_m*U).*Jlam)*i_U;
                Rrf = Rrf -(const_CSA*B(l))^2*(1/30)*(sigma_const)*Am_m'*U*((i_U*Am_m*U).*Jlam)*i_U;
            end            
        end        
        %% Учёт корреляции диполь-диполей
        if (D_D_corr_flag==1)&&(read_from_file_flag==0)
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
            for idx = 1:length(angles)
                spin_i = i_idx(idx);  j = j_idx(idx);   spin_k = k_idx(idx);  m = m_idx(idx);                         
                %first pair
                A_cs2_1 = Iup{j}*Idn{spin_i}+Idn{j}*Iup{spin_i}-4*Iz{j}*Iz{spin_i};
                A_up_1 = Iz{j}*Iup{spin_i}+Iup{j}*Iz{spin_i};
                A_dn_1 = Iz{j}*Idn{spin_i}+Idn{j}*Iz{spin_i};
                A4_1 = Iup{j}*Iup{spin_i};
                A5_1 = Idn{j}*Idn{spin_i};
                % dimentions
                A_cs2_m_1=kron(A_cs2_1, eye_dim) - kron(eye_dim, A_cs2_1');
                A_up_m_1=kron(A_up_1, eye_dim) - kron(eye_dim, A_up_1');
                A_dn_m_1=kron(A_dn_1, eye_dim) - kron(eye_dim, A_dn_1');
                A_4_m_1=kron(A4_1, eye_dim) - kron(eye_dim, A4_1');
                A_5_m_1=kron(A5_1, eye_dim) - kron(eye_dim, A5_1');
                %second pair
                A_cs2_2 = Iup{m}*Idn{spin_k}+Idn{m}*Iup{spin_k}-4*Iz{m}*Iz{spin_k};
                A_up_2 = Iz{m}*Iup{spin_k}+Iup{m}*Iz{spin_k};
                A_dn_2 = Iz{m}*Idn{spin_k}+Idn{m}*Iz{spin_k};
                A4_2 = Iup{m}*Iup{spin_k};
                A5_2 = Idn{m}*Idn{spin_k};
                % dimentions
                A_cs2_m_2=kron(A_cs2_2, eye_dim) - kron(eye_dim, A_cs2_2');
                A_up_m_2=kron(A_up_2, eye_dim) - kron(eye_dim, A_up_2');
                A_dn_m_2=kron(A_dn_2, eye_dim) - kron(eye_dim, A_dn_2');
                A_4_m_2=kron(A4_2, eye_dim) - kron(eye_dim, A4_2');
                A_5_m_2=kron(A5_2, eye_dim) - kron(eye_dim, A5_2');
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
                if OPTIMIZATION.use_sparse
                    Rrf = Rrf -const_rel*(1/80)*A_cs2_m_1'*U*(Jlam_sparse *(i_U*A_cs2_m_2*U))*i_U;
                    Rrf = Rrf -const_rel*(3/40)*A_up_m_1'*U*(Jlam_sparse *(i_U*A_up_m_2*U))*i_U;
                    Rrf = Rrf -const_rel*(3/40)*A_dn_m_1'*U*(Jlam_sparse *(i_U*A_dn_m_2*U))*i_U;
                    Rrf = Rrf -const_rel*(3/40)*A_4_m_1'*U*(Jlam_sparse *(i_U*A_4_m_2*U))*i_U;
                    Rrf = Rrf -const_rel*(3/40)*A_5_m_1'*U*(Jlam_sparse *(i_U*A_5_m_2*U))*i_U;
    
                    Rrf = Rrf -const_rel*(1/80)*A_cs2_m_2'*U*(Jlam_sparse *(i_U*A_cs2_m_1*U))*i_U;
                    Rrf = Rrf -const_rel*(3/40)*A_up_m_2'*U*(Jlam_sparse *(i_U*A_up_m_1*U))*i_U;
                    Rrf = Rrf -const_rel*(3/40)*A_dn_m_2'*U*(Jlam_sparse *(i_U*A_dn_m_1*U))*i_U;
                    Rrf = Rrf -const_rel*(3/40)*A_4_m_2'*U*(Jlam_sparse *(i_U*A_4_m_1*U))*i_U;
                    Rrf = Rrf -const_rel*(3/40)*A_5_m_2'*U*(Jlam_sparse *(i_U*A_5_m_1*U))*i_U;
                else
                    Rrf = Rrf -const_rel*(1/80)*A_cs2_m_1'*U*((i_U*A_cs2_m_2*U).*Jlam)*i_U;
                    Rrf = Rrf -const_rel*(3/40)*A_up_m_1'*U*((i_U*A_up_m_2*U).*Jlam)*i_U;
                    Rrf = Rrf -const_rel*(3/40)*A_dn_m_1'*U*((i_U*A_dn_m_2*U).*Jlam)*i_U;
                    Rrf = Rrf -const_rel*(3/40)*A_4_m_1'*U*((i_U*A_4_m_2*U).*Jlam)*i_U;
                    Rrf = Rrf -const_rel*(3/40)*A_5_m_1'*U*((i_U*A_5_m_2*U).*Jlam)*i_U;
    
                    Rrf = Rrf -const_rel*(1/80)*A_cs2_m_2'*U*((i_U*A_cs2_m_1*U).*Jlam)*i_U;
                    Rrf = Rrf -const_rel*(3/40)*A_up_m_2'*U*((i_U*A_up_m_1*U).*Jlam)*i_U;
                    Rrf = Rrf -const_rel*(3/40)*A_dn_m_2'*U*((i_U*A_dn_m_1*U).*Jlam)*i_U;
                    Rrf = Rrf -const_rel*(3/40)*A_4_m_2'*U*((i_U*A_4_m_1*U).*Jlam)*i_U;
                    Rrf = Rrf -const_rel*(3/40)*A_5_m_2'*U*((i_U*A_5_m_1*U).*Jlam)*i_U;
                end                                          
            end
        end
        %% CSA corr
        if (CSA_D_corr_flag==1)&&(read_from_file_flag==0)
            angles = {[1, 2], [1, 3], [1, 4], [1, 6], [2, 4], [2, 5], [2, 6]};%NH-NN
            r_CSA_mas=[r12 r13 r14 r16 r24 r25 r26];        
            phi_mas=[phi_12 phi_13 phi_14 phi_16 phi_24 phi_25 phi_26];
            i_idx=zeros(1, length(angles));
            j_idx=zeros(1, length(angles));        
            for o=1:length(angles)
                i_idx(o)=angles{o}(1);
                j_idx(o)=angles{o}(2);            
            end
            for idx = 1:length(angles)
                spin_i = i_idx(idx);  j = j_idx(idx);                    
                %dipole
                A_cs2_1 = Iup{j}*Idn{spin_i}+Idn{j}*Iup{spin_i}-4*Iz{j}*Iz{spin_i};
                A_up_1 = Iz{j}*Iup{spin_i}+Iup{j}*Iz{spin_i};
                A_dn_1 = Iz{j}*Idn{spin_i}+Idn{j}*Iz{spin_i};
                A4_1 = Iup{j}*Iup{spin_i};
                A5_1 = Idn{j}*Idn{spin_i};
                % dimentions
                A_cs2_m_1=kron(A_cs2_1, eye_dim) - kron(eye_dim, A_cs2_1');
                A_up_m_1=kron(A_up_1, eye_dim) - kron(eye_dim, A_up_1');
                A_dn_m_1=kron(A_dn_1, eye_dim) - kron(eye_dim, A_dn_1');
                A_4_m_1=kron(A4_1, eye_dim) - kron(eye_dim, A4_1');
                A_5_m_1=kron(A5_1, eye_dim) - kron(eye_dim, A5_1');
                %CSA
                Ap = Iup{1}+Iup{2};
                Am = Idn{1}+Idn{2};
                Az = Iz{1}+Iz{2};
                Ap_m=kron(Ap, eye_dim) - kron(eye_dim, Ap');
                Am_m=kron(Am, eye_dim) - kron(eye_dim, Am');
                Az_m=kron(Az, eye_dim) - kron(eye_dim, Az');
                % correlation constant
                r = r_CSA_mas(idx);            
                phi = phi_mas(idx); 
                const_CSA=1e3*gn*beta/h;
                sigma_corr=2*sigmaZZ-sigmaXX-sigmaYY-3*(sigmaXX-sigmaYY)*cos(2*(phi-psi));

                if (idx<2)
                    const_rel=-sigma_corr*const_CSA*B(l)*1e3 *gn^2 * beta^2 / (h*r^3);                    
                else
                    const_rel=-sigma_corr*const_CSA*B(l)*1e3 *gn*g * beta^2 / (h*r^3);                    
                end                
                % Вклад в релаксационный оператор
                if OPTIMIZATION.use_sparse
                    Rrf = Rrf -const_rel*(1/60)*A_cs2_m_1'*U*(Jlam_sparse *(i_U*Az_m*U))*i_U;
                    Rrf = Rrf -const_rel*(1/40)*A_up_m_1'*U*(Jlam_sparse *(i_U*Ap_m*U))*i_U;
                    Rrf = Rrf -const_rel*(1/40)*A_dn_m_1'*U*(Jlam_sparse *(i_U*Am_m*U))*i_U;
    
                    Rrf = Rrf -const_rel*(1/60)*Az_m'*U*(Jlam_sparse *(i_U*A_cs2_m_1*U))*i_U;
                    Rrf = Rrf -const_rel*(1/40)*Ap_m'*U*(Jlam_sparse *(i_U*A_up_m_1*U))*i_U;
                    Rrf = Rrf -const_rel*(1/40)*Am_m'*U*(Jlam_sparse *(i_U*A_dn_m_1*U))*i_U;
                else
                    Rrf = Rrf -const_rel*(1/60)*A_cs2_m_1'*U*((i_U*Az_m*U).*Jlam)*i_U;
                    Rrf = Rrf -const_rel*(1/40)*A_up_m_1'*U*((i_U*Ap_m*U).*Jlam)*i_U;
                    Rrf = Rrf -const_rel*(1/40)*A_dn_m_1'*U*((i_U*Am_m*U).*Jlam)*i_U;
    
                    Rrf = Rrf -const_rel*(1/60)*Az_m'*U*((i_U*A_cs2_m_1*U).*Jlam)*i_U;
                    Rrf = Rrf -const_rel*(1/40)*Ap_m'*U*((i_U*A_up_m_1*U).*Jlam)*i_U;
                    Rrf = Rrf -const_rel*(1/40)*Am_m'*U*((i_U*A_dn_m_1*U).*Jlam)*i_U; 
                end
                           
            end
        end
        %R_DD_CSA_corr = Rrf;
        %save(['6_spin_sym_Rel_mat_DD_CSA_3.3_5.3/' num2str(tau(p)*1e9) '_' num2str(log(B(l))/log(10), 3) '.mat'], 'R_DD_CSA_corr');
         
        % if (read_from_file_flag==1)
        %     load(['6_spin_sym_Rel_mat_DD_CSA_3.3_5.3/' num2str(tau(p)*1e9) '_' num2str(log(B(l))/log(10), 3) '.mat']);
        %     Rrf = R_DD_CSA_corr;
        % end
        %% flip C_N
        if (flip_flag==1)
            A_flip_1=dJ2*0.5*(Iup{1}*(Idn{3}-Idn{6})+Idn{1}*(Iup{3}-Iup{6})+Iz{1}*(Iz{3}-Iz{6}));
            A_flip_2=dJ1*0.5*(Iup{2}*(Idn{3}-Idn{6})+Idn{2}*(Iup{3}-Iup{6})+Iz{2}*(Iz{3}-Iz{6}));
            A_flip_3=1e3 * g * beta * B(l) * d_sig / h .* (Iz{3}-Iz{6});
            A_flip_36 = A_flip_1+A_flip_2+A_flip_3;
            A_flip_36_m=kron(A_flip_36, eye_dim) - kron(eye_dim, A_flip_36'); 
            
            A_flip_1=dJ1*0.5*(Iup{1}*(Idn{4}-Idn{5})+Idn{1}*(Iup{4}-Iup{5})+Iz{1}*(Iz{4}-Iz{5}));
            A_flip_2=dJ2*0.5*(Iup{2}*(Idn{4}-Idn{5})+Idn{2}*(Iup{4}-Iup{5})+Iz{2}*(Iz{4}-Iz{5}));
            A_flip_3=1e3 * g * beta * B(l) * d_sig / h .* (Iz{4}-Iz{5});
            A_flip_45 = A_flip_1+A_flip_2+A_flip_3;
            A_flip_45_m=kron(A_flip_45, eye_dim) - kron(eye_dim, A_flip_45');      
            
            Jlam=zeros(dim^2, dim^2);
            for k1 = 1:dim^2
                for k2 = 1:dim^2            	    
            	    Jlam(k1,k2)=1/(1/tau_flip+1i*(lam(k1,k1)-lam(k2,k2))); %0->inf 2xtimes slower^M
                end
            end
            Rrf = Rrf -A_flip_36_m'*U*((U\A_flip_36_m*U).*Jlam)*i_U;   
            Rrf = Rrf -A_flip_45_m'*U*((U\A_flip_45_m*U).*Jlam)*i_U;  
        end
        %% Оператор эволюции и среднее время жизни синглета        
        diff_M = -1i*(kron(H_total, eye_dim)-kron(eye_dim, conj(H_total)))+gather(Rrf);        
        % Преобразование начальной матрицы плотности
        rv0 = reshape(ro0, [dim^2, 1]);        
        % Вычисление времени жизни через преобразование Лапласа
        s = 1e-5;
        I0 = eye_dim2;
        %PS = Iz{4};
        rv_s = (I0 * s - diff_M) \ rv0;
        rho_s = reshape(rv_s, [dim, dim]);
        prob_S_s = trace(PS * rho_s);
        tau_S(l) = prob_S_s - 1/(4 * s);
        %tau_S(l) = prob_S_s;
        if (T1_azo_flag==1)
            tau_S(l) = prob_S_s;
        end
        
        rv_s = (I0 * s - diff_M) \ rv_s;
        rho_s = reshape(rv_s, [dim, dim]);
        prob_S_s = trace(PS * rho_s);
        ttau_S(l) = prob_S_s - 1/(4 * s^2);
        %ttau_S(l) = prob_S_s;
        if (T1_azo_flag==1)
            ttau_S(l) = prob_S_s;
        end
        
        disp(l);
        disp(['field = ' num2str(B(l))]);
        disp(['Ts = ' num2str(real(ttau_S(l)/tau_S(l)))]);
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
            initial_guess = [0.1, 10, 0.5, 100];
            lb = [0, 0, 0, 0]; % Нижние границы
            ub = [1, inf, 1, inf];   % Верхние границы
            params_fit = lsqcurvefit(exp_model, initial_guess, t_mas, real(kin_S), lb, ub);
            p1_kin(l)=params_fit(1);
            T1_kin_1(l)=params_fit(2);
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
    
    %hold on;
    %plot(B, real(ttau_S ./ tau_S), 'DisplayName', ['\tau_c = ' num2str(tau(p)*1e9) ' ns'], 'LineWidth', 2);
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
elapsed_time = toc;
fprintf('\n=== Выполнение завершено ===\n');
fprintf('Общее время: %.4f секунд (%.2f минут) (%.2f часов)\n', elapsed_time, elapsed_time/60, elapsed_time/3600);
