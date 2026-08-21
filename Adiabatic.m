clear;
% Всё считаем в герцах
%% Константы
beta = 5.05; % 1e-24
g = 5.585;
gn=-0.567;
h = 1.054; % 1e-27
%% Параметры расчёта
NB_half=1600;
NB = NB_half*2;
%B = linspace(5.3, 3.3, NB);
%B = 1 * 10.^(B);
Find_S=zeros(NB, 1);
t_mas=zeros(NB, 1);
v_mas=zeros(NB, 1);
B=zeros(NB, 1);
v0=1;
v_max=2;
tv=0.2;
% Количество спинов
n_spins = 4;
dim = 2^n_spins;
%% Параметры системы
% Химические (ppm)
sigmaXX=-789e-6;
sigmaYY=-146e-6;
sigmaZZ=136e-6;
sigma1 = 509.94e-6;
sigma2 = 509.94e-6;
J12 = 16 * 2 * pi; % Гц переводим в цикл. частотуdJ12 = 0; % Гц
if (n_spins==4)
    sigma3 = 7.925e-6; sigma4 = 7.591e-6;
    sigma_mas=[sigma1 sigma2 sigma3 sigma4];
    % two chained protons on one side
    J13 = (-0.42) * 2 * pi; J14 = (0.16) * 2 * pi;
    J23 = (1.67) * 2 * pi; J24 = (0.2) * 2 * pi;
    J34 = 7.96 * 2 * pi;
    J_mas = [J12 J13 J14 J23 J24 J34];
end
if (n_spins==5)
    sigma3 = 7.925e-6; sigma4 = 7.591e-6; sigma5 = 7.557e-6;
    sigma_mas=[sigma1 sigma2 sigma3 sigma4 sigma5];
    % one side
    J13 = (-0.42) * 2 * pi; J14 = (0.16) * 2 * pi; J15 = (-0.32) * 2 * pi;
    J23 = (1.67) * 2 * pi; J24 = (0.2) * 2 * pi; J25 = (0.31) * 2 * pi;

    J34 = 7.96 * 2 * pi; J35 = 1.22 * 2 * pi; 
    J45 = 7.37 * 2 * pi; 

    J_mas = [J12 J13 J14 J15 J23 J24 J25 J34 J35 J45];
end
if (n_spins==6)   
    %sym
    sigma3 = 7.925e-6; sigma4 = 7.925e-6; sigma5 = 7.925e-6; sigma6 = 7.925e-6;
    J13 = (1.67) * 2 * pi; J14 = (-0.42) * 2 * pi; J15 = (-0.42) * 2 * pi; J16 = (1.67) * 2 * pi;
    J23 = (-0.42) * 2 * pi; J24 = (1.67) * 2 * pi; J25 = (1.67) * 2 * pi; J26 = (-0.42) * 2 * pi;

    J34 = 0 * 2 * pi; J35 = 0 * 2 * pi; J36 = 2.11 * 2 * pi;
    J45 = 2.11 * 2 * pi; J46 = 0 * 2 * pi;
    J56 = 0 * 2 * pi;
    %side
    sigma3 = 7.925e-6; sigma4 = 7.925e-6; sigma5 = 7.591e-6; sigma6 = 7.591e-6;
    J13 = (-0.42) * 2 * pi; J14 = (-0.42) * 2 * pi; J15 = (0.16) * 2 * pi; J16 = (0.16) * 2 * pi;
    J23 = (1.67) * 2 * pi; J24 = (1.67) * 2 * pi; J25 = (0.2) * 2 * pi; J26 = (0.2) * 2 * pi;

    J34 = 2.11 * 2 * pi; J35 = 7.96 * 2 * pi; J36 = 0.58 * 2 * pi;
    J45 = 0.58 * 2 * pi; J46 = 7.96 * 2 * pi;
    J56 = 1.53 * 2 * pi;
    %chain
%     sigma3 = 7.591e-6; sigma4 = 7.925e-6; sigma5 = 7.925e-6; sigma6 = 7.591e-6;
%     J13 = (0.2) * 2 * pi; J14 = (1.67) * 2 * pi; J15 = (-0.42) * 2 * pi; J16 = (0.16) * 2 * pi;
%     J23 = (0.16) * 2 * pi; J24 = (-0.42) * 2 * pi; J25 = (1.67) * 2 * pi; J26 = (0.2) * 2 * pi;
% 
%     J34 = 7.96 * 2 * pi; J35 = 0 * 2 * pi; J36 = 0 * 2 * pi;
%     J45 = 0 * 2 * pi; J46 = 0 * 2 * pi;
%     J56 = 7.96 * 2 * pi;

    sigma_mas=[sigma1 sigma2 sigma3 sigma4 sigma5 sigma6];    
    J_mas = [J12 J13 J14 J15 J16 J23 J24 J25 J26 J34 J35 J36 J45 J46 J56];
end
if (n_spins==7)   
    %side
    sigma3 = 7.925e-6; sigma4 = 7.925e-6; sigma5 = 7.591e-6; sigma6 = 7.591e-6; sigma7 = 7.557e-6;
    
    J13 = (-0.42) * 2 * pi; J14 = (-0.42) * 2 * pi; J15 = (0.16) * 2 * pi; J16 = (0.16) * 2 * pi;J17 = (-0.32) * 2 * pi;
    J23 = (1.67) * 2 * pi; J24 = (1.67) * 2 * pi; J25 = (0.2) * 2 * pi; J26 = (0.2) * 2 * pi; J27 = (0.31) * 2 * pi;

    J34 = 2.11 * 2 * pi; J35 = 7.96 * 2 * pi; J36 = 0.58 * 2 * pi;J37 = 1.22 * 2 * pi;
    J45 = 0.58 * 2 * pi; J46 = 7.96 * 2 * pi; J47 = 1.22 * 2 * pi;
    J56 = 1.53 * 2 * pi; J57 = 7.37 * 2 * pi;
    J67 = 7.37 * 2 * pi;
    sigma_mas=[sigma1 sigma2 sigma3 sigma4 sigma5 sigma6 sigma7];    
    J_mas = [J12 J13 J14 J15 J16 J17 J23 J24 J25 J26 J27 J34 J35 J36 J37 J45 J46 J47 J56 J57 J67];
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
PS = (eye(dim, dim)-4*Iz{1}*Iz{2}-2*(Iup{1}*Idn{2}+Idn{1}*Iup{2}))/4;
rho = ((eye(dim, dim)-4*Iz{1}*Iz{2}-2*(Iup{1}*Idn{2}+Idn{1}*Iup{2}))/4)/(2^(n_spins-2));
B_even = load('Field profile even steps 85 gear 410mA.dat');
B_odd = zeros(2150, 1);
for i=1:2150
    B_odd(i)=(B_even(i, 2)+B_even(i+1, 2))/2;
end
%% Основные циклы по времени корреляции и магнитному полю
for l = 1:NB
    if (mod(l, 2)==0)
        B(l)=B_even(fix(l/2)+1, 2)*10;
    else
        B(l)=B_odd(fix((l-1)/2)+1)*10;
    end
    if (l>1)
        if (l<(NB_half+1))      
            v=(v0-v_max)*exp(-t_mas(l-1)/tv)+v_max;
        else        
            v=(v0-v_max)*exp(-(2*t_mas(NB_half)-t_mas(l-1))/tv)+v_max; 
        end
    end
    if (l==1)  
        v=v0;
        t_mas(l)=1/v;
    else
        t_mas(l)=t_mas(l-1)+1/v;
    end
    
    
    %% Гамильтониан Зеемана        
    H_zeeman = zeros(dim, dim);
    H_zeeman = H_zeeman - 1e3 * gn * beta * B(l) * (1 - sigma_mas(1)) / h .* (Iz{1}+Iz{2});
    for k = 3:n_spins                       
        H_zeeman = H_zeeman - 1e3 * g * beta * B(l) * (1 - sigma_mas(k)) / h .* Iz{k};
    end        
    %% Гамильтониан скалярного взаимодействия
    H_J = zeros(dim, dim);                

    % Все пары диполь-дипольного взаимодействия
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
    rho=expm(-1i*H_total/v)*rho*expm(1i*H_total/v);
    Find_S(l) = trace(PS * rho);   
    v_mas(l)=v;
end
plot(B, Find_S);
%plot(t_mas, B);
%plot(t_mas, Find_S);
%plot(t_mas, B);
grid on;
set(gca, 'XScale', 'log');
xlabel('Магнитное поле, Гс');
ylabel('Ps');


