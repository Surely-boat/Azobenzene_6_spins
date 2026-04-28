clear;
% Всё считаем в герцах
%% Константы
beta = 5.05; % 1e-24
g = 5.585;
gn=-0.567;
h = 1.054; % 1e-27
%% Параметры расчёта
NB = 900;
B = linspace(2.3, 5.3, NB);
B = 1 * 10.^(B);
find_ac_flag=1;
ac_tres=5;
hist_flag=0;
% Количество спинов
n_spins = 12;
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
if (n_spins==12)   
    %side
    sigma3 = 7.925e-6; sigma4 = 7.925e-6; sigma5 = 7.591e-6; sigma6 = 7.591e-6; sigma7 = 7.557e-6;
    sigma8 = 7.925e-6; sigma9 = 7.925e-6; sigma10 = 7.591e-6; sigma11 = 7.591e-6; sigma12 = 7.557e-6;
    
    J13 = (-0.42) * 2 * pi; J14 = (-0.42) * 2 * pi; J15 = (0.16) * 2 * pi; J16 = (0.16) * 2 * pi;J17 = (-0.32) * 2 * pi;
    J23 = (1.67) * 2 * pi; J24 = (1.67) * 2 * pi; J25 = (0.2) * 2 * pi; J26 = (0.2) * 2 * pi; J27 = (0.31) * 2 * pi;

    J34 = 2.11 * 2 * pi; J35 = 7.96 * 2 * pi; J36 = 0.58 * 2 * pi;J37 = 1.22 * 2 * pi;
    J45 = 0.58 * 2 * pi; J46 = 7.96 * 2 * pi; J47 = 1.22 * 2 * pi;
    J56 = 1.53 * 2 * pi; J57 = 7.37 * 2 * pi;
    J67 = 7.37 * 2 * pi;

    J18 = (1.67) * 2 * pi; J19 = (1.67) * 2 * pi; J1_10 = (0.2) * 2 * pi; J1_11 = (0.2) * 2 * pi; J1_12 = (0.31) * 2 * pi;
    J28 = (-0.42) * 2 * pi; J29 = (-0.42) * 2 * pi; J2_10 = (0.16) * 2 * pi; J2_11 = (0.16) * 2 * pi;J2_12 = (-0.32) * 2 * pi;
    
    J38=0; J39=0; J3_10=0; J3_11=0; J3_12=0;
    J48=0; J49=0; J4_10=0; J4_11=0; J4_12=0;
    J58=0; J59=0; J5_10=0; J5_11=0; J5_12=0;
    J68=0; J69=0; J6_10=0; J6_11=0; J6_12=0;
    J78=0; J79=0; J7_10=0; J7_11=0; J7_12=0;

    J89 = 2.11 * 2 * pi; J8_10 = 7.96 * 2 * pi; J8_11 = 0.58 * 2 * pi;J8_12 = 1.22 * 2 * pi;
    J9_10 = 0.58 * 2 * pi; J9_11 = 7.96 * 2 * pi; J9_12 = 1.22 * 2 * pi;
    J10_11 = 1.53 * 2 * pi; J10_12 = 7.37 * 2 * pi;
    J11_12 = 7.37 * 2 * pi;

    sigma_mas=[sigma1 sigma2 sigma3 sigma4 sigma5 sigma6 sigma7 sigma8 sigma9 sigma10 sigma11 sigma12];    
    J_mas = [J12 J13 J14 J15 J16 J17 J18 J19 J1_10 J1_11 J1_12 ...
        J23 J24 J25 J26 J27 J28 J29 J2_10 J2_11 J2_12 ...
        J34 J35 J36 J37 J38 J39 J3_10 J3_11 J3_12 ...
        J45 J46 J47 J48 J49 J4_10 J4_11 J4_12 J56 J57 J58 J59 J5_10 J5_11 J5_12 ...
        J67 J68 J69 J6_10 J6_11 J6_12 J78 J79 J7_10 J7_11 J7_12 ...
        J89 J8_10 J8_11 J8_12 J9_10 J9_11 J9_12 J10_11 J10_12 J11_12];
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

%% Гамильтониан скалярного взаимодействия
H_J = zeros(dim, dim);
all_pairs = nchoosek(1:n_spins, 2);    
i_idx = all_pairs(:, 1);
j_idx = all_pairs(:, 2);         
for idx = 1:size(all_pairs, 1)
    i = i_idx(idx);
    j = j_idx(idx);           
    H_J = H_J + J_mas(idx)*(Iz{i}*Iz{j} + 0.5*(Iup{i}*Idn{j} + Idn{i}*Iup{j}));
end
%% Основные циклы по времени корреляции и магнитному полю
Energy=zeros(NB, dim);
for l = 1:NB
    %% Гамильтониан Зеемана        
    H_zeeman = zeros(dim, dim);
    H_zeeman = H_zeeman - 1e3 * gn * beta * B(l) * (1 - sigma_mas(1)) / h .* (Iz{1}+Iz{2});
    for k = 3:n_spins                       
        H_zeeman = H_zeeman - 1e3 * g * beta * B(l) * (1 - sigma_mas(k)) / h .* Iz{k};
    end      
    %% Полный гамильтониан и диагонализация
    H_total = H_zeeman + H_J;        
    % Диагонализация
    [V, D] = eig(H_total);
    for k=1:dim
        Energy(l, k)=D(k,k);
    end    

    timestamp = datestr(now, 'yyyy_mm_dd__HH_MM_SS');       
    fileID = fopen(['data/finish_' num2str(l) '_' timestamp '.txt'],'w');    
    fprintf(fileID,'help');
    fclose(fileID);
end
hold on;
B_ac=[];
if (find_ac_flag==1)
    for i=1:dim
        for j=(i+1):dim
            if ((Energy(1, i)-Energy(1, j))*(Energy(NB, i)-Energy(NB, j)))>0
                to_analize=abs(Energy(:, i)-Energy(:, j));
                tf = islocalmin(abs(to_analize), 'MinProminence', 0.3);
                indices = find(tf);
                values = to_analize(tf);
                %[values, indices] = findpeaks(to_analize);
                for k=1:length(values)
                    if(abs(values(k))<ac_tres)
                        if (hist_flag==0)
                            %scatter([B(indices(k)) B(indices(k))], [Energy(indices(k), i) Energy(indices(k), j)], 75, '*');
                            B_ac = [B_ac, B(indices(k))];                            
                        end
                        
                    end
                end                
            end
            
        end
    end
end

timestamp = datestr(now, 'yyyy_mm_dd__HH_MM_SS');
fileID = fopen(['data/LAC_'  timestamp '.txt'],'w');
fprintf(fileID,'%6.6f\n', B_ac');
fclose(fileID);
if (hist_flag==0)
    for l = 1:dim
        %plot(B, Energy(:, l)); 
    %     if(l==40)
    %         plot(B, Energy(:, l));    
    %     end
    end
    %ylim([-200 200])
    %ylim([-5e8 5e8])
else    
    %data1=get(findobj(open('6_nuc_side_Ts_DD_CSA_60.fig'), 'Type','line'), {'XData','YData'});

%     [counts, edges] = histcounts(B_ac, 'BinWidth', 100);
%     centers = (edges(1:end-1) + edges(2:end)) / 2;
%     scale_factor = 100;
%     % Увеличиваем частоты
%     counts = counts * scale_factor;   
%     hold on;
%     b = bar(centers, counts, 'Displayname', 'кол-во антиперес. N*50');
%     b.FaceAlpha = 0.3;
    %histogram(log(B_ac)/log(10), 'BinWidth', 0.1);
    %plot(data1{1,1},data1{1,2}, 'DisplayName', '6_nuc_side_Ts_DD_CSA_60 ps', 'LineWidth', 2)

end
%plot(data{1,1},data{1,2}, 'k', 'DisplayName', '6_nuc_side_Ts_DD_CSA_60 ps', 'LineWidth', 2)
%legend
%xlabel('Магнитное поле, Гс');
%ylabel('Частота, с^{-1}');
%set(gca, 'XScale', 'log');


%B_ac=(J12^2-J34^2)^(1/2)/(1e3*g*beta*abs(sigma3-sigma4)/h);
%xline(B_ac, ':k', 'a-c','LineWidth',1);
%legend
%grid on;
