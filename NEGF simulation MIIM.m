
clear all

%ctes
hbar=1.06e-34;
q=1.6e-19;
Kb=8.617e-5;
T=300;
E_T=Kb*T

%effective mass
%m=0.07*9.11e-31;
m = 0.1 * 9.11e-31%masse eff moyenne pour les oxydes

%position vector
Dx= 1e-10; %position mesh (m)   %Dense grid required the m_eff approximation breaks down outside the parabolic regime of H
N=250; %position points number
XX=(1:N).*Dx;%position vector
L=Dx*N; %length

%potential versus position
U=zeros(N,1);

%step


Phi_Au = 5.10;
Phi_Ag = 4.26;

X_TiO2  = 3.90;

X_Al2O3 = 2.30; %val Th 1.35


h_TiO2  = Phi_Au - X_TiO2;
h_Al2O3 = Phi_Au - X_Al2O3;
%h_TiO2 = 0.4
%h_Al2O3= 0.3

w_TiO2  = 1.2e-9;%2.5
w_Al2O3 = 1.3e-9;%3

pts_TiO2  = round(w_TiO2 / Dx);
pts_Al2O3 = round(w_Al2O3 / Dx);

%Centering
Total_Pts_Oxydes = pts_TiO2 + pts_Al2O3;
i_start = round(N/2 - Total_Pts_Oxydes/2);

idx_T_end   = i_start + pts_TiO2 - 1;
idx_A_start = idx_T_end + 1;
idx_A_end   = idx_A_start + pts_Al2O3 - 1;



Pente_Vide = linspace(Phi_Au, Phi_Ag, Total_Pts_Oxydes)';


U = zeros(N,1);

U(i_start : idx_T_end) = Pente_Vide(1 : pts_TiO2) - X_TiO2;

U(idx_A_start : idx_A_end) = Pente_Vide(pts_TiO2+1 : end) - X_Al2O3;






%energy vector
dE=10e-3; %energy mesh (eV)
ES=(min(U)-0.1:dE:max(U)+0.5);

dES=ES(2)-ES(1);%pas
NE=length(ES);

%effective mass Hamiltonian
t=hbar.^2/(2*m*Dx.^2*q);
H=(2*t*diag(ones(1,N)))-(t*diag(ones(1,N-1),1))-(t*diag(ones(1,N-1),-1));








% Green functions
sigS=zeros(N); sigD=zeros(N);
DOS=[]; Trans=[]; J_E=[];

V_bias = 0.1;

for k=1:NE % loop over the energies
    E = ES(k);

    ck_L = 1 - ((E - U(1))/(2*t));


    if abs(ck_L) <= 1
        ka_L = acos(ck_L);
    else
        ka_L = -1j * acosh(abs(ck_L));
    end
    sigS(1,1) = -t * exp(1j * ka_L);
    gam1 = 1j * (sigS - sigS');

    ck_R = 1 - ((E - U(N))/(2*t));

    if abs(ck_R) <= 1
        ka_R = acos(ck_R);
    else
        ka_R = -1j * acosh(abs(ck_R));
    end
    sigD(N,N) = -t * exp(1j * ka_R);
    gam2 = 1j * (sigD - sigD');

    % --- 3. FONCTION DE GREEN ---
    GS = inv((E * eye(N)) - H - diag(U) - sigS - sigD);


    % LDOS
    SpecFunc = 1j * (GS - GS');
    DOS = [DOS, diag(real(SpecFunc))];

    % Transmission
    Trans(k) = real(trace(gam1 * GS * gam2 * GS'));

    %  J(E)
    mu_L = 0;
    mu_R = V_bias;

    fL = 1 / (1 + exp((E - mu_L)/E_T));
    fR = 1 / (1 + exp((E - mu_R)/E_T));

    J_E(k) = Trans(k) * (fL - fR);
end


figure(2); clf;
plot(ES, Trans, 'LineWidth', 2);
xlabel("Energy (eV)");
ylabel("Transmission T(E)");
title("Transmission Probability (Linear Scale)");
grid on;


figure(5); clf;
plot(ES, real(J_E), 'LineWidth', 2);
xlabel("Energy (eV)");
ylabel("Spectral Current J(E)");
title(['Spectral Current flow for V = ' num2str(V_bias) ' V']);
grid on;

%%%% FIGURE 6 : LDOS  %%%%
figure(6); clf;
set(gcf, 'Position', [100, 100, 1200, 600]);

XX_nm = XX * 1e9;

LDOS_Log = log10(abs(DOS') + 1e-20);

imagesc(XX_nm, ES, LDOS_Log);
shading interp
colormap('jet');
set(gca, 'YDir', 'normal');
colorbar;
caxis([-6, 0]);

hold on
plot(XX_nm, U, 'w', 'LineWidth', 3);
hold off

xlabel("Position (nm)");
ylabel("Energy (eV)");
title("Local Density of States & Band Diagram");

CD_nm = max(XX_nm) / 2;
xlim([CD_nm - 10, CD_nm + 10]);
















%%%%%%%%%%%%% I-V and RR calculation


eps_r_TiO2  = 50;
eps_r_Al2O3 = 9;


denom = (eps_r_TiO2 * w_Al2O3) + (eps_r_Al2O3 * w_TiO2);
frac_V_TiO2  = (eps_r_Al2O3 * w_TiO2) / denom;
frac_V_Al2O3 = (eps_r_TiO2 * w_Al2O3) / denom;

V_sweep_max = 2.5; % Volts
V_points = 41;
V_list = linspace(-V_sweep_max, V_sweep_max, V_points);
I_list = zeros(size(V_list));

U_eq = U;

fprintf('Lancement du sweep I-V sur %d points...\n', V_points);

%  (Voltage Sweep)
for iv = 1:length(V_list)
    V_app = V_list(iv);


    U_bias = zeros(N,1);

    drop_TiO2 = -V_app * frac_V_TiO2;
    slope_TiO2 = drop_TiO2 / pts_TiO2;
    U_bias(i_start : idx_T_end) = (1:pts_TiO2)' * slope_TiO2;

    drop_Al2O3 = -V_app * frac_V_Al2O3;
    slope_Al2O3 = drop_Al2O3 / pts_Al2O3;
    offset_Al2O3 = drop_TiO2;
    U_bias(idx_A_start : idx_A_end) = offset_Al2O3 + (1:pts_Al2O3)' * slope_Al2O3;


    U_bias(idx_A_end+1 : end) = -V_app;

    U_total = U_eq + U_bias;


    Current_Density_E = zeros(NE, 1);

    for k = 1:NE
        E = ES(k);


        ck_L = 1 - ((E - U_total(1))/(2*t));
        if abs(ck_L) <= 1, ka_L = acos(ck_L); else, ka_L = -1j*acosh(abs(ck_L)); end
        sigS_val = -t * exp(1j * ka_L);

        ck_R = 1 - ((E - U_total(N))/(2*t));
        if abs(ck_R) <= 1, ka_R = acos(ck_R); else, ka_R = -1j*acosh(abs(ck_R)); end
        sigD_val = -t * exp(1j * ka_R);

        Sigma_L = sparse(1, 1, sigS_val, N, N);
        Sigma_R = sparse(N, N, sigD_val, N, N);

        Gamma_L = 1j * (Sigma_L - Sigma_L');
        Gamma_R = 1j * (Sigma_R - Sigma_R');

        % Fonction  Green rtr
        G_R = inv((E * eye(N)) - H - diag(U_total) - Sigma_L - Sigma_R);

        T_E = real(trace(Gamma_L * G_R * Gamma_R * G_R'));

        % Distrib Fermi
        mu_L = 0;
        mu_R = -V_app;
        fL = 1 / (1 + exp((E - mu_L)/E_T));
        fR = 1 / (1 + exp((E - mu_R)/E_T));
        Current_Density_E(k) = T_E * (fL - fR);
    end


    coeff = 2 * q / 4.135e-15;
    I_list(iv) = coeff * trapz(ES, Current_Density_E);

    if mod(iv, 5) == 0, fprintf('V = %.2f V terminé.\n', V_app); end
end


figure(10); clf;
semilogy(V_list, abs(I_list), 'b-o', 'LineWidth', 2, 'MarkerSize', 4);
xlabel('Bias Voltage (V)');
ylabel('|Current| (A)');
title('I-V Characteristic of MIIM Structure');
grid on;


V_target = V_sweep_max;
[~, idx_pos] = min(abs(V_list - V_target));
[~, idx_neg] = min(abs(V_list - (-V_target)));

I_pos = I_list(idx_pos);
I_neg = I_list(idx_neg);
RR = abs(I_pos) / abs(I_neg);

text(V_sweep_max/2, min(abs(I_list))*10, ...
    sprintf('RR @ %.1fV = %.2f', V_target, RR), ...
    'FontSize', 12, 'BackgroundColor', 'white', 'EdgeColor', 'black');

fprintf('\n=== RÉSULTATS ===\n');
fprintf('Courant à +%.1f V : %e A\n', V_target, I_pos);
fprintf('Courant à -%.1f V : %e A\n', V_target, I_neg);
fprintf('Rectification Ratio : %.2f\n', RR);
