// ============================================================
// Lista 4 - Exercício 3 - Controle e Servomecanismos
// Tanques comunicantes: controlabilidade, observabilidade,
// realimentação de estados, observador e funções de Lyapunov.
// ============================================================
clear; clc; close;

// --- Parâmetros Físicos do Sistema ---
C1 = 1.5;   // Área do Tanque 1 (m^2)
C2 = 0.8;   // Área do Tanque 2 (m^2)
R1 = 2.0;   // Resistência da válvula 1 (s/m^2)
R2 = 3.0;   // Resistência da válvula 2 (s/m^2)

// --- Matrizes em Espaço de Estados ---
A = [-1/(R1*C1), 1/(R1*C1); ...
      1/(R1*C2), -(1/(R1*C2) + 1/(R2*C2))];
B = [1/C1; 0];
C = [0, 1];
D = 0;

sys = syslin('c', A, B, C, D);

// --- a. CONTROLABILIDADE E OBSERVABILIDADE ---
Pc = [B, A*B];
Po = [C; C*A];

det_Pc = det(Pc);
det_Po = det(Po);
rank_Pc = rank(Pc);
rank_Po = rank(Po);

mprintf("--- a. Controlabilidade e Observabilidade ---\n");
mprintf("Det(Pc) = %.4f, Posto(Pc) = %d\n", det_Pc, rank_Pc);
mprintf("Det(Po) = %.4f, Posto(Po) = %d\n", det_Po, rank_Po);

if rank_Pc == 2 then
    mprintf("O sistema é COMPLETAMENTE CONTROLÁVEL.\n");
end
if rank_Po == 2 then
    mprintf("O sistema é COMPLETAMENTE OBSERVÁVEL.\n");
end

// --- b. REALIMENTAÇÃO DE ESTADOS ---
// Polos alocados: zeta = 0.6, wn = 1.0 rad/s -> s = -0.6 +- j0.8
zeta = 0.6;
wn = 1.0;
p1 = -zeta*wn + %i*wn*sqrt(1 - zeta^2);
p2 = -zeta*wn - %i*wn*sqrt(1 - zeta^2);
poles_cl = [p1, p2];

// Alocação do ganho K
K = ppol(A, B, poles_cl);

// Ganho de pré-filtro Kr para erro estático nulo
A_cl = A - B*K;
dc_gain = C * inv(-A_cl) * B;
Kr = 1 / dc_gain;

mprintf("\n--- b. Realimentação de Estados ---\n");
mprintf("K  = [%.4f, %.4f]\n", K(1), K(2));
mprintf("Kr = %.4f\n", Kr);

// --- c. OBSERVADOR DE ESTADOS E REALIMENTAÇÃO DINÂMICA ---
// Polos do observador duas vezes mais à esquerda: s = -1.2 +- j1.6
p_obs1 = -1.2 + %i*1.6;
p_obs2 = -1.2 - %i*1.6;
poles_obs = [p_obs1, p_obs2];

// Alocação do ganho L via dualidade
L = ppol(A', C', poles_obs)';

mprintf("\n--- c. Observador de Estados ---\n");
mprintf("L = [%.4f; %.4f]\n", L(1), L(2));

// Matriz estendida do sistema combinado (Princípio da Separação)
A_comb = [A - B*K, B*K; zeros(2,2), A - L*C];
B_comb = [B*Kr; zeros(2,1)];
C_comb = [C, zeros(1,2)];
D_comb = 0;

sys_comb = syslin('c', A_comb, B_comb, C_comb, D_comb);

// --- d. SIMULAÇÃO E COMPARAÇÃO DAS RESPOSTAS ---
t = 0:0.01:15;
u_step = ones(t);

// Realimentação direta de estados
sys_sf = syslin('c', A_cl, B*Kr, C, D);
y_sf = csim('step', t, sys_sf);

// Realimentação dinâmica de saída (com erro de estimação inicial e0 = [-0.2; -0.1])
x0_comb = [0; 0; -0.2; -0.1];
y_dof = csim(u_step, t, sys_comb, x0_comb);

// Gráfico comparativo com fundo branco
f = gcf();
f.background = color("white");

plot(t, y_sf, 'r-', 'LineWidth', 2);
plot(t, y_dof, 'b--', 'LineWidth', 2);
plot(t, u_step, 'k:', 'LineWidth', 1.5);

a = gca();
a.background = color("white");
a.box = "on";
a.font_size = 2;
a.x_label.font_size = 3;
a.y_label.font_size = 3;
a.title.font_size = 4;

xlabel("Tempo (s)");
ylabel("Nível do Tanque 2 h_2(m)");
title("Comparação da Resposta ao Degrau - Exercício 3");
legend(["Realimentação de Estados Direta"; "Realimentação Dinâmica de Saída (Observador)"; "Referência r(t)"], 4);
xgrid(1);

// Salvamento automático da figura
xs2png(f, "resposta_tanques.png");

// --- e. FUNÇÕES DE LYAPUNOV ---
Q1 = eye(2,2);
P1 = lyap(A_cl, -Q1, 'c');

Q_comb = eye(4,4);
P_comb = lyap(A_comb, -Q_comb, 'c');

mprintf("\n--- e. Matrizes de Lyapunov P ---\n");
mprintf("P1 (Realimentação de Estados):\n");
disp(P1);
mprintf("Autovalores de P1 (devem ser > 0):\n");
disp(spec(P1));

mprintf("\nP_comb (Realimentação Dinâmica de Saída):\n");
disp(P_comb);
mprintf("Autovalores de P_comb (devem ser > 0):\n");
disp(spec(P_comb));