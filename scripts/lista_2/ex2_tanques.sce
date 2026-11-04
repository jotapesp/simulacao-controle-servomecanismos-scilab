// =====================================================================
// Exercicio 2)
// Sistema de dois tanques comunicantes.
//
// Precisa rodar ex2_blocos.ssp antes desse script
// =====================================================================

// parametros do enunciado
C1 = 1.5;
C2 = 0.8;
R1 = 2;
R2 = 3;

// o diagrama é aberto com o comando abaixo e simulado pela interface:
// xcos

// ---------------------------------------------------------------------
// Niveis dos dois tanques
// ---------------------------------------------------------------------
scf();
plot(h1.time, h1.values, 'r', h2.time, h2.values, 'b');
xlabel("t(s)", "fontsize", 3);
ylabel("nível(m)", "fontsize", 3);
title("Níveis dos tanques", "fontsize", 3);
legend(["h_1(t)", "h_2(t)"], 4);
xgrid();

// ---------------------------------------------------------------------
// Diferenca entre os niveis
// ---------------------------------------------------------------------
scf();
plot(h1.time, h1.values - h2.values, 'g');
xlabel("t(s)", "fontsize", 3);
ylabel("h_1 - h_2(m)", "fontsize", 3);
title("Diferença de níveis entre os tanques", "fontsize", 3);
legend("h_1(t) - h_2(t)", 4);
xgrid();

// ---------------------------------------------------------------------
// Vazoes nas valvulas
// ---------------------------------------------------------------------
scf();
plot(q1.time, q1.values, 'r', q2.time, q2.values, 'b');
xlabel("t [s]", "fontsize", 3);
ylabel("vazão [m³/s]", "fontsize", 3);
title("Vazões nas válvulas", "fontsize", 3);
legend(["q_1(t)", "q_2(t)"], 4);
xgrid();

// ---------------------------------------------------------------------
// Verificacao da relacao  q - q2 = C1*h1' + C2*h2'
// As derivadas vem direto do diagrama (entradas dos integradores),
// nao sao calculadas numericamente a partir dos niveis.
// ---------------------------------------------------------------------
scf();
plot(q2.time, 0.1 - q2.values, "k-", ..
     x1p.time, C1*x1p.values + C2*x2p.values, "r--");
xlabel("t [s]", "fontsize", 3);
ylabel("vazão [m³/s]", "fontsize", 3);
title("Verificação: q - q_2 = C_1 h_1'' + C_2 h_2''", "fontsize", 3);
legend(["q - q_2", "C_1·h_1'' + C_2·h_2''"], 1);
xgrid;
