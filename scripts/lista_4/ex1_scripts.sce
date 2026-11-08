
clear;
clc;
clf();

// =====================================================
// EXERCICIO 1 - LISTA 4
// Sistema massa-mola-amortecedor com base movel
// =====================================================

// 1. Parametros do sistema
m = 400;       // Massa (kg)
b = 200;       // Amortecimento (N.s/m)
k = 100;       // Constante elastica (N/m)
A = 0.05;      // Amplitude da base (m)

// Frequencias de excitacao (rad/s)
w1 = 10;
w2 = 20;

// Coeficiente de amortecimento e frequencia natural
wn = sqrt(k/m);
zeta = b/(2*sqrt(m*k));

disp("Coeficiente de amortecimento:");
disp(zeta);

disp("Frequencia natural (rad/s):");
disp(wn);

disp("Frequencia natural (Hz):");
disp(wn/(2*%pi));

// Verificacao da condicao do enunciado
if zeta > 0.4 & zeta < 0.8 then
    disp("Condicao de amortecimento satisfeita.");
else
    disp("Condicao de amortecimento NAO satisfeita.");
end

// 2. Funcao de transferencia Y(s)/U(s)
// Y(s)/U(s) = (b*s + k)/(m*s^2 + b*s + k)

s = poly(0, "s");

num = b*s + k;
den = m*s^2 + b*s + k;

G = syslin("c", num, den);

disp("Funcao de transferencia:");
disp(G);

// 3. Ganhos nas frequencias solicitadas
H1 = horner(num, %i*w1)/horner(den, %i*w1);
H2 = horner(num, %i*w2)/horner(den, %i*w2);

ganho1 = abs(H1);
ganho2 = abs(H2);

fase1 = atan(imag(H1), real(H1))*180/%pi;
fase2 = atan(imag(H2), real(H2))*180/%pi;

amplitude1 = A*ganho1;
amplitude2 = A*ganho2;

disp(" ");
disp("===== RESULTADOS EM REGIME PERMANENTE =====");

disp("Ganho para w = 10 rad/s:");
disp(ganho1);

disp("Fase para w = 10 rad/s (graus):");
disp(fase1);

disp("Amplitude de saida para w = 10 rad/s (m):");
disp(amplitude1);

disp("Ganho para w = 20 rad/s:");
disp(ganho2);

disp("Fase para w = 20 rad/s (graus):");
disp(fase2);

disp("Amplitude de saida para w = 20 rad/s (m):");
disp(amplitude2);

// 4. Diagramas de Bode

scf(1);
clf();
bode(G);
xtitle("Diagramas de Bode - Sistema mecanico", ...
       "Frequencia angular (rad/s)", ...
       "Magnitude (dB) / Fase (graus)");

// 5. Simulacao temporal com deff
// Excitacoes senoidais da base

deff("y = entrada10(t)", ...
     "y = 0.05*sin(10*t)");

deff("y = entrada20(t)", ...
     "y = 0.05*sin(20*t)");

// Vetor de tempo
dt = 0.001;
t = 0:dt:5;

// Entradas
u1 = entrada10(t);
u2 = entrada20(t);

// Respostas completas, incluindo o transitório
y1 = csim(u1, t, G);
y2 = csim(u2, t, G);

// 6. Graficos das respostas temporais

scf(2);
clf();

// Excitacao de 10 rad/s
subplot(2,1,1);
plot(t, u1, "b--");
plot(t, y1, "r");
xgrid();
xtitle("Resposta para w = 10 rad/s", ...
       "Tempo (s)", "Deslocamento (m)");
legend("Entrada da base", "Deslocamento da massa");

// Excitacao de 20 rad/s
subplot(2,1,2);
plot(t, u2, "b--");
plot(t, y2, "r");
xgrid();
xtitle("Resposta para w = 20 rad/s", ...
       "Tempo (s)", "Deslocamento (m)");
legend("Entrada da base", "Deslocamento da massa");

// 7. Comparacao das amplitudes simuladas
// Considera os ultimos 20% do intervalo como trecho final.

n = length(t);
i0 = round(0.8*n);

amp_sim1 = (max(y1(i0:n)) - min(y1(i0:n)))/2;
amp_sim2 = (max(y2(i0:n)) - min(y2(i0:n)))/2;

disp(" ");
disp("===== COMPARACAO DAS AMPLITUDES =====");

disp("Amplitude teorica para 10 rad/s (m):");
disp(amplitude1);

disp("Amplitude simulada no trecho final para 10 rad/s (m):");
disp(amp_sim1);

disp("Amplitude teorica para 20 rad/s (m):");
disp(amplitude2);

disp("Amplitude simulada no trecho final para 20 rad/s (m):");
disp(amp_sim2);
