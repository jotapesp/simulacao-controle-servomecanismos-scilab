// Definicao da variavel s
s = poly(0, 's');

// Planta e Compensador
G = 2 / (0.5*s + 1);
C = (s + 4.5) / s;

// Malha Fechada: T(s) = (C*G) / (1 + C*G)
L = C * G;
T = syslin('c', L / (1 + L));

// Simulacao da resposta ao degrau
t = 0:0.01:5;
y = csim('step', t, T);

// Grafico
figure(1);
plot(t, y, 'b-', 'LineWidth', 2);
xgrid();
xtitle("Resposta ao Degrau - Sistema em Malha Fechada", "Tempo (s)", "Saida y(t)");