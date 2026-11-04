// ============================================================
// Lista 4 - Exercício 2 - Controle e Servomecanismos
// Planta g(s) = 1/(s(s+1)(s+2)), controlador c(s) = 6,
// realimentação unitária.
// ============================================================
clear; clc;

s = poly(0, 's');
g = syslin('c', 1 / (s*(s+1)*(s+2)));   // planta
c = 6;                                   // controlador proporcional
L = c * g;                               // FTMA: L(s) = 6/(s(s+1)(s+2))

// --- Margens (confirmação numérica da estabilidade marginal) ---
[gm, fg]  = g_margin(L);   // fg em Hz
[pm, fpm] = p_margin(L);
mprintf("Margem de ganho = %.4f dB em w = %.4f rad/s\n", gm, 2*%pi*fg);
mprintf("Margem de fase  = %.4f graus\n", pm);

// --- a. DIAGRAMA DE NYQUIST ---
// Na figura: os números sobre a curva são frequências em Hz;
//  azul tracejado  = semieixo real negativo (fase -180°)
//  azul grosso     = segmento da origem até onde a curva cruza -180°
//  vermelho tracej.= círculo unitário (achatado pela escala do eixo y)
//  vermelho grosso = segmento da origem até onde a curva cruza |L| = 1
// Para c = 6, os dois cruzamentos caem no ponto -1 (w = sqrt(2) rad/s).
scf(0); clf();
show_margins(L, 'nyquist');
xtitle("Diagrama de Nyquist - c(s) = 6 (Estabilidade Marginal)");

// --- b. DIAGRAMA DE BLACK ---
scf(1); clf();
black(L);
xtitle("Diagrama de Black/Nichols - c(s) = 6 (Estabilidade Marginal)");

// --- c. SIMULAÇÃO PARA u(t) = sen(sqrt(2) t) ---
// Malha fechada: H(s) = L/(1+L) = 6/(s^3 + 3s^2 + 2s + 6)
//                     = 6/((s+3)(s^2+2))  -> polos em -3 e +-j*sqrt(2)
H = syslin('c', 6 / (s^3 + 3*s^2 + 2*s + 6));
disp("Polos de malha fechada:", roots(H.den));

w  = sqrt(2);          // rad/s
f  = w / (2*%pi);      // = 0,225 Hz  (período 2*pi/w = 4,44 s)
Tlim = %pi / w;        // limite de Nyquist: T < 1/(2f) = pi/sqrt(2) = 2,22 s
mprintf("f = %.4f Hz, T_limite (Nyquist) = %.4f s\n", f, Tlim);

// Tempo de amostragem adotado: fs = 20 f (20 amostras por período,
// 10x a taxa de Nyquist 2f). ws = 20*sqrt(2) = 28,3 rad/s,
// fs = 4,5 Hz, T = 2*pi/ws = 0,222 s. Também T < 1/3 s, constante de
// tempo do polo mais rápido (s = -3).
T  = 2*%pi / (20*w);
tf = 60;
t  = 0:T:tf;
u  = sin(w*t);
y  = csim(u, t, H);

// Envoltória teórica da ressonância: |y| ~ (3*sqrt(2))/(2*sqrt(11)) * t
A = 3*sqrt(2) / (2*sqrt(11));   // ≈ 0,64

scf(2); clf();
plot(t, u, 'b', t, y, 'r', t, A*t, 'k--', t, -A*t, 'k--');
xgrid();
xtitle("Resposta em malha fechada a u(t) = sen(sqrt(2) t), T = 0,222 s", "t (s)", "amplitude");
legend(["u(t)"; "y(t)"; "envoltória +-0,64 t"], 2);

// --- Efeito do tempo de amostragem (teorema da amostragem) ---
Ts = [T, Tlim, 3];
tit = ["T = 0,222 s (fs = 20f, adequado)"; "T = pi/sqrt(2) = 2,22 s (limite de Nyquist)"; "T = 3 s (aliasing)"];
scf(3); clf();
for k = 1:3
    tk = 0:Ts(k):tf;
    uk = sin(w*tk);
    yk = csim(uk, tk, H);
    subplot(3, 1, k);
    plot(tk, uk, 'b.-', tk, yk, 'r.-');
    xgrid();
    xtitle(tit(k), "t (s)", "amplitude");
end
