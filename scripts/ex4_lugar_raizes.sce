// =====================================================================
// Exercicio 4.b)
// Lugar das Raizes para
//     g(s) = 1 / (s^6 + 48s^5 + 817s^4 + 5720s^3 + 15075s^2 + 27000s + 16875)
// seguindo os procedimentos completos da aula 14, o que foi feito no scilab 
// segue cálculos matemáticos e no final plotagem com evans():
// =====================================================================
clear; clc;

s = %s;
// =====================================================================
// i. Polos e zeros
// =====================================================================
coef_g = [16875, 27000, 15075, 5720, 817, 48, 1];   // ordem CRESCENTE
deng   = poly(coef_g, 's', 'coeff');
g      = syslin('c', 1, deng);

disp(roots(deng));      // -1 ; -1+-2i ; -15 (triplo)

// n = 6 polos, m = 0 zeros  ->  6 ramos, todos tendendo ao infinito
n = 6;  m = 0;

// =====================================================================
// ii. LR sobre o eixo real
//     Polos reais: -1 e -15 (triplo). Polos complexos nao contam
//     (vem aos pares). Numero impar de polos+zeros a direita:
//        s > -1        : 0  -> par    -> NAO pertence
//        -15 < s < -1  : 1  -> impar  -> PERTENCE
//        s < -15       : 4  -> par    -> NAO pertence
//     Logo o LR no eixo real e' o segmento [-15, -1].
// =====================================================================

// =====================================================================
// iii. Assintotas
// =====================================================================
sigma_a = -coef_g($-1)/coef_g($) / (n-m);
theta_a = (2*(1:(n-m)) - 1)*180/(n-m);

disp(sigma_a); // -8
disp(theta_a); // 30  90  150  210  270  330

// =====================================================================
// iv. Pontos de ramificacao e seus ganhos
//     Como n_g = 1, a eq. caracteristica é deng(s) + K = 0, logo
//     K = -deng(s)  e  dK/ds = 0  <=>  deng'(s) = 0 (d_deng).
// =====================================================================
d_deng = derivat(deng); // deng'(s)
cand   = roots(d_deng);  //candidatos

disp(cand);  // -15 ; -15 ; -7.8097 ; -1.0952 +- 1.1669i

K_15 = -horner(deng, -15);  // 0 -> é o proprio polo, descarta
K_7  = -horner(deng, -7.8096726);  // 127514.13  -> VALIDO
K_C  = -horner(deng, -1.0951637 + %i*1.1668519); // 2081 - j8190 -> complexo, descarta

disp(K_15);  disp(K_7);  disp(K_C);

// Observacao: K_7 = 127514 esta ACIMA do K limite de estabilidade
// (ver item vi). A quebra existe no tracado mas fora da faixa estável.

// =====================================================================
// vi. Angulos de partida
//     theta_p = 180 - soma(angulos dos demais polos ate p) + soma(zeros)
// =====================================================================

// --- do polo complexo p = -1 + j2 ---
p      = -1 + 2*%i;
outros = [-1; -1-2*%i; -15; -15; -15];              // os 5 demais polos
ang_p  = sum(atan(imag(p-outros), real(p-outros))) * 180/%pi;
theta_p = pmodulo(180 - ang_p + 180, 360) - 180;

disp(ang_p);   // 204.39
disp(theta_p);  // -24.39 graus

// --- do polo triplo em -15 (multiplicidade = 3) ---
// theta_k = [180 + 360(k-1) - soma(angulos)] / 3 ,  k = 1,2,3
pt   = -15;   mult = 3;
o15  = [-1; -1+2*%i; -1-2*%i];  // exclui as 3 copias de -15
ang_15   = sum(atan(imag(pt-o15), real(pt-o15))) * 180/%pi;
theta_15 = pmodulo((180 + 360*(0:mult-1) - ang_15)/mult + 180, 360) - 180;

disp(theta_15);         // 0 ; 120 ; -120 graus

// =====================================================================
// v. Cruzamento do eixo imaginario - Routh-Hurwitz
// =====================================================================
k    = poly(0, 'k');
rou  = routh_t(g, k);
col1 = rou(:,1);

// Pega a primeira coluna importa apenas:
numeradores_routh = col1.num;

k_5 = roots(numeradores_routh(5));  // -1.05e6
k_6 = roots(numeradores_routh(6));  // -7.47e6 ; +37662.88
k_7 = roots(numeradores_routh(7));  // -7.47e6 ; +37662.88 ; -16875

disp(k_5);  disp(k_6);  disp(k_7);

// Intersecao das condicoes:  -16875 < k < 37662.88
// (o denominador de cada linha e' a linha anterior, entao analisando de
//  cima para baixo basta olhar os numeradores)
K_lim = 37662.88;

// polinomio auxiliar:
// A linha 6 (linha do s^1) é a que zera, entao a auxiliar sai da
// linha 5 (linha do s^2):  A(s) = d1*s^2 + d2
poli_aux_a = horner(rou(5,1), K_lim);  // primeiro coef da linha 5: 11076.566
poli_aux_b = horner(rou(5,2), K_lim);  // segundo coef da linha 5: 54537.880
polinomio_auxiliar_routh = poly([poli_aux_b, 0, poli_aux_a], 's', 'coeff');
cruza_imaginario = roots(polinomio_auxiliar_routh);

disp(polinomio_auxiliar_routh);
disp(cruza_imaginario);  // +- 2.218945 i

// =====================================================================
// Desenho LR
// =====================================================================
scf();  evans(g,  50000);   // k_max = 50000: mostra o cruzamento do eixo imaginario
scf();  evans(g, 150000);   // k_max = 150000: mostra o ponto de ramificacao em -7.81
