// =====================================================================
// Exercicio 4(b) e 4(d) - Lista de Exercicios 2
// Controle e Servomecanismos - UFMS - Prof. Victor L. Yoshimura
//
// Projeto do compensador via lugar das raizes, nas duas plantas.
// Estrutura adotada: c(s) = kc (s - z) / ( s (s - p) )
//   o termo 1/s e' o integrador, exigido pelo erro de posicao nulo
//   o par (z,p) e' a parcela de avanco, que fornece a deficiencia angular
// =====================================================================
clear; clc;

s = %s;

// ---------------------------------------------------------------------
// Plantas
// ---------------------------------------------------------------------
coef_g = [16875, 27000, 15075, 5720, 817, 48, 1];
deng   = poly(coef_g, 's', 'coeff');
g      = syslin('c', 1, deng);

num_r  = [2.9630e-04, -5.9259e-05, 7.9014e-06];   // ordem crescente
den_r  = [5, 7, 3, 1];
numg_r = poly(num_r, 's', 'coeff');
deng_r = poly(den_r, 's', 'coeff');
gr     = syslin('c', numg_r, deng_r);

// ---------------------------------------------------------------------
// 1. Raiz dominante desejada
//    sigma = 3.2 e beta = 35 graus, a mesma do item (a)
// ---------------------------------------------------------------------
sd = -3.2 + %i*2.24;

// ---------------------------------------------------------------------
// 2. Deficiencia angular
//    malha aberta ja conhecida = g(s)/s  (planta + integrador)
//    deficit = 180 - arg[ g(sd)/sd ]
// ---------------------------------------------------------------------
ar     = horner(g,  sd);
ar2    = horner(gr, sd);

ang_ar  = atan(imag(ar),  real(ar))  * 180/%pi;      // -97.9271
ang_ar2 = atan(imag(ar2), real(ar2)) * 180/%pi;      // -90.6975
ang_sd  = atan(imag(sd),  real(sd))  * 180/%pi;      // 145.0080

deficit   = 180 - (ang_ar  - ang_sd);
deficit   = pmodulo(deficit + 180, 360) - 180;       // 62.9351
deficit_r = 180 - (ang_ar2 - ang_sd);
deficit_r = pmodulo(deficit_r + 180, 360) - 180;     // 55.7054

disp(deficit);  disp(deficit_r);

// o pmodulo reduz o angulo para (-180, 180]. Sem ele o resultado sai
// 422.9351, que e' o mesmo angulo mas nao serve para ler como "quanto falta"

// ---------------------------------------------------------------------
// 3. Escolha do zero
//    O zero do compensador vira zero da malha fechada, entao e' ele que
//    se escolhe; o polo p e' o preco angular. Adotou-se z = -3.5,
//    proximo da parte real de sd.
// ---------------------------------------------------------------------
z   = -3.5;
ar3 = horner((s - z), sd);
angd_1 = atan(imag(ar3), real(ar3)) * 180/%pi;       // 82.3719

// ---------------------------------------------------------------------
// 4. Angulo que sobra para o polo, e o polo
//    arg(sd - z) - arg(sd - p) = deficit
//    tan(ang_p) = Im(sd) / ( Re(sd) - p )
// ---------------------------------------------------------------------
ang_p  = angd_1 - deficit;                           // 19.4367
ang_pr = angd_1 - deficit_r;                         // 26.6664

disp(ang_p);  disp(ang_pr);

p  = real(sd) - imag(sd)/tan(ang_p  * %pi/180);      // -9.5478
pr = real(sd) - imag(sd)/tan(ang_pr * %pi/180);      // -7.6603

disp(p);  disp(pr);

// ---------------------------------------------------------------------
// 5. Compensadores (sem o ganho) e ganho pelo criterio de magnitude
//    | kc c(sd) g(sd)/sd | = 1
// ---------------------------------------------------------------------
c  = syslin('c', (s - z), (s - p));
cr = syslin('c', (s - z), (s - pr));

val   = horner(c,  sd) * horner(g/s,  sd);
val_r = horner(cr, sd) * horner(gr/s, sd);

kc   = 1/abs(val);                                   // 669058.92
kc_r = 1/abs(val_r);                                 // 492148.17

disp(kc);  disp(kc_r);

// compensadores completos
comp   = syslin('c', kc  *(s - z), s*(s - p));
comp_r = syslin('c', kc_r*(s - z), s*(s - pr));

disp(comp);  disp(comp_r);

// ---------------------------------------------------------------------
// 6. VERIFICACAO - malha fechada
//    As condicoes de angulo e de magnitude posicionam UM par de polos
//    em sd. Elas nao dizem nada sobre os demais, entao e' obrigatorio
//    conferir todos.
// ---------------------------------------------------------------------
FTMF   = (comp   * g ) /. 1;
FTMF_r = (comp_r * gr) /. 1;

disp(roots(FTMF.den));
// -17.462 ; -14.574 +- 3.135i ; -6.845 ; -3.200 +- 2.240i ; +1.154 +- 2.109i

disp(roots(FTMF_r.den));
// -6.420 ; -3.200 +- 2.240i ; +1.080 +- 2.011i

disp(max(real(roots(FTMF.den))));      // +1.154
disp(max(real(roots(FTMF_r.den))));    // +1.080

// O par dominante ficou exatamente em sd, como projetado, mas em ambos
// os casos sobra um par no semiplano DIREITO: as malhas fechadas sao
// instaveis. O ganho exigido esta muito acima do limite de estabilidade
// do sistema compensado (75764 na completa e 59042 na reduzida).

// ---------------------------------------------------------------------
// 7. Limite de estabilidade do sistema compensado
//    A equacao caracteristica e'  s(s-p)*deng + kc*(s-z) = 0 .
//    Varrendo kc e olhando a maior parte real das raizes, ve-se onde
//    o sistema cruza para o semiplano direito.
// ---------------------------------------------------------------------
disp("planta completa:");
for kk = [50000 70000 75000 76000 100000 669059]
    disp([kk, max(real(roots( s*(s-p)*deng + kk*(s-z) )))]);
end
// o cruzamento ocorre entre 75000 e 76000  ->  K_lim = 75764

disp("planta reduzida:");
for kk = [40000 55000 59000 60000 100000 492148]
    disp([kk, max(real(roots( s*(s-pr)*deng_r + kk*(s-z)*numg_r )))]);
end
// o cruzamento ocorre entre 59000 e 60000  ->  K_lim = 59042

// Conclusao: os ganhos exigidos pelo criterio de magnitude (669059 e
// 492148) estao cerca de 9 vezes acima dos limites de estabilidade.
