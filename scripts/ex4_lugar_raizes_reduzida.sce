// =====================================================================
// Exercicio 4.c) e 4.d) 
// (c) Reducao de ordem da planta pelo met odo polos dominantes
// (d) Lugar das Raizes da planta reduzida e analise da regiao Omega
// =====================================================================
clear; clc;

s = %s;

// =====================================================================
// Planta completa
// =====================================================================
coef_g = [16875, 27000, 15075, 5720, 817, 48, 1];
deng   = poly(coef_g, 's', 'coeff');
g      = syslin('c', 1, deng);

disp(roots(deng));      // -1 ; -1+-2i ; -15 (triplo)

// =====================================================================
// (c) REDUCAO DE ORDEM
// Buscou-se eliminar os polos com parte real cinco
// vezes acima da dos polos mantidos.
//    polos mantidos : -1 e -1+-j2   (Re = -1)
//    polo eliminado : -15 (triplo)  (15 >= 5 x 1)
// =====================================================================

// seguindo os procedimentos das aulas
dg    = derivat(g);
dg2   = derivat(dg);

// derivadas de g em s = 0 (lado direito do sistema linear)
g_0   = horner(g,   0);  //  5.9259e-05
dg_0  = horner(dg,  0);  // -9.4815e-05
dg2_0 = horner(dg2, 0);  //  1.9753e-04

disp([g_0; dg_0; dg2_0]);

// polos mantidos
p1_dom_planta = -1 + %i*2;
p2_dom_planta = -1 - %i*2;
p3_dom_planta = -1;

// matriz:  linha i, coluna j  =  -(i-1)! / p_j^i
mat_pi = [ -1/p1_dom_planta,   -1/p2_dom_planta,   -1/p3_dom_planta   ;
           -1/p1_dom_planta^2, -1/p2_dom_planta^2, -1/p3_dom_planta^2 ;
           -2/p1_dom_planta^3, -2/p2_dom_planta^3, -2/p3_dom_planta^3 ];
mat_0  = [g_0; dg_0; dg2_0];

disp(mat_pi);

// novos residuos bi
br = mat_pi \ mat_0; // aqui nao precisou usar pdiv pois deu certo de cara assim
disp(br); // -4.1481e-05 +- 1.8765e-05i  ;  9.0864e-05

// gr(s) = somatorio br_j/(s - p_j) 
gr = br(1)/(s-p1_dom_planta) + br(2)/(s-p2_dom_planta) + br(3)/(s-p3_dom_planta);
gr = syslin('c', poly(real(coeff(gr.num)),'s','coeff'), ...
                 poly(real(coeff(gr.den)),'s','coeff'));
disp(gr);
//        7.9014e-06 s^2 - 5.9259e-05 s + 2.9630e-04
//   gr = ------------------------------------------
//                  s^3 + 3s^2 + 7s + 5

// verificacoes
disp(horner(gr, 0) - g_0);    // ~0
disp(horner(derivat(gr), 0) - dg_0);   // ~0
disp(horner(derivat(derivat(gr)), 0) - dg2_0);// ~0
disp(roots(gr.num));  // 3.7499 +- 4.8412i

// comparacao com a planta completa
t = 0:0.01:8;
scf();
plot(t, csim('step', t, g),  'b');
plot(t, csim('step', t, gr), 'r--');
legend(['planta completa (6a ordem)'; 'reduzida (3a ordem)'], 4);
xgrid(); xtitle('4c - reducao de ordem', 't [s]', 'y(t)');

// =====================================================================
// (d) LUGAR DAS RAIZES DA PLANTA REDUZIDA
// Seguindo mesmo procedimento
// =====================================================================

// i. polos e zeros 
// 3 polos (-1, -1+-j2), 2 zeros (3.7499 +- j4.8412, no semiplano DIREITO)
// n - m = 1 -> um unico ramo tende ao infinito

// ii. LR sobre o eixo real
// unico polo real: -1 ; nenhum zero real
//   s > -1 : 0 a direita -> par   -> nao pertence
//   s < -1 : 1 a direita -> impar -> PERTENCE
// LR real = (-inf, -1]

// iii. assintotas
n_r = 3;  m_r = 2;
sigma_ar = real(sum(roots(gr.den)) - sum(roots(gr.num))) / (n_r - m_r);
theta_ar = (2*(1:(n_r-m_r)) - 1)*180/(n_r-m_r);
disp(sigma_ar);   // -10.4999
disp(theta_ar);   // 180

// iv. pontos de ramificacao
pol_quebra    = derivat(gr.num)*gr.den - derivat(gr.den)*gr.num;
candidatos_sb = roots(pol_quebra);
disp(candidatos_sb);    // -1.0950 +- 1.1580i  ;  8.5950 +- 6.6508i

// testando  K real e positivo
for i = 1:size(candidatos_sb,'*')
    Ki = -1/horner(gr, candidatos_sb(i));
    mprintf("sb = %9.4f %+9.4fi   K = %14.2f %+14.2fi\n", ...
            real(candidatos_sb(i)), imag(candidatos_sb(i)), real(Ki), imag(Ki));
end
// todos dao K COMPLEXO -> NAO HA ponto de ramificacao neste LR

// vi. angulo de partida de p = -1 + j2
//  theta_p = 180 - soma(angulos dos demais polos) + soma(angulos dos zeros)
pr = p1_dom_planta;
outros = [p2_dom_planta; p3_dom_planta];
zeros_gr = roots(gr.num);
ang_pol = sum(atan(imag(pr-outros),   real(pr-outros)))   * 180/%pi;   // 180.00
ang_zer = sum(atan(imag(pr-zeros_gr), real(pr-zeros_gr))) * 180/%pi;   // -24.34
theta_pr = pmodulo(180 - ang_pol + ang_zer + 180, 360) - 180;
disp(theta_pr); // -24.34 graus

// v. cruzamento do eixo imaginario/estabilidade
k = poly(0, 'k');
rou_r= routh_t(gr, k);
col1_r = rou_r(:,1);
numeradores_routh_r = col1_r.num;

kr_2 = roots(numeradores_routh_r(2));   // -379687.50
kr_3 = roots(numeradores_routh_r(3));   // -931076.47 ; +36701.47
kr_4 = roots(numeradores_routh_r(4));   // -931076.47 ; +36701.47 ; -16875
disp(kr_2); disp(kr_3); disp(kr_4);

kr_lim = kr_3(2);                       // 36701.47
disp(kr_lim);

// polinomio auxiliar: a linha 3 (do s^1) e' a que zera, entao a
// auxiliar sai da linha 2 (do s^2):  A(s) = a2*s^2 + a0
poli_aux_ra = horner(rou_r(2,1), kr_lim);
poli_aux_rb = horner(rou_r(2,2), kr_lim);
polinomio_auxiliar_routh_r = poly([poli_aux_rb, 0, poli_aux_ra], 's', 'coeff');
cruza_imaginario_gr = roots(polinomio_auxiliar_routh_r);
disp(cruza_imaginario_gr);              // +- 2.1966106 i

// =====================================================================
// (d) ANALISE DA REGIAO OMEGA: por que o compensador P nao serve
// =====================================================================

// ganho que leva o polo real a fronteira sigma = 3
Kr_sigma = -1/horner(gr, -3);
disp(Kr_sigma);  // 29347.77

// polos de malha fechada nesse ganho
disp(roots(gr.den + Kr_sigma*gr.num));
//   -3.0000  ;  -0.1159 +- 2.1335i  (o par complexo domina e está fora)

// grafico: LR sobre a regiao Omega
b40 = 40*%pi/180;   sg = 3;   Lp = 15;

scf();
xv = [-sg, -Lp, -Lp, -sg];
yv = [ sg*tan(b40), Lp*tan(b40), -Lp*tan(b40), -sg*tan(b40)];
xfpoly(xv, yv);   e = gce();   e.background = color(220,235,250);

plot([0 -Lp], [0  Lp*tan(b40)], 'k--');
plot([0 -Lp], [0 -Lp*tan(b40)], 'k--');
plot([-sg -sg], [-Lp*tan(b40) Lp*tan(b40)], 'k--');
plot([0 0], [-Lp*tan(b40) Lp*tan(b40)], 'black');

evans(gr, 40000);
