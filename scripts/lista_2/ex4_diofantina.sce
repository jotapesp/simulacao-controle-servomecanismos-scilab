// =====================================================================
// Exercicio 4.a) 
//
// Projeto de compensador via equacao diofantina para
//     g(s) = 1 / (s^6 + 48s^5 + 817s^4 + 5720s^3 + 15075s^2 + 27000s + 16875)
// com: erro estatico de posicao nulo, dominancia de 2a ordem,
//      beta < 40 graus e sigma > 3.
// =====================================================================
clear; clc;

s = %s;

// =====================================================================
// 1. Planta
// =====================================================================
coef_g = [16875, 27000, 15075, 5720, 817, 48, 1];
deng   = poly(coef_g, 's', 'coeff');
g      = syslin('c', 1, deng);

disp(roots(deng));      // -1, -1+-2i e -15 (triplo)

// =====================================================================
// 2. Planta aumentada
// A planta é tipo 0 (deng(0) = 16875), entao o erro de posicao so
// zera com um integrador no compensador. 
// Procedimento da Aula 13:
// trocar g por g~ = g/s (passa o polo na origem do compensador para a
// planta fictícia e resolve a diofantina para c~ = sc)
// Logo d_g~ = s*deng, de grau ng = 7 (e nao 6).
// =====================================================================
dengt = s*deng;

disp(degree(dengt)); // 7

// =====================================================================
// 3. Polos desejados de malha fechada
// Dominantes: escolhido sigma = 3.2 (> 3) e beta = 35 (deg < 40)
//   Im = sigma*tan(beta) = 2.2407 ; zeta = cos(beta) = 0.819
// Rapidos (10): duplos em -16, -17, -18, -19, -20
//   critério tomado foi que sejam pelo menos 5x "maiores" que os dominantes
// =====================================================================
sig = 3.2;
bet = 35*%pi/180;
wd  = sig*tan(bet);

p1 = -sig + wd*%i;
p2 = -sig - wd*%i;
polos = [p1, p2, -16,-16, -17,-17, -18,-18, -19,-19, -20,-20];

q = poly(polos, 's', 'roots');

disp(wd);               // 2.2406641
disp(cos(bet));         // 0.8191520  (zeta)
disp(degree(q));        // 12

// =====================================================================
// 4. Dimensoes e matriz de Sylvester  (desenvolvimento teorico)
//
// Diofantina:  n_c~*n_g + d_c~*d_g~ = q, com n_g = 1
// S_D (sylvester denominador) => (nc+ng+1)x(nc+1) ; 
// S_N (numerador) => (nc+ng+1)x(mc+1)
// linhas = nc + ng + 1 = degree(q) + 1
// S quadrada  <=>  mc = ng - 1 = 6
// c = n_c~/(s*d_c~) propria  =>  nc + 1 >= mc  =>  nc = 5
// => S_D 13x6 , S_N 13x7 , S 13x13
// =====================================================================
nc = 5;  mc = 6;  L = nc + 7 + 1;

alfa = coeff(dengt);

S = zeros(L, L);
for j = 1:nc+1, S(j:j+7, j) = alfa'; end
for j = 1:mc+1, S(j, nc+1+j) = 1;    end

disp(S);
disp(det(S));
disp(rcond(S));  // ~1e-17 : mal condicionada, S\q e' inseguro

// =====================================================================
// 5. Solucao da diofantina por divisao polinomial (pdiv pois estava dando
// erro " mal condicionada")
//
// Com n_g = 1 a diofantina vira  q = d_c~*d_g~ + n_c~ , deg(n_c~) < deg(d_g~)
// pdiv devolve [resto, quociente].
// =====================================================================
[numc, denc] = pdiv(q, dengt);

disp(denc); // d_c~
disp(numc);  // n_c~

// =====================================================================
// 6. Compensador e malha fechada
// =====================================================================
c    = syslin('c', numc, s*denc);
FTMF = (c*g) /. 1;

disp(roots(denc)); // compensador estavel (nao  garantido pelo metodo)
disp(roots(numc));      // zeros de malha fechada
disp(horner(FTMF, 0));  // 1 -> ganho unitario -> erro de posicao nulo

// =====================================================================
// 7. Resposta ao degrau
// Comparada com a 2a ordem pura dos polos dominantes.
// =====================================================================
t     = 0:0.001:2.5; // 2.5s escolhido pelo ts
y_dio = csim('step', t, FTMF);

wn = sqrt(sig^2 + wd^2);
zt = sig/wn;
s2 = syslin('c', wn^2, s^2 + 2*zt*wn*s + wn^2);
y_2a = csim('step', t, s2);

disp(y_dio($));     // 1.0003  (valor final)
disp(max(y_dio));  // 1.8239  (sobressinal 82.4%)
disp(t(find(y_dio == max(y_dio))));  // 0.324 s (instante de pico)
disp(max(y_2a)); // 1.0112  (previsto por zeta = 0.819)

scf();
plot(t, y_dio, 'b');
plot(t, y_2a, 'k--');
legend(['diofantina (malha fechada completa)'; '2a ordem dos polos dominantes']);
xgrid();
xtitle('4a - resposta ao degrau', 't [s]', 'y(t)');

// =====================================================================
// Zona de desempenho garantido (Omega) no plano s
// aqui plota os polos e zeros na região omega
// =====================================================================
b40 = 40*%pi/180;
sg  = 3;
Lp  = 25;

scf();
xv = [-sg, -Lp, -Lp, -sg];
yv = [ sg*tan(b40), Lp*tan(b40), -Lp*tan(b40), -sg*tan(b40)];
xfpoly(xv, yv);
e = gce();  e.background = color(220,235,250);

plot([0 -Lp], [0  Lp*tan(b40)], 'k--');   // raio beta = +40 deg
plot([0 -Lp], [0 -Lp*tan(b40)], 'k--');   // raio beta = -40 deg
plot([-sg -sg], [-Lp*tan(b40) Lp*tan(b40)], 'k--');   // reta sigma = 3

plot(real(polos), imag(polos), 'bx');  // 12 polos de malha fechada
z = roots(numc);
plot(real(z), imag(z), 'ro');  // 6 zeros de malha fechada
plot([-sig -sig], [wd -wd], 'k*'); // par dominante

legend(['Omega'; 'beta = 40 deg'; 'beta = -40 deg'; 'sigma = 3'; ...
        'polos MF'; 'zeros MF'; 'par dominante']);
xgrid();
xtitle('Plano s - zona de desempenho garantido', 'Re', 'Im');
