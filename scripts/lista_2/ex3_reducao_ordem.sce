// =====================================================================
// Exercicio 3)
//
// (a) simulacao da planta com u = t^i, i = 0,1,2
// (b) reducao do modelo a 5a, 4a e 3a ordens (polos dominantes, eq. 9.15)
// (c) simulacao das reducoes e comparacao com o modelo completo
// =====================================================================
clear; clc;

// ---------------------------------------------------------------------
// Planta
// ---------------------------------------------------------------------
s   = %s;
num = s + 6;
den = s^6 + 11*s^5 + 52*s^4 + 142*s^3 + 241*s^2 + 231*s + 90;
G   = syslin('c', num, den);

// ---------------------------------------------------------------------
// (a) Entradas u = t^i, i = 0,1,2
// ---------------------------------------------------------------------
t  = 0:0.01:10;
u1 = t;
u2 = t.^2;

y0 = csim('step', t, G);
y1 = csim(u1, t, G);
y2 = csim(u2, t, G);

// graficos da planta completa (item a)
scf(10); clf;
plot(t, y0, "k-");
xtitle("Resposta ao degrau (u = 1)", "t [s]", "y(t)"); xgrid;

scf(11); clf;
plot(t, y1, "k-");
xtitle("Resposta a rampa (u = t)", "t [s]", "y(t)"); xgrid;

scf(12); clf;
plot(t, y2, "k-");
xtitle("Resposta a parabola (u = t^2)", "t [s]", "y(t)"); xgrid;

// ---------------------------------------------------------------------
// Derivadas de g(s) em s = 0  (lado direito do sistema linear)
// ---------------------------------------------------------------------
G1 = derivat(G);
G2 = derivat(G1);
G3 = derivat(G2);
G4 = derivat(G3);

G_0  = horner(G , 0);
G1_0 = horner(G1, 0);
G2_0 = horner(G2, 0);
G3_0 = horner(G3, 0);
G4_0 = horner(G4, 0);

// ---------------------------------------------------------------------
// Polos da planta, em ordem de dominancia
//   -1, -1+2i, -1-2i  (Re = -1)   >   -2   >   -3, -3
// ---------------------------------------------------------------------
p1 = -1;  p2 = -2;  p3 = -1+2*%i;  p4 = -1-2*%i;  p5 = -3;  p6 = -3;

// =====================================================================
// (b) Reducoes
// =====================================================================

// ---------------------- 5a ordem: descarta um -3 ---------------------
polos5 = [p1, p2, p3, p4, p5];
n  = 5;
A5 = zeros(n,n);
for i = 1:n
    for j = 1:n
        A5(i,j) = -factorial(i-1)/polos5(j)^i;
    end
end
B5 = [G_0; G1_0; G2_0; G3_0; G4_0];
b5 = A5\B5;

G_5g = 0;
for i = 1:n
    G_5g = G_5g + b5(i)/(s - polos5(i));
end
nr5 = poly(clean(real(coeff(G_5g.num))), 's', 'coeff');
dr5 = poly(clean(real(coeff(G_5g.den))), 's', 'coeff');
G_5 = syslin('c', nr5, dr5);

// ---------------------- 4a ordem: descarta o outro -3 ----------------
polos4 = [p1, p2, p3, p4];
n  = 4;
A4 = zeros(n,n);
for i = 1:n
    for j = 1:n
        A4(i,j) = -factorial(i-1)/polos4(j)^i;
    end
end
B4 = [G_0; G1_0; G2_0; G3_0];
b4 = A4\B4;

G_4g = 0;
for i = 1:n
    G_4g = G_4g + b4(i)/(s - polos4(i));
end
nr4 = poly(clean(real(coeff(G_4g.num))), 's', 'coeff');
dr4 = poly(clean(real(coeff(G_4g.den))), 's', 'coeff');
G_4 = syslin('c', nr4, dr4);

// ---------------------- 3a ordem: descarta o -2 ----------------------
polos3 = [p1, p3, p4];
n  = 3;
A3 = zeros(n,n);
for i = 1:n
    for j = 1:n
        A3(i,j) = -factorial(i-1)/polos3(j)^i;
    end
end
B3 = [G_0; G1_0; G2_0];
b3 = A3\B3;

G_3g = 0;
for i = 1:n
    G_3g = G_3g + b3(i)/(s - polos3(i));
end
nr3 = poly(clean(real(coeff(G_3g.num))), 's', 'coeff');
dr3 = poly(clean(real(coeff(G_3g.den))), 's', 'coeff');
G_3 = syslin('c', nr3, dr3);

// ---------------------- resultados da reducao ------------------------
disp(b5, "residuos b (5a ordem):");   disp(G_5);
disp(b4, "residuos b (4a ordem):");   disp(G_4);
disp(b3, "residuos b (3a ordem):");   disp(G_3);

// verificacao: gr(0) deve ser igual a g(0) = 1/15 nas tres reducoes
mprintf("g(0) = %.7f | g5(0) = %.7f | g4(0) = %.7f | g3(0) = %.7f\n", ...
        G_0, horner(G_5,0), horner(G_4,0), horner(G_3,0));

// =====================================================================
// (c) Simulacao das reducoes e comparacao
// =====================================================================
y0_5g = csim('step', t, G_5);   y1_5g = csim(u1, t, G_5);   y2_5g = csim(u2, t, G_5);
y0_4g = csim('step', t, G_4);   y1_4g = csim(u1, t, G_4);   y2_4g = csim(u2, t, G_4);
y0_3g = csim('step', t, G_3);   y1_3g = csim(u1, t, G_3);   y2_3g = csim(u2, t, G_3);

leg = ["completo (6a)", "5a ordem", "4a ordem", "3a ordem"];

scf(0); clf;
plot(t, y0, "k-", t, y0_5g, "r--", t, y0_4g, "b-.", t, y0_3g, "g:");
xtitle("Resposta ao degrau (u = 1)", "t [s]", "y(t)"); legend(leg); xgrid;

scf(1); clf;
plot(t, y1, "k-", t, y1_5g, "r--", t, y1_4g, "b-.", t, y1_3g, "g:");
xtitle("Resposta a rampa (u = t)", "t [s]", "y(t)"); legend(leg); xgrid;

scf(2); clf;
plot(t, y2, "k-", t, y2_5g, "r--", t, y2_4g, "b-.", t, y2_3g, "g:");
xtitle("Resposta a parabola (u = t^2)", "t [s]", "y(t)"); legend(leg); xgrid;

// erro maximo de cada reducao em relacao ao modelo completo
mprintf("\nerro maximo |y - yr|\n");
mprintf("  degrau   : 5a %.3e | 4a %.3e | 3a %.3e\n", ...
        max(abs(y0-y0_5g)), max(abs(y0-y0_4g)), max(abs(y0-y0_3g)));
mprintf("  rampa    : 5a %.3e | 4a %.3e | 3a %.3e\n", ...
        max(abs(y1-y1_5g)), max(abs(y1-y1_4g)), max(abs(y1-y1_3g)));
mprintf("  parabola : 5a %.3e | 4a %.3e | 3a %.3e\n", ...
        max(abs(y2-y2_5g)), max(abs(y2-y2_4g)), max(abs(y2-y2_3g)));
