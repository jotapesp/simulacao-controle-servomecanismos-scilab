// =====================================================================
// Exercício 4. Zona de desempenho garantido (Omega) no plano s
// =====================================================================
b40 = 40*%pi/180;
sg  = 3;
Lp  = 10;

scf();
xv = [-sg, -Lp, -Lp, -sg];
yv = [ sg*tan(b40), Lp*tan(b40), -Lp*tan(b40), -sg*tan(b40)];
xfpoly(xv, yv);
e = gce();  e.background = color(220,235,250);

plot([0 -Lp], [0  Lp*tan(b40)], 'k--');   // raio beta = +40 deg
plot([0 -Lp], [0 -Lp*tan(b40)], 'k--');   // raio beta = -40 deg
plot([-sg -sg], [-Lp*tan(b40) Lp*tan(b40)], 'k--');   // reta sigma = 3

legend(['Omega'; 'beta = 40 deg'; 'beta = -40 deg'; 'sigma = 3']);
xgrid();
xtitle('Plano s - zona de desempenho garantido', 'Re', 'Im');
