```scilab
--> s=%s

 s = [polynomial] of s

  s

--> coef_g = [16875, 27000, 15075, 5720, 817, 48, 1];

--> deng   = poly(coef_g, 's', 'coeff');

--> g = syslin('c', 1, deng);

--> g

 g = [rational] of s

                            1                            
   ----------------------------------------------------  
   16875 +27000s +15075s^2 +5720s^3 +817s^4 +48s^5 +s^6  

--> gr

Undefined variable: gr

--> // 7.9014e-06 s^2 - 5.9259e-05 s + 2.9630e-04

--> num_r = [2.9630e-04, - 5.9259e-05, 7.9014e-06]

 num_r = [1x3 double]

   0.0002963  -0.0000593   0.0000079

--> // s^3 + 3s^2 + 7s + 5

--> den_r = [5, 7, 3, 1]

 den_r = [1x4 double]

   5.   7.   3.   1.

--> numg_r=poly(num_r, 's', 'coeff')

 numg_r = [polynomial] of s

  0.0002963 -0.0000593s +0.0000079s^2

--> deng_r=poly(den_r, 's', 'coeff')

 deng_r = [polynomial] of s

  5 +7s +3s^2 +s^3

--> gr = syslin('c', numg_r, deng_r);

--> gr

 gr = [rational] of s

   0.0002963 -0.0000593s +0.0000079s^2  
   -----------------------------------  
            5 +7s +3s^2 +s^3            

--> sd = -3.2 + %i*2.24;

--> sd

 sd = 

  -3.2 + 2.24i

--> ar=horner(g, sd)

 ar = 

  -0.0000024 - 0.0000172i

--> ar2=horner(gr, sd)

 ar2 = 

  -0.0000002 - 0.0000175i

--> ar = atan(0.0000172/0.0000024)

 ar = 

   1.4321566

--> a=0.0000172/0.0000024

 a = 

   7.1666667

--> arctan(a)

Undefined variable: arctan

--> ang_ar = atan(imag(ar), real(ar)) * 180/%pi;

--> ang_ar

 ang_ar = 

   0.

--> ang_a = atan(imag(a), real(a)) * 180/%pi;

--> ang_a

 ang_a = 

   0.

--> ar = horner(g, sd);

--> ang_ar = atan(imag(ar), real(ar)) * 180/%pi;

--> ang_ar

 ang_ar = 

  -97.927139

--> disp(180 - (ang - 90));

Undefined variable: ang

--> disp(180 - (ang_ar - 90));

   367.92714

--> disp(180 - (ang_ar - 145));

   422.92714

--> roots(deng)

 ans = [6x1 double]

  -15.000073 + 0.i       
  -14.999963 + 0.0000633i
  -14.999963 - 0.0000633i
  -1.        + 2.i       
  -1.        - 2.i       
  -1.        + 0.i       

--> ar2 = horner(gr, sd);

--> ar2

 ar2 = 

  -0.0000002 - 0.0000175i

--> ang_ar2 = atan(imag(ar2), real(ar2)) * 180/%pi;

--> ang_ar2

 ang_ar2 = 

  -90.697465

--> ang_sd=atan(imag(sd), real(sd)) * 180/%pi;

--> ang_sd

 ang_sd = 

   145.00798

--> deficit = ang_sd + ang_ar

 deficit = 

   47.080840

--> deficit = -ang_sd+ang_ar

 deficit = 

  -242.93512

--> deficit = 180 - deficit

 deficit = 

   422.93512

--> defict_360 = 422-360

 defict_360 = 

   62.

--> deficit=defict_360

 deficit = 

   62.

--> // para a planta reduzida:

--> deficit_r = -ang_sd+ang_ar2

 deficit_r = 

  -235.70545

--> deficit_r = 180 - deficit_r

 deficit_r = 

   415.70545

--> deficit_r = deficit_r - 360

 deficit_r = 

   55.705445

--> gr

 gr = [rational] of s

   0.0002963 -0.0000593s +0.0000079s^2  
   -----------------------------------  
            5 +7s +3s^2 +s^3            

--> z=-3.5

 z = 

  -3.5

--> arg(sd-z)

Undefined variable: arg

--> ar3=horner((s-z), sd)

 ar3 = 

   0.3 + 2.24i

--> angd_1=atan(imag(ar3), real(ar3)) * 180/%pi;

--> angd_1

 angd_1 = 

   82.371850

--> p = 9.6

 p = 

   9.6

--> pr = 8.52

 pr = 

   8.52

--> p = -p

 p = 

  -9.6

--> pr = -pr

 pr = 

  -8.52

--> c= syslin('c', (s-z), (s-p))

 c = [rational] of s

   3.5 +s  
   ------  
   9.6 +s  

--> cr = yslin('c', (s-z), (s-pr))

Undefined variable: yslin

--> cr = syslin('c', (s-z), (s-pr))

 cr = [rational] of s

   3.5 +s   
   -------  
   8.52 +s  

--> val = horner(c, sd)*horner(g, sd)

 val = 

   0.0000048 - 0.0000033i

--> abs(val)

 ans = 

   0.0000058

--> 1/abs(val)

 ans = 

   172538.09

--> kc = 1/abs(val)

 kc = 

   172538.09

--> val_r = horner(cr, sd)*horner(gr, sd)

 val_r = 

   0.0000059 - 0.0000036i

--> kc_r = 1/(abs(val_r))

 kc_r = 

   145715.37

--> val = horner(c, sd)*horner(s*g, sd)

 val = 

  -0.0000078 + 0.0000213i

--> val_r = horner(cr, sd)*horner(gr/s, sd)

 val_r = 

  -0.0000018 - 0.0000001i

--> val = horner(c, sd)*horner(g/s, sd)

 val = 

  -0.0000015 - 3.799D-09i

--> kc = 1/abs(val)

 kc = 

   673950.67

--> kc_r = 1/(abs(val_r))

 kc_r = 

   569178.47
```