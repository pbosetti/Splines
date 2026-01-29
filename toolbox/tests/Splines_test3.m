%%-------------------------------------------------------------------------%%
%%                                                                         %%
%%  Copyright (C) 2020                                                     %%
%%                                                                         %%
%%         , __                 , __                                       %%
%%        /|/  \               /|/  \                                      %%
%%         | __/ _   ,_         | __/ _   ,_                               %%
%%         |   \|/  /  |  |   | |   \|/  /  |  |   |                       %%
%%         |(__/|__/   |_/ \_/|/|(__/|__/   |_/ \_/|/                      %%
%%                           /|                   /|                       %%
%%                           \|                   \|                       %%
%%                                                                         %%
%%      Enrico Bertolazzi                                                  %%
%%      Dipartimento di Ingegneria Industriale                             %%
%%      Università degli Studi di Trento                                   %%
%%      email: enrico.bertolazzi@unitn.it                                  %%
%%                                                                         %%
%%-------------------------------------------------------------------------%%

% Ottieni la dimensione dello schermo [left bottom width height]
scrn = get(0, 'ScreenSize');

% Calcola dimensioni all'80%
width  = round(scrn(3) * 0.8);
height = round(scrn(4) * 0.8);

% Centra la figura
left = round((scrn(3) - width) / 2);
bottom = round((scrn(4) - height) / 2);

% Crea la figura
figure('Position', [left bottom width height]);

X = -2:0.01:2;
Y = -2:0.01:2;
[XX,YY] = ndgrid(X,Y);
ZZ  = peaks(XX,YY);
bl  = Spline2D('bilinear',X,Y,ZZ);
bc  = Spline2D('bicubic',X,Y,ZZ);
ba  = Spline2D('bicubic_akima',X,Y,ZZ);
bb  = Spline2D('bicubic_vanleer',X,Y,ZZ);
bp  = Spline2D('bicubic_pchip',X,Y,ZZ);
bq  = Spline2D('biquintic',X,Y,ZZ);
bqa = Spline2D('biquintic_akima',X,Y,ZZ);
bqb = Spline2D('biquintic_vanleer',X,Y,ZZ);
bqp = Spline2D('biquintic_pchip',X,Y,ZZ);

surf(XX,YY,ZZ), view(145,-2), set(gca,'Fontsize',16);

X = -2:0.1:2;
Y = -2:0.1:2;
[XX,YY] = ndgrid(X,Y);
ZZ = peaks(XX,YY);

Z1 = bl.eval(XX,YY);
Z2 = bc.eval(XX,YY);
Z3 = ba.eval(XX,YY);
Z4 = bb.eval(XX,YY);
Z5 = bp.eval(XX,YY);
Z6 = bq.eval(XX,YY);
Z7 = bqa.eval(XX,YY);
Z8 = bqb.eval(XX,YY);
Z9 = bqp.eval(XX,YY);

if any(isnan(Z1(:)))
  fprintf('Found NaN on bilinear\n');
end
if any(isinf(Z2(:)))
  fprintf('Found Inf on bicubic\n');
end
if any(isnan(Z3(:)))
  fprintf('Found NaN on bicubic akima\n');
end
if any(isinf(Z4(:)))
  fprintf('Found Inf on bicubic vanleer\n');
end
if any(isnan(Z5(:)))
  fprintf('Found NaN on bicubic pchip\n');
end
if any(isinf(Z6(:)))
  fprintf('Found Inf on biquintic\n');
end
if any(isnan(Z7(:)))
  fprintf('Found NaN on biquintic akima\n');
end
if any(isinf(Z8(:)))
  fprintf('Found Inf on biquintic vanleer\n');
end
if any(isinf(Z9(:)))
  fprintf('Found Inf on biquintic pchip\n');
end

subplot(3,3,1);
surf(XX,YY,Z1,'Linestyle',':'), view(145,40), set(gca,'Fontsize',16);
title('bilinear');

subplot(3,3,2);
surf(XX,YY,Z1,'Linestyle',':'), view(145,40), set(gca,'Fontsize',16);
title('bicubic');

subplot(3,3,3);
surf(XX,YY,Z1,'Linestyle',':'), view(145,40), set(gca,'Fontsize',16);
title('bicubic akima');

subplot(3,3,4);
surf(XX,YY,Z1,'Linestyle',':'), view(145,40), set(gca,'Fontsize',16);
title('bicubic vanleer');

subplot(3,3,5);
surf(XX,YY,Z1,'Linestyle',':'), view(145,40), set(gca,'Fontsize',16);
title('bicubic pchip');

subplot(3,3,6);
surf(XX,YY,Z1,'Linestyle',':'), view(145,40), set(gca,'Fontsize',16);
title('biquintic');

subplot(3,3,7);
surf(XX,YY,Z1,'Linestyle',':'), view(145,40), set(gca,'Fontsize',16);
title('biquintic akima');

subplot(3,3,8);
surf(XX,YY,Z1,'Linestyle',':'), view(145,40), set(gca,'Fontsize',16);
title('biquintic vanleer');

subplot(3,3,9);
surf(XX,YY,Z1,'Linestyle',':'), view(145,40), set(gca,'Fontsize',16);
title('biquintic pchip');
