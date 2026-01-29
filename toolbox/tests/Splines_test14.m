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

close all;

% Ottieni la dimensione dello schermo [left bottom width height]
scrn = get(0, 'ScreenSize');

% Calcola dimensioni all'80%
width  = round(scrn(3) * 0.8);
height = round(scrn(4) * 0.8);

% Centra la figura
left = round((scrn(3) - width) / 2);
bottom = round((scrn(4) - height) / 2);

type = {
    'bilinear',
    'bicubic',
    'bicubic_akima',
    'bicubic_bessel',
    'bicubic_pchip',
    'biquintic',
    'biquintic_akima',
    'biquintic_bessel',
    'biquintic_pchip'
    };

for k=1:9

  % Crea la figura
  figure('Position', [left bottom width height]);

  set(gca,'Fontsize',16);

  disableDefaultInteractivity(gca)

  S = Spline2D(type{k});
  S.build('OLDAeroMap.json');
  %S.build('test_surf.json');
  
  x_min = S.x_min();
  x_max = S.x_max();
  y_min = S.y_min();
  y_max = S.y_max();

  X = x_min:0.05:x_max;
  Y = y_min:0.01:y_max;
  [XX,YY] = ndgrid(X,Y);
  
  for kk=1:6

    switch (kk)
    case 1
        deriv = 'f(x, y)';
        ZZ = S.eval(XX, YY);
    case 2
        deriv = '∂f/∂x';
        ZZ = S.eval_Dx(XX, YY);
    case 3
        deriv = '∂f/∂y';
        ZZ = S.eval_Dy(XX, YY);
    case 4
        deriv = '∂²f/∂x∂y';
        ZZ = S.eval_Dxy(XX, YY);
    case 5
        deriv = '∂²f/∂x²';
        ZZ = S.eval_Dxx(XX, YY);
    case 6
        deriv = '∂²f/∂y²';
        ZZ = S.eval_Dyy(XX, YY);
    end
    
    if any(isnan(ZZ(:)))
      fprintf('Trovati NaN sulla superfice');
    end
    if any(isinf(ZZ(:)))
      fprintf('Trovati Inf sulla superfice');
    end

    subplot(2,3,kk);

    surf(XX,YY,ZZ,'Linestyle',':');

    axis tight

    %zlim([-1,70]);

    view(60,60);

    title(sprintf('%s %s',type{k},deriv));
  end
end

