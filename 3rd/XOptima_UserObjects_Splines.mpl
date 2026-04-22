#==============================================================================
#==============================================================================
#!
#! \page SPLINES Interface with Splines C++ libhrary
#! \brief This module provides ....
#!
#! bla bla
#!
#!
#!
#!
#!
#!
export register_DataSpline::static := proc($)
  local methods, methods1;

  #   ____        _ _            ____       _
  #  / ___| _ __ | (_)_ __   ___/ ___|  ___| |_
  #  \___ \| '_ \| | | '_ \ / _ \___ \ / _ \ __|
  #   ___) | |_) | | | | | |  __/___) |  __/ |_
  #  |____/| .__/|_|_|_| |_|\___|____/ \___|\__|
  #        |_|
  methods := [
    table([
      "name"        = "eval", # metodo da mappare
      "alias"       = ["evaluate"],
      "nargs"       = 1, # numero di argomenti (variabile)
      "npars"       = 1, # numero di parametri (da mettere in coda)
      "derivatives" = 3  # derivate--> eval_D, eval_DD, eval_DDD
    ]),

    table([
      "name"        = "eval_D",
      "alias"       = ["evaluate_D"],
      "nargs"       = 1,
      "npars"       = 1,
      "derivatives" = 0
    ]),
    table([
      "name"        = "eval_DD",
      "alias"       = ["evaluate_DD"],
      "nargs"       = 1,
      "npars"       = 1,
      "derivatives" = 0
    ]),
    table([
      "name"        = "eval_DDD",
      "alias"       = ["evaluate_DDD"],
      "nargs"       = 1,
      "npars"       = 1,
      "derivatives" = 0
    ]),
    #table([ "name" = "eval_DDDD",  "alias" = ["evaluate_DDDD"] , "nargs" = 1, "npars" = 1, "derivatives" = 0 ]),
    #table([ "name" = "eval_DDDDD", "alias" = ["evaluate_DDDDD"], "nargs" = 1, "npars" = 1, "derivatives" = 0 ]),

    table([
      "name"        = "eval2",
      "nargs"       = 1,
      "npars"       = 2,
      "derivatives" = 3
    ]),

    table([ "name" = "x_min", "alias" = ["xMin"], "nargs" = 0, "npars" = 0, "derivatives" = 0 ]),
    table([ "name" = "y_min", "alias" = ["yMin"], "nargs" = 0, "npars" = 0, "derivatives" = 0 ]),
    table([ "name" = "x_max", "alias" = ["xMax"], "nargs" = 0, "npars" = 1, "derivatives" = 0 ]),
    table([ "name" = "y_max", "alias" = ["yMax"], "nargs" = 0, "npars" = 1, "derivatives" = 0 ])
  ];

  XOptima:-register_object_class( [
    "class"          = "PINS#SplineSet",
    "methods"        = methods,
    "parameters"     = [],
    "is_mesh_object" = false,
    "header"         = "#include <PINS_Core/PINS_Core.hh>",
    "namespace"      = "SplinesLoad"
  ] );

  #   ____        _ _
  #  / ___| _ __ | (_)_ __   ___
  #  \___ \| '_ \| | | '_ \ / _ \
  #   ___) | |_) | | | | | |  __/
  #  |____/| .__/|_|_|_| |_|\___|
  #        |_|
  methods := [
    table([ "name" = "x_begin", "alias" = ["xBegin"], "nargs" = 0, "npars" = 0, "derivatives" = 0 ]),
    table([ "name" = "y_begin", "alias" = ["yBegin"], "nargs" = 0, "npars" = 0, "derivatives" = 0 ]),
    table([ "name" = "x_end",   "alias" = ["xEnd"],   "nargs" = 0, "npars" = 0, "derivatives" = 0 ]),
    table([ "name" = "y_end",   "alias" = ["yEnd"],   "nargs" = 0, "npars" = 0, "derivatives" = 0 ]),
    table([ "name" = "x_min",   "alias" = ["xMin"],   "nargs" = 0, "npars" = 0, "derivatives" = 0 ]),
    table([ "name" = "x_max",   "alias" = ["xMax"],   "nargs" = 0, "npars" = 0, "derivatives" = 0 ]),
    table([ "name" = "y_min",   "alias" = ["yMin"],   "nargs" = 0, "npars" = 0, "derivatives" = 0 ]),
    table([ "name" = "y_max",   "alias" = ["yMax"],   "nargs" = 0, "npars" = 0, "derivatives" = 0 ]),

    table([
      "name"        = "eval", # metodo da mappare
      "alias"       = ["evaluate"],
      "nargs"       = 1, # numero di argomenti (variabile)
      "npars"       = 0, # numero di parametri (da mettere in coda)
      "derivatives" = 5  # derivate--> eval_D, eval_DD, eval_DDD
    ]),
    table([
      "name"        = "eval_D", # metodo da mappare
      "alias"       = ["evaluate_D"],
      "nargs"       = 1, # numero di argomenti (variabile)
      "npars"       = 0, # numero di parametri (da mettere in coda)
      "derivatives" = 0  # derivate--> eval_D, eval_DD, eval_DDD
    ]),
    table([
      "name"        = "eval_DD", # metodo da mappare
      "alias"       = ["evaluate_DD"],
      "nargs"       = 1, # numero di argomenti (variabile)
      "npars"       = 0, # numero di parametri (da mettere in coda)
      "derivatives" = 0  # derivate--> eval_D, eval_DD, eval_DDD
    ]),
    table([
      "name"        = "eval_DDD", # metodo da mappare
      "alias"       = ["evaluate_DDD"],
      "nargs"       = 1, # numero di argomenti (variabile)
      "npars"       = 0, # numero di parametri (da mettere in coda)
      "derivatives" = 0  # derivate--> eval_D, eval_DD, eval_DDD
    ]),
    table([
      "name"        = "eval_DDDD", # metodo da mappare
      "alias"       = ["evaluate_DDDD"],
      "nargs"       = 1, # numero di argomenti (variabile)
      "npars"       = 0, # numero di parametri (da mettere in coda)
      "derivatives" = 0  # derivate--> eval_D, eval_DD, eval_DDD
    ]),
    table([
      "name"        = "eval_DDDDD", # metodo da mappare
      "alias"       = ["evaluate_DDDDD"],
      "nargs"       = 1, # numero di argomenti (variabile)
      "npars"       = 0, # numero di parametri (da mettere in coda)
      "derivatives" = 0  # derivate--> eval_D, eval_DD, eval_DDD
    ])
  ];

  XOptima:-register_object_class([
    "class"          = "PINS#Spline1D", # cubic[extrapolate,natural,parabolic,not_a_knot], akima, bessel, pchip, linear, constant, quintic[cubic,pchip,akima,bessel]
    "methods"        = methods,
    "parameters"     = [],
    "is_mesh_object" = false,
    "header"         = "#include <PINS_Core/PINS_Core.hh>",
    "namespace"      = "SplinesLoad"
  ]);

  XOptima:-register_object_class([
    "class"          = "PINS#CubicSpline",
    "methods"        = methods,
    "parameters"     = [],
    "is_mesh_object" = false,
    "header"         = "#include <PINS_Core/PINS_Core.hh>",
    "namespace"      = "SplinesLoad"
  ]);

  XOptima:-register_object_class([
    "class"          = "PINS#AkimaSpline",
    "methods"        = methods,
    "parameters"     = [],
    "is_mesh_object" = false,
    "header"         = "#include <PINS_Core/PINS_Core.hh>",
    "namespace"      = "SplinesLoad"
  ]);

  XOptima:-register_object_class([
    "class"          = "PINS#BesselSpline",
    "methods"        = methods,
    "parameters"     = [],
    "is_mesh_object" = false,
    "header"         = "#include <PINS_Core/PINS_Core.hh>",
    "namespace"      = "SplinesLoad"
  ]);

  XOptima:-register_object_class([
    "class"          = "PINS#PchipSpline",
    "methods"        = methods,
    "parameters"     = [],
    "is_mesh_object" = false,
    "header"         = "#include <PINS_Core/PINS_Core.hh>",
    "namespace"      = "SplinesLoad"
  ]);

  XOptima:-register_object_class([
    "class"          = "PINS#LinearSpline",
    "methods"        = methods,
    "parameters"     = [],
    "is_mesh_object" = false,
    "header"         = "#include <PINS_Core/PINS_Core.hh>",
    "namespace"      = "SplinesLoad"
  ]);

  XOptima:-register_object_class([
    "class"          = "PINS#ConstantSpline",
    "methods"        = methods,
    "parameters"     = [],
    "is_mesh_object" = false,
    "header"         = "#include <PINS_Core/PINS_Core.hh>",
    "namespace"      = "SplinesLoad"
  ]);

  #    ___        _       _   _
  #   / _ \ _   _(_)_ __ | |_(_) ___
  #  | | | | | | | | '_ \| __| |/ __|
  #  | |_| | |_| | | | | | |_| | (__
  #   \__\_\\__,_|_|_| |_|\__|_|\___|

  XOptima:-register_object_class([
    "class"          = "PINS#QuinticSpline",
    "methods"        = methods,
    "parameters"     = [],
    "is_mesh_object" = false,
    "header"         = "#include <PINS_Core/PINS_Core.hh>",
    "namespace"      = "SplinesLoad"
  ]);

  #   ____        _ _            ____  ____
  #  / ___| _ __ | (_)_ __   ___|___ \|  _ \
  #  \___ \| '_ \| | | '_ \ / _ \ __) | | | |
  #   ___) | |_) | | | | | |  __// __/| |_| |
  #  |____/| .__/|_|_|_| |_|\___|_____|____/
  #        |_|
  methods := [
    table([
      "name"        = "eval", # metodo da mappare
      "alias"       = ["evaluate"],
      "nargs"       = 2,      # numero di argomenti (variabile)
      "npars"       = 0,      # numero di parametri (da mettere in coda)
      "derivatives" = 2       # derivate--> eval_D, eval_DD, eval_DDD
    ]),
    table([
      "name"        = "eval_D_1",
      "alias"       = ["evaluate_D_1"],
      "nargs"       = 2,
      "npars"       = 0,
      "derivatives" = 0
    ]),
    table([
      "name"        = "eval_D_2",
      "alias"       = ["evaluate_D_2"],
      "nargs"       = 2,
      "npars"       = 0,
      "derivatives" = 0
    ]),
    table([
      "name"        = "eval_D_1_1",
      "alias"       = ["evaluate_D_1_1"],
      "nargs"       = 2,
      "npars"       = 0,
      "derivatives" = 0
    ]),
    table([
      "name"        = "eval_D_1_2",
      "alias"       = ["evaluate_D_1_2"],
      "nargs"       = 2,
      "npars"       = 0,
      "derivatives" = 0
    ]),
    table([
      "name"        = "eval_D_2_2",
      "alias"       = ["evaluate_D_2_2"],
      "nargs"       = 2,
      "npars"       = 0,
      "derivatives" = 0
    ]),
    table([ "name" = "x_min", "alias" = ["xMin"], "nargs" = 0, "npars" = 0, "derivatives" = 0 ]),
    table([ "name" = "x_max", "alias" = ["xMax"], "nargs" = 0, "npars" = 0, "derivatives" = 0 ]),
    table([ "name" = "y_min", "alias" = ["yMin"], "nargs" = 0, "npars" = 0, "derivatives" = 0 ]),
    table([ "name" = "y_max", "alias" = ["yMax"], "nargs" = 0, "npars" = 0, "derivatives" = 0 ]),
    table([ "name" = "z_min", "alias" = ["zMin"], "nargs" = 0, "npars" = 0, "derivatives" = 0 ]),
    table([ "name" = "z_max", "alias" = ["zMax"], "nargs" = 0, "npars" = 0, "derivatives" = 0 ])
  ];

  XOptima:-register_object_class([
    "class"          = "PINS#Spline2D",
    "methods"        = methods,
    "parameters"     = [],
    "is_mesh_object" = false,
    "header"         = "#include <PINS_Core/PINS_Core.hh>",
    "namespace"      = "SplinesLoad"
  ]);

  XOptima:-register_object_class([
    "class"          = "PINS#BilinearSpline",
    "methods"        = methods,
    "parameters"     = [],
    "is_mesh_object" = false,
    "header"         = "#include <PINS_Core/PINS_Core.hh>",
    "namespace"      = "SplinesLoad"
  ]);

  XOptima:-register_object_class([
    "class"          = "PINS#BiCubicSpline",
    "methods"        = methods,
    "parameters"     = [],
    "is_mesh_object" = false,
    "header"         = "#include <PINS_Core/PINS_Core.hh>",
    "namespace"      = "SplinesLoad"
  ]);

  XOptima:-register_object_class([
    "class"          = "PINS#Akima2Dspline",
    "methods"        = methods,
    "parameters"     = [],
    "is_mesh_object" = false,
    "header"         = "#include <PINS_Core/PINS_Core.hh>",
    "namespace"      = "SplinesLoad"
  ]);

  XOptima:-register_object_class([
    "class"          = "PINS#BiQuinticSplineBase",
    "methods"        = methods,
    "parameters"     = [],
    "is_mesh_object" = false,
    "header"         = "#include <PINS_Core/PINS_Core.hh>",
    "namespace"      = "SplinesLoad"
  ]);

  #   ____  _                _
  #  | __ )| | ___ _ __   __| |
  #  |  _ \| |/ _ \ '_ \ / _` |
  #  | |_) | |  __/ | | | (_| |
  #  |____/|_|\___|_| |_|\__,_|
  #  

  XOptima:-register_object_class([
    "class"          = "PINS#Spline1Dblend",
    "methods"        =  [
      table([
        "name"        = "x_begin", # metodo da mappare
        "alias"       = ["xBegin"],
        "nargs"       = 1,      # numero di argomenti (variabile)
        "npars"       = 0,      # numero di parametri (da mettere in coda)
        "derivatives" = 0       # derivate--> eval_D, eval_DD, eval_DDD
      ]),
      table([
        "name"        = "x_end",
        "alias"       = ["xEnd"],
        "nargs"       = 1,
        "npars"       = 0,
        "derivatives" = 0
      ]),
      table([
        "name"        = "y_begin",
        "alias"       = ["yBegin"],
        "nargs"       = 1,
        "npars"       = 0,
        "derivatives" = 0 
      ]),
      table([
        "name"        = "y_end",
        "alias"       = ["yEnd"],
        "nargs"       = 1,
        "npars"       = 0,
        "derivatives" = 0
      ]),
      table([
        "name"        = "eval",
        "alias"       = ["evaluate"],
        "nargs"       = 1,
        "npars"       = 1,
        "derivatives" = 2 
      ]),
      table([
        "name"        = "eval_D",
        "alias"       = [],
        "nargs"       = 1,
        "npars"       = 1,
        "derivatives" = 0
      ]),
      table([
        "name"        = "eval_DD",
        "alias"       = [],
        "nargs"       = 1,
        "npars"       = 1,
        "derivatives" = 0
      ]),
      table([
        "name"        = "eval_DDD",
        "alias"       = [],
        "nargs"       = 1,
        "npars"       = 1,
        "derivatives" = 0
      ]),
      table([
        "name"        = "eval_DDDD",
        "alias"       = [],
        "nargs"       = 1,
        "npars"       = 1,
        "derivatives" = 0
      ]),
      table([
        "name"        = "eval_DDDDD",
        "alias"       = [],
        "nargs"       = 1,
        "npars"       = 1,
        "derivatives" = 0
      ])
    ],
    "parameters"     = [],
    "is_mesh_object" = false,
    "header"         = "#include <PINS_Core/PINS_Core.hh>",
    "namespace"      = "SplinesLoad"
  ]);

  XOptima:-register_object_class([
    "class"          = "PINS#Spline2Dblend",
    "methods"        =  [
      table([
        "name"        = "eval",
        "alias"       = ["evaluate"],
        "nargs"       = 2,
        "npars"       = 1,
        "derivatives" = 2 
      ]),
      table([
        "name"        = "Dx",
        "alias"       = [],
        "nargs"       = 2,
        "npars"       = 1,
        "derivatives" = 0
      ]),
      table([
        "name"        = "Dy",
        "alias"       = [],
        "nargs"       = 2,
        "npars"       = 1,
        "derivatives" = 0
      ]),
      table([
        "name"        = "Dxx",
        "alias"       = [],
        "nargs"       = 2,
        "npars"       = 1,
        "derivatives" = 0
      ]),
      table([
        "name"        = "Dxy",
        "alias"       = [],
        "nargs"       = 2,
        "npars"       = 1,
        "derivatives" = 0
      ]),
      table([
        "name"        = "Dyy",
        "alias"       = [],
        "nargs"       = 2,
        "npars"       = 1,
        "derivatives" = 0
      ])
    ],
    "parameters"     = [],
    "is_mesh_object" = false,
    "header"         = "#include <PINS_Core/PINS_Core.hh>",
    "namespace"      = "SplinesLoad"
  ]);

end proc:
