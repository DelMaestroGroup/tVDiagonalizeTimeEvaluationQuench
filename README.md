# tVDiagonalizeTimeEvaluationQuench

##Example

* `julia tV_main_quench_q0R1PH1.jl  --save-states --V0=0.0 --V=0.765 --time-min=0.0 --time-max=1.0 --time-step=0.1 --site-max=1 --ee 4  16 8;

                                                     (and)
                                                     
julia tV_main_quench_q0R1PH1.jl  --load-states --ftime-min=0.0 --ftime-max=1.0  --V0=0.0 --V=0.765 --time-min=0.3 --time-max=0.7 --time-step=0.1 --site-max=1 --out output.dat --ee 4  16 8;


                                                     (or)


julia tV_main_quench_q0R1PH1.jl --states-file StatesOut.dat --save-states --V0=0.0 --V=0.765 --time-min=0.0 --time-max=1.0 --time-step=0.1 --site-max=1 --ee 4  16 8;

                                                     (and)
                                                     
julia tV_main_quench_q0R1PH1.jl --states-file StatesOut.dat --load-states --V0=0.0 --V=0.765 --time-min=0.3 --time-max=0.7 --time-step=0.1 --site-max=1 --out output.dat --ee 4  16 8;


                                                     (or)


julia tV_main_quench_q0R1PH1.jl  --V0=0.0 --V=0.765 --time-min=0.3 --time-max=0.7 --time-step=0.1 --site-max=1 --out output.dat --ee 4  16 8;


- - - -
