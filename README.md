# tVDiagonalizeTimeEvaluationQuench

#Examples
- - - -
* julia tV_main_quench_q0R1PH1.jl --states-file StatesOut.dat --save-states --V0=0.0 --V=0.765 --time-min=0.0 --time-max=1.0 --time-step=0.1 --site-max=1 --ee 4  16 8

Saves the time evolved state (t=0.0 to 1.0 with Δt=0.1) to "StatesOut.dat"

* julia tV_main_quench_q0R1PH1.jl --states-file StatesOut.dat --load-states --V0=0.0 --V=0.765 --time-min=0.3 --time-max=0.7 --time-step=0.1 --site-max=1 --out output.dat --ee 4  16 8

Load the time evolved state (t=0.3 to 0.7 with Δt=0.1) from "StatesOut.dat"
- - - -

* julia tV_main_quench_q0R1PH1.jl --states-file StatesOut.dat --save-states --V0=0.0 --V=0.765 --time-min=0.0 --time-max=1.0 --time-step=0.1 --site-max=1 --ee 4  16 8

Saves the time evolved state (t=0.0 to 1.0 with Δt=0.1) to "StatesOut.dat"

* julia tV_main_quench_q0R1PH1.jl --states-file StatesOut.dat --load-states --V0=0.0 --V=0.765 --time-min=0.3 --time-max=0.7 --time-step=0.1 --site-max=1 --out output.dat --ee 4  16 8

Load the time evolved state (t=0.3 to 0.7 with Δt=0.1) from "StatesOut.dat"
- - - -

* julia tV_main_quench_q0R1PH1.jl  --save-states --V0=0.0 --V=0.765 --time-min=0.0 --time-max=1.0 --time-step=0.1 --site-max=1 --ee 4  16 8
File name for the time evolved state is generated using the input parameters.

* julia tV_main_quench_q0R1PH1.jl  --load-states --ftime-min=0.0 --ftime-max=1.0  --V0=0.0 --V=0.765 --time-min=0.3 --time-max=0.7 --time-step=0.1 --site-max=1 --out output.dat --ee 4  16 8
Uses --ftime-min=0.0 and --ftime-max=1.0  to generate the input file nameu
- - - -

* julia tV_main_quench_q0R1PH1.jl  --V0=0.0 --V=0.765 --time-min=0.3 --time-max=0.7 --time-step=0.1 --site-max=1 --out output.dat --ee 4  16 8
The time evolved state is not saved

