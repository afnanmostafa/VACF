### Files:

* nVACF.m = MATLAB function for plotting VACF (normalized)

* vacf_plotter.m = MATLAB solo script for plotting VACF (normalized)

P.S. please be noted LAMMPS does not do a time series analysis for **__compute vacf__**, rather it plots the vacf according to **"Each atom’s contribution to the VACF is its current velocity vector dotted into its initial velocity vector at the time the compute was specified."** So, V(t).V(0) is repeated for several steps of time (t=0,1,2,......N). No lag analysis is done.

#### lag analysis:

AVG(SUM((V(0).V(0), V(1).V(0), V(2).V(1), V(3).V(2), ......, V(N).V(N-1)))) for lag = 1


AVG(SUM((V(0).V(0), V(1).V(0), V(2).V(0), V(3).V(1), V(4).V(2), ......, V(N).V(N-2)))) for lag = 2


AVG(SUM((V(0).V(0), V(1).V(0), V(2).V(0), V(3).V(0), V(4).V(0), V(5).V(0), V(6).V(1), V(7).V(2), ......, V(N).V(N-5)))) for lag = 5

This can continue up to **N+1** (considering 0 lag) lags where N = # of time samples. Then AVG(ALL LAGS) = take average of all lags been used ==> final vacf

codes in this directory do this cumbersome boring memory-intensive calculation.

LAMMPS does something like this:

((V(0).V(0), V(1).V(0), V(2).V(0), V(3).V(0), ......, V(N).V(0)) ==> final vacf for 0th, 1st, 2nd, ...., Nth step 

ALSO, all the above calculations must consider atomic contribution ==> _i.e._, V(2).V(0) would have to do this for M atoms in the system. 