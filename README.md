# traj4c
MATLAB implementation of a four times continuously differentiable trajectory using sin² jerk and optional constant-acceleration and constant-velocity phases for minimal travel time.

## Parameters
* t - current simulation time in seconds
* t_0 - trajectory start time in seconds; for t <= t_0, all outputs are 0
* s_0 - starting position in meters
* s_end - end position in meters
* v_lim - max. allowable velocity in meters/second; expects a positiv value
* a_lim -  max. allowable acceleration in meters/second²; expects a positiv value
* j_lim -  max. allowable jerk in meters/second³; expects a positiv value

## Return Values
* s - current position in meters
* v - current velocity in meters/second
* a - current acceleration in meters/second²
* j - current jerk in meters/second³
* active - <code>True</code> if currently between start and end position, <code>False</code> otherwise
* done - <code>True</code> if end position is reached, <code>False</code> otherwise
