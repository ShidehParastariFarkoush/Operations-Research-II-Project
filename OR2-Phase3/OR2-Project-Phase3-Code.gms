
$title Port Optimization Model

Sets
    i   Ships           /1*8/
    j   berths           /1*4/
    j_limited(j)     /1, 2, 4/;

alias(i,k);

set shipPairs(i,k);
shipPairs(i,k) = ord(i) < ord(k);

Scalar
    M  "Big M constant" /80/;

Parameters
    ca(j)      "Dock capacity"   /1 57, 2 40, 3 52, 4 30/
    c(i)       "Containers on each ship"  /1 6, 2 18, 3 20, 4 5, 5 11, 6 19, 7 20, 8 13/
    max_l(j)   "Max length at each dock"    /1 118, 2 236, 3 120, 4 138/
    l(i)       "Length of each ship"      /1 114, 2 126, 3 115, 4 102, 5 102, 6 131, 7 119, 8 138/
    p(i)       "Priority of each ship"   /1 0.963797, 2 0.773138, 3 0.441551, 4 0.981939, 5 0.190926, 6 0.348681, 7 0.970487, 8 0.445607/
    epsilon(j) "Lag time at each dock"   /1 6.922094, 2 9.727825, 3 9.564813, 4 5.748308/
    anch(j)    "number of anches on each berth" /1 2, 2 2, 3 3, 4 2/
    t(i,j)     "Unloading time matrix"
    A(i)       "Arrival time of each ship"  /1 4.76553, 2 1.00459, 3 3.653, 4 6.05418, 5 1.86919, 6 2.11799, 7 4.62295, 8 6.60118/
    C_ec(j)    "Cost for extra capacity at dock j"  /1 5.87154, 2 5.880447, 3 8.864247, 4 6.64032/
    E(j)       "Max extra allowable capacity at dock j" /1 29, 2 15, 3 49, 4 17/;



Table t(i,j)
       1             2             3             4
    1  9.306945068   8.320459113   9.041338719   8.655545623
    2  8.499993353   9.905633818   9.993113985   8.089112765
    3  9.720322075   9.206381222   8.763211972   8.567236436
    4  9.349929694   8.913662302   9.371722971   9.32369264
    5  8.265956289   9.535675628   9.964826498   9.938776321
    6  9.226653641   8.088521266   8.008110288   8.267945054
    7  9.882004543   8.605721124   8.732291203   9.796392489
    8  8.62872761    9.097964368   8.872061915   8.129988352;


Binary Variables
    x(i,j)     "Assignment of ship i to berth j"
    y(i,k)   "Ordering of ships i and q at berth j (i ? q)"
    delta1(i)  "if W(i) >= 10"
    delta2(i)  "if W(i) >= 40";

Positive Variables
    W(i)       "Waiting time for ship i"
    d1(i)      "= W(i) or 0 (depends om delta1(i))"
    d2(i)      "= W(i) or 0 (depends om delta2(i))";

Integer Variables
    u(j)      "Extra capacity used at berth j";

Variable z;

Equation obj;
obj..
    sum((i,j), 5 * t(i,j) * x(i,j))
  + sum(j, C_ec(j) * u(j))
  + sum(i, 25 * W(i) - 15 * d1(i) + 150 * delta1(i) - 5 * d2(i) + 200 * delta2(i))
  =e= z;


Equations
    one_berth_per_ship
    exact_three_ships_berth3
    max_ships_limited(j)
    length_limit(i,j)
    berth_capacity(j)
    berth_extra_capacity_limit(j)
    priority_con_1(i,k)
    priority_con_2(i,k)
    Sequencing_1(i,k,j)
    Sequencing_2(i,k,j)
    waiting_upper_con(i)
    delta1_lower_con(i)
    delta2_lower_con(i)
    d1_conditional_con_1(i)
    d2_conditional_con_1(i)
    d1_conditional_con_2(i)
    d2_conditional_con_2(i);



* Each ship assigned to exactly one berth
one_berth_per_ship(i).. sum(j, x(i,j)) =e= 1;

* Number of ships assigned to each berth
* For j = 3
exact_three_ships_berth3..
    sum(i, x(i,'3')) =e= 3;

max_ships_limited(j)$(j_limited(j))..
    sum(i, x(i,j)) =l= anch(j);

* Ship length must not exceed berth capacity
length_limit(i,j)..
    l(i) * x(i,j) =l= max_l(j);

* Total containers at each berth within capacity + extra
berth_capacity(j)..
    sum(i, c(i) * x(i,j)) =l= ca(j) + u(j);

* Extra capacity must be within allowed extension
berth_extra_capacity_limit(j)..
    u(j) =l= E(j);

* Ships pririty
priority_con_1(i,k)$(shipPairs(i,k))..
    p(i) - p(k) =l= y(i,k);

priority_con_2(i,k)$(shipPairs(i,k))..
    p(k) - p(i) =l= 1 - y(i,k);


Sequencing_1(i,k,j)$(ord(i) < ord(k))..
    W(k) + A(k) =g= W(i) + A(i) + t(i,j) + epsilon(j) - M * (3 - y(i,k) - x(i,j) - x(k,j));

Sequencing_2(i,k,j)$(ord(i) < ord(k))..
    W(i) + A(i) =g= W(k) + A(k) + t(k,j) + epsilon(j) - M * (y(i,k) + 2 - x(i,j) - x(k,j));

waiting_upper_con(i)..
    W(i) =l= 80 * (1 - p(i));

delta1_lower_con(i)..
    W(i) - 10 =l= M * delta1(i);

delta2_lower_con(i)..
    W(i) - 40 =l= M * delta2(i);

d1_conditional_con_1(i)..
    d1(i) =l= W(i);

d2_conditional_con_1(i)..
    d2(i) =l= W(i);

d1_conditional_con_2(i)..
    d1(i) =l= M * delta1(i);

d2_conditional_con_2(i)..
    d2(i) =l= M * delta2(i);


Model port_model /all/;
Solve port_model using mip minimizing z;
Display z.l;
Display W.l;

Parameter assignment(i,j) "Final assignment of ships to berths";

assignment(i,j)$(x.l(i,j) > 0.5) = x.l(i,j);

display assignment;






















