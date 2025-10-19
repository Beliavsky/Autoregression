# OpenAI Codex Fortran experiment

On Windows, I used OpenAI Codex CLI to generate the code here, using the instructions in prompt.txt. The `Approval Mode` was `Full Access`.
The code was required to run with gfortran and ifx with debugging options and to pass `fortitude check`.
The code is commented so that documentation can be created using [FORD](https://github.com/Fortran-FOSS-Programmers/ford). Codex took about 15 minutes to create a working program. Here is sample truncated output.

```
AR Monte Carlo experiment configuration:
  Number of observations: 1000
  Number of simulations: 10
  Maximum AR order fitted: 6
  Maximum ACF lag compared: 20
  True AR coefficients:  0.100  0.100 -0.100
  Innovation standard deviation:  1.000
Normal RNG diagnostics:
  Mean:        0.000444
  Variance:    1.001254
  Kurtosis:    2.991856
Empirical versus theoretical ACF (first simulation):
  Lag  0 : empirical =   1.000000, theoretical =   1.000000
  Lag  1 : empirical =   0.118298, theoretical =   0.100000
  Lag  2 : empirical =   0.127890, theoretical =   0.100000
  Lag  3 : empirical =  -0.011970, theoretical =  -0.080000
  Lag  4 : empirical =   0.024423, theoretical =  -0.008000
  Lag  5 : empirical =  -0.038231, theoretical =  -0.018800
  Lag  6 : empirical =   0.005852, theoretical =   0.005320
  Lag  7 : empirical =  -0.018721, theoretical =  -0.000548
  Lag  8 : empirical =  -0.007431, theoretical =   0.002357
  Lag  9 : empirical =   0.022611, theoretical =  -0.000351
  Lag 10 : empirical =  -0.010081, theoretical =   0.000255
  Maximum absolute difference across lags:   0.068030
Simulation 1
  AR order: 0
    AR(0) model (white noise)
    Innovation variance:     1.011054
    AICC:                  10.993043
    BIC:                   10.993043
  AR order: 1
    True coefficients:        0.100000
    Estimated coefficients:   0.118298
    Innovation variance:     0.996904
    AICC:                  -1.096326
    BIC:                    3.807421
  AR order: 2
    True coefficients:        0.100000   0.100000
    Estimated coefficients:   0.104634   0.115512
    Innovation variance:     0.983603
    AICC:                 -12.521209
    BIC:                   -2.717735
  AR order: 3
    True coefficients:        0.100000   0.100000  -0.100000
    Estimated coefficients:   0.109266   0.119709  -0.040106
    Innovation variance:     0.982021
    AICC:                 -12.118906
    BIC:                    2.580263
  AR order: 4
    True coefficients:        0.100000   0.100000  -0.100000   0.000000
    Estimated coefficients:   0.109892   0.117840  -0.041812   0.015614
    Innovation variance:     0.981781
    AICC:                 -10.346636
    BIC:                    9.244184
  AR order: 5
    True coefficients:        0.100000   0.100000  -0.100000   0.000000   0.000000
    Estimated coefficients:   0.110471   0.116289  -0.037442   0.019689  -0.037078
    Innovation variance:     0.980431
    AICC:                  -9.702207
    BIC:                   14.776207
  AR order: 6
    True coefficients:        0.100000   0.100000  -0.100000   0.000000   0.000000   0.000000
    Estimated coefficients:   0.110802   0.116113  -0.037108   0.018651  -0.038064   0.008926
    Innovation variance:     0.980353
    AICC:                  -7.757650
    BIC:                   21.604290
  Best order by AICC: 2
  Best order by BIC : 2

<snip>

Selection summary over simulations:
  Order 0: AICC chosen  0 times, BIC chosen  0 times
  Order 1: AICC chosen  0 times, BIC chosen  3 times
  Order 2: AICC chosen  1 times, BIC chosen  3 times
  Order 3: AICC chosen  8 times, BIC chosen  4 times
  Order 4: AICC chosen  0 times, BIC chosen  0 times
  Order 5: AICC chosen  0 times, BIC chosen  0 times
  Order 6: AICC chosen  1 times, BIC chosen  0 times
Total CPU time (seconds):     0.0156
```

