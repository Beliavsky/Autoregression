program xar
  !!* Monte Carlo study of AR model selection for an AR(3) process.
  ! The program simulates multiple series, checks diagnostics, and compares
  ! information criteria across candidate autoregressive orders.
  use ar_toolkit, only: dp, generate_normal, simulate_ar_series, autocorrelations, &
       theoretical_acf_ar, mean_value, variance_value, kurtosis_value, yule_walker_fit, &
       compute_information_criteria
  implicit none

  integer, parameter :: nobs = 1000  !! Number of observations per series
  integer, parameter :: nsim = 10  !! Number of simulated series
  integer, parameter :: max_order = 6  !! Maximum AR order considered in fitting
  integer, parameter :: maxlag = 20  !! Maximum lag for ACF comparison
  integer, parameter :: n_rng_check = 200000  !! Sample size for RNG diagnostics
  real(dp), parameter :: phi_true(3) = [0.1_dp, 0.1_dp, -0.1_dp]  !! True AR(3) coefficients
  real(dp), parameter :: sigma_eps = 1.0_dp  !! Innovation standard deviation

  real(dp), allocatable :: series(:)  !! Simulated series storage
  real(dp), allocatable :: acf_emp(:)  !! Empirical autocorrelation values
  real(dp), allocatable :: acf_theory(:)  !! Theoretical autocorrelation values
  real(dp), allocatable :: phi_hat(:)  !! Estimated AR coefficients for a fit
  real(dp), allocatable :: phi_true_pad(:)  !! True coefficients truncated to model order
  real(dp), allocatable :: aicc_values(:)  !! AICC values for candidate orders
  real(dp), allocatable :: bic_values(:)  !! BIC values for candidate orders
  integer, allocatable :: order_counts_aicc(:)  !! Tally of AICC-selected orders
  integer, allocatable :: order_counts_bic(:)  !! Tally of BIC-selected orders
  real(dp), allocatable :: rng_vals(:)  !! Storage for RNG diagnostics

  real(dp) :: rng_mean, rng_var, rng_kurt  !! RNG diagnostic statistics
  real(dp) :: start_time, end_time  !! CPU time measurements
  real(dp) :: max_diff  !! Maximum absolute discrepancy between ACFs
  real(dp) :: aicc, bic, sigma2_hat  !! Information criteria and innovation variance
  real(dp) :: huge_value  !! Helper constant for guard values
  integer :: isim, order, lag, best_aicc_order, best_bic_order  !! Loop counters and selected orders
  integer :: info  !! Status code from AR fitting

  call cpu_time(start_time)
  huge_value = huge(1.0_dp)

  allocate(series(nobs))
  allocate(acf_emp(0:maxlag))
  allocate(acf_theory(0:maxlag))
  allocate(aicc_values(0:max_order))
  allocate(bic_values(0:max_order))
  allocate(order_counts_aicc(0:max_order))
  allocate(order_counts_bic(0:max_order))
  order_counts_aicc = 0
  order_counts_bic = 0

  print '(a)', 'AR Monte Carlo experiment configuration:'
  print '(a,i0)', '  Number of observations: ', nobs
  print '(a,i0)', '  Number of simulations: ', nsim
  print '(a,i0)', '  Maximum AR order fitted: ', max_order
  print '(a,i0)', '  Maximum ACF lag compared: ', maxlag
  print '(a,*(f6.3,1x))', '  True AR coefficients: ', phi_true
  print '(a,f6.3)', '  Innovation standard deviation: ', sigma_eps

  allocate(rng_vals(n_rng_check))
  call generate_normal(n_rng_check, rng_vals)
  rng_mean = mean_value(rng_vals)
  rng_var = variance_value(rng_vals)
  rng_kurt = kurtosis_value(rng_vals)
  print '(a)', 'Normal RNG diagnostics:'
  print '(a,f10.6)', '  Mean:      ', rng_mean
  print '(a,f10.6)', '  Variance:  ', rng_var
  print '(a,f10.6)', '  Kurtosis:  ', rng_kurt
  deallocate(rng_vals)

  call theoretical_acf_ar(phi_true, sigma_eps**2, maxlag, acf_theory)

  do isim = 1, nsim
    call simulate_ar_series(phi_true, sigma_eps, 500, series)
    call autocorrelations(series, maxlag, acf_emp)

    if (isim == 1) then
      print '(a)', 'Empirical versus theoretical ACF (first simulation):'
      do lag = 0, maxlag
        print '("  Lag ", i2, " : empirical = ", f10.6, ", theoretical = ", f10.6)', &
             & lag, acf_emp(lag), acf_theory(lag)
        if (lag == 10) exit
      end do
      max_diff = maxval(abs(acf_emp - acf_theory))
      print '(a,f10.6)', '  Maximum absolute difference across lags: ', max_diff
    end if

    print '(a,i0)', 'Simulation ', isim
    do order = 0, max_order
      call yule_walker_fit(series, order, phi_hat, sigma2_hat, info)
      if (info /= 0) then
        aicc = huge_value
        bic = huge_value
      else
        call compute_information_criteria(nobs, order, sigma2_hat, aicc, bic)
      end if
      aicc_values(order) = aicc
      bic_values(order) = bic

      print '(a,i0)', '  AR order: ', order
      if (order > 0) then
        allocate(phi_true_pad(order))
        phi_true_pad = 0.0_dp
        phi_true_pad(1:min(order, size(phi_true))) = phi_true(1:min(order, size(phi_true)))
        print '(a,*(f10.6,1x))', '    True coefficients:      ', phi_true_pad
        if (info == 0) then
          print '(a,*(f10.6,1x))', '    Estimated coefficients: ', phi_hat
        else
          print '(a)', '    Estimated coefficients: (fit failed)'
        end if
        deallocate(phi_true_pad)
      else
        print '(a)', '    AR(0) model (white noise)'
      end if
      if (info == 0) then
        print '(a,f12.6)', '    Innovation variance: ', sigma2_hat
      else
        print '(a,i0)', '    Fit status code: ', info
      end if
      print '(a,f12.6)', '    AICC:               ', aicc_values(order)
      print '(a,f12.6)', '    BIC:                ', bic_values(order)

      if (allocated(phi_hat)) deallocate(phi_hat)
    end do

    best_aicc_order = 0
    best_bic_order = 0
    do order = 0, max_order
      if (order == 0 .or. aicc_values(order) < aicc_values(best_aicc_order)) best_aicc_order = order
      if (order == 0 .or. bic_values(order) < bic_values(best_bic_order)) best_bic_order = order
    end do
    order_counts_aicc(best_aicc_order) = order_counts_aicc(best_aicc_order) + 1
    order_counts_bic(best_bic_order) = order_counts_bic(best_bic_order) + 1

    print '(a,i0)', '  Best order by AICC: ', best_aicc_order
    print '(a,i0)', '  Best order by BIC : ', best_bic_order
  end do

  print '(a)', 'Selection summary over simulations:'
  do order = 0, max_order
    print '("  Order ", i1, ": AICC chosen ", i2, " times, BIC chosen ", i2, " times")', &
         & order, order_counts_aicc(order), order_counts_bic(order)
  end do

  call cpu_time(end_time)
  print '(a,f10.4)', 'Total CPU time (seconds): ', end_time - start_time

  deallocate(series)
  deallocate(acf_emp)
  deallocate(acf_theory)
  deallocate(aicc_values)
  deallocate(bic_values)
  deallocate(order_counts_aicc)
  deallocate(order_counts_bic)
end program xar
