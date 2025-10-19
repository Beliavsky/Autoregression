module ar_toolkit
  !!* Utilities for simulating and analysing autoregressive models.
  ! The module gathers reusable procedures for Monte Carlo experiments and
  ! model diagnostics.
  implicit none
  private
  integer, parameter, public :: dp = selected_real_kind(15, 307)  !! Double precision kind parameter
  public :: mean_value, variance_value, kurtosis_value
  public :: generate_normal, simulate_ar_series
  public :: autocovariances, autocorrelations, theoretical_acf_ar
  public :: yule_walker_fit, compute_information_criteria
contains
  pure function mean_value(x) result(val)
    !!* Compute the arithmetic mean of the supplied sample.
    ! The mean is returned as zero when the input is empty.
    real(dp), intent(in) :: x(:)  !! Sample data
    real(dp) :: val  !! Sample mean
    integer :: n

    n = size(x)
    if (n == 0) then
      val = 0.0_dp
    else
      val = sum(x) / real(n, dp)
    end if
  end function mean_value

  pure function variance_value(x) result(val)
    !!* Compute the population variance of a sample.
    ! The variance uses a divisor equal to the sample size.
    real(dp), intent(in) :: x(:)  !! Sample data
    real(dp) :: val  !! Population variance
    real(dp) :: mu
    integer :: n

    n = size(x)
    if (n == 0) then
      val = 0.0_dp
      return
    end if

    mu = mean_value(x)
    val = sum((x - mu)**2) / real(n, dp)
  end function variance_value

  pure function kurtosis_value(x) result(val)
    !!* Compute the population kurtosis of a sample.
    ! The kurtosis is defined as the normalised fourth central moment.
    real(dp), intent(in) :: x(:)  !! Sample data
    real(dp) :: val  !! Population kurtosis
    real(dp) :: mu, denom
    integer :: n

    n = size(x)
    if (n == 0) then
      val = 0.0_dp
      return
    end if

    mu = mean_value(x)
    denom = sum((x - mu)**2) / real(n, dp)
    if (denom <= 0.0_dp) then
      val = 0.0_dp
    else
      val = (sum((x - mu)**4) / real(n, dp)) / (denom * denom)
    end if
  end function kurtosis_value

  subroutine generate_normal(n, sample)
    !!* Generate independent standard normal variates using the Box-Muller transform.
    ! The procedure fills the supplied array with pseudo-random draws.
    integer, intent(in) :: n  !! Number of variates required
    real(dp), intent(out) :: sample(n)  !! Output variates
    integer :: i
    real(dp) :: u(2), r, theta

    i = 1
    do while (i <= n)
      call random_number(u)
      if (u(1) <= 0.0_dp) cycle
      r = sqrt(-2.0_dp * log(u(1)))
      theta = 2.0_dp * acos(-1.0_dp) * u(2)
      sample(i) = r * cos(theta)
      if (i + 1 <= n) then
        sample(i + 1) = r * sin(theta)
      end if
      i = i + 2
    end do
  end subroutine generate_normal

  subroutine simulate_ar_series(phi, sigma, burn_in, series)
    !!* Simulate an autoregressive process of arbitrary order.
    ! A burn-in period is discarded to reduce dependence on starting values.
    real(dp), intent(in) :: phi(:)  !! AR coefficients ordered by increasing lag
    real(dp), intent(in) :: sigma  !! Innovation standard deviation
    integer, intent(in), optional :: burn_in  !! Length of discarded burn-in
    real(dp), intent(out) :: series(:)  !! Simulated time series
    integer :: n, p, b, total_length, t, j
    real(dp), allocatable :: innovations(:)
    real(dp), allocatable :: work(:)
    real(dp) :: value

    p = size(phi)
    n = size(series)
    if (present(burn_in)) then
      b = burn_in
    else
      b = 250
    end if
    total_length = n + b

    allocate(innovations(total_length))
    call generate_normal(total_length, innovations)
    innovations = sigma * innovations

    allocate(work(total_length))
    work = 0.0_dp

    do t = 1, total_length
      value = innovations(t)
      do j = 1, p
        if (t - j > 0) value = value + phi(j) * work(t - j)
      end do
      work(t) = value
    end do

    series = work(b + 1:total_length)

    deallocate(innovations)
    deallocate(work)
  end subroutine simulate_ar_series

  pure subroutine autocovariances(x, maxlag, acov)
    !!* Compute sample autocovariances up to a specified lag.
    ! The divisor equals the sample size to match population conventions.
    real(dp), intent(in) :: x(:)  !! Input series
    integer, intent(in) :: maxlag  !! Maximum lag required
    real(dp), intent(out) :: acov(0:maxlag)  !! Estimated autocovariances
    integer :: n, k, t
    real(dp) :: mu

    n = size(x)
    acov = 0.0_dp
    if (n == 0) return

    mu = mean_value(x)
    do k = 0, maxlag
      do t = k + 1, n
        acov(k) = acov(k) + (x(t) - mu) * (x(t - k) - mu)
      end do
      acov(k) = acov(k) / real(n, dp)
    end do
  end subroutine autocovariances

  pure subroutine autocorrelations(x, maxlag, acf)
    !!* Compute sample autocorrelations up to a specified lag.
    ! The values are normalised by the lag-zero autocovariance.
    real(dp), intent(in) :: x(:)  !! Input series
    integer, intent(in) :: maxlag  !! Maximum lag required
    real(dp), intent(out) :: acf(0:maxlag)  !! Estimated autocorrelations
    real(dp) :: acov(0:maxlag)
    integer :: k

    call autocovariances(x, maxlag, acov)
    if (acov(0) <= 0.0_dp) then
      acf = 0.0_dp
      return
    end if

    do k = 0, maxlag
      acf(k) = acov(k) / acov(0)
    end do
  end subroutine autocorrelations

  pure subroutine theoretical_acf_ar(phi, sigma2, maxlag, acf)
    !!* Compute the theoretical autocorrelation function of an AR process.
    ! The routine solves the Yule-Walker equations for the supplied coefficients.
    real(dp), intent(in) :: phi(:)  !! AR coefficients ordered by lag
    real(dp), intent(in) :: sigma2  !! Innovation variance
    integer, intent(in) :: maxlag  !! Maximum lag required
    real(dp), intent(out) :: acf(0:maxlag)  !! Theoretical autocorrelations
    integer :: p, dim, i, j, k, idx
    real(dp) :: matrix(size(phi) + 1, size(phi) + 1)
    real(dp) :: rhs(size(phi) + 1)
    real(dp) :: gamma_vec(size(phi) + 1)
    real(dp) :: ext_gamma(0:maxlag)

    p = size(phi)
    dim = p + 1
    matrix = 0.0_dp
    rhs = 0.0_dp

    do i = 1, dim
      matrix(i, i) = 1.0_dp
    end do

    do j = 1, p
      matrix(1, j + 1) = matrix(1, j + 1) - phi(j)
    end do
    rhs(1) = sigma2

    do k = 1, p
      do j = 1, p
        idx = abs(k - j) + 1
        matrix(k + 1, idx) = matrix(k + 1, idx) - phi(j)
      end do
    end do

    call solve_linear_system(matrix, rhs, gamma_vec)

    do k = 0, p
      ext_gamma(k) = gamma_vec(k + 1)
    end do
    if (maxlag > p) then
      do k = p + 1, maxlag
        ext_gamma(k) = 0.0_dp
        do j = 1, p
          ext_gamma(k) = ext_gamma(k) + phi(j) * ext_gamma(k - j)
        end do
      end do
    end if

    if (ext_gamma(0) <= 0.0_dp) then
      acf = 0.0_dp
    else
      do k = 0, maxlag
        acf(k) = ext_gamma(k) / ext_gamma(0)
      end do
    end if
  end subroutine theoretical_acf_ar

  pure subroutine solve_linear_system(a, b, x)
    !!* Solve a dense linear system using Gaussian elimination with partial pivoting.
    ! The routine modifies its inputs and is intended for small systems.
    real(dp), intent(inout) :: a(:, :)  !! Coefficient matrix
    real(dp), intent(inout) :: b(:)  !! Right-hand side vector
    real(dp), intent(out) :: x(:)  !! Solution vector
    integer :: n, i, j, k, pivot_row
    real(dp) :: factor, max_val, temp

    n = size(b)
    x = 0.0_dp

    do k = 1, n
      pivot_row = k
      max_val = abs(a(k, k))
      do i = k + 1, n
        if (abs(a(i, k)) > max_val) then
          max_val = abs(a(i, k))
          pivot_row = i
        end if
      end do

      if (max_val == 0.0_dp) cycle

      if (pivot_row /= k) then
        call swap_rows(a, b, k, pivot_row)
      end if

      do i = k + 1, n
        if (a(k, k) == 0.0_dp) cycle
        factor = a(i, k) / a(k, k)
        do j = k, n
          a(i, j) = a(i, j) - factor * a(k, j)
        end do
        b(i) = b(i) - factor * b(k)
      end do
    end do

    do i = n, 1, -1
      temp = b(i)
      do j = i + 1, n
        temp = temp - a(i, j) * x(j)
      end do
      if (a(i, i) == 0.0_dp) then
        x(i) = 0.0_dp
      else
        x(i) = temp / a(i, i)
      end if
    end do
  end subroutine solve_linear_system

  pure subroutine swap_rows(a, b, irow, jrow)
    !!* Swap two rows in the matrix and right-hand side vector.
    ! The helper is used by the Gaussian elimination routine.
    real(dp), intent(inout) :: a(:, :)  !! Matrix to pivot
    real(dp), intent(inout) :: b(:)  !! Right-hand side vector
    integer, intent(in) :: irow  !! First row index
    integer, intent(in) :: jrow  !! Second row index
    real(dp) :: temp_row(size(a, 2))
    real(dp) :: temp_val

    if (irow == jrow) return
    temp_row = a(irow, :)
    a(irow, :) = a(jrow, :)
    a(jrow, :) = temp_row

    temp_val = b(irow)
    b(irow) = b(jrow)
    b(jrow) = temp_val
  end subroutine swap_rows

  subroutine yule_walker_fit(x, order, phi, sigma2, info)
    !!* Estimate an autoregressive model using the Yule-Walker equations.
    ! Reflection coefficients are obtained via the Durbin-Levinson recursion.
    real(dp), intent(in) :: x(:)  !! Input time series
    integer, intent(in) :: order  !! Desired AR order
    real(dp), allocatable, intent(out) :: phi(:)  !! Estimated AR coefficients
    real(dp), intent(out) :: sigma2  !! Innovation variance estimate
    integer, intent(out) :: info  !! Status flag (0 indicates success)
    integer :: n, j, m
    real(dp) :: acov0, denom, kappa
    real(dp), allocatable :: gamma_hat(:)
    real(dp), allocatable :: phi_work(:)
    real(dp), allocatable :: phi_prev(:)

    n = size(x)
    info = 0
    if (order < 0) then
      info = 1
      sigma2 = 0.0_dp
      allocate(phi(0))
      return
    end if
    if (order >= n) then
      info = 2
      sigma2 = 0.0_dp
      allocate(phi(0))
      return
    end if

    allocate(phi(order))
    if (order == 0) then
      sigma2 = variance_value(x)
      if (sigma2 <= 0.0_dp) sigma2 = tiny(sigma2) + 1.0e-12_dp
      return
    end if

    allocate(gamma_hat(0:order))
    call autocovariances(x, order, gamma_hat)
    acov0 = gamma_hat(0)
    if (acov0 <= 0.0_dp) then
      info = 3
      sigma2 = 0.0_dp
      phi = 0.0_dp
      deallocate(gamma_hat)
      return
    end if

    allocate(phi_work(order))
    allocate(phi_prev(order))
    phi_work = 0.0_dp
    sigma2 = acov0

    do m = 1, order
      denom = sigma2
      if (denom <= 0.0_dp) then
        info = 4
        exit
      end if
      kappa = gamma_hat(m)
      do j = 1, m - 1
        kappa = kappa - phi_work(j) * gamma_hat(m - j)
      end do
      kappa = kappa / denom

      phi_prev(1:order) = phi_work(1:order)
      phi_work(m) = kappa
      do j = 1, m - 1
        phi_work(j) = phi_prev(j) - kappa * phi_prev(m - j)
      end do
      sigma2 = sigma2 * (1.0_dp - kappa * kappa)
    end do

    if (sigma2 <= 0.0_dp) then
      sigma2 = tiny(sigma2) + 1.0e-12_dp
    end if

    phi = phi_work

    deallocate(gamma_hat)
    deallocate(phi_work)
    deallocate(phi_prev)
  end subroutine yule_walker_fit

  pure subroutine compute_information_criteria(nobs, order, sigma2, aicc, bic)
    !!* Compute information criteria for an AR model.
    ! The formulas assume a zero-mean process estimated by Yule-Walker.
    integer, intent(in) :: nobs  !! Sample size
    integer, intent(in) :: order  !! AR order
    real(dp), intent(in) :: sigma2  !! Innovation variance estimate
    real(dp), intent(out) :: aicc  !! Corrected Akaike information criterion
    real(dp), intent(out) :: bic  !! Bayesian information criterion
    real(dp) :: k_par, n_real, aic

    n_real = real(nobs, dp)
    k_par = real(order, dp)
    if (sigma2 <= 0.0_dp) then
      aicc = huge(1.0_dp)
      bic = huge(1.0_dp)
      return
    end if

    aic = n_real * log(sigma2) + 2.0_dp * k_par
    if (nobs - order - 1 > 0) then
      aicc = aic + (2.0_dp * k_par * (k_par + 1.0_dp)) / (n_real - k_par - 1.0_dp)
    else
      aicc = huge(1.0_dp)
    end if
    bic = n_real * log(sigma2) + k_par * log(n_real)
  end subroutine compute_information_criteria

end module ar_toolkit
