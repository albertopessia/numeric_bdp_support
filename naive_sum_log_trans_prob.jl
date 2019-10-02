using KahanSummation;
using SpecialFunctions;

# support functions aware of floating-point arithmetic
function log1mexp(x)
  if x > -log(2)
    log(-expm1(x))
  else
    log1p(-exp(x))
  end
end

function logexpm1(x)
  if x >= 37
    x
  elseif x >= 19
    x - exp(-x)
  else
    log(expm1(x))
  end
end

function log1pexp(x)
  if x <= -37
    exp(x)
  elseif x <= 18
    log1p(exp(x))
  elseif x <= 33.3
    x + exp(-x)
  else
    x
  end
end

# these functions are the log-addends in the transition probability summation
# depending on parameter values, they are designed to reduce floating-point
# rounding error as much as possible
function log_addend_equal_rates_v1(h, i, j, θ, ω)
  log(i) + lgamma(i + j - h) - lgamma(1 + h) - lgamma(i + 1 - h) -
  lgamma(j + 1 - h) + (i + j - 2 * h) * ω + h * log((1 - θ) / (1 + θ))
end

function log_addend_equal_rates_v2(h, i, j, θ, ω)
  log(i) + lgamma(i + j - h) - lgamma(1 + h) - lgamma(i + 1 - h) -
  lgamma(j + 1 - h) + (i + j - 2 * h) * ω + h * log((θ - 1) / (θ + 1))
end

function log_addend_birth_rate_greater_v1(h, i, j, θ, ω)
  log(i) + lgamma(i + j - h) - lgamma(1 + h) - lgamma(i + 1 - h) -
  lgamma(j + 1 - h) + i * ω +
  (i + j - 2 * h) * (log1mexp(θ) - log1mexp(ω + θ)) +
  h * log(-expm1(θ - ω) / expm1(θ + ω))
end

function log_addend_birth_rate_greater_v2(h, i, j, θ, ω)
  log(i) + lgamma(i + j - h) - lgamma(1 + h) - lgamma(i + 1 - h) -
  lgamma(j + 1 - h) + i * ω +
  (i + j - 2 * h) * (log1mexp(θ) - log1mexp(ω + θ)) +
  h * log(expm1(θ - ω) / expm1(θ + ω))
end

function log_addend_death_rate_greater_v1(h, i, j, θ, ω)
  log(i) + lgamma(i + j - h) - lgamma(1 + h) - lgamma(i + 1 - h) -
  lgamma(j + 1 - h) + j * ω +
  (i + j - 2 * h) * (log1mexp(θ) - log1mexp(ω + θ)) +
  h * log(-expm1(θ - ω) / expm1(θ + ω))
end

function log_addend_death_rate_greater_v2(h, i, j, θ, ω)
  log(i) + lgamma(i + j - h) - lgamma(1 + h) - lgamma(i + 1 - h) -
  lgamma(j + 1 - h) + j * ω +
  (i + j - 2 * h) * (log1mexp(θ) - log1mexp(ω + θ)) +
  h * log(expm1(θ - ω) / expm1(θ + ω))
end

# naive summation algorithm with compensated summation
function log_trans_prob(
  i::F,
  j::F,
  t::F,
  λ::F,
  μ::F
)::F where {
  F <: AbstractFloat
}
  # anything lower than ϵ (subnormal numbers) are treated as zero
  ϵ = floatmin(F)

  # a = max(i, j)
  # b = min(i, j)
  a, b = (j <= i) ? (i, j) : (j, i)

  if (t < ϵ) || ((λ < ϵ) && (μ < ϵ))
    # time is zero or both parameters are zero
    if i == j
      # probability is 1 => log(1) = 0
      0
    else
      # probability is 0 => by convention we set log(0) = -Inf
      -Inf
    end
  elseif λ ≈ μ
    # rates are equal for our floating-point machine
    ξ = 1 / λ
    θ = λ * t
    ω = log(θ / (1 + θ))

    if j > 0
      if t ≈ ξ
        # in this case we have a simple closed form solution

        # log(2) is automatically converted to a Float64 therefore we use
        # log(F(2)) for type stability
        lgamma(i + j) - lgamma(i) - lgamma(j + 1) - (i + j) * log(F(2))
      elseif t < ξ
        # addends are all positive and we can safely evaluate each value
        x = [log_addend_equal_rates_v1(h, i, j, θ, ω) for h = 0:b]

        # we have logarithms but we must sum exp(x)
        # in this case each and every exp(x[u]) is very close to zero (their sum
        # must be in [0, 1] by definition of probability)

        # use Kahan summation algorithm and other numerical tricks to reduce
        # rounding error
        sort!(x, rev=true)

        # log(exp(x[1]) + exp(x[2]) + ... + exp(x[n])) =
        # = x[1] + log(1 + exp(x[2] - x[1]) + ... + exp(x[n] - x[1]))
        x[1] + log1p(sum_kbn(exp.(x[2:end] .- x[1])))
      else
        # addends are alternating in sign and we have a huge problem
        # compute the positive and negative values separately
        x = [log_addend_equal_rates_v2(h, i, j, θ, ω) for h = 0:2:b]
        y = [log_addend_equal_rates_v2(h, i, j, θ, ω) for h = 1:2:b]

        # we sort them and try to subtract similar magnitude numbers
        sort!(x, rev=true)
        sort!(y, rev=true)

        # we cannot use exp(x) and exp(y) directly because they might have a
        # magnitude that is considered Inf by our machine
        # we normalize them by the absolute maximum
        M = max(x[1], y[1])

        # when b is odd, length(x) and length(y) are the same
        # when b is even, length(x) = length(y) + 1
        z = if (b % 2) != 0
          exp.(x .- M) .- exp.(y .- M)
        else
          [exp.(x[1:(end - 1)] .- M) .- exp.(y .- M); exp(x[end] - M)]
        end

        s = sum_kbn(z)

        # the sum must be positive because otherwise the probability would
        # be negative. If `s` is negative it means that numerical error is too
        # big
        if s > 0
          M + log(s)
        else
          -Inf
        end
      end
    else
      # j = 0
      if t ≈ ξ
        - i * log(F(2))
      else
        i * log(ω)
      end
    end
  elseif μ < ϵ
    # pure birth process
    if j >= i
      θ = λ * t

      lgamma(j) + (j - i) * logexpm1(θ) -
      (lgamma(i) + lgamma(j - i + 1) + j * θ)
    else
      # final population size cannot be lower than the initial population size
      -Inf
    end
  elseif λ < ϵ
    # pure death process
    if j <= i
      θ = μ * t

      lgamma(i + 1) + (i - j) * logexpm1(θ) -
      (lgamma(j + 1) + lgamma(i - j + 1) + i * θ)
    else
      # final population size cannot be higher than the initial population size
      -Inf
    end
  else
    ξ = log(λ / μ) / (λ - μ)

    if j > 0
      if t ≈ ξ
        # in this case we have a simple closed form solution
        lgamma(i + j) + i * log(μ) + j * log(λ) -
        (lgamma(i) + lgamma(j + 1) + (i + j) * log(λ + μ))
      elseif λ > μ
        θ = (μ - λ) * t
        ω = log(μ / λ)

        if t < ξ
          # addends are all positive and we can safely evaluate each value
          x = [log_addend_birth_rate_greater_v1(h, i, j, θ, ω) for h = 0:b]

          # we have logarithms but we must sum exp(x)
          # in this case each and every exp(x[u]) is very close to zero (their
          # sum must be in [0, 1] by definition of probability)

          # use Kahan summation algorithm and other numerical tricks to reduce
          # rounding error
          sort!(x, rev=true)

          # log(exp(x[1]) + exp(x[2]) + ... + exp(x[n])) =
          # = x[1] + log(1 + exp(x[2] - x[1]) + ... + exp(x[n] - x[1]))
          x[1] + log1p(sum_kbn(exp.(x[2:end] .- x[1])))
        else
          # addends are alternating in sign and we have a huge problem
          # compute the positive and negative values separately
          x = [log_addend_birth_rate_greater_v2(h, i, j, θ, ω) for h = 0:2:b]
          y = [log_addend_birth_rate_greater_v2(h, i, j, θ, ω) for h = 1:2:b]

          # we sort them and try to subtract similar magnitude numbers
          sort!(x, rev=true)
          sort!(y, rev=true)

          # we cannot use exp(x) and exp(y) directly because they might have a
          # magnitude that is considered Inf by our machine
          # we normalize them by the absolute maximum
          M = max(x[1], y[1])

          # when b is odd, length(x) and length(y) are the same
          # when b is even, length(x) = length(y) + 1
          z = if (b % 2) != 0
            exp.(x .- M) .- exp.(y .- M)
          else
            [exp.(x[1:(end - 1)] .- M) .- exp.(y .- M); exp(x[end] - M)]
          end

          s = sum_kbn(z)

          # the sum must be positive because otherwise the probability would
          # be negative. If `s` is negative it means that numerical error is too
          # big
          if s > 0
            M + log(s)
          else
            -Inf
          end
        end
      else
        θ = (λ - μ) * t
        ω = log(λ / μ)

        if t < ξ
          # addends are all positive and we can safely evaluate each value
          x = [log_addend_death_rate_greater_v1(h, i, j, θ, ω) for h = 0:b]

          # we have logarithms but we must sum exp(x)
          # in this case each and every exp(x[u]) is very close to zero (their
          # sum must be in [0, 1] by definition of probability)

          # use Kahan summation algorithm and other numerical tricks to reduce
          # rounding error
          sort!(x, rev=true)

          # log(exp(x[1]) + exp(x[2]) + ... + exp(x[n])) =
          # = x[1] + log(1 + exp(x[2] - x[1]) + ... + exp(x[n] - x[1]))
          x[1] + log1p(sum_kbn(exp.(x[2:end] .- x[1])))
        else
          # addends are alternating in sign and we have a huge problem
          # compute the positive and negative values separately
          x = [log_addend_death_rate_greater_v2(h, i, j, θ, ω) for h = 0:2:b]
          y = [log_addend_death_rate_greater_v2(h, i, j, θ, ω) for h = 1:2:b]

          # we sort them and try to subtract similar magnitude numbers
          sort!(x, rev=true)
          sort!(y, rev=true)

          # we cannot use exp(x) and exp(y) directly because they might have a
          # magnitude that is considered Inf by our machine
          # we normalize them by the absolute maximum
          M = max(x[1], y[1])

          # when b is odd, length(x) and length(y) are the same
          # when b is even, length(x) = length(y) + 1
          z = if (b % 2) != 0
            exp.(x .- M) .- exp.(y .- M)
          else
            [exp.(x[1:(end - 1)] .- M) .- exp.(y .- M); exp(x[end] - M)]
          end

          s = sum_kbn(z)

          # the sum must be positive because otherwise the probability would
          # be negative. If `s` is negative it means that numerical error is too
          # big
          if s > 0
            M + log(s)
          else
            -Inf
          end
        end
      end
    else
      # j = 0
      if t ≈ ξ
        i * log(μ / (λ + μ))
      else
        θ = (λ - μ) * t
        ω = log(λ / μ)
        i * log(expm1(θ) / expm1(ω + θ))
      end
    end
  end
end
