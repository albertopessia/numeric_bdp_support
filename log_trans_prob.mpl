# Maple 2018.2 file to be called with
# read "log_trans_prob.mpl"

# we assume that b > 0!

with(LinearAlgebra):

trans_prob_unequal := proc(a, b, t, lambda, mu)
  local h, theta, psi, phi, alpha, beta, gamma:

  theta := lambda - mu:
  psi := lambda + mu:

  phi := (exp(theta * t) - 1) / (lambda * exp(theta * t) - mu):
  alpha := mu * phi:
  beta := lambda * phi:
  gamma := 1 - psi * phi:

  binomial(a + b - 1, a - 1) * alpha^a * beta^b + sum(binomial(a, h) * binomial(a + b - h - 1, a - 1) * alpha^(a - h) * beta^(b - h) * gamma^h, h = 1 .. min(a, b))
end proc:

trans_prob_equal := proc(a, b, t, lambda)
  local h, theta, alpha, gamma:

  theta := lambda * t:

  alpha := theta / (1 + theta):
  gamma := (1 - theta) / (1 + theta):

  binomial(a + b - 1, a - 1) * alpha^(a + b) + sum(binomial(a, h) * binomial(a + b - h - 1, a - 1) * alpha^(a + b - 2 * h) * gamma^h, h = 1 .. min(a, b))
end proc:

log_trans_prob := proc(a, b, t, lambda, mu)
  local s, lp:

  s := 0:
  lp := 0:

  if (lambda <> mu) then
    s := trans_prob_unequal(a, b, t, lambda, mu):
  else
    s := trans_prob_equal(a, b, t, lambda):
  end if:

  if evalf[20](s) > 0 then
    lp := ln(s):
  else
    lp := ln(normal(s)):
  end if:

  lp
end proc:

log_trans_prob_float := proc(a, b, t, lambda, mu)
  evalf[100](log_trans_prob(a, b, t, lambda, mu))
end proc:

# generate data for Figures 1 and 3
N := 1000:
A := 25:
B := 35:
T := 2:
L := 1:
M := 3:

values := Array([seq(0, i = 1 .. N)]):

for n to N do
  values[n] := log_trans_prob_float(A, B, T, L, n * M / N):
end do:

values := Transpose(Matrix(convert(values, list))):

Export("log_trans_prob_values_univariate.csv", values):

# generate data for Figure 4 (beware: it takes a long time)
N := 100:
A := 200:
B := 100:
T := 1:
L := 10:
M := 10:

values := Array([seq(0, i = 1 .. (N * N))]):

c0 := 1:
for c1 to N do
  for c2 to N do
    values[c0] := log_trans_prob_float(A, B, T, c1 * L / N, c2 * M / N):
    c0 := c0 + 1:
  end do:
end do:

values := Transpose(Matrix(convert(values, list))):

Export("log_trans_prob_values_bivariate.csv", values):
