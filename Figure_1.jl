include("naive_sum_log_trans_prob.jl");
using DelimitedFiles;
using Plots;

N = 1_000;
F = Float64;
I = F(25);
J = F(35);
T = F(2);
L = F(1);
M = F(3);

# values obtained with Maple 2018.2
correct_values = readdlm("log_trans_prob_values_univariate.csv", ',', BigFloat);

μset = [M * F(n) / N for n = 1:N];
logp = [log_trans_prob(I, J, T, L, F(μ)) for μ = μset];

relative_error = abs.(one(F) .- logp ./ correct_values);

# we cannot plot infinite values
relative_error[isinf.(relative_error)] .= NaN;

pl = plot(μset, relative_error, seriestype=:line, colour=:black,
          yscale=:log10, xlab="\\mu", ylab="Relative error", legend=false,
          size=(1600, 800));

savefig(pl, "Figure_1.svg")
