using SimpleBirthDeathProcess;
using DelimitedFiles;
using Plots;

N = 1_000;
F = Float64;
I = 25;
J = 35;
T = 2;
L = F(1);
M = F(3);

# values obtained with Maple 2018.2
correct_values = readdlm("log_trans_prob_values_univariate.csv", ',', BigFloat);

μset = [M * F(n) / N for n = 1:N];
logp = [trans_prob(I, J, T, [L; F(μ)]) for μ = μset];

relative_error = abs.(1 .- logp ./ correct_values);

pl = plot(μset, relative_error, seriestype=:line, colour=:black,
          yscale=:log10, xlab="\\mu", ylab="Relative error",
          legend=false, size=(1920, 1080));

savefig(pl, "Figure_3.svg")
