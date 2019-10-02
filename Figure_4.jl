using SimpleBirthDeathProcess;
using DelimitedFiles;
using Plots;

N = 200;
F = Float64;
I = 200;
J = 100;
T = 1;
L = F(10);
M = F(10);

# values obtained with Maple 2018.2
correct_values = readdlm("log_trans_prob_values_bivariate.csv", ',', BigFloat);
correct_values = reshape(correct_values, N, :);

λset = [L * F(n) / N for n = 1:N];
μset = [M * F(n) / N for n = 1:N];
logp = map((λ, μ) -> trans_prob(I, J, T, [λ; μ]),
           repeat(reshape(λset, 1, :), N, 1),
           repeat(μset, 1, N));

relative_error = abs.(1 .- logp ./ correct_values);

colgr = ColorGradient([colorant"#ffffff", colorant"#525252"], [0.0, 1.0]);
pl = plot(λset, μset, log10.(relative_error), seriestype=:contour,
          levels=[-20; -16; -14; -13; -12], fill=true, contour_labels=false,
          legend=true, linecolour=:black, aspect_ratio=:equal, fillcolour=colgr,
          xlim=(0, 10), ylim=(0, 10), xlab="\\lambda", ylab="\\mu");

savefig(pl, "Figure_4.svg")
