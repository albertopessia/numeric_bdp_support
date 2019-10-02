# Warning: long computations ahead!
# Run this script in parallel on a high-performance computing server
using SimpleBirthDeathProcess;
using DelimitedFiles;
using Statistics;

function run_simulations(N, t, S, parameter_set)
  tot_sim = 5 * N
  τ = t / S

  for (id, i, η) = parameter_set
    state = zeros(Int64, S + 1, tot_sim)

    n = 1
    while n <= tot_sim
      y = rand_discrete(i, S, τ, η)

      if y.state[2] > 0
        copyto!(state, 1 + (n - 1) * (S + 1), y.state, 1, S + 1)
        n += 1
      end
    end

    file_name = string("simulation_", lpad(id, 2, '0'), ".csv")
    fp = open(file_name, "w")
    write(fp, string("# i = ", i, ", t = ", t, ", λ = ", η[1], ", μ = ", η[2],
                     ", θ = ", η[1] - η[2], "\n"))
    writedlm(fp, state, ',')
    close(fp)
  end

  nothing
end

function find_mle(N, t, parameter_set)
  λ = zeros(Float64, 4, N)
  μ = zeros(Float64, 4, N)

  for u = [1; 3; 5]
    for (id, i, η) = parameter_set
      fill!(λ, zero(Float64))
      fill!(μ, zero(Float64))

      file_name = string("simulation_", lpad(id, 2, '0'), ".csv")
      state = readdlm(file_name, ',', Int64, '\n', comments=true)
      x = zeros(Int64, size(state, 1), u)

      for n = 1:N
        lower_bound = 1 + (n - 1) * u
        upper_bound = lower_bound + u - 1
        copyto!(x, state[:, lower_bound:upper_bound])

        for (iter, idx) = enumerate((1:8:9, 1:4:9, 1:2:9, 1:1:9))
          k = length(idx)
          S = k - 1
          τ = t / S
          z = x[idx, :]
          extinct = findall(z[2, :] .== 0)

          for q = extinct
            o = zeros(Int64, k)

            while o[2] == 0
              copyto!(o, rand_discrete(i, S, τ, η).state)
            end

            copyto!(z, 1 + (q - 1) * k, o, 1, k)
          end

          y = observation_discrete_time_equal(τ, z)
          ρ = mle(y)

          λ[iter, n] = ρ[1]
          μ[iter, n] = ρ[2]
        end
      end

      for (iter, S) = enumerate([1; 2; 4; 8])
        str = string("mle_", lpad(id, 2, '0'), "_i_", lpad(i, 4, '0'),
                     "_u_", u, "_S_", S, ".csv")
        fp = open(str, "w")
        write(fp, "lambda,mu\n")
        writedlm(fp, hcat(λ[iter, :], μ[iter, :]), ',')
        close(fp)
      end
    end
  end

  nothing
end

function evaluate_estimates(parameter_set)
  fp = open("simulation_results.csv", "w")
  write(fp, string("u,i,S,lambda,bias_lambda,rmse_lambda,mu,bias_mu,rmse_mu,",
                   "theta,bias_theta,rmse_theta\n"))

  for u = [1; 3; 5]
    for (id, i, η) = parameter_set
      θ = η[1] - η[2]

      for (iter, S) = enumerate([1; 2; 4; 8])
        file_name = string("mle_", lpad(id, 2, '0'), "_i_", lpad(i, 4, '0'),
                           "_u_", u, "_S_", S, ".csv")
        params = readdlm(file_name, ',', Float64, '\n', header=true)

        N = size(params[1], 1)
        λ = params[1][:, 1]
        μ = params[1][:, 2]

        bias = λ .- η[1]
        bias_lambda = mean(bias)
        mse_lambda = mean(bias.^2)

        bias = μ .- η[2]
        bias_mu = mean(bias)
        mse_mu = mean(bias.^2)

        bias = (λ .- μ) .- θ
        bias_theta = mean(bias)
        mse_theta = mean(bias.^2)

        str = string(u, ",", i, ",", S, ",", η[1], ",", bias_lambda, ",",
                     sqrt(mse_lambda), ",", η[2], ",", bias_mu, ",",
                     sqrt(mse_mu), ",", θ, ",", bias_theta, ",",
                     sqrt(mse_theta), "\n")
        write(fp, str)
      end
    end
  end

  close(fp)

  nothing
end

N = 1_000_000;
t = 10;
S = 8;

brg_i_0010_1 = (
# E[N(t)] = 2 * i, SD[N(t)] = 1.25 * i
(1, 10, [0.30541797643422590; 0.236103258378231370]),
# E[N(t)] = 2 * i, SD[N(t)] = 1.5 * i
(2, 10, [0.42455264809296650; 0.355237930036971960]),
# E[N(t)] = 2 * i, SD[N(t)] = 2.0 * i
(3, 10, [0.72780453958794260; 0.658489821531948000]));

brg_i_0010_2 = (
# E[N(t)] = 4 * i, SD[N(t)] = 1.25 * i
(4, 10, [0.15956825719140408; 0.020938821079415016]),
# E[N(t)] = 4 * i, SD[N(t)] = 1.5 * i
(5, 10, [0.19927981441098427; 0.060650378298995215]),
# E[N(t)] = 4 * i, SD[N(t)] = 2.0 * i
(6, 10, [0.30036377824264300; 0.161734342130653900]));

brg_i_0100_1 = (
# E[N(t)] = 2 * i, SD[N(t)] = 1.25 * i
(7, 100, [2.7422635330902840; 2.6729488150342893]),
# E[N(t)] = 2 * i, SD[N(t)] = 1.5 * i
(8, 100, [3.9336102496776895; 3.864295531621695]),
# E[N(t)] = 2 * i, SD[N(t)] = 2.0 * i
(9, 100, [6.9661291646274510; 6.896814446571456]));

brg_i_0100_2 = (
# E[N(t)] = 4 * i, SD[N(t)] = 1.25 * i
(10, 100, [0.97185010941009; 0.8332206732981009]),
# E[N(t)] = 4 * i, SD[N(t)] = 1.5 * i
(11, 100, [1.368965681605892; 1.2303362454939029]),
# E[N(t)] = 4 * i, SD[N(t)] = 2.0 * i
(12, 100, [2.379805319922479; 2.24117588381049]));

brg_i_1000_1 = (
# E[N(t)] = 2 * i, SD[N(t)] = 1.25 * i
(13, 1000, [27.110719099650860; 27.041404381594866]),
# E[N(t)] = 2 * i, SD[N(t)] = 1.5 * i
(14, 1000, [39.024186265524920; 38.954871547468926]),
# E[N(t)] = 2 * i, SD[N(t)] = 2.0 * i
(15, 1000, [69.349375415022520; 69.280060696966530]));

brg_i_1000_2 = (
# E[N(t)] = 4 * i, SD[N(t)] = 1.25 * i
(16, 1000, [9.0946686315969500; 8.9560391954849600]),
# E[N(t)] = 4 * i, SD[N(t)] = 1.5 * i
(17, 1000, [13.065824353554970; 12.927194917442980]),
# E[N(t)] = 4 * i, SD[N(t)] = 2.0 * i
(18, 1000, [23.174220736720837; 23.035591300608850]));

drg_i_0010_1 = (
# E[N(t)] = i / 2, SD[N(t)] = i / 4
(19, 10, [0.0519860385419959; 0.12130075659799043]),
# E[N(t)] = i / 2, SD[N(t)] = i / 2
(20, 10, [0.3119162312519754; 0.38123094930796990]),
# E[N(t)] = i / 2, SD[N(t)] = i
(21, 10, [1.3516370020918933; 1.42095172014788780]));

drg_i_0010_2 = (
# E[N(t)] = i / 4, SD[N(t)] = i / 4
(22, 10, [0.1617343421306539; 0.30036377824264300]),
# E[N(t)] = i / 4, SD[N(t)] = i / 2
(23, 10, [0.8548815226905992; 0.99351095880258830]),
# E[N(t)] = i / 4, SD[N(t)] = i
(24, 10, [3.6274702449303806; 3.76609968104236970]));

drg_i_0100_1 = (
# E[N(t)] = i / 2, SD[N(t)] = i / 4
(25, 100, [0.8317766166719344; 0.9010913347279289]),
# E[N(t)] = i / 2, SD[N(t)] = i / 2
(26, 100, [3.4310785437717293; 3.5003932618277240]),
# E[N(t)] = i / 2, SD[N(t)] = i
(27, 100, [13.828286252170908; 13.897600970226904]));

drg_i_0100_2 = (
# E[N(t)] = i / 4, SD[N(t)] = i / 4
(28, 100, [2.2411758838104900; 2.3798053199224790]),
# E[N(t)] = i / 4, SD[N(t)] = i / 2
(29, 100, [9.1726476894099440; 9.3112771255219330]),
# E[N(t)] = i / 4, SD[N(t)] = i
(30, 100, [36.898534911807760; 37.037164347919740]));

drg_i_1000_1 = (
# E[N(t)] = i / 2, SD[N(t)] = i / 4
(31, 1000, [8.6296823979713190; 8.6989971160273140]),
# E[N(t)] = i / 2, SD[N(t)] = i / 2
(32, 1000, [34.622701668969270; 34.692016387025260]),
# E[N(t)] = i / 2, SD[N(t)] = i
(33, 1000, [138.59477875296108; 138.66409347101705]));

drg_i_1000_2 = (
# E[N(t)] = i / 4, SD[N(t)] = i / 4
(34, 1000, [23.035591300608850; 23.1742207367208370]),
# E[N(t)] = i / 4, SD[N(t)] = i / 2
(35, 1000, [92.350309356603380; 92.488938792715370]),
# E[N(t)] = i / 4, SD[N(t)] = i
(36, 1000, [369.60918158058150; 369.74781101669350]));
