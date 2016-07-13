using SQLite

db = SQLite.DB("data.sqlite")


using Convex
include("quantum.jl")
include("RainsProbWithPT.jl")

ps = collect(0.4:0.1:1)
δs = collect(0.05:0.05:0.95)
# ϵ = 0.05
ϵ = false

entries = length(ps)*length(δs)
stateName = "Rerun Double Werner"
n_copies = 2

times = []

firstRun = true
maxRetries = 2
retries = 0
foundNaN = false

function generateState(parameter, n_copies)
  return copies(wernerState(parameter),  n_copies)
end

function saveData(F, p_succ, δ_max, δ_min, stateName, p, eprF)
  q = string("INSERT INTO `RainsProb` (`fidelity`, `p_succ`, `delta_max`, `delta_min`, `state`, `p`, `eprFidelity`)
    VALUES (" , F , ", ", p_succ ,", ", δ_max ,", ", δ_min ,", '", stateName, "', ", p ,", ", eprF, ")")
  SQLite.query(db, q)
end

function hasBeenComputed(δ_max, δ_min, stateName, p)
  q = string("SELECT COUNT(*) FROM `RainsProb` WHERE `delta_max` = ", δ_max ," AND `delta_min` = ", δ_min ," AND p = ", p ," AND `state` = '", stateName, "'")
  result = SQLite.query(db, q)
  count = result.data[1][1].value

  return count > 0
end

function hms(timeleft)
  timeleft = round(Int, timeleft)
  hours = floor(Int, timeleft/3600)
  timeleft = timeleft - hours * 3600

  minutes = floor(Int, timeleft/60)
  timeleft = timeleft - minutes * 60

  return string(hours, "h ", minutes, "m ", timeleft, "s")
end

# sometimes the solved gives an error for unknown reasons
# therefore, retry if an error occurs
while firstRun == true || (foundNaN == true && retries < maxRetries)
  entriesDone = 0
  firstRun = false
  foundNaN = false
  retries = retries + 1

  for p in ps
    for δ in δs
      # start timer
      starttime = time()

      if ϵ ≢ false
        δ_max = δ + ϵ
        δ_min = δ - ϵ
      else
        δ_max = δ
        δ_min = 0
      end

      # check if it is already computed
      computed = hasBeenComputed(δ_max, δ_min, stateName, p)

      # only compute data if no entry is available
      if !computed
        x = generateState(p, n_copies)
      	state = sortAB(x, 2, n_copies);
        println("computing (δ_max, δ_min, p, state) = (", δ_max, ", ", δ_min, ", ", p, ", ", stateName, ")")

        (F, p_succ) = RainsProb(state, 2^n_copies, 2^n_copies, 2, δ_min, δ_max, !firstRun)
        eprF = eprFidelity(generateState(p, 1));

        println("found (F, p_succ) = " , (F, p_succ))

        if isnan(F)
          foundNaN = true
        else
          saveData(F, p_succ, δ_max, δ_min, stateName, p, real(eprF))

          stoptime = time()
          times = [times; stoptime - starttime]
          timeleft = (entries - entriesDone) * mean(times)
          entriesDone = entriesDone + 1
          println("Time to go: ", hms(timeleft), ". computation time: ", hms(stoptime-starttime), " | ", entriesDone, "/" , entries , "\n")
        end
      else
        println("already computed (δ_max, δ_min, p, state) = (", δ_max, ", ", δ_min, ", ", p, ", ", stateName, "). \n")
        entriesDone = entriesDone + 1
      end
    end
  end
end
