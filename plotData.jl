include("BEPQuantum.jl")

using BEPQuantum
using SQLite
using Gadfly
using DataFrames
using Compose

# settings, change these for different plots
ps = 0.1:0.1:0.9
states = ["Double Ronald", "Double Werner", "Double Ronald 2"]

for p in ps
  for state in states
    # some states have input fidelity ≤ 0.5, so they can be discarded
    if state == "Double Werner" && p < 1/3
      continue
    end
    if state == "Double Ronald 2" && p < 0.5
      continue
    end
    # Connect to database
    db = SQLite.DB("data.sqlite")

    # calculate the performance of the protocols
    if state == "Double Werner"
      stateFn = wernerState
    elseif state == "Double Ronald"
      stateFn = ronaldState
    elseif state == "Double Ronald 2"
      stateFn = ronald2State
    end
    F, p_succ = deutsch(stateFn(p))
    deutschLayer = layer(y = [F], x = [p_succ], Geom.point, Theme(default_color=colorant"cyan"))
    F, p_succ = bennett(stateFn(p))
    bennettLayer = layer(y = [F], x = [p_succ], Geom.point, Theme(default_color=colorant"orange"))

    # get data from DB
    q = string("SELECT * FROM `RainsProb` WHERE
      (`state` = '", state, "') AND
      `p` = '", p, "' AND `delta_min` = '0'")
      println(q)
    table = SQLite.query(db, q)
    # get rid of the nullables introduced by SQLite
    for c in 1:ncol(table)
      table[c] = table[c].values
    end
    PPTdf = copy(table)

    q = string("SELECT * FROM `RainsProb` WHERE
      (`state` = '1ext ", state, "') AND
      `p` = '", p, "' AND `delta_min` = '0'")
      println(q)
    table = SQLite.query(db, q)
    # get rid of the nullables introduced by SQLite
    for c in 1:ncol(table)
      table[c] = table[c].values
    end
    kextdf = copy(table)

    q = string("SELECT * FROM `RainsProb` WHERE
      (`state` = '1ext ", state, " Only Sym') AND
      `p` = '", p, "' AND `delta_min` = '0'")
      println(q)
    table = SQLite.query(db, q)
    # get rid of the nullables introduced by SQLite
    for c in 1:ncol(table)
      table[c] = table[c].values
    end
    onlySymdf = copy(table)


    function crossShape(x, y)
      lineA =  line([(1 + x,  1 + y), (-1 + x, -1 + y)])
      lineB =  line([(1 + x, -1 + y), (-1 + x,  1 + y)])
      cross = compose(context(), lineA, lineB, stroke("black"))
      return
    end

    PPTLayer = layer(PPTdf,
      x=:p_succ,
      y=:fidelity,
      Geom.line,
      Theme(default_color=colorant"blue"))

    layers = typeof(PPTLayer)[]
    names = AbstractString[]
    colors = AbstractString[]
    if nrow(PPTdf) > 0
      push!(layers, PPTLayer)
      push!(names, "PPT")
      push!(colors, "blue")
    end

    if nrow(kextdf) > 0
      kextLayer = layer(kextdf,
        x=:p_succ,
        y=:fidelity,
        Geom.point,
        Theme(default_color=colorant"red"))
      push!(layers, kextLayer)
      push!(names, "1-ext")
      push!(colors, "red")
    end

    if nrow(onlySymdf) > 0
      onlySymLayer = layer(onlySymdf,
      x=:p_succ,
      y=:fidelity,
      Geom.point,
      Theme(default_color=colorant"green"))
      push!(layers, onlySymLayer)
      push!(names, "1-ext only sym")
      push!(colors, "green")
    end

    unshift!(layers, deutschLayer)
    unshift!(layers, bennettLayer)
    unshift!(names, "Deutsch protocol")
    unshift!(names, "Bennett protocol")
    unshift!(colors, "cyan")
    unshift!(colors, "orange")

    if state == "Double Ronald 2"
      F, p_succ = epl(p)
      eplLayer = layer(y = [F], x = [p_succ], Geom.point, Theme(default_color=colorant"purple"))
      unshift!(layers, eplLayer)
      unshift!(names, "EPL protocol")
      unshift!(colors, "purple")
    end

    length(layers)
    length(names)
    names
    length(colors)
    pl = plot(
      Theme(panel_fill=colorant"white"),
      layers...,
      Stat.xticks(ticks=[0:0.1:1]),
      Stat.yticks(ticks=[0.5:0.1:1]),
      Guide.xlabel("Probability of success"),
      Guide.ylabel("Fidelity"),
      Guide.title(string(state, " state (p = ", p, ")")),
      Guide.manual_color_key("Legend",
        names,
        colors))

    draw(PNG(string("plots/",state,p,".png"), 600px, 300px), pl)
    draw(PDF(string("pdfplots/",state,p,".pdf"), 4inch, 3inch), pl)
  end
end
