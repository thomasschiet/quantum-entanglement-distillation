using SQLite
using Gadfly
using DataFrames

db = SQLite.DB("data.sqlite")

c_id = 1
c_fidelty = 2
c_p_succ = 3
c_delta_min = 4
c_delta_max = 5
c_state = 6
c_p = 7
c_eprFidelity = 8
function getColumn(table, columnId)
  map((F) -> F.value, collect(table.data[columnId]))
end

q = "SELECT * FROM `RainsProb` WHERE
  `state` = 'Double Roland' AND
  `p` = '0.4'"
table = SQLite.query(db, q)

plot(Coord.cartesian(xmin = 0, xmax = 1, ymin = 0, ymax = 1),
  Geom.point,
  Guide.XLabel("p"),
  Guide.YLabel("F"),
  Stat.xticks(ticks=[0:0.1:1]),
  Stat.yticks(ticks=[0:0.1:1]),
  x = getColumn(table, p),
  y = getColumn(table, c_fidelty))
