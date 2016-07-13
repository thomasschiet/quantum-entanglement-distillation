# Distillation project
In this project `generateData.jl` generates data and stores this in the SQLite database, specifically in the file `data.sqlite`. This data can be read using any SQLite browser (I use [sqlitebrowser](http://sqlitebrowser.org/), because it can automatically graph my queries).

## PPT relaxations
The following files are important for the optimization using PPT relaxations.

### RainsProb.jl
The `RainsProb` function solves the optimization program as described in the notes, but also includes an upper bound `δ_max` on the success probability `p_succ`.

### generateData.jl
To input a specific state, the function `generateState` has to be adjusted. It will input two copies of the state returned by `generateState` into the `RainsProb` program.

After it has successfully optimized, an entry is created in the `RainsProb` table with the optimal value and the parameters used.

### formatData.js
The file `formatData.js` stores the data from the SQLite database into `.dat` files such that they can be plotted using PGFPlots. It can be run using [node.js](http://www.nodejs.org/), simply run the command `node formatData.js`.



## ϵ-net
By running the file `epsilonNet.jl`, an epsilon net will be created and the program is run for just a couple of random points. It only outputs the identity matrix at the moment so something is going wrong. The optimization is described in `lib/constructEpsnet.jl`, together with the methods used to construct the epsilon net.
