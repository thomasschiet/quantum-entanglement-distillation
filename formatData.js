var sqlite3 = require('sqlite3');
var db = new sqlite3.Database('data.sqlite');

var fs = require('fs');

p=5
db.serialize(function() {
  var resultingfile = 'p_succ\tfidelity\tepr\n';
  db.each("SELECT * FROM RainsProb WHERE `state` = 'Double Ronald' AND p = 0." + p + " AND delta_min = 0 ORDER BY p_succ ASC", function(err, row) {
    resultingfile = resultingfile + (row.p_succ + "\t" + row.fidelity + '\t' + row.eprFidelity +'\n');
    fs.writeFile('data/ronald/p0_'+p+'.dat', resultingfile, function(err) { if (err) console.log(err);});
  })
});

db.close();
