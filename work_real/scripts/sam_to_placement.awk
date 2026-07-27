BEGIN { OFS = "\t" }

and($2, 4) { next }

{
  strand = and($2, 16) ? "-" : "+"
  ctg = $3
  pos = $4 - 1
  n = 0
  span = 0
  cig = $6

  for (i = 1; i <= length(cig); i++) {
    c = substr(cig, i, 1)
    if (c >= "0" && c <= "9") {
      n = n * 10 + c
    } else {
      if (c ~ /[MDN=X]/) span += n
      n = 0
    }
  }

  hap = (ctg ~ /_MATERNAL$/) ? "mat" : ((ctg ~ /_PATERNAL$/) ? "pat" : "NA")
  if (include_mapq) {
    print $1, hap, ctg, pos, pos + span, strand, $5
  } else {
    print $1, hap, ctg, pos, pos + span, strand
  }
}
