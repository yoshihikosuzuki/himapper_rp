BEGIN { OFS = "\t" }

{
  primary = 0
  for (i = 13; i <= NF; i++) {
    if ($i == "tp:A:P") {
      primary = 1
      break
    }
  }
  if (!primary) next

  hap = ($6 ~ /_MATERNAL$/) ? "mat" : (($6 ~ /_PATERNAL$/) ? "pat" : "NA")
  print $1, hap, $6, $8, $9, $5
}
