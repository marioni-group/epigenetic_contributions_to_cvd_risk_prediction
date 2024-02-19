csv = read.csv('/Users/shirin/Documents/Edinburgh/Manuscripts/Troponin/New_journal/Concordace/concordance_cTnI_corrected.csv')
full_annots = readRDS("/Volumes/marioni-lab/Ola/Lab/EpiScores/Annotations/EpiAnnots_prepped.RDS")
joined = merge(csv, full_annots, by.x = "protein", by.y = "ID")
joined$diff = joined$c_troponin - joined$c_null_t

final = data.frame(
  joined$protein,
  joined$Gene,
  joined$c_null_t,
  joined$c_troponin,
  joined$diff,
  joined$Group,
  joined$String_annots
)

write.csv(final, '/Users/shirin/Documents/Edinburgh/Manuscripts/Troponin/New_journal/Concordace/final.csv')
