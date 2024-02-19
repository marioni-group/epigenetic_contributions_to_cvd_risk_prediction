####################################### CVD Event definitions #######################################

# For this study, the primary composite CVD outcome was any event included 
# in the national ASSIGN risk score definition of CVD,13 
# including any International Classification of Diseases, 
# 10th Revision codes I20–25, G45, I60–69, and death from CVD (I00–I99), 
# as well, and OPCS-4 (Office of Population Censuses and Surveys:
# Classification of Interventions and Procedures version 4) procedure codes
# L29.5, L31.1, K40–46, K49, and K75 (procedures comprising carotid endarterectomy,
# carotid angioplasty, coronary artery bypass graft, and percutaneous transluminal coronary angioplasty). 
# CVD death included deaths coded from underlying causes I00 to I99; all others were classified as non-CVD deaths. 
# Other clinical outcomes (fatal or nonfatal) included coronary heart disease (I00–I25), 
# myocardial infarction (I21, I22), ischemic stroke (I63, I64, G45), any malignancy (C00–C97), 
# and hospitalization for heart failure (I50, I42.0, I42.6, I42.7, I42.9, I11.0).

# Defining CVD
CVD = c("I20", "I21", "I22", "I23", "I24", "I25", "G45", 
        "I60", "I61", "I62", "I63", "I64", "I65", "I66", "I67", "I68", "I69",
        "I00", "I01", "I02", "I03", "I04", "I05", "I06", "I07", "I08", "I09", "I10",
        "I11", "I12", "I13", "I14", "I15", "I16", "I17", "I18", "I19", 
        "I26", "I27", "I28", "I29", "I30", "I31", "I32", "I33", "I34", "I35", "I36", "I37", "I38", "I39",
        "I40", "I41", "I42", "I43", "I44", "I45", "I46", "I47", "I48", "I49",
        "I50", "I51", "I52", "I53", "I54", "I55", "I56", "I57", "I58", "I59",
        "I70", "I71", "I72", "I73", "I74", "I75", "I76", "I77", "I78", "I79",
        "I80", "I81", "I82", "I83", "I84", "I85", "I86", "I87", "I88", "I89",
        "I90", "I91", "I92", "I93", "I94", "I95", "I96", "I97", "I98", "I99",
        "L29.5", "L31.1", 
        "K40", "K41", "K42", "K43", "K44", "K45", "K46", "K49", "K75", 
        "I42.0", "I42.6", "I42.7", "I42.9", 
        "I11.0")

CVD_death = c("I00", "I01", "I02", "I03", "I04", "I05", "I06", "I07", "I08", "I09",
              "I10", "I11", "I12", "I13", "I14", "I15", "I16", "I17", "I18", "I19",
              "I20", "I21", "I22", "I23", "I24", "I25", "I26", "I27", "I28", "I29",
              "I30", "I31", "I32", "I33", "I34", "I35", "I36", "I37", "I38", "I39",
              "I40", "I41", "I42", "I43", "I44", "I45", "I46", "I47", "I48", "I49",
              "I50", "I51", "I52", "I53", "I54", "I55", "I56", "I57", "I58", "I59",
              "I60", "I61", "I62", "I63", "I64", "I65", "I66", "I67", "I68", "I69",
              "I70", "I71", "I72", "I73", "I74", "I75", "I76", "I77", "I78", "I79",
              "I80", "I81", "I82", "I83", "I84", "I85", "I86", "I87", "I88", "I89",
              "I90", "I91", "I92", "I93", "I94", "I95", "I96", "I97", "I98", "I99")

CVD_without_death = c("G45", "L29.5", "L31.1", "K40", "K41", "K42", "K43", "K44", "K45",
                      "K46", "K49", "K75", "I42.0", "I42.6", "I42.7", "I42.9", "I11.0")


CVD_by_def = c("I20", "I21", "I22", "I23", "I24", "I25", "G45", 
               "I60", "I61", "I62", "I63", "I64", "I65", "I66", "I67", "I68", "I69",
               "L29.5", "L31.1", 
               "K40", "K41", "K42", "K43", "K44", "K45", "K46", "K49", "K75",
               "I21", "I22",
               "I50", "I42.0", "I42.6", "I42.7", "I42.9", "I42", "I11.0", "I11",
               "I00", "I01", "I02", "I03", "I04", "I05", "I06", "I07", "I08", "I09",
               "I10", "I11", "I12", "I13", "I14", "I15", "I16", "I17", "I18", "I19",
               "I20", "I21", "I22", "I23", "I24", "I25")


CVD_non_fatal = c("I20", "I21", "I22", "I23", "I24", "I25", "G45", 
                  "I60", "I61", "I62", "I63", "I64", "I65", "I66", "I67", "I68", "I69",
                  "L29.5", "L31.1", 
                  "K40", "K41", "K42", "K43", "K44", "K45", "K46", "K49", "K75")


CHD_stroke_MI = c("I00", "I01", "I02", "I03", "I04", "I05", "I06", "I07", "I08", "I09",
                  "I10", "I11", "I12", "I13", "I14", "I15", "I16", "I17", "I18", "I19",
                  "I20", "I21", "I22", "I23", "I24", "I25", 
                  "I63", "I64", "G45" )

CHD_stroke_MI_heart_failure = c("I00", "I01", "I02", "I03", "I04", "I05", "I06", "I07", "I08", "I09",
                                "I10", "I11", "I12", "I13", "I14", "I15", "I16", "I17", "I18", "I19",
                                "I20", "I21", "I22", "I23", "I24", "I25", 
                                "I63", "I64", "G45",
                                "I50", "I42.0", "I42.6", "I42.7", "I42.9", "I11.0", "I42")

CHD = c("I00", "I01", "I02", "I03", "I04", "I05", "I06", "I07", "I08", "I09",
        "I10", "I11", "I12", "I13", "I14", "I15", "I16", "I17", "I18", "I19",
        "I20", "I21", "I22", "I23", "I24", "I25")

MI = c("I21", "I22")
stroke = c("I63", "I64", "G45")
heart_failure = c("I50", "I42.0", "I42.6", "I42.7", "I42.9", "I11.0")