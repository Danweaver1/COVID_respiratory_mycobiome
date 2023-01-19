#### This code was adapted from atable resources by Armin Strobel eg.:
# https://journal.r-project.org/archive/2019/RJ-2019-001/RJ-2019-001.pdf

### Initial setup ####
#load libraries
library(atable)
library(multgee)
library(readxl)

#get data
#set directory 
setwd("")
#input clinical data from excel file
data.raw <- read_excel("", sheet = "samples")

#### create tables of full raw data #####
#cut to variables want to show
data.trim <- data.raw[c(3,6,10:18,20:34,36:53,55,61,63)]
## Reformatting 
#replace numbers/acronyms where necessary 
library(stringr)
#cut$LifeSupport <- cut$LifeSupport %>% str_replace("None","none") %>% str_replace("MV", "mechanical") %>% str_replace("CV","mechanical & noninvasive") %>% str_replace("NIV","noninvasive")

#remove any duplicated patient data 
data.unique <- unique(data.trim)

colnames(data.unique)
#Set more appropriate classes:
data.reformatted <- within(data.unique, {
  CAPANoCAPA = factor(CAPANoCAPA, levels = c(1,0), labels = c("yes","no"))
  patient = ordered(seq(1,length(patient)))
  #DaysfromICUadmintoCAPA  - continuous 
  BALGM1ODI = factor(BALGM1ODI, levels = c(1,0), labels = c("positive","negative"))  
  BAL_LFD = factor(BAL_LFD, levels = c(1,0), labels = c("positive","negative"))   
  PosBALPCR = factor(PosBALPCR, levels = c(1,0), labels = c("positive","negative"))   
  PosBALculture = factor(PosBALculture, levels = c(1,0), labels = c("positive","negative"))
  PosBronchAspculture = factor(PosBronchAspculture, levels = c(1,0), labels = c("positive","negative"))
  PosTAAspCult = factor(PosTAAspCult, levels = c(1,0), labels = c("positive","negative"))
  PosTAGM = factor(PosTAGM, levels = c(1,0), labels = c("positive","negative"))
  PosPCRfromTA = factor(PosPCRfromTA, levels = c(1,0), labels = c("positive","negative"))
  SerumGM05 = factor(SerumGM05, levels = c(1,0), labels = c("positive","negative"))
  Institution =factor(Institution, levels = c(1,8,9), labels = c("Graz", "Genoa", "Rennes"))
  Sex = factor(Sex, levels = c(0,1), labels = c("male", "female"))
  #Age  - continuous 
  Ethnicity = factor(Ethnicity, levels = c(5,6,7), labels = c("caucasian", "unknown" , "other"), order=TRUE)
  #BMInew  - continuous
  BMI30 = factor(BMI30, levels = c(1,0), labels = c("yes","no"))
  UnderlHemOnc = factor(UnderlHemOnc, levels = c(1,0), labels = c("yes","no"))
  UnderlSOT = factor(UnderlSOT, levels = c(1,0), labels = c("yes","no"))
  UnderlCardiovasc = factor(UnderlCardiovasc, levels = c(1,0), labels = c("yes","no"))
  UnderlPulm = factor(UnderlPulm, levels = c(1,0), labels = c("yes","no"))
  UnderlDM = factor(UnderlDM, levels = c(1,0), labels = c("yes","no"))
  Smoking = factor(Smoking, levels = c(1,0), labels = c("yes","no"))
  Corticos = factor(Corticos, levels = c(1,0), labels = c("yes","no"))
  Tocilizumab = factor(Tocilizumab, levels = c(1,0), labels = c("yes","no"))
  Azithro = factor(Azithro, levels = c(1,0), labels = c("yes","no"))
  MechanVnet = factor(MechanVnet, levels = c(1,0), labels = c("yes","no"))
  NIVent = factor(NIVent, levels = c(1,0), labels = c("yes","no"))
  ECMOnew = factor(ECMOnew, levels = c(1,0), labels = c("yes","no"))
  SystemicAFbeforeatICUAdmin = factor(SystemicAFbeforeatICUAdmin, levels = c(1,0), labels = c("yes","no"))
  AFinitforCAPA = factor(AFinitforCAPA, levels = c(1,0), labels = c("yes","no"))
  WhichAF = factor(WhichAF, levels = c("PosaIsa","Posa Isa","Isav","IsavCasp","LAmb","Vori"), labels = c("PosaIsa","Posa Isa","Isav","IsavCasp","LAmb","Vori"))
  Vori = factor(Vori, levels = c(1,0), labels = c("yes","no"))
  Isa = factor(Isa, levels = c(1,0), labels = c("yes","no"))
  LipAmb = factor(LipAmb, levels = c(1,0), labels = c("yes","no"))
  DAmB = factor(DAmB, levels = c(1,0), labels = c("yes","no"))
  Posa = factor(Posa, levels = c(1,0), labels = c("yes","no"))
  Echinocandin = factor(Echinocandin, levels = c(1,0), labels = c("yes","no"))
  Combination = factor(Combination, levels = c(1,0), labels = c("yes","no"))
  AFTreatmentOutcopme = factor(AFTreatmentOutcopme, levels = c("yes","no","partial"))
  PrimaryTreatment = factor(PrimaryTreatment, levels = c(1,0), labels = c("yes","no"))
  #DurationICUDays - continuous
  PataliveDay28or32 = factor(PataliveDay28or32, levels = c(1,0), labels = c("yes","no"))
  SurvivalEOF = factor(SurvivalEOF, levels = c(1,0), labels = c("yes","no"))
  #BALGMODI  - continuous
  LifeSupport = factor(LifeSupport, levels = c("MV","ECMO","NIV","CV","None"), labels = c("mechanical","ECMO","noninvasive", "mechanical & noninvasive","none"),order=TRUE)
  })


#get variable names in required format for defining columns when creating table
var_names <- toString(shQuote(colnames(data.reformatted), type = "cmd"))
cat(var_names, "\n")

#First, create a table that contains demographic and clinical characteristics for each group. The target
#variables are sex , age and baselinescore ; the variable trt acts as the grouping variable:
#change output format with:format_to - Possible values are 'Latex', 'Word', 'Raw', 'HTML', 'Console', 'markdown', 'md'.
#Word, can be further processed with e.g. flextable of flextable Gohel (2018).
the_table <- atable::atable(data.reformatted,
                            target_cols = c( "CAPANoCAPA", "DaysfromICUadmintoCAPA", "BALGM1ODI", "BAL_LFD",
                                            "PosBALPCR", "PosBALculture", "PosBronchAspculture", "PosTAAspCult",
                                            "PosTAGM", "PosPCRfromTA", "SerumGM05", "Institution", "Sex", "Age",
                                            "Ethnicity", "BMI30", "UnderlHemOnc", "UnderlSOT",
                                            "UnderlCardiovasc", "UnderlPulm", "UnderlDM", "Smoking",
                                            "Corticos", "Tocilizumab", "Azithro", "MechanVnet",
                                            "NIVent", "ECMOnew", "SystemicAFbeforeatICUAdmin", "AFinitforCAPA",
                                            "WhichAF", "Vori", "Isa", "LipAmb", "DAmB", "Posa", "Echinocandin",
                                            "Combination", "AFTreatmentOutcopme", "PrimaryTreatment", "DurationICUDays",
                                            "PataliveDay28or32", "SurvivalEOF", "LifeSupport" ), 
                            format_to = "Word")


# print in Word with packages flextable and officer
library(flextable)
MyFTable <- flextable::regulartable(data = the_table)
# left aligned first column:
#MyFTable <- flextable::align(MyFTable, align = "left", j = 1)
#save on disc. Not run here:
doc <- officer::read_docx()
doc <- flextable::body_add_flextable(doc, value = MyFTable)
print(doc, target = "basic_table_full_data.docx")



####create tables for final dataset (ie. patients whose samples passed QC) ##################
##subset clinical data based on mycobiome samples which passed QC
### Read files in ####
otu_mat <- as.matrix(read.table("out", sep =",", header=TRUE))
#Define row names from otu column "X"
row.names(otu_mat) <- otu_mat[,1]
#Remove the column from the matrix
otu_mat <- otu_mat[,-1]
#subset clinical data for samples available 
#alternatively, can subset later using subset_samples(physeq, <sample data column> =="Yes")
data.passqc <- subset(data.raw, Sample %in% colnames(otu_mat) )
#cut to variables want to show
data.trim <- data.passqc[c(3,6,10:18,20:34,36:53,55,61,63)]
## Reformatting 
#replace numbers/acronyms where necessary 
library(stringr)
#cut$LifeSupport <- cut$LifeSupport %>% str_replace("None","none") %>% str_replace("MV", "mechanical") %>% str_replace("CV","mechanical & noninvasive") %>% str_replace("NIV","noninvasive")

#remove any duplicated patient data 
data.unique <- unique(data.trim)

colnames(data.unique)
#Set more appropriate classes:
data.reformatted <- within(data.unique, {
  CAPANoCAPA = factor(CAPANoCAPA, levels = c(1,0), labels = c("yes","no"))
  patient = ordered(seq(1,length(patient)))
  #DaysfromICUadmintoCAPA  - continuous 
  BALGM1ODI = factor(BALGM1ODI, levels = c(1,0), labels = c("positive","negative"))  
  BAL_LFD = factor(BAL_LFD, levels = c(1,0), labels = c("positive","negative"))   
  PosBALPCR = factor(PosBALPCR, levels = c(1,0), labels = c("positive","negative"))   
  PosBALculture = factor(PosBALculture, levels = c(1,0), labels = c("positive","negative"))
  PosBronchAspculture = factor(PosBronchAspculture, levels = c(1,0), labels = c("positive","negative"))
  PosTAAspCult = factor(PosTAAspCult, levels = c(1,0), labels = c("positive","negative"))
  PosTAGM = factor(PosTAGM, levels = c(1,0), labels = c("positive","negative"))
  PosPCRfromTA = factor(PosPCRfromTA, levels = c(1,0), labels = c("positive","negative"))
  SerumGM05 = factor(SerumGM05, levels = c(1,0), labels = c("positive","negative"))
  Institution =factor(Institution, levels = c(1,8,9), labels = c("Graz", "Genoa", "Rennes"))
  Sex = factor(Sex, levels = c(0,1), labels = c("male", "female"))
  #Age  - continuous 
  Ethnicity = factor(Ethnicity, levels = c(5,6,7), labels = c("caucasian", "unknown" , "other"), order=TRUE)
  #BMInew  - continuous
  BMI30 = factor(BMI30, levels = c(1,0), labels = c("yes","no"))
  UnderlHemOnc = factor(UnderlHemOnc, levels = c(1,0), labels = c("yes","no"))
  UnderlSOT = factor(UnderlSOT, levels = c(1,0), labels = c("yes","no"))
  UnderlCardiovasc = factor(UnderlCardiovasc, levels = c(1,0), labels = c("yes","no"))
  UnderlPulm = factor(UnderlPulm, levels = c(1,0), labels = c("yes","no"))
  UnderlDM = factor(UnderlDM, levels = c(1,0), labels = c("yes","no"))
  Smoking = factor(Smoking, levels = c(1,0), labels = c("yes","no"))
  Corticos = factor(Corticos, levels = c(1,0), labels = c("yes","no"))
  Tocilizumab = factor(Tocilizumab, levels = c(1,0), labels = c("yes","no"))
  Azithro = factor(Azithro, levels = c(1,0), labels = c("yes","no"))
  MechanVnet = factor(MechanVnet, levels = c(1,0), labels = c("yes","no"))
  NIVent = factor(NIVent, levels = c(1,0), labels = c("yes","no"))
  ECMOnew = factor(ECMOnew, levels = c(1,0), labels = c("yes","no"))
  SystemicAFbeforeatICUAdmin = factor(SystemicAFbeforeatICUAdmin, levels = c(1,0), labels = c("yes","no"))
  AFinitforCAPA = factor(AFinitforCAPA, levels = c(1,0), labels = c("yes","no"))
  WhichAF = factor(WhichAF, levels = c("PosaIsa","Posa Isa","Isav","IsavCasp","LAmb","Vori"), labels = c("PosaIsa","Posa Isa","Isav","IsavCasp","LAmb","Vori"))
  Vori = factor(Vori, levels = c(1,0), labels = c("yes","no"))
  Isa = factor(Isa, levels = c(1,0), labels = c("yes","no"))
  LipAmb = factor(LipAmb, levels = c(1,0), labels = c("yes","no"))
  DAmB = factor(DAmB, levels = c(1,0), labels = c("yes","no"))
  Posa = factor(Posa, levels = c(1,0), labels = c("yes","no"))
  Echinocandin = factor(Echinocandin, levels = c(1,0), labels = c("yes","no"))
  Combination = factor(Combination, levels = c(1,0), labels = c("yes","no"))
  AFTreatmentOutcopme = factor(AFTreatmentOutcopme, levels = c("yes","no","partial"))
  PrimaryTreatment = factor(PrimaryTreatment, levels = c(1,0), labels = c("yes","no"))
  #DurationICUDays - continuous
  PataliveDay28or32 = factor(PataliveDay28or32, levels = c(1,0), labels = c("yes","no"))
  SurvivalEOF = factor(SurvivalEOF, levels = c(1,0), labels = c("yes","no"))
  #BALGMODI  - continuous
  LifeSupport = factor(LifeSupport, levels = c("MV","ECMO","NIV","CV","None"), labels = c("mechanical","ECMO","noninvasive", "mechanical & noninvasive","none"),order=TRUE)
})


#get variable names in required format for defining columns when creating table
var_names <- toString(shQuote(colnames(data.reformatted), type = "cmd"))
cat(var_names, "\n")

#First, create a table that contains demographic and clinical characteristics for each group. The target
#variables are sex , age and baselinescore ; the variable trt acts as the grouping variable:
#change output format with:format_to - Possible values are 'Latex', 'Word', 'Raw', 'HTML', 'Console', 'markdown', 'md'.
#Word, can be further processed with e.g. flextable of flextable Gohel (2018).
the_table <- atable::atable(data.reformatted,
                            target_cols = c( "CAPANoCAPA", "DaysfromICUadmintoCAPA", "BALGM1ODI", "BAL_LFD",
                                             "PosBALPCR", "PosBALculture", "PosBronchAspculture", "PosTAAspCult",
                                             "PosTAGM", "PosPCRfromTA", "SerumGM05", "Institution", "Sex", "Age",
                                             "Ethnicity", "BMI30", "UnderlHemOnc", "UnderlSOT",
                                             "UnderlCardiovasc", "UnderlPulm", "UnderlDM", "Smoking",
                                             "Corticos", "Tocilizumab", "Azithro", "MechanVnet",
                                             "NIVent", "ECMOnew", "SystemicAFbeforeatICUAdmin", "AFinitforCAPA",
                                             "WhichAF", "Vori", "Isa", "LipAmb", "DAmB", "Posa", "Echinocandin",
                                             "Combination", "AFTreatmentOutcopme", "PrimaryTreatment", "DurationICUDays",
                                             "PataliveDay28or32", "SurvivalEOF", "LifeSupport" ), 
                            format_to = "Word")

the_table_grouped <- atable::atable(data.reformatted,
                                    target_cols = c( "DaysfromICUadmintoCAPA", "BALGM1ODI", "BAL_LFD",
                                                     "PosBALPCR", "PosBALculture", "PosBronchAspculture", "PosTAAspCult",
                                                     "PosTAGM", "PosPCRfromTA", "SerumGM05", "Institution", "Sex", "Age",
                                                     "Ethnicity", "BMI30", "UnderlHemOnc", "UnderlSOT",
                                                     "UnderlCardiovasc", "UnderlPulm", "UnderlDM", "Smoking",
                                                     "Corticos", "Tocilizumab", "Azithro", "MechanVnet",
                                                     "NIVent", "ECMOnew", "SystemicAFbeforeatICUAdmin", "AFinitforCAPA",
                                                     "WhichAF", "Vori", "Isa", "LipAmb", "DAmB", "Posa", "Echinocandin",
                                                     "Combination", "AFTreatmentOutcopme", "PrimaryTreatment", "DurationICUDays",
                                                     "PataliveDay28or32", "SurvivalEOF", "LifeSupport" ), 
                                    group_col = "CAPANoCAPA", format_to = "Word")

# print in Word with packages flextable and officer
library(flextable)
MyFTable <- flextable::regulartable(data = the_table)
# left aligned first column:
#MyFTable <- flextable::align(MyFTable, align = "left", j = 1)
#save on disc. Not run here:
doc <- officer::read_docx()
doc <- flextable::body_add_flextable(doc, value = MyFTable)
print(doc, target = "basic_table_final_data.docx")

MyFTable_grouped <- flextable::regulartable(data = the_table_grouped)
# left aligned first column:
#MyFTable <- flextable::align(MyFTable, align = "left", j = 1)
#save on disc. Not run here:
doc <- officer::read_docx()
doc <- flextable::body_add_flextable(doc, value = MyFTable_grouped)
print(doc, target = "basic_table_final_data_grouped_CAPA.docx")


