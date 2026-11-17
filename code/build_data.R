studies = read.csv("../data/raw/ctgov/studies.csv")
ratios = read.csv("../data/raw/ctgov/hazard_ratios.csv")
arms = read.csv("../data/raw/ctgov/arms.csv")

hemonc = "../data/raw/hemonc_2026-08-28/"
drugs = read.csv(paste0(hemonc, "tables/drugs.csv"))
concepts = read.csv(paste0(hemonc, "concept_stage.csv"))
synonyms = read.csv(paste0(hemonc, "concept_synonym_stage.csv"))

pfs.re = "progression[- ]?free|\\bPFS\\b"
os.re = "overall survival|\\bOS\\b"
other.re = paste0("(disease|event|relapse|recurrence|invasive disease)",
                  "[- ]?free survival|\\b(DFS|EFS|RFS|iDFS)\\b")
not.pfs.re = paste0("PFS ?2\\b|PFS-2|second[- ](treatment|progression|PFS)",
                    "|to second|next[- ]line|on subsequent|intracranial",
                    "|\\bCNS\\b|brain|GRFS|GVHD")

found = function(pattern, x) grepl(pattern, x, ignore.case = TRUE)
norm.title = function(x) gsub("[^a-z0-9]", "", tolower(x))

# PFS and OS hazard ratios with a two-sided 95% CI between two groups
tt = ratios$outcome_title
is.os = found(os.re, tt) & !found(pfs.re, tt) & !found(other.re, tt)
is.pfs = !is.os & found(pfs.re, tt) & !found(not.pfs.re, tt)
two = lengths(strsplit(ratios$group_ids, "|", fixed = TRUE)) == 2
keep = (is.os | is.pfs) & two & ratios$ci_pct %in% 95 &
  ratios$ci_sides != "ONE_SIDED" & !is.na(ratios$ci_lower) &
  !is.na(ratios$ci_upper)
hr = ratios[keep, ]
hr$ep = ifelse(is.os[keep], "os", "pfs")
hr$row = 1:nrow(hr)

arm.titles = strsplit(hr$group_titles, " || ", fixed = TRUE)
hr$arm1 = sapply(arm.titles, function(x) x[1])
hr$arm2 = sapply(arm.titles, function(x) x[2])
hr$key = paste(pmin(norm.title(hr$arm1), norm.title(hr$arm2)),
               pmax(norm.title(hr$arm1), norm.title(hr$arm2)))
# some studies post several HRs for the same endpoint and arms (eg subgroups)
# we count them and keep the 1st one
hr$n.hr = ave(hr$row, hr$nct_id, hr$ep, hr$key, FUN = length)
hr = hr[!duplicated(hr[, c("nct_id", "ep", "key")]), ]

# a PFS and an OS hazard ratio on the same two groups, the first one posted
pairs = merge(hr[hr$ep == "pfs", ], hr[hr$ep == "os", ],
              by = c("nct_id", "key"), suffixes = c(".pfs", ".os"))
pairs = pairs[order(pairs$row.pfs), ]
n.pairs = table(pairs$nct_id)
pairs = pairs[!duplicated(pairs$nct_id), ]
pairs = merge(pairs, studies, by = "nct_id")

trials = data.frame(nct.id = pairs$nct_id,
                    title = pairs$brief_title,
                    conditions = pairs$conditions,
                    start.year = as.numeric(substr(pairs$start_date, 1, 4)),
                    allocation = pairs$allocation,
                    phase = pairs$phases,
                    n.pairs = as.numeric(n.pairs[pairs$nct_id]),
                    pfs.n.hr = pairs$n.hr.pfs, os.n.hr = pairs$n.hr.os,
                    arm1 = pairs$arm1.os, arm2 = pairs$arm2.os,
                    pfs.outcome = pairs$outcome_title.pfs,
                    pfs.hr = pairs$hr.pfs, pfs.lower = pairs$ci_lower.pfs,
                    pfs.upper = pairs$ci_upper.pfs,
                    os.outcome = pairs$outcome_title.os,
                    os.hr = pairs$hr.os, os.lower = pairs$ci_lower.os,
                    os.upper = pairs$ci_upper.os)

groups = arms[arms$source != "protocol_arm", ]
groups = groups[!duplicated(paste(groups$nct_id, groups$label)), ]
desc1 = groups$description[match(paste(trials$nct.id, trials$arm1),
                                 paste(groups$nct_id, groups$label))]
desc2 = groups$description[match(paste(trials$nct.id, trials$arm2),
                                 paste(groups$nct_id, groups$label))]
desc1[is.na(desc1)] = ""
desc2[is.na(desc2)] = ""

# drug names and synonyms (codes such as MK-3475) from HemOncKB,
# plus short forms that some arm titles use
components = concepts[concepts$concept_class_id == "Component", ]
syn = merge(synonyms, components, by.x = "synonym_concept_code",
            by.y = "concept_code")
lexicon = data.frame(name = c(syn$synonym_name, drugs$drug),
                     drug = c(syn$concept_name, drugs$drug))
lexicon = lexicon[nchar(lexicon$name) >= 5, ]
short = data.frame(name = c("pembro", "nivo", "atezo", "durva", "ipi"),
                   drug = c("Pembrolizumab", "Nivolumab", "Atezolizumab",
                            "Durvalumab", "Ipilimumab"))
lexicon = rbind(lexicon, short)
lexicon$key = paste0(" ", gsub("[^a-z0-9]+", " ", tolower(lexicon$name)), " ")
ici.drugs = drugs$drug[grepl("PD-1|PD-L1|CTLA-4", drugs$main_class)]
lexicon$ici = lexicon$drug %in% ici.drugs

# "placebo matched to X" does not mean the arm gets X
clean = function(x) {
  x = tolower(x)
  x = gsub(paste0("placebo[ -]+(matched to|matched with|matching|to match",
                  "|for|of|to) +[a-z0-9-]+"), " ", x)
  x = gsub("[a-z0-9-]+[- ]matching +placebo", " ", x)
  paste0(" ", gsub("[^a-z0-9]+", " ", x), " ")
}

# known: the text names a drug, or a placebo, supportive care or observation
match.arm = function(x) {
  x = clean(x)
  drug = rep(FALSE, length(x))
  ici = rep(FALSE, length(x))
  for (k in 1:nrow(lexicon)) {
    hit = grepl(lexicon$key[k], x, fixed = TRUE)
    drug = drug | hit
    if (lexicon$ici[k]) ici = ici | hit
  }
  other = grepl(" placebo | pbo | bsc | supportive care | observation ", x)
  data.frame(known = drug | other, ici = ici)
}

title1 = match.arm(trials$arm1)
text1 = match.arm(desc1)
title2 = match.arm(trials$arm2)
text2 = match.arm(desc2)

# an arm is read from its title, and from its description only if the title
# says nothing. An ICI seen only in a description can be crossover after
# progression, so those comparisons stay undetermined.
known = (title1$known | text1$known) & (title2$known | text2$known)
from.desc = (!title1$known & text1$ici) | (!title2$known & text2$ici)
trials$ici.status = ifelse(!known | from.desc, "undetermined",
                           ifelse(title1$ici & title2$ici, "ICI in both arms",
                                  ifelse(title1$ici | title2$ici,
                                         "ICI in one arm", "no ICI")))

tumors = c(hematologic = paste0("lymphoma|leuk|myeloma|myelodysplastic",
                                "|hodgkin|lymphocytic|myeloid|myelofibrosis",
                                "|macroglobulin|waldenstr|amyloidosis",
                                "|plasma cell|gvhd|graft"),
           lung = "lung|NSCLC|SCLC",
           breast = "breast",
           prostate = "prostat",
           colorectal = "colorectal|colon|rectal",
           gynecologic = "ovar|cervical|endometri|uterine|fallopian|peritoneal",
           gastroesophageal = "gastric|stomach|esophag|gastroesoph",
           "liver or biliary" = "hepatocellular|liver|biliary|cholangio",
           "renal or urothelial" = "renal|kidney|urothelial|bladder",
           melanoma = "melanoma",
           "head and neck" = paste0("head and neck|nasopharyn|oropharyn",
                                    "|hypopharyn|laryn|oral cavity"),
           pancreas = "pancrea")

# first matching category wins, from the conditions and then from the title
trials$tumor = "other"
for (k in 1:length(tumors)) {
  hit = trials$tumor == "other" & found(tumors[k], trials$conditions)
  trials$tumor[hit] = names(tumors)[k]
}
for (k in 1:length(tumors)) {
  hit = trials$tumor == "other" & found(tumors[k], trials$title)
  trials$tumor[hit] = names(tumors)[k]
}

write.csv(trials, "../data/processed/trials_first_pass.csv", row.names = FALSE)
