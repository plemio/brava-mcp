"""What SQL cannot know about this table.

Shipping a database removes the ceiling on what can be asked, and moves the
whole risk to interpretation: a syntactically perfect query can still be wrong
about the science. Every item here is a mistake that reads as a valid result,
and several were made during this server's own development.

Carried by `schema()` so the traps arrive with the columns, not in a README
nobody fetched.
"""

BETA_PROVENANCE = (
    "beta and se ALWAYS come from the inverse-variance-weighted Burden "
    "meta-analysis, including on rows where you read p_skato or p_skat. There is "
    "no SKAT-O effect size. Reporting beta as 'the SKAT-O effect' is wrong."
)

SYNONYMOUS_IS_A_CONTROL = (
    "mask='synonymous' is a CALIBRATION CONTROL, not a biological category. A "
    "significant synonymous result means residual test inflation. Exclude it from "
    "any 'which genes are associated' query unless you are auditing calibration."
)

ANCESTRY_OVERLAP = (
    "ancestry='All' is the cross-ancestry meta-analysis and 'non_EUR' pools four "
    "of the five superpopulations. Both overlap the individual strata, so "
    "counting or averaging across every ancestry double-counts the same people. "
    "Aggregate over EUR, AFR, AMR, EAS, SAS only."
)

HET_ONLY_ON_META = (
    "p_het is the cross-cohort heterogeneity test and is only populated on "
    "ancestry='All'. It is null elsewhere by construction, not by omission."
)

UNDERFLOW = (
    "p = 0 is a real value, not missing: SAIGE's p-value underflows below "
    "5e-324, and upstream floors it rather than nulling it. Those are the MOST "
    "significant results in the table. Order by p ASC puts them first, which is "
    "correct; treating 0 as 'no result' inverts the answer."
)

NULLS_ARE_FAILURES = (
    "A null p or beta is a degenerate SAIGE stratum, a computational failure, "
    "not evidence of no association. Filter them out; do not read them as null "
    "results."
)

ROWS_ARE_NOT_FINDINGS = (
    "Every gene-trait pair appears up to 12 times, once per mask x MAF cutoff. "
    "count(*) counts tests, not findings. Use count(DISTINCT gene) or "
    "count(DISTINCT (gene, trait)) unless you mean tests."
)

POWER_NOT_BIOLOGY = (
    "Strata differ in size by an order of magnitude (EUR ~555k against AFR ~38k "
    "for LDL cholesterol). A result significant in one stratum and not another is "
    "usually a power difference, not an ancestry-specific effect: a difference in "
    "significance is not a significant difference. Compare betas and confidence "
    "intervals before claiming specificity."
)

MULTIPLE_TESTING = (
    "Thresholds from the flagship paper: gene x mask Bonferroni 1.39e-7, "
    "gene-level Cauchy 2.5e-6, suggestive 1e-4, variant-level 1.82e-8. p<0.05 is "
    "meaningless on a table of 65 million tests."
)

AMBIGUOUS_SYMBOLS = (
    "25 gene symbols are shared by more than one Ensembl gene. Join or filter on "
    "ensg when identity matters; filtering on symbol can silently merge two genes "
    "or hide the one that was never tested."
)

ALL = [
    BETA_PROVENANCE,
    SYNONYMOUS_IS_A_CONTROL,
    ANCESTRY_OVERLAP,
    HET_ONLY_ON_META,
    UNDERFLOW,
    NULLS_ARE_FAILURES,
    ROWS_ARE_NOT_FINDINGS,
    POWER_NOT_BIOLOGY,
    MULTIPLE_TESTING,
    AMBIGUOUS_SYMBOLS,
]

RECIPES = [
    ("what does one gene do",
     "SELECT trait, mask, maf, p_skato, beta FROM results "
     "WHERE gene='PCSK9' AND ancestry='All' AND mask<>'synonymous' "
     "AND p_skato IS NOT NULL ORDER BY p_skato LIMIT 20"),
    ("top genes for one trait, one row per gene",
     "SELECT gene, min(p_skato) p, arg_min(beta, p_skato) beta FROM results "
     "WHERE trait_id='T2Diab' AND ancestry='All' AND mask<>'synonymous' "
     "GROUP BY gene ORDER BY p LIMIT 20"),
    ("most pleiotropic genes",
     "SELECT gene, count(DISTINCT trait) traits FROM results "
     "WHERE ancestry='All' AND p_skato < 1.39e-7 AND mask<>'synonymous' "
     "GROUP BY gene ORDER BY traits DESC LIMIT 20"),
    ("screen candidates, keeping the ones that clear nothing",
     "SELECT gene, min(p_skato) p FROM results WHERE ancestry='All' "
     "AND trait_id='LDLC' AND gene IN ('PCSK9','LDLR','ACAN') GROUP BY gene ORDER BY p"),
    ("significant in one ancestry and not another (read as a lead, see the power trap)",
     "SELECT a.gene, a.trait, a.p_skato FROM results a "
     "WHERE a.ancestry='AFR' AND a.p_skato < 2.5e-6 AND NOT EXISTS ("
     "SELECT 1 FROM results e WHERE e.ancestry='EUR' AND e.gene_idx=a.gene_idx "
     "AND e.pheno=a.pheno AND e.p_skato < 2.5e-6) ORDER BY a.p_skato LIMIT 20"),
    ("who contributed to a trait's analysis",
     "SELECT biobank, country, sum(n) n, sum(cases) cases FROM biobank_sizes "
     "WHERE trait_id='T2Diab' GROUP BY biobank, country ORDER BY n DESC"),
]
