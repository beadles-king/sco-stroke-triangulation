# =============================================================================
# Reproducibility repository - Wang et al., transportability-aware triangulation
# of metabolic stroke evidence across the Shanghai Cooperation Organisation.
# Paths are resolved by config.R; see README.md.
# =============================================================================
if (!exists("RAWDAT")) {
  if (file.exists("config.R")) source("config.R") else
    stop("Run this script from the code/ directory (the folder containing config.R).")
}

# ===========================================================================
# Independent verification of the geographical coherence metrics in Table 1.
#
# Derives land-border adjacency from Natural Earth polygons, computes connected
# components for each of the eleven blocs, and compares the result against the
# tabulated values. Runs top to bottom with no editing required.
#
# Requires: sf, rnaturalearth, rnaturalearthdata, igraph
#   install.packages(c("sf","rnaturalearth","rnaturalearthdata","igraph"))
# ===========================================================================

library(sf)
library(rnaturalearth)
library(igraph)

sf_use_s2(TRUE)

## -------------------------------------------------------- membership ------
## ISO 3166-1 alpha-3. Counts must match n_members in Table 1.

membership <- list(
  
  `Nordic Council` = c("DNK","FIN","ISL","NOR","SWE"),
  
  `GCC` = c("BHR","KWT","OMN","QAT","SAU","ARE"),
  
  `EAEU` = c("ARM","BLR","KAZ","KGZ","RUS"),
  
  `EU-27` = c("AUT","BEL","BGR","HRV","CYP","CZE","DNK","EST","FIN","FRA","DEU",
              "GRC","HUN","IRL","ITA","LVA","LTU","LUX","MLT","NLD","POL","PRT",
              "ROU","SVK","SVN","ESP","SWE"),
  
  # Timor-Leste acceded on Oct 26, 2025, before the Jan 1, 2026 census date.
  `ASEAN` = c("BRN","KHM","IDN","LAO","MYS","MMR","PHL","SGP","THA","TLS","VNM"),
  
  # Full members only. Turkmenistan is an associate member and is excluded;
  # Ukraine's withdrawal completed in 2018 and Georgia's in 2009.
  `CIS` = c("ARM","AZE","BLR","KAZ","KGZ","MDA","RUS","TJK","UZB"),
  
  # Syria's membership was restored in 2023.
  `Arab League` = c("DZA","BHR","COM","DJI","EGY","IRQ","JOR","KWT","LBN","LBY",
                    "MRT","MAR","OMN","PSE","QAT","SAU","SOM","SDN","SYR","TUN",
                    "ARE","YEM"),
  
  `SAARC` = c("AFG","BGD","BTN","IND","MDV","NPL","PAK","LKA"),
  
  `African Union` = c("DZA","AGO","BEN","BWA","BFA","BDI","CPV","CMR","CAF","TCD",
                      "COM","COG","COD","CIV","DJI","EGY","GNQ","ERI","SWZ","ETH",
                      "GAB","GMB","GHA","GIN","GNB","KEN","LSO","LBR","LBY","MDG",
                      "MWI","MLI","MRT","MUS","MAR","MOZ","NAM","NER","NGA","RWA",
                      "STP","SEN","SYC","SLE","SOM","ZAF","SSD","SDN","TGO","TUN",
                      "UGA","TZA","ZMB","ZWE"),
  
  `SCO` = c("BLR","CHN","IND","IRN","KAZ","KGZ","PAK","RUS","TJK","UZB"),
  
  `Commonwealth` = c("ATG","AUS","BHS","BGD","BRB","BLZ","BWA","BRN","CMR","CAN",
                     "CYP","DMA","SWZ","FJI","GAB","GMB","GHA","GRD","GUY","IND",
                     "JAM","KEN","KIR","LSO","MWI","MYS","MDV","MLT","MUS","MOZ",
                     "NAM","NRU","NZL","NGA","PAK","PNG","RWA","KNA","LCA","VCT",
                     "WSM","SYC","SLE","SGP","SLB","ZAF","LKA","TGO","TON","TTO",
                     "TUV","UGA","GBR","TZA","VUT","ZMB"),
  
  `BRICS` = c("BRA","RUS","IND","CHN","ZAF","EGY","ETH","IRN","ARE","IDN"),
  
  `OIC` = c("AFG","ALB","DZA","AZE","BHR","BGD","BEN","BRN","BFA","CMR","TCD",
            "COM","CIV","DJI","EGY","GAB","GMB","GIN","GNB","GUY","IDN","IRN",
            "IRQ","JOR","KAZ","KWT","KGZ","LBN","LBY","MYS","MDV","MLI","MRT",
            "MAR","MOZ","NER","NGA","OMN","PAK","PSE","QAT","SAU","SEN","SLE",
            "SOM","SDN","SUR","SYR","TJK","TGO","TUN","TUR","TKM","UGA","ARE",
            "UZB","YEM")
)

expected <- data.frame(
  bloc = c("Nordic Council","GCC","EAEU","CIS","EU-27","ASEAN","SAARC",
           "Arab League","African Union","SCO","Commonwealth","BRICS","OIC"),
  n_members    = c(5, 6, 5,  9, 27, 11, 8, 22, 54, 10, 56, 10, 57),
  n_components = c(3, 2, 2,  2,  5,  3, 3,  4,  7,  1, 39,  8, 11),
  largest_n    = c(3, 5, 4,  8, 22,  9, 6, 18, 48, 10, 12,  3, 43),
  stringsAsFactors = FALSE
)

stopifnot(all(lengths(membership) == expected$n_members[
  match(names(membership), expected$bloc)]))

## -------------------------------------------------------- geometry --------

world <- ne_countries(scale = "medium", type = "countries", returnclass = "sf")

nmcol <- if ("admin" %in% names(world)) "admin" else "name"

# Natural Earth codes iso_a3 as "-99" for several features, in two distinct
# situations that must be handled differently:
#
#  (a) recognised sovereign states whose code is simply absent from the field,
#      notably France and Norway. iso_a3_eh or adm0_a3 recovers these.
#
#  (b) self-declared or disputed territories carried as records separate from
#      the state they are attributed to, notably Somaliland and Western Sahara.
#      These have no code of their own to recover. Dropping them truncates the
#      parent country: without Somaliland, Djibouti sits ~611 km from what
#      remains of Somalia and their shared border disappears; without Western
#      Sahara, Morocco and Mauritania cease to be neighbours.
#
# We assign (b) to the parent location used by the GBD, which models none of
# these territories separately, so that the polygon set matches the location
# hierarchy used everywhere else in this study. Western Sahara carries a valid
# code of its own (ESH) rather than "-99", so it is caught by name and by code.
# Reassigning it restores the Morocco-Mauritania border; this changes no
# component count, because Mauritania is already connected through Mali,
# Senegal, and Algeria, but it removes an otherwise unexplained near miss.

pick <- function(x, y) ifelse(is.na(x) | x == "-99", y, x)
world$ISO3 <- pick(world$iso_a3,
                   pick(if ("iso_a3_eh" %in% names(world)) world$iso_a3_eh else NA,
                        if ("adm0_a3"   %in% names(world)) world$adm0_a3   else NA))

parent <- c(
  "Somaliland"      = "SOM",
  "W. Sahara"       = "MAR",
  "Western Sahara"  = "MAR",
  "N. Cyprus"       = "CYP",
  "Northern Cyprus" = "CYP",
  "Kosovo"          = "SRB"
)
hit <- match(world[[nmcol]], names(parent))
use <- !is.na(hit) & (is.na(world$ISO3) |
                        world$ISO3 %in% c("-99", "SOL", "SAH", "ESH", "CYN", "KOS"))
world$ISO3[use] <- unname(parent[hit[use]])

world$ISO3[world$ISO3 == "PSX"] <- "PSE"   # Natural Earth code for Palestine

unresolved <- world[[nmcol]][is.na(world$ISO3) | world$ISO3 == "-99"]
if (length(unresolved)) {
  cat("Features dropped for want of a country code:",
      paste(sort(unresolved), collapse = "; "), "\n")
}
world <- world[!is.na(world$ISO3) & world$ISO3 != "-99", ]
world <- st_make_valid(world)

# Aggregate rather than deduplicate: a country may legitimately be carried as
# several records, and keeping only the first truncates it.
dups <- sort(names(which(table(world$ISO3) > 1)))
if (length(dups)) {
  cat("ISO codes carried as multiple records (aggregated):",
      paste(dups, collapse = " "), "\n")
}
ids  <- sort(unique(world$ISO3))
geom <- do.call(c, lapply(ids, function(k)
  st_cast(st_union(st_geometry(world)[world$ISO3 == k]), "MULTIPOLYGON")))
world <- st_sf(ISO3 = ids, geometry = geom)

## -------------------------------------------------------- components ------
## Contact is a non-empty polygon intersection. No distance tolerance is used:
## a tolerance wide enough to matter is also wide enough to join pairs that
## must stay separate, the tightest constraint being the Johor Strait between
## Malaysia and Singapore at roughly 600 m.
##
## An earlier version failed on Djibouti-Somalia and the cause was not the
## predicate but the polygon set: see the aggregation step above. Any pair that
## is close but not intersecting is now reported rather than silently dropped,
## so a repeat of that failure is visible instead of being absorbed by a
## tolerance.

NEAR_MISS_M <- 5000   # report, do not join

component_stats <- function(iso, bloc_name) {
  
  present <- iso[iso %in% world$ISO3]
  absent  <- setdiff(iso, world$ISO3)
  
  sub <- world[match(present, world$ISO3), ]
  tch <- st_intersects(sub, sparse = TRUE)
  tch <- lapply(seq_along(tch), function(i) setdiff(tch[[i]], i))
  
  # flag near misses: pairs within NEAR_MISS_M but not intersecting
  nearby <- st_is_within_distance(sub, sub, dist = NEAR_MISS_M, sparse = TRUE)
  nm <- character(0)
  for (i in seq_along(nearby)) {
    j <- setdiff(nearby[[i]], c(i, tch[[i]]))
    j <- j[j > i]
    if (length(j)) nm <- c(nm, paste0(present[i], "-", present[j]))
  }
  
  edges <- do.call(rbind, lapply(seq_along(tch), function(i) {
    if (length(tch[[i]]) == 0) return(NULL)
    cbind(present[i], present[tch[[i]]])
  }))
  
  # Absent members are island or micro-states missing from the 1:50m polygon
  # set; they can only ever be isolated vertices, so add them back as such.
  g <- make_empty_graph(directed = FALSE) + vertices(iso)
  if (!is.null(edges)) g <- add_edges(g, as.vector(t(edges)))
  g <- simplify(g)
  
  cm <- components(g)
  
  list(near_miss = nm,
       row = data.frame(
         bloc         = bloc_name,
         n_members    = length(iso),
         n_components = cm$no,
         largest_n    = max(cm$csize),
         largest_pct  = round(100 * max(cm$csize) / length(iso)),
         stringsAsFactors = FALSE),
       absent = absent,
       sizes  = sort(cm$csize, decreasing = TRUE))
}

res <- Map(component_stats, membership, names(membership))

observed <- do.call(rbind, lapply(res, `[[`, "row"))
rownames(observed) <- NULL

## -------------------------------------------------------- comparison ------

check <- merge(observed, expected, by = c("bloc", "n_members"),
               suffixes = c("_observed", "_tabulated"))
check$components_agree <- check$n_components_observed == check$n_components_tabulated
check$largest_agree    <- check$largest_n_observed    == check$largest_n_tabulated

cat("\n=== Comparison against Table 1 ===\n")
print(check[order(-check$largest_pct),
            c("bloc","n_members","n_components_observed","n_components_tabulated",
              "largest_n_observed","largest_n_tabulated",
              "components_agree","largest_agree")],
      row.names = FALSE)

cat("\n=== Component size profiles ===\n")
for (b in names(res)) {
  cat(sprintf("%-16s %s\n", b, paste(res[[b]]$sizes, collapse = " ")))
}

cat("\n=== Members absent from the polygon set (added back as isolates) ===\n")
for (b in names(res)) {
  a <- res[[b]]$absent
  if (length(a)) cat(sprintf("%-16s %s\n", b, paste(a, collapse = " ")))
}

if (all(check$components_agree) && all(check$largest_agree)) {
  cat("\nAll eleven blocs reproduce. Table 1 Panel B verified.\n")
} else {
  cat("\nDISCREPANCY. Resolve before submission — see the checks below.\n")
}

## -------------------------------------------------------- sanity ---------
## Pairs that MUST NOT be adjacent (bridge, causeway, or strait, not a land
## border) and pairs that MUST be adjacent. A failure here means the geometry,
## not the bloc lists, is at fault.

# Note: GBR-IRL is NOT in this list. The United Kingdom shares a land border
# with Ireland through Northern Ireland, so adjacency there is correct. No bloc
# analysed contains both (Ireland left the Commonwealth in 1949), so it does not
# affect any component count.
must_not <- list(c("MYS","SGP"), c("DNK","SWE"), c("BHR","SAU"),
                 c("ITA","GRC"), c("EGY","SAU"), c("AUS","NZL"),
                 c("IDN","PHL"), c("CYP","TUR"), c("MDA","RUS"),
                 c("COM","MOZ"))
must     <- list(c("CHN","RUS"), c("CHN","PAK"), c("IND","PAK"),
                 c("IRN","PAK"), c("BLR","RUS"), c("KGZ","TJK"),
                 c("MYS","IDN"), c("DEU","DNK"), c("FIN","SWE"),
                 c("DJI","SOM"), c("PSE","JOR"), c("PSE","EGY"),
                 c("GUY","SUR"), c("IDN","TLS"), c("MAR","MRT"),
                 c("AZE","RUS"))

adjacent <- function(a, b) {
  if (!all(c(a, b) %in% world$ISO3)) return(NA)
  as.logical(lengths(st_intersects(world[world$ISO3 == a, ],
                                   world[world$ISO3 == b, ])) > 0)
}

cat("\n=== Near misses (within", NEAR_MISS_M, "m but not joined) ===\n")
for (b in names(res)) {
  if (length(res[[b]]$near_miss))
    cat(sprintf("%-16s %s\n", b, paste(res[[b]]$near_miss, collapse = " ")))
}

cat("\n=== Sanity checks ===\n")
for (p in must)     cat(sprintf("  adjacent %s-%s : %s (expect TRUE)\n",
                                p[1], p[2], adjacent(p[1], p[2])))
for (p in must_not) cat(sprintf("  adjacent %s-%s : %s (expect FALSE)\n",
                                p[1], p[2], adjacent(p[1], p[2])))

## Note on Malaysia-Indonesia: they DO share a land border on Borneo, so the
## pair appears in `must`, not `must_not`. Singapore is the maritime case.

## -------------------------------------------------------- notes ----------
# If a bloc disagrees by exactly one component, inspect its edge list directly:
#
#   b   <- "EU-27"
#   iso <- membership[[b]]
#   sub <- world[match(iso[iso %in% world$ISO3], world$ISO3), ]
#   tch <- st_touches(sub, sparse = TRUE)
#   setNames(lapply(tch, function(i) sub$ISO3[i]), sub$ISO3)
#
# Known Natural Earth quirks worth checking if something is off:
#   - France's polygon may include French Guiana, creating a South American
#     land border. Irrelevant for EU-27 (no South American members) but check
#     if you add other blocs.
#   - The Netherlands and the United Kingdom may carry Caribbean territories.
#   - Russia's Kaliningrad exclave borders Poland and Lithuania; this matters
#     only if Russia and an EU member are ever in the same bloc list.
#   - Very short borders (Zambia-Botswana at Kazungula) are sometimes rendered
#     as a point rather than a line, which st_touches still accepts but some
#     alternative predicates do not.