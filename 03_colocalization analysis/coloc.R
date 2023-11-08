#!/usr/bin/Rscript

args <- commandArgs(TRUE)
gwas_data <- args[1]
xQTL_data <- args[2]

suppressMessages(library(tidyverse))
suppressMessages(library(data.table))


############ Read table from GWAS data ############
gwas <- read.table(paste(gwas_data,"coloc.snp","txt",sep="."),header=F,sep="\t")
colnames(gwas) <- c("snp","chr","position","ref","alt","maf","gene","beta","varbeta","pvalue")

############ Read table from xQTL data ############
xQTL = fread(paste(xQTL_data,"coloc.snp","txt",sep=".")) %>% 
       rename(snp  = `V1`,
	          chr = `V2`,
			  position = `V3`,
			  ref = `V4`,
			  alt = `V5`,
			  maf = `V6`,
			  probes = `V7`,
			  beta=`V8`,
			  varbeta=`V9`,
			  pvalue = `V10`
			  ) %>%
		select(snp, chr, position, ref, alt, maf, probes, beta, varbeta, pvalue)

############ Extract lead GWAS snps ############
lead_gSNP = gwas %>%
            group_by(gene) %>%
	        arrange(pvalue) %>%
		    distinct(gene, .keep_all = TRUE) %>%
			rename_at(
				    vars(starts_with("be")),funs(str_replace(.,"beta","beta_eqtl"))) %>%
			rename_at(
					vars(starts_with("var")),funs(str_replace(.,"varbeta","varbeta_eqtl"))) %>%
			rename_at(
				    vars(starts_with("pva")),funs(str_replace(.,"pvalue","pvalue_eqtl")))

############ Extract lead QTL snps ############
lead_xSNP = xQTL %>%
            group_by(probes) %>%
			arrange(pvalue) %>%
		    distinct(probes, .keep_all = TRUE)

############ Defining Colocalization Data Pairs ############
coloc_pair_list = apply(
            lead_gSNP %>%
			mutate(gene2 = gene) %>%
			column_to_rownames("gene2") %>%
			mutate_at(
			      .vars = vars(c("position", ends_with("_eqtl"))),
				  .funs = as.numeric
			),
			MARGIN = 1,
			FUN = function(x) {
				  result                   = as.list(x)
				  result[["maf"]]          = as.numeric(result[["maf"]])
				  result[["position"]]     = as.numeric(result[["position"]])
				  result[["beta_eqtl"]]    = as.numeric(result[["beta_eqtl"]])
				  result[["varbeta_eqtl"]] = as.numeric(result[["varbeta_eqtl"]])
				  result[["pvalue_eqtl"]]  = as.numeric(result[["pvalue_eqtl"]])
				  return(result)
			},
				  simplify = FALSE
)

############ Finding lead xQTL snps that overlapped with lead_gSNP ############
overlapped_xSNP = lapply(
    X = coloc_pair_list,
	    FUN = function(x) {
			 xQTL %>% filter(snp %in% x[["snp"]])
		}
)

############ Calculate xQTL sites that have highest LD value with lead_gSNP ############
lead_probe = mapply(
    FUN = function(coloc_list, meqtl_overlap) {
        if(length(meqtl_overlap[["snp"]]) == 0) {
            return(list(
                probe        = NA,
                beta_meqtl   = NA,
                varbeta      = NA,
                pvalue_meqtl = NA
            ))
        }

        if(length(unique(meqtl_overlap[["probes"]])) == 1) {
            return(
                list(
                    probe         = meqtl_overlap[["probes"]][1],
                    beta_meqtl    = meqtl_overlap[["beta"]][1],
                    varbeta_meqtl = meqtl_overlap[["varbeta"]][1],
                    pvalue_meqtl  = meqtl_overlap[["pvalue"]][1]
                )
            )
        } else {
            lead_esnp = coloc_list[["snp"]]
            lead_mesnps_query = unique(lead_xSNP$snp[lead_xSNP$probes %in% meqtl_overlap$probes])
            ld_with_lead_esnp = sapply(
                X = lead_mesnps_query,
                FUN = function(x) {
                    if (length(x) == 0) {
                        return(0)
                    } else {
                        tryCatch({
                            ld_matrix = LDlinkR::LDpair(
                                var1  = x,
                                var2  = lead_esnp,
                                pop   = "CEU",
                                token = "7ff48dd58066"
                            )
                          return(ld_matrix[["r2"]][1])
                        },
                        error = function(error) {
                            print(paste(x, lead_esnp, error, sep = "\t"))
                            return(0)
                        })
                    }
                },
                simplify = "array"
            )
            max_ld = max(ld_with_lead_esnp)
            if(max_ld == 0 | is.na(max_ld)) {
                return(list(
                    probe        = NA,
                    beta_meqtl   = NA,
                    varbeta      = NA,
                    pvalue_meqtl = NA
                ))
            } else {
  max_ld_snp = names(ld_with_lead_esnp)[which(ld_with_lead_esnp == max_ld)]
  max_ld_probe   = (lead_xSNP %>% filter(snp %in% max_ld_snp) %>% arrange(pvalue))$probes[1]
  max_ld_beta    = (lead_xSNP %>% filter(snp %in% max_ld_snp) %>% arrange(pvalue))$beta[1]
  max_ld_varbeta = (lead_xSNP %>% filter(snp %in% max_ld_snp) %>% arrange(pvalue))$varbeta[1]
  max_ld_pvalue  = (lead_xSNP %>% filter(snp %in% max_ld_snp) %>% arrange(pvalue))$pvalue[1]
                return(
                    list(
                        probe        = max_ld_probe,
                        beta_meqtl   = max_ld_beta,
                        varbeta      = max_ld_varbeta,
                        pvalue_meqtl = max_ld_pvalue
                    )
                )
            }
        }
    },
  coloc_pair_list,
  overlapped_xSNP,
  SIMPLIFY = FALSE
)

eSNP_by_gene = sapply(
    X = unique(gwas[["gene"]]),
    FUN = function(egene) {
        gwas %>% filter(gene == egene)
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

xSNP_by_gene = lapply(
    X = eSNP_by_gene,
    FUN = function(esnps) {
        xQTL %>% filter(snp %in% esnps[["snp"]])
    }
)

coloc_result = mapply(
    FUN = function(df_1, df_2, n_1, n_2, type_1 = "cc", type_2 = "quant") {
        if (nrow(df_1) == 0 | nrow(df_2) == 0) {
            return(
                list(
                    n_snps = 0,
                    PP3    = 0,
                    PP4    = 0
                )
            )
        }
        df_1 = df_1 %>%
            rename(MAF = maf, p = pvalue) %>%
            arrange(p) %>%
            distinct(snp, .keep_all = TRUE) %>%
            select(snp, position, p, beta, varbeta, MAF) %>%
            filter(!is.na(MAF)) %>%
            as.list()
        df_1[["N"]] = n_1
        df_1[["type"]] = type_1

        df_2 = df_2 %>%
            rename(MAF = maf, p = pvalue) %>%
            arrange(p) %>%
            distinct(snp, .keep_all = TRUE) %>%
            select(snp, position, p, beta, varbeta, MAF) %>%
            filter(!is.na(MAF)) %>%
            as.list()
        df_2[["N"]] = n_2
        df_2[["type"]] = type_2

        if (length(df_1[["snp"]])== 0 | length(df_2[["snp"]]) == 0) {
            return(
                list(
                    n_snps = 0,
                    PP3    = 0,
                    PP4    = 0
                )
            )
        }

        if (is.null(coloc::check_dataset(df_1)) & is.null(coloc::check_dataset(df_2))) {
            invisible(capture.output({
                coloc_result = coloc::coloc.abf(
                    dataset1 = df_1,
                    dataset2 = df_2,
					p1 = 1e-04,
					p2 = 1e-04,
					p12 = 1e-04
                )
            }))
            return(
                list(
                    n_snps = coloc_result[["summary"]][["nsnps"]],
                    PP3    = coloc_result[["summary"]][["PP.H3.abf"]],
                    PP4    = coloc_result[["summary"]][["PP.H4.abf"]]
                )
            )
        } else {
            return(
                list(
                    n_snps = 0,
                    PP3    = 0,
                    PP4    = 0
                )
            )
        }
    },
    df_1      = eSNP_by_gene,
    df_2      = xSNP_by_gene,
    n_1       = 200000,
    n_2       = 664,
    SIMPLIFY  = FALSE,
    USE.NAMES = TRUE
)

final_coloc_result_list = mapply(
    FUN = function(coloc_pairs_eqtl, coloc_pairs_meqtl, coloc_result) {
        return(c(
            coloc_pairs_eqtl,
            coloc_pairs_meqtl,
            coloc_result
        ))
    },
    coloc_pairs_eqtl  = coloc_pair_list,
    coloc_pairs_meqtl = lead_probe,
    coloc_result      = coloc_result,
    SIMPLIFY          = FALSE
)

final_coloc_result_table = do.call(rbind, final_coloc_result_list) %>%
    as.data.frame() %>%
		   mutate_at(
			        .vars = vars(c("snp", "chr", "ref", "alt", "gene", "probe")),
					.funs = as.character
					) %>%
		   mutate_at(
				    .vars = vars(-c("snp", "chr", "ref", "alt", "gene", "probe")),
					.funs = as.numeric
				    ) %>%
		  mutate_all(
					.funs = function(x) {
					ifelse(is.na(x) | x == "NA", NA, x)
										}
					) %>%
		  mutate(PP4 = ifelse(is.na(probe), 0, PP4)) %>%
		  arrange(desc(PP4))

write.table(final_coloc_result_table,file = paste(xQTL_data,gwas_data,"coloc","snp.csv",sep="."),col.names = T,row.names = F,sep = "\t",append = F,quote = F)


