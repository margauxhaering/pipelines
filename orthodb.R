# this script retrieve symbols and human orthologs symbols from OrthoDB from a list of gene names 


library(httr)
library(jsonlite)
library(dplyr)

ortho_df <- read.csv(file.choose(), stringsAsFactors = FALSE, sep = ';')
ortho_df <- ortho_df %>%
  distinct(Name, .keep_all = TRUE)


get_symbol <- function(id) { #function to get the symbol of a gene in orthoDB
  details_url <- paste0("https://data.orthodb.org/v12/ogdetails?id=", id)
  details <- fromJSON(content(GET(details_url), "text"))
  symbol <- details$data$xrefs$name[details$data$xrefs$type == 'Gene_Name']
  if (length(symbol) == 0) return('')
  
  return(symbol)
}

OrthoDB_human <- function(gene_name){
  tryCatch({
    search_url <- paste0("https://data.orthodb.org/v12/genesearch?query=", gene_name) #look for the gene
    response <- GET(search_url)
    gene_data <- fromJSON(content(response, "text"))
    human_taxon_id <- "9606_0"
    gene_id <- gene_data[["gene"]][["gene_id"]][["param"]]
    symbol <- get_symbol(gene_id) #get the symbol of the gene
    if(any(gene_data$orthologs_in_model_organisms$organism$id == human_taxon_id) ==F){  #if no human orthologs :
      symbols ='' 
    }else{
      nb <- which(gene_data$orthologs_in_model_organisms$organism$id == human_taxon_id) #else get the line where human orthologs are
      ortho_ids <- gene_data$orthologs_in_model_organisms$genes[[nb]]$gene_id$param #get the ids
      symbols <- sapply(ortho_ids, get_symbol, USE.NAMES = FALSE) #retrieve the symbols
      symbols <- paste(symbols, collapse = ", ")
    }
    return(data.frame('gene_name' = gene_name, 'symbol' = symbol, 'human_ortho' = symbols))
  }, error = function(e) {
    message("Error occurred for gene: ", gene_name) #if an error, retrieve the problematic gene
    message("Error details: ", e$message)
    return(NULL)
  })
}
result_df <- do.call(rbind, lapply(ortho_df$Name, OrthoDB_human))


write.csv(result_df, "Desktop/orthologs2.tsv", row.names = F)


ortho <- read.csv(file.choose(), stringsAsFactors = FALSE, sep = '\t')
aliases <- read.csv(file.choose(), stringsAsFactors = FALSE, sep = '\t')
df_merged <- merge(aliases, ortho, by = "Name", all.x = TRUE)


gene  = 'LOC111318824'

OrthoDB_human(gene)



#### pour l'histoire lol

OrthoDB_human <- function(gene_name){
  tryCatch({
    search_url <- paste0("https://data.orthodb.org/v12/genesearch?query=", gene_name) #look for the gene
    response <- GET(search_url)
    gene_data <- fromJSON(content(response, "text"))
    human_taxon_id <- "9606"
    metazoa_taxon_id <- "33208"
    if(any(gene_data$orthologs_in_model_organisms$organism$id == '9606_0') ==F){  #if no results
      human_ortholog1 ='' 
    }else{
      nb <- which(gene_data$orthologs_in_model_organisms$organism$id == '9606_0')
      human_ortholog1 <- gene_data$orthologs_in_model_organisms$genes[[nb]]$gene_id$id
      human_ortholog1 <- paste( human_ortholog1, collapse = ", ")
    }
    orthodb_id <- gene_data$orthologs_in_model_organisms$lca$lca_cluster_id[gene_data$orthologs_in_model_organisms$lca$lca_tax_id == metazoa_taxon_id][1] #get the orthodb_id
    orthologs_url <- paste0("https://data.orthodb.org/v12/tab?id=", orthodb_id, "&species=", human_taxon_id) #look for orthologs for homo sapiens only
    response <- GET(orthologs_url)
    if(content(response, "text") ==""){  #if no results
      human_ortholog ='' 
    }else{
      orthologs_data <- read.delim(text = content(response, "text"), header = TRUE, sep = "\t", stringsAsFactors = FALSE) #put results in a table
      human_ortholog <- orthologs_data$pub_gene_id #get the gene names
      human_ortholog <- paste(human_ortholog, collapse = ", ") #paste them together
    }
    human_orthologs <- paste(human_ortholog1, human_ortholog, sep=', ', collapse = ", ")
    human_orthologs <- sapply(strsplit(human_orthologs, ", "), function(x) paste(unique(trimws(x)), collapse = ", "))
    human_orthologs <- sapply(strsplit(human_orthologs, ","), function(x) {
      filtered_genes <- x[!grepl("^\\d+$", trimws(x))]  # Remove purely numeric values
      paste(filtered_genes, collapse = ",")  
    })
    return(human_orthologs)
  }, error = function(e) {
    message("Error occurred for gene: ", gene_name) #if an error, retrieve the problematic gene
    message("Error details: ", e$message)
    return(NULL)
  })
}


