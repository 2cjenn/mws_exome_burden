
library(data.table)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(tictoc)
library(jsonlite)
library(httr2)

authenticate <- request("https://eggwhite.eu.auth0.com/oauth/token") |>
  req_body_json(
    list(
      "client_id" = "",
      "client_secret" = "",
      "audience" = "predictions.ai.bio-prodict.nl",
      "grant_type" = "client_credentials")) |>
  req_perform() |> 
  resp_body_json() 

access_token <- authenticate$access_token

# Load variants file
helix_uniques <- data.table::fread("/shared/MWS_Regeneron_Data/MWS_regeneron/Jennifer/annotation/helix/MWS_variants.txt", header=TRUE)

mws_variants <- helix_uniques %>%
  mutate(variant = paste0(Feature, "/", residue, "/", acid2)) %>%
  pull(variant)

chunks <- split(mws_variants, ceiling(seq_along(mws_variants)/5000))
n <- length(chunks)
print(n)
pages <- vector("list", n)

for(i in seq(1, n)) {
  print(paste0("chunk ", i))
  chunk <- chunks[[i]]
  variants <- paste(chunk, collapse="\n")
  
  annotations <- request("https://api.helixlabs.ai/batch") |>
    req_headers(
      Authorization = paste0("Bearer ", access_token)
    ) |>
    req_body_form("variants" = variants) |>
    req_perform() |>
    resp_body_string() |> 
    fromJSON(simplifyVector=TRUE, flatten=TRUE)
  
  if(is.data.frame(annotations)){
    pages[[i]] <- annotations %>%
      filter(!is.na(transcript_id)) %>%
      distinct()
  }
  Sys.sleep(2) # courtesy nap
}

total <- rbind_pages(pages) %>% 
  distinct()

helix <- total %>%
  left_join(
    helix_uniques %>% 
      select(Feature, HGVSp, acid1, residue, acid2) %>% 
      mutate(residue = as.numeric(residue)
             ), 
    join_by(transcript_id == Feature, residue_number == residue, variant_type == acid2))

data.table::fwrite(helix, "/shared/MWS_Regeneron_Data/MWS_regeneron/Jennifer/annotation/helix/MWS_anno.txt")