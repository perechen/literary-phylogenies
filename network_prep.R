library(tidyverse)
library(tidytext)
library(topicmodels)


## load data
hathi_lda <- readRDS("data/lda_model_22k_texts.rds") # LDA model
meta <- readRDS("data/titles_prepared.rds") # metadata

## get the topic-in-doc table
docprobs <- tidy(hathi_lda,"gamma")

# rm(hathi_lda)
# gc()

## construct "wide" document-term matrix
## average each chunk probabilities 
dtm <- docprobs %>% 
  ## !!!! handles doc ids and filters those notpresent in metadata
  mutate(document = str_remove(document, "_[0-9]*$")) %>% 
  filter(document %in% meta$docid) %>% 
  ## averages topics across chunks
  group_by(document, topic) %>% 
  summarize(gamma=mean(gamma)) %>% 
  ## makes table `wide`
  pivot_wider(names_from = topic,values_from = gamma)

## save ids for later
ids <- dtm$document

## scale + manhatann (delta for speed, for now)
delta_dist=scale(dtm[,-1]) %>% dist(method="manhattan")

## normalize by number of topics
delta_dist <- as.matrix(delta_dist)/200

rownames(delta_dist) <- ids
colnames(delta_dist) <- ids

saveRDS(delta_dist,"distances_delta_lda.rds")

delta_dist <- readRDS("data/distances_delta_lda.rds")

## join metadata to 2k ids
year_names=tibble(docid = rownames(delta_dist)) %>% 
  left_join(meta) %>% 
  select(docid, true_comp_year) %>% 
  rename(year=true_comp_year)

source("scripts/networks_cutoffs_etc.R")
## this function takes a matrix of distances and id+year metadata 
## for each text X in year T
## looks at ALL distances that X has, looks at their distribution
## selects only neighbors of text X (percentile cutoff base determined by alpha, default=0.05)
## returns the average time distance between these neighbors and text X
## apply function to each row of the matrix
temporal_diff=apply(delta_dist,1,function(x) time_distance_neighbours(x, 
                                                                       years=year_names,
                                                                      alpha = 0.005))  


median(temporal_diff[2,]) # 0.065
mean(temporal_diff[2,]) # 0.066

median(temporal_diff[1,]) # 15
mean(temporal_diff[1,]) # 16.4

## plot histogram of median temporal distances
tibble(median_temp_diff = temporal_diff[1,]) %>%
  mutate(mean_time = mean(median_temp_diff),
         median_time = median(median_temp_diff)) %>% 
  ggplot(aes(median_temp_diff)) +
  geom_histogram(bins=50,fill="red") +
  geom_vline(aes(xintercept=mean_time),size=2) +
  geom_vline(aes(xintercept=median_time),linetype=2,size=2) + labs(x = "Median temporal difference between neighboring texts (0.005 quantile)",
                                                                   y="Count") + theme_minimal()



### Let's say our window is 12 years (mean temporal distance between neighbors)

time_window = 12 # how longs link persist
time_step = 1 # how long is the 'rolling' window
start_year = 1870 # star at 
end_year = 1960 # end at (optimally it is max(year) - time_window)
cutoff = 0.001 # cutoff for determining 'meaningful' connections between books


## metadata
novels <- tibble(docid=rownames(delta_dist)) %>% 
  left_join(meta)


while(start_year <= end_year) {
  
  message(paste0("Now at ", start_year, " + ", time_window))
  
  # first select novels in the time frame
  novels_at_i=novels %>%
    filter(true_comp_year >= start_year & true_comp_year < start_year+time_window)
  
  # skip if there is no book in a particular year
  if(!start_year %in% novels_at_i$true_comp_year) {
    start_year = start_year + time_step
    next
  }
  
  ## current ids
  ids = novels_at_i$docid
  
  ## distances in time window
  dist_at_t <- delta_dist[novels_at_i$docid,novels_at_i$docid]
  
  
  ids_with_year = novels_at_i %>% select(docid,true_comp_year)
  
  
  ## consider only distances from a book in a particular year to the future (determined by window)
  lookahead_novels = ids_with_year %>% filter(true_comp_year == start_year)
  lookahead_dist = dist_at_t[lookahead_novels$docid,]
  
  ## handle one book cases
  if(nrow(lookahead_novels) == 1) {
    lookahead_dist = matrix(lookahead_dist, nrow=1)
    colnames(lookahead_dist) = colnames(dist_at_t)
  }
  
  # determine the distance that is considered "meaningful" based on the percentile cutoff
  median_cutoff = apply(lookahead_dist, 1, percentile_cutoff,p=cutoff) %>% median() %>% round(4)
  
  # function gets the nearest neighbors of each text based on cutoff
  nn=apply(lookahead_dist, 1, nearest_neighbours,cutoff=median_cutoff)
  
  
  ## construct edge list from the results, removes duplicates
  network_at_i = bind_rows(nn) %>% 
    unnest(cols=neighbors) %>%
    mutate(time = start_year) %>% 
    group_by(target_text,neighbors) %>% 
    mutate(edge_id = paste(sort(unique(c(target_text,neighbors))), collapse=" ")) %>% # get connections sorted and pasted (we don't care of the direction of the link and don't want to count the same thing twice)
    group_by(edge_id) %>% 
    distinct(edge_id,.keep_all = T) %>% # remove repeating connections
    ungroup() %>% 
    select(edge_id, time) %>% 
    separate(edge_id, c("target", "source"), sep = " ")
  
  ## append to file
  write_tsv(network_at_i,file="data/edges_list.csv",append = T)
  
  start_year = start_year + time_step
  
}



####### ADDENDUM #########

## nodes to authors

authors = meta %>% select(docid, author) %>% mutate(author=str_remove_all(author," |,"))
edges = read_tsv("data/edges_list_005.csv",col_names = c("target","source", "time"))

edges
edges_aut =edges %>%
  left_join(authors, by=c("target"="docid")) %>%
  rename(author_source=author) %>%
  mutate(author_source=case_when(is.na(author_source) ~ target,
                                                                                                                                    T ~ author_source)) %>% left_join(authors,by=c("source"="docid")) %>% 
  rename(author_target=author) %>% 
  mutate(author_target=case_when(is.na(author_target) ~ source,
                                 T ~ author_target)) %>% select(author_source,author_target,time)

## remove duplicated connections
edges_aut_upd = edges_aut %>%
  filter(author_target != author_source) %>% # remove self-loops
  rowwise() %>% 
  mutate(edge_id = paste(sort(unique(c(author_source,author_target,time))), collapse=" ")) %>% # get connections sorted and pasted (we don't care of the direction of the link and don't want to count the same thing twice)
  group_by(edge_id) %>% 
  distinct(edge_id,.keep_all = T) %>% # remove repeating connections
  select(-edge_id)

edges_aut_upd<-edges_aut_upd %>% ungroup %>% select(-edge_id)

write_tsv(edges_aut_upd, "data/edges_001_authors.tsv", col_names=F)
