library(tidyverse)
library(vegan)
otu_tab11 = (data_list[[1]]) # sample by taxa
otu_tab = as.data.frame(t(otu_tab))
otu_tab =  tibble(Group =  rownames(otu_tab), otu_tab)

otu_tab <- otu_tab %>%
  select(Group, starts_with("sp")) %>%
  pivot_longer(-Group) %>%
  group_by(Group) %>%
  mutate(total = sum(value)) %>%
  filter(total > 1800) %>%
  group_by(name) %>%
  mutate(total = sum(value)) %>%
  filter(total != 0) %>%
  ungroup() %>%
  select(-total)

otu_tab_data <- otu_tab %>%
  pivot_wider(names_from="name", values_from="value") %>%
  select(-Group) %>%
  as.data.frame()
raremax = min(rowSums(otu_tab_data))
otu_tab_rarefy <- rarefy(otu_tab_data, raremax) %>%
  as_tibble(rownames = "Group") %>%
  select(Group, vegan =value)  %>%
  pivot_wider()


rare_data <- otu_tab_rarefy %>%
  pivot_wider(names_from="Samples", values_from="value") %>%
  as.data.frame()


 





days_wanted <- c(0:9, 141:150)
shared <- read_tsv("Power_Project/Datasets/mice.shared.txt") %>%
  select(Group, starts_with("Otu")) %>%
  mutate(day = str_replace(Group, ".*D", "")) %>%
  filter(day %in% days_wanted) %>%
  select(-day) %>%
  pivot_longer(-Group) %>%
  group_by(Group) %>%
  mutate(total = sum(value)) %>%
  filter(total > 1800) %>%
  group_by(name) %>%
  mutate(total = sum(value)) %>%
  filter(total != 0) %>%
  ungroup() %>%
  select(-total)


View(shared)
rand <- shared %>%
  uncount(value) %>%
  mutate(rand_name = sample(name)) %>%
  select(-name) %>%
  count(Group, rand_name)

View(shared)

shared_group_count <- shared %>%
  group_by(Group) %>%
  summarize(n = sum(value))

rand_group_count <- rand %>%
  group_by(Group) %>%
  summarize(n = sum(n))

inner_join(shared_group_count, rand_group_count,  by="Group")



shared_otu_count <- shared %>%
  group_by(name) %>%
  summarize(n = sum(value))

rand_otu_count <- rand %>%
  group_by(rand_name) %>%
  summarize(n = sum(n))

inner_join(shared_otu_count, rand_otu_count,  by=c("name" = "rand_name" ))


rand_group_count %>%
  ggplot(aes(x=n)) + geom_histogram()

View(rand)
rand_df <- rand %>%
  pivot_wider(names_from="rand_name", values_from="n", values_fill = 0) %>%
  as.data.frame()

rownames(rand_df) <- rand_df$Group
rand_df <- rand_df[, -1]

rand_matrix <- as.matrix(rand_df)

norare_dist_matrix <- vegdist(rand_matrix, method="bray")
rare_dist_matrix <- avgdist(rand_matrix, dmethod="bray", sample=1826)

norare_dist_tibble <- norare_dist_matrix %>%
  as.matrix() %>%
  as_tibble(rownames="sample") %>%
  pivot_longer(-sample) %>%
  filter(name < sample)

rare_dist_tibble <- rare_dist_matrix %>%
  as.matrix() %>%
  as_tibble(rownames="sample") %>%
  pivot_longer(-sample) %>%
  filter(name < sample)

comparison <- inner_join(norare_dist_tibble, rare_dist_tibble, by=c("sample", "name")) %>%
  select(sample, name, norare=value.x, rare=value.y) %>%
  inner_join(., rand_group_count, by=c("sample" = "Group")) %>%
  inner_join(., rand_group_count, by=c("name" = "Group")) %>%
  mutate(n_diff = abs(n.x-n.y)) %>%
  select(-n.x, -n.y)

comparison %>%
  ggplot(aes(x=norare, y=rare, color=n_diff)) +
  geom_point(size=0.25, alpha=0.25) +
  geom_smooth()

comparison %>%
  pivot_longer(cols=c("norare", "rare"), names_to="type", values_to="dist") %>%
  ggplot(aes(x=n_diff,  y=dist)) +
  geom_point(size=0.25, alpha=0.25) +
  facet_wrap(~type, nrow=2)