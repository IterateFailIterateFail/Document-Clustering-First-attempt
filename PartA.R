rm(list=ls())
dev.off()

library(tm) 
library(SnowballC)
library(topicmodels)
library(wordcloud)
library(cluster)
library(igraph) 
library(lsa)
library(stringr)
library(wordcloud)
library(tidyverse)
library(quanteda)



#Definite cropus loading and cleaning

source <- DirSource("./docs")

docs <- VCorpus(source)

#might help structuring?



docs <- tm_map(docs, removePunctuation)
docs <- tm_map(docs,content_transformer(tolower))
docs <- tm_map(docs, removeNumbers)


#segt usaful variables
k <- 6
docs_length <- length(docs)
graph_dir <- "graphs/"
group_colors <- c("orange","red","green3","blue","cyan","magenta","yellow","gray","orchid",
                  "lightslateblue","lightsalmon","green1","burlywood","darkkhaki","ghostwhite")
col_pallete <- palette(group_colors)
#Remove stopwords 
add_stop <- c("can","one","will","use","see","way","may","said","look","like","also","figur","ibi","issu",
              "ask","go","ie","us","go","plot")

docs <- tm_map(docs, removeWords, stopwords("english"))
docs <- tm_map(docs, stripWhitespace)

#make document term matrix after stemming, and limmiting word length to between 4-15
docs <- tm_map(docs,stemDocument)
control <- list(
  wordLengths=c(1,15),
  bounds = list(global = c(4,docs_length))
)
docs <- tm_map(docs, removeWords, add_stop)
docs <- tm_map(docs, stripWhitespace)
#regualr term frezuney matirx
dtm <- DocumentTermMatrix(docs,control = control)

freq <- colSums(as.matrix(dtm))

#create sort order (asc)
ord <- order(freq,decreasing=TRUE)

#inspect most frequently occurring terms
sing_df<-data.frame(term=names(freq),occurrences=freq)
ggplot(subset(sing_df, occurrences>150), aes(reorder(term,occurrences), occurrences))+
  geom_bar(stat="identity") + theme(axis.text.x=element_text(angle=90, hjust=1))+
  labs(title = "Most Frequent single term",
       x = "",
       y = "Occurrences")
unlink(paste(graph_dir,"sing_freq.png",sep = ""))
ggsave(paste(graph_dir,"sing_freq.png",sep = ""))



#bigram checks
BigramTokenizer <-  function(x) unlist(lapply(ngrams(words(x), 2), paste, collapse = " "), use.names = FALSE)
dtmbi <- DocumentTermMatrix(docs, control = list(tokenize = BigramTokenizer))
freqbi <- colSums(as.matrix(dtmbi))
ordbi <- order(freqbi,decreasing=TRUE)

#inspect most frequently occurring bigrams
bi_df<-data.frame(term=names(freqbi),occurrences=freqbi)
ggplot(subset(bi_df, occurrences>30), aes(reorder(term,occurrences), occurrences))+
  geom_bar(stat="identity") + theme(axis.text.x=element_text(angle=90, hjust=1))+
  labs(title = "Most Frequent Bigram",
       x = "",
       y = "Occurrences")
unlink(paste(graph_dir,"Bi_freq.png",sep = ""))
ggsave(paste(graph_dir,"Bi_freq.png",sep = ""))




#simple word cloud
set.seed(4234)
unlink(paste(graph_dir,"simple_word_clpud.png",sep = ""))
png(filename = paste(graph_dir,"simple_word_clpud.png",sep = ""), width = 700, height = 700)
wordcloud(names(freq),freq,max.words=30,colors=brewer.pal(6,"Dark2"))
dev.off()

#for cluseting get cosine similarity metrics for TDM-LSA and convert to matrixs. Distance is questionalbe as 

#Wierd atmept at LSA based clustering
tdm_matrix <- as.matrix(TermDocumentMatrix(docs,control = control))
tdm_matrix_lsa <- lw_tf(tdm_matrix) * gw_idf(tdm_matrix)

cosineSim <- function(x){
  as.dist(x%*%t(x)/(sqrt(rowSums(x^2) %*% t(rowSums(x^2)))))
}

lsaSpace <- lsa(tdm_matrix_lsa, dimcalc_share()) # create LSA space
cs <- cosineSim(lsaSpace$dk)
cd <- 1-cs
cd_matrix <- as.matrix(cd)

#m<-as.matrix(dtm)
# cs <- cosineSim(m)
# cd <- 1-cs
# cd_matrix <- as.matrix(cd)

grouping_df <- as.data.frame(cbind(colnames(tdm_matrix)))

set.seed(41244)
hier_cosine <- hclust(cd,method="ward.D")

#plot, use hang to ensure that labels fall below tree
unlink(paste(graph_dir,"Hierechy.png",sep = ""))
png(filename = paste(graph_dir,"Hierechy.png",sep = ""), width = 700, height = 700)
plot(hier_cosine, hang=-1,xlab = "")
rect.hclust(hier_cosine,k)
dev.off()
#cut into k subtrees.
hclusters_cosine <- cutree(hier_cosine,k)

#make clusplot
unlink(paste(graph_dir,"clusplot_h.png",sep = ""))
png(filename = paste(graph_dir,"clusplot_h.png",sep = ""), width = 700, height = 700)
clusplot(cd_matrix,hclusters_cosine, color=T, shade=T, labels=3, lines=0)
dev.off()

unlink(paste(graph_dir,"shil_h.png",sep = ""))
png(filename = paste(graph_dir,"shil_h.png",sep = ""), width = 700, height = 700)
plot(silhouette(hclusters_cosine,cd_matrix))
dev.off()

#bind to dataframe
grouping_df <- cbind(grouping_df,hclusters_cosine)
names(grouping_df)[ncol(grouping_df)] <- "H_cluster"
#makes it easier to sort
grouping_df$H_cluster <-  sprintf("%02d",as.integer(grouping_df$H_cluster,"","")) 



set.seed(941)
#k means clusering!
kfit_cs <- kmeans(cd, k, nstart=docs_length)
kfit_clusters <- kfit_cs$cluster

#cluster plot of kmeans
unlink(paste(graph_dir,"clusplot_k.png",sep = ""))
png(filename = paste(graph_dir,"clusplot_k.png",sep = ""), width = 700, height = 700)
clusplot(cd_matrix, kfit_clusters, color=T, shade=T, labels=3, lines=0)
dev.off()

#shiloute plot
unlink(paste(graph_dir,"shil_k.png",sep = ""))
png(filename = paste(graph_dir,"shil_k.png",sep = ""), width = 700, height = 700)
plot(silhouette(kfit_clusters,cd_matrix))
dev.off()

#bind to data frame
grouping_df <- cbind(grouping_df,kfit_clusters)
names(grouping_df)[ncol(grouping_df)] <- "K_cluster"
grouping_df$K_cluster <-  sprintf("%02d",as.integer(grouping_df$K_cluster,"","")) 

#WSS diffrence in clusters
wss <- 2:(length(docs)-1)
for (i in 2:(length(docs)-1)) wss[i] <- sum(kmeans(cd,centers=i,nstart=25)$withinss)
unlink(paste(graph_dir,"WSS.png",sep = ""))
png(filename = paste(graph_dir,"WSS.png",sep = ""), width = 700, height = 700)
plot(2:(length(docs)-1), wss[2:(length(docs)-1)], type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares") 
dev.off()

#graph netwrok
#number chosen is arbitary
cs[cs < max(cs)/(3)] <- 0
cs <- round(cs,3)
g <- graph.adjacency(as.matrix(cs), weighted=T, mode = "undirected")

# set name to just neumebr
V(g)$label <- str_extract(V(g)$name,"[0-9]+")

#set asthetics choices
V(g)$label.color <- "black"
E(g)$color <- "grey"
E(g)$width <- E(g)$weight*7

#find communitie
comm_fg <- fastgreedy.community(g)
V(g)$color <- kfit_clusters
# 
#plot(g, layout=layout_nicely,palette = col_pallete)
grouping_df <- cbind(grouping_df,comm_fg$membership)
names(grouping_df)[ncol(grouping_df)] <- "Com_member"
grouping_df$Com_member <-  sprintf("%02d",as.integer(grouping_df$Com_member,"","")) 

unlink(paste(graph_dir,"shil_g.png",sep = ""))
png(filename = paste(graph_dir,"shil_g.png",sep = ""), width = 700, height = 700)
plot(silhouette(comm_fg$membership,cd_matrix))
dev.off()


#topic modelling

burnin <- 1000
iter <- 2000

thin <- 500

nstart <- 5

seed <- list(2303,321,323,1231,731)
best <- TRUE


ldaOut <- LDA(dtm,k, method="Gibbs", control=
                list(nstart=nstart, seed = seed, best=best, burnin = burnin, iter = iter, thin=thin))
ldaOut_topics <-as.matrix(topics(ldaOut))

topics <- terms(ldaOut,10)
#bind to df
grouping_df <- cbind(grouping_df,ldaOut_topics)
names(grouping_df)[ncol(grouping_df)] <- "LDA_topics"
grouping_df$LDA_topics <-  sprintf("%02d",as.integer(grouping_df$LDA_topics,"","")) 

unlink(paste(graph_dir,"shil_l.png",sep = ""))
png(filename = paste(graph_dir,"shil_l.png",sep = ""), width = 700, height = 700)
plot(silhouette(ldaOut_topics,cd_matrix))
dev.off()


rownames(grouping_df) <- NULL
#save topics
write.csv(topics,file=paste("topic_terms.csv"),row.names = F)

#set corpus for futreu use
corpus<-corpus(docs)
docvars(corpus)<-cbind(docvars(corpus),grouping_df)

#this is useful for legends
number_text <- NULL

for(i in 1:k)number_text <- c(number_text,as.character(i))

groupings <- names(grouping_df)[2:length(names(grouping_df))]

#for each grouping
for(group in groupings){
  #get 10 most common terms
  sing_weight <-tokens(corpus) %>% 
    dfm() %>% 
    dfm_select(max_nchar = 20) %>% 
    textstat_frequency(n= 10,groups = group)
  #and get 10 most frequent bigrams
  bi_weight <-tokens(corpus) %>% 
    tokens_ngrams(n = 2) %>% 
    dfm() %>% 
    dfm_select(max_nchar = 20) %>% 
    textstat_frequency(n= 10,groups = group)
  
  
  #graohing regular terms
  plot<-ggplot(data = sing_weight, aes(x = nrow(sing_weight):1, y = frequency,fill = group)) +
    geom_col() +
    facet_wrap(~ group, scales = "free") +
    coord_flip() +
    scale_x_continuous(breaks = nrow(sing_weight):1,
                       labels = sing_weight$feature) +
    labs(x = NULL, y = "Relative frequency",title = paste(group,"Single terms"))+
    scale_fill_manual(values = group_colors,name = "")
  
  print(plot)
  unlink(paste(graph_dir,group,"_sing_freq.png",sep = ""))
  ggsave(paste(graph_dir,group,"_sing_freq.png",sep = ""),width = 7,height = 7,limitsize = F)
  
  
  #Graphing Bigrams
  plot<-ggplot(data = bi_weight, aes(x = nrow(bi_weight):1, y = frequency,fill = group)) +
    geom_col() +
    facet_wrap(~ group, scales = "free") +
    coord_flip() +
    scale_x_continuous(breaks = nrow(bi_weight):1,
                       labels = bi_weight$feature) +
    labs(x = NULL, y = "Relative frequency",title = paste(group,"Bi-grams"))+
    scale_fill_manual(values = group_colors,name = "")
  
  print(plot)
  unlink(paste(graph_dir,group,"_bi_freq.png",sep = ""))
  ggsave(paste(graph_dir,group,"_bi_freq.png",sep = ""),width = 12,height = 10,limitsize = F)
  
  #grpah network graphs
  unlink(paste(graph_dir,group,"_graph.png",sep = ""))
  png(filename = paste(graph_dir,group,"_graph.png",sep = ""), width = 800, height = 800)
  
  set.seed(13131)
  layout1 <- layout_nicely(g)
  V(g)$color <- as.numeric(unlist(grouping_df[group]))
  plot(g, layout=layout1,palette = col_pallete,main = group)
  
  legend("bottomleft", inset=.02, title="Groups",
         legend = number_text, fill=group_colors, horiz=F, cex=0.8)
  dev.off()
  
}

print_table <- grouping_df %>% 
  select(-c(Com_member)) %>% 
  rename(Doc = V1)

write.csv(print_table,file=paste("H_cluster.csv"),row.names = F)
