a <- readRDS('./db/tfclass.rds')


# make tree ---------------------------------------------------------------
# from igraph
# family => class -> superclass 


require(igraph)
require(tidyverse)
nodes <- a$family%>%mutate(class='family')%>%rbind(a$class%>%select(-about)%>%mutate(class='class'))%>%
  rbind(a$superclass%>%select(-about) %>%mutate(class='superclass'))%>%
  rbind(data.frame(id='-',name='root',class='root'))
edges <- data.frame(from='-',to=a$superclass%>%pull(id)%>%unique,stringsAsFactors = F)%>%
  rbind(do.call(rbind,lapply( a$superclass$id%>%unique, function(i) data.frame(from=i,to=grep(paste0('^',i),a$class$id,value = T)))))%>%
  rbind(do.call(rbind,lapply(a$class$id%>%unique, function(i) data.frame(from=i,to=grep(paste0('^',i),a$family$id,value = T)))))

g <- graph_from_data_frame(edges,vertices = nodes)%>% as_tbl_graph()

require(ggraph)
require(tidygraph)



ggraph(g,'dendrogram',circular=T)+
  geom_edge_diagonal(colour="grey") +
  geom_node_point(aes(color=class))+
  geom_node_text(aes(label=name),size=3)+
  coord_fixed()+theme_void() 

ggsave('tfclass.png')
