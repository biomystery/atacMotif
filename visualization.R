a <- readRDS('./db/tfclass.rds')


# make tree ---------------------------------------------------------------
# from igraph
# family => class -> superclass 


require(igraph)
require(tidyverse)
nodes <- a$family%>%mutate(class='family')%>%
  rbind(a$class%>%select(-about)%>%mutate(class='class'))%>%
  rbind(data.frame(id='2.7',name='C2HC',class='class',stringsAsFactors = F))%>%
  rbind(a$superclass%>%select(-about) %>%mutate(class='superclass'))%>%
  rbind(data.frame(id='-',name='root',class='root'))
edges <- data.frame(from='-',to=nodes%>%filter(class=="superclass")%>%pull(id),stringsAsFactors = F)%>%
  rbind(do.call(rbind,lapply( nodes%>%filter(class=="superclass")%>%pull(id),
                              function(i) 
                                data.frame(from=i,to=grep(paste0('^',i,'.'),
                                                          nodes%>%filter(class=="class")%>%pull(id),
                                                          value = T)))))%>%
  rbind(do.call(rbind,lapply(nodes%>%filter(class=="class")%>%pull(id), 
                             function(i) data.frame(from=i,to=grep(paste0('^',i,'.'),
                                                                   nodes%>%filter(class=="family")%>%pull(id),value = T)))))

g <- graph_from_data_frame(edges,vertices = nodes)%>% as_tbl_graph()

require(ggraph)
require(tidygraph)



ggraph(g,'dendrogram',circular=T)+
  geom_edge_diagonal(colour="grey") +
  geom_node_point(aes(color=class))+
  geom_node_text(aes(label=name),size=3)+
  coord_fixed()+theme_void() 

ggsave('tfclass_circle.png')

ggraph(g,'dendrogram',circular=F)+
  geom_edge_diagonal(colour="grey") +
  geom_node_point(aes(color=class))+
  geom_node_text(aes(label=name,y=y-.02),size=2.5,angle=45,hjust=1)+
  theme_void() +
  expand_limits(y=c(-1,3))
ggsave('tfclass.png',width = 16,height =6)
