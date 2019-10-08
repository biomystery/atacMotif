a <- readRDS('./db/tfclass.rds')
a$class <- rbind(a$class,data.frame(stringsAsFactors = F,id="2.7",name='C2HC',about="The known transcription factors exhibiting a DNA-binding domain of intertwined zinc fingers belong to the DM group of proteins, named after Drosophila Doublesex and C. elegans Meb-3. Intertwined are one zinc finger of C2HC and one of HCC2 type, each binding one zinc ion. Together, they follow the consensus pattern C--C--h--------H---C----c-c-c, where upper and lower case letter refer to either zinc finger, and dashes symbolize anyother amino acid. The most C-terminal cysteine is part of a short alpha-helical structure, which extends into a carboxy-terminal tail that stably folds into an alpha helix upon binding to DNA. This recognition helix binds in the minor groove of the DNA. In contrast to other minor-groove binders, it does not induce a sharp bending of the DNA. (PMID 16382145)"))
saveRDS(a,'./db/tfclass.rds')

# make tree ---------------------------------------------------------------
# from igraph
# family => class -> superclass 


require(igraph)
require(tidyverse)
nodes <- a$family%>%mutate(class='family')%>%
  rbind(a$class%>%select(-about)%>%mutate(class='class'))%>%
#  rbind(data.frame(id='2.7',name='C2HC',class='class',stringsAsFactors = F))%>%
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

require(ggraph)
require(tidygraph)
g <- graph_from_data_frame(edges,vertices = nodes)%>% as_tbl_graph()


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
