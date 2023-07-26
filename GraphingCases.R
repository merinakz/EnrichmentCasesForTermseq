rm(list = ls())
library(ggplot2)

data <- read.csv("~/Desktop/Chlcase.csv")
data$Value <- as.numeric(data$Value)
data$Strand <- as.character(data$Strand)
data$Case <- as.character(data$Case)
data$Value <-  round(data$Value, 1)
ggplot(data,  aes(fill=Case, y=Value, x=Strand)) + 
  geom_bar(position="stack", stat="identity") +
  geom_text(size = 4, position = position_stack(vjust = 0.5), colour = "black" , aes(label = Value)) +
  scale_fill_manual(values = c("#66CCFF", "#FECC66", "#CCCCCC", "#FC6766", "#66FFCC")) +
  ggtitle("Chloramphenicol",) + 
  ylab("Percentage") + 
  xlab("Strand") +
  theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=25), axis.title.y = element_text(size=25),
        plot.title = element_text(size = 40))+
  coord_flip()
