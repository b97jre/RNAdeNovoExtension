library(ggplot2)
datadir = "/Users/johanreimegard/Vetenskap/Data/deNovoExtension/paper2"

readDistFile = "Figure2.data.csv"


readDistribution = read.csv(paste(datadir,readDistFile, sep = "/")) 

readDistribution$UniqueFraction = readDistribution$Unique/readDistribution$Total 
readDistribution$NonUniqueFraction = readDistribution$Not_Unique/readDistribution$Total 
readDistribution$properlyPaired = (readDistribution$Unique+readDistribution$Not_Unique)/readDistribution$Total 
readDistribution$x = paste(readDistribution$Species,readDistribution$Assembly,readDistribution$Extended, sep ="+")
readDistribution$xLabel = paste(readDistribution$Species,readDistribution$Assembly, sep ="+")

Before = readDistribution[readDistribution$Extended=="Before", ]
After = readDistribution[readDistribution$Extended=="After", ]
After$diff = After$properlyPaired-Before$properlyPaired
After$x = paste(After$Species,After$Assembly, sep ="+")
Diff = After

ggplot(After, aes(x = x, y = diff, fill=Extended))+ geom_bar(stat ="identity")+
theme(axis.text.x = element_text(angle = 90, hjust = 1))+
scale_y_continuous(name="Difference in properly paired reads")
ggsave(paste(datadir, "ReadDistributionDiff.pdf", sep = "/"))


#Figure 2

ggplot(readDistribution, aes(x = x, y = properlyPaired, fill=Extended))+
 geom_bar(stat ="identity", guide = guide_legend(reverse=TRUE),labels=readDistribution$xLabel) +
  scale_x_discrete(breaks=readDistribution$x, labels=readDistribution$xLabel)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(paste(datadir, "ReadDistribution.pdf", sep = "/"))

#



