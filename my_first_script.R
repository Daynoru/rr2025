setwd("D:/uni/RR3/")

library(openxlsx)

test_table <- read.xlsx(xlsxFile = "Test_table.xlsx")

str(test_table)

unique(test_table$Species)

hist(test_table$PO.activity)
hist(test_table$GST.activity)
hist(test_table$CAT.activity)

mean(test_table$PO.activity)
median(test_table$PO.activity)
mean(test_table$GST.activity)
median(test_table$GST.activity)
mean(test_table$CAT.activity, na.rm=TRUE)
median(test_table$CAT.activity, na.rm=TRUE)

shapiro.test(test_table$PO.activity)
shapiro.test(test_table[test_table$
                          Species=="E. coli",]$PO.activity)
shapiro.test(test_table[test_table$
                          Species=="E. verrucosus",]$PO.activity)

wilcox.test(test_table$PO.activity ~ 
              test_table$Species)

#RR5
setwd("D:\\uni\\RR5_data")
library(ggplot2)
library(openxlsx)

tbl <- read.xlsx("Test_table_2.xlsx", sheet = 2)

str(tbl)

hist(tbl$PO.activity)
hist(tbl$Hemocyte.count)

unique(tbl$Species)
library(ggpubr)

ggplot(data=tbl,
       aes(x=Group, y=GST.activity)) +
 expand_limits(y=0) + 
  geom_boxplot()

ggplot(data=tbl, 
       aes(x=Group, y=PO.activity)) +
  geom_jitter(width = 0.1)+
  geom_boxplot(outliers = FALSE) 

library(scales)

ggplot(data=tbl,
       aes(x=Group, y=CAT.activity)) +
  expand_limits(y=0) +
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  scale_y_continuous(
    labels=comma_format(decimal.mark=","))

ggplot(data=tbl, aes(x=Group, y=CAT.activity)) +
  expand_limits(y=0) +
  geom_boxplot() + 
  facet_wrap(~Species)


ggplot(data=tbl, aes(x=Group, y=GST.activity, fill=Species)) +
  expand_limits(y=0) + #y=0 включаем
  geom_boxplot(show.legend = FALSE) + #не показывать легенду к цвету
  scale_fill_manual(values=c("#D2AA6D", "forestgreen")) + #указываем цвета (RGB или имя)
  facet_wrap(~Species)  #панели по видам

ggplot(data=tbl, aes(x=Group, y=GST.activity, fill=Species)) +
  expand_limits(y=0) + #y=0 включаем
  geom_boxplot(show.legend = FALSE) + #не показывать легенду к цвету
  scale_fill_manual(values=c("#D2AA6D", "forestgreen")) + #указываем цвета (RGB или имя)
  facet_wrap(~Species) +  #панели по видам
  ylab("GST activity, a.u.") + 
  xlab("") + #название оси Y
  theme_bw(base_size = 16) + #увеличим размер шрифта + белый фон
  theme(strip.text = element_text(face="italic")) #курсив

plot.PO <- ggplot(data=tbl, aes(x=Group, y=PO.activity)) +
  expand_limits(y=0) + #y=0 включаем
  geom_boxplot(show.legend = FALSE) + #боксплоты
  facet_wrap(~Species) #панели по видам
plot.PO1 <- plot.PO + 
  ylab("PO activity, a.u.") + xlab("") + #название оси Y
  theme_bw(base_size = 16) + #увеличим размер шрифта + белый фон
  theme(strip.text = element_text(face="italic")) #курсив
plot.PO2 <- plot.PO1 +
  aes(fill=Species) + # добавляем заливку по виду
  scale_fill_manual(values=c("#D2AA6D", "forestgreen")) #указываем цвета (RGB или имя)
plot.PO2 #чтобы отобразить

plot.PO1 + geom_pwc(method = "wilcox_test", label="p.adj")
ggsave("PO_with_stats.png", device=png, width=20, height=12, units="cm")

ggsave(filename="GST.png", device=png, width=16, height=12, units="cm", dpi=300)

yeast.growth <- read.xlsx("253_2023_12863_MOESM2_ESM.xlsx", sheet = 3)
ggplot(yeast.growth, aes(x=Condition, y=RGR2)) + #данные для графика
  geom_boxplot(color = "forestgreen", fill = "pink") + #рисуем боксплот
  geom_pwc(method="wilcox_test", label = " {p.adj.signif}", 
            p.adjust.method = "holm", ref.group = 1) +
geom_jitter(width = 0.1, colour = "forestgreen") 
    


ggsave("results1.png", device=png, width=16, height=12, units="cm")
