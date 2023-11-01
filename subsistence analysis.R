
#subsistence analysis
library(survival)
library(survminer)
library(ggplot2)
library(dplyr)
## DFS cutoff40

table<-read.table("C://Users//11042//Desktop//b.csv",header=TRUE,sep=',')
table2<-read.table("C://Users//11042//Desktop//f.csv",header=TRUE,sep=',')
table<-table %>% select(renameV1:order) %>% filter(renameV1 %in% table2[,1])

#table<-table %>% select(renameV1:order) %>% filter(Risk.22.2Group=="low")
fit <- survfit(Surv(DFS.M.used, Event_DFS.used) ~ Risk.22.2Group, data = table)
grid.draw.ggsurvplot <- function(x){survminer:::print.ggsurvplot(x, newpage = FALSE)}
res<-ggsurvplot(fit, 
                data = table,
                pval = TRUE, 
                risk.table = TRUE,
                tables.height=0.3,
                xlab = "Time in months",
                legend.title="Risk",
                xlim = c(0, 58),
                legend.labs = c("High", "Low"),
                break.time.by = 10,
                font.title = 16,
                font.x =  14,
                font.y = 14,
                font.tickslab = 12,
                legend = c(0.1, 0.5))
ggsave("survival_DFS.png", plot = res, dpi=300, width = 10, height = 7, units = "in")

b<-ggsurvplot(fit, 
              data = table,
              pval = TRUE, 
              risk.table = TRUE,
              tables.height=0.3,
              xlab = "Time in months",
              legend.title="Chemo",
              xlim = c(59.9, 65),
              legend.labs = c("No", "Yes"),
              break.x.by = 60,
              break.y.by = 0.05,
              font.title = 16,
              font.x =  14,
              font.y = 14,
              font.tickslab = 12,
              legend = c(0.1, 0.5),
              risk.table.col = "strata", # Change risk table color by groups
              ggtheme = theme_bw())
ggsave("b.png", plot = b, dpi=300, width = 20, height = 10, units = "in")

a<-ggsurvplot(fit, 
              data = table,
              pval = TRUE, 
              risk.table = TRUE,
              tables.height=0.3,
              xlab = "Time in months",
              legend.title="Chemo",
              xlim = c(82, 89),
              legend.labs = c("No", "Yes"),
              break.x.by = 82,
              break.y.by = 0.05,
              font.title = 16,
              font.x =  14,
              font.y = 14,
              font.tickslab = 12,
              legend = c(0.1, 0.5),
              risk.table.col = "strata", # Change risk table color by groups
              ggtheme = theme_bw())
ggsave("a.png", plot = a, dpi=300, width = 20, height = 10, units = "in")

yes<-table %>% select(renameV1:order) %>% filter(Risk.22.2Group=="high")
result.quantile.DFS.M. <- quantile(yes$DFS.M.)
round(result.quantile.DFS.M.,0)
no<-table %>% select(renameV1:order) %>% filter(Risk.22.2Group=="low")
result.quantile.DFS.M. <- quantile(no$DFS.M.)
round(result.quantile.DFS.M.,0)
result.quantile.DFS.M. <- quantile(table$DFS.M.)
round(result.quantile.DFS.M.,0)

res.cox <- coxph(Surv(DFS.M.used, Event_DFS.used) ~ Risk.22.2Group, data = table)
c<-summary(res.cox)
d<-c$conf.int
round(d[,2],4) 
round(1/d[,3],4)
round(1/d[,4],4)

summary(res.cox)

## DMFS cutoff40


#table2<-read.csv("1059_sample_DMFS/over12/over12_cutoff40/1059_samples_date_over12_cutoff45.csv",header=TRUE,sep=',')
#table2<-table2 %>% select(NO.ori:age.1) %>% filter(Chemo=="1" | Chemo=="0")
fit_new <- survfit(Surv(DMFS.M., Event_DMFS) ~ Risk.22.2Group, data = table)
grid.draw.ggsurvplot <- function(x){survminer:::print.ggsurvplot(x, newpage = FALSE)}
res2<-ggsurvplot(fit, 
                 data = table,
                 pval = TRUE, 
                 risk.table = TRUE,
                 tables.height=0.3,
                 xlab = "Time in months",
                 legend.title="Risk",
                 xlim = c(0, 58),
                 legend.labs = c("High", "Low"),
                 break.time.by = 10,
                 font.title = 16,
                 font.x =  14,
                 font.y = 14,
                 font.tickslab = 12,
                 legend = c(0.1, 0.5))
ggsave("survival_DMFS.png", plot = res2, dpi=300, width = 10, height = 7, units = "in")

b<-ggsurvplot(fit_new, 
              data = table,
              pval = TRUE, 
              risk.table = TRUE,
              tables.height=0.3,
              xlab = "Time in months",
              legend.title="Risk",
              xlim = c(59.9, 65),
              legend.labs = c("High", "Low"),
              break.x.by = 60,
              break.y.by = 0.05,
              font.title = 16,
              font.x =  14,
              font.y = 14,
              font.tickslab = 12,
              legend = c(0.1, 0.5),
              risk.table.col = "strata", # Change risk table color by groups
              ggtheme = theme_bw())
ggsave("b.png", plot = b, dpi=300, width = 20, height = 10, units = "in")

a<-ggsurvplot(fit_new, 
              data = table2,
              pval = TRUE, 
              risk.table = TRUE,
              tables.height=0.3,
              xlab = "Time in months",
              legend.title="Chemo",
              xlim = c(80, 89),
              legend.labs = c("No", "Yes"),
              break.x.by = 82,
              break.y.by = 0.05,
              font.title = 16,
              font.x =  14,
              font.y = 14,
              font.tickslab = 12,
              legend = c(0.1, 0.5),
              risk.table.col = "strata", # Change risk table color by groups
              ggtheme = theme_bw())
ggsave("a.png", plot = a, dpi=300, width = 20, height = 10, units = "in")

yes<-table %>% select(renameV1:order) %>% filter(Risk.22.2Group=="high")
result.quantile.DMFS.M. <- quantile(yes$DMFS.M.)
round(result.quantile.DMFS.M.,0)
no<-table %>% select(renameV1:order) %>% filter(Risk.22.2Group== "low")
result.quantile.DMFS.M. <- quantile(no$DMFS.M.)
round(result.quantile.DMFS.M.,0)
result.quantile.DMFS.M. <- quantile(table$DMFS.M.)
round(result.quantile.DMFS.M.,0)

res.cox <- coxph(Surv(DMFS.M., Event_DMFS) ~ Risk.22.2Group, data = table)
c<-summary(res.cox)
d<-c$conf.int
round(d[,2],4) 
round(1/d[,3],4)
round(1/d[,4],4)

summary(res.cox)
