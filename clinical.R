setwd("E:/work/Project/wes")
library(data.table)
library(ggplot2)

pe_data <- fread("normal.csv", stringsAsFactors = FALSE)
pe_data <- pe_data[,!c("签收时间", "申请单号")]
pe_data <- dcast(pe_data[, 1:5], 姓名 + 性别 + 年龄 ~ 项目名称, value.var = "结果", fun.aggregate = list)
pe_data <- pe_data[!sapply(`APO-B/APO-A1`, function(x)identical(x, character(0))),]#remove character(0)
ca <- pe_data[!sapply(糖类抗原CA125, function(x)identical(x, character(0))),]


#bind case and sample normal####
norm <- fread("norm_col.csv", stringsAsFactors = FALSE)
patient <- fread("patient.csv", stringsAsFactors = FALSE)
patient_ca <- fread("patient_CA125.csv", stringsAsFactors = FALSE)
patient <- patient[patient_ca, on = "样本编号"]

set.seed(2)
norm_sample <- norm[,.SD[sample(.N, 40)],by = 性别]#smaple 40 people

col_name <- fread("col_name.csv", stringsAsFactors = FALSE)
co <- colnames(col_name)

mix <- rbind(patient[,..co], norm_sample[,..co])
mix$case <- as.integer(mix$样本编号 > 150)#mark case and control

#fwrite(mix, "mix.csv")

#HBV
small3 <- norm[`乙肝e抗体（HBe-Ab）` == "阳性"&`乙肝核心抗体（HBc-Ab）` == "阳性"&`乙肝表面抗原（HBs-Ag)` == "阳性",]$样本编号
big3 <- norm[`乙肝e抗原(HBe-Ag）` == "阳性"&`乙肝核心抗体（HBc-Ab）` == "阳性"&`乙肝表面抗原（HBs-Ag)` == "阳性",]


#statistics####
mix <- read.csv("mix.csv", stringsAsFactors = FALSE)
mix[, 6:35] <- lapply(mix[, 6:35], as.numeric)

#normality test
shap <- apply(mix[, 6:35], 2, shapiro.test)
normality <- setDT(shap)[2]#keep p-value #shap class is transfered to data.table
normal_distribution <- c("LY", "MO", "RBC", "HGB", "PLT", "CREA")

#split data by distribution
gauss <- mix[normal_distribution]
gauss$case <- mix$case 
non <- mix[!names(mix) %in% normal_distribution]

#t-test
stat_t <- data.frame(unlist(t.test(LY ~ case, data = gauss, var.equal = TRUE)))
for (n in c(2:6)){
    t <- t.test(gauss[, n] ~ gauss$case, var.equal = TRUE)
    t <- data.frame(unlist(t))
    stat_t <- cbind(stat_t, t)
}

colnames(stat_t) <- colnames(gauss[, 1:6])
stat_t <- data.frame(t(stat_t))

write.csv(stat_t, "t_test_result.csv", row.names = TRUE)

#U-test
stat_u <- data.frame(unlist(wilcox.test(WBC ~ case, data = non, exact=FALSE)))
for (n in c(7:29)){
    u <- wilcox.test(non[, n] ~ non$case, exact=FALSE)
    u <- data.frame(unlist(u))
    stat_u <- cbind(stat_u, u)
}

colnames(stat_u) <- colnames(non[, 6:29])
stat_u <- data.frame(t(stat_u))

write.csv(stat_u, "U_test_result.csv", row.names = TRUE)


v <- cbind(stat_t$p.value, stat_u$p.value)

#output summary table####

#get mean +se
mixp <- setDT(mix)
mixp <- melt(mixp, id.vars = c("样本编号", "姓名", "性别", "年龄", "case"))
mixp <- mixp[!is.na(value),]
mean_case <- mixp[variable %in% normal_distribution, .(mean = round(mean(value), digits = 2), 
                                                       se = round(sd(value)/sqrt(.N), digits = 2)),
                  by = c("case", "variable")]

mean_case$count <- paste(mean_case$mean, "±", mean_case$se)
mean_output <- dcast(mean_case[, c(1, 2, 5)], variable ~ case, value.var = "count")
fwrite(mean_output, "mean_output.csv")

#get IQR
iqr <- mixp[!variable %in% normal_distribution, .(Q1 = round(quantile(value, 0.25), digits = 2),
                                                  Q2 = round(quantile(value, 0.50), digits = 2),
                                                  Q3 = round(quantile(value, 0.75), digits = 2)),
            by = .(case, variable)]

iqr$count <- paste0(iqr$Q2, " (", iqr$Q1, ", ", iqr$Q3, ")")
iqr_output <- dcast(iqr[, c(1, 2, 6)], variable ~ case, value.var = "count")

#bind
all_output <- rbind(mean_output, iqr_output)
all_output <- all_output[match(unique(mixp$variable), variable),]
fwrite(all_output, "all_output.csv")
