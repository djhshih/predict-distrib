# Assess how well marginal methylation distributions can predict methylation
# classes in ovarian cancer

library(io)
library(glmnet)

pheno0 <- qread("data/ccoc-ts_sample-info_stage3.tsv");
params <- qread("out/ccoc-ts_methy_beta_pmfmix_params.rds");

out.fn <- filename("ccoc-ts", path="out", tag=c("methy", "beta"));
pdf.fn <- insert(out.fn, ext="pdf");

W <- params$W;
rownames(W) <- sub("-02", "", rownames(W));
rownames(W) <- sub("R$", "", rownames(W));

idx <- match(rownames(W), pheno0$sample_id);
pheno <- pheno0[idx, ];
pheno$cluster[pheno$cluster == "mixture"] <- "hypomethylated";
pheno$cluster[pheno$sample_type == "normal"] <- "normal";

table(pheno$cluster, useNA="always")

pheno[is.na(pheno$cluster), ]

idx <- !is.na(pheno$cluster);
pheno <- pheno[idx, ];
W <- W[idx, ];
y <- factor(pheno$cluster);


dim(W)
dim(pheno)
table(y)

stopifnot(pheno$sample_id == rownames(W))

fit <- cv.glmnet(W, y, family="multinomial");
fit


scores <- predict(fit, W);
scores <- scores[, , 1];

yhat <- colnames(scores)[apply(scores, 1, which.max)];

cmat <- table(y, yhat);

qwrite(cmat, insert(out.fn, tag=c("elnet", "confusion-mat"), ext="mtx"))

fisher.test(cmat, simulate=TRUE)

confusion_matrix_subset <- function(scores, y, classes) {
	idx <- y %in% classes;
	y.cc <- factor(y[idx], levels=classes);
	yhat.cc <- factor(classes[apply(scores[idx, classes], 1, which.max)], levels=classes);
	table(y.cc, yhat.cc)
}

refs <- c("clear_cell", "high_grade_serous", "mucinous", "endometrioid", "normal");
names(refs) <- refs;
query <- c("hypomethylated");

cmats <- lapply(refs,
	function(ref) {
		confusion_matrix_subset(scores, y, c(ref, query))
	}
);

lapply(cmats, prop.table, margin=1)

# ---

library(randomForest)

fit.rf <- randomForest(W, y);
# high estimated out-of-bag error rate!
fit.rf

# model memorizes training data?
yhat <- predict(fit.rf, W);
table(y, yhat)

# ----

library(mmalign)
pca_plot(t(W), pheno, aes(color=cluster))
pca_plot(t(W), pheno, aes(color=cluster), dims=3:4)

mu <- unlist(lapply(params$theta, function(th) th$mu));
lambda <- unlist(lapply(params$theta, function(th) th$lambda));

v <- seq(0.005, 0.995, by=0.01);


lF <- exp(params$lF);
dm <- melt(lF, varnames=c("component", "bin"));
dm <- left_join(dm, data.frame(bin=1:ncol(lF), x=v));
dm$component <- factor(dm$component);

qdraw(
	ggplot(dm, aes(x=x, y=value, colour=component)) +
		theme_classic() +
		geom_line() +
		xlab("methylation beta value") + ylab("density") +
		theme(legend.position = "none")
	,
	file = insert(pdf.fn, "mixture")
)


library(dplyr)
library(reshape2)

pmfs <- W %*% exp(params$lF);
colnames(pmfs) <- 1:ncol(pmfs);

d <- melt(pmfs, varnames=c("sample_id", "bin"));
d <- left_join(d, data.frame(bin=1:ncol(pmfs), x=v));
d <- left_join(d, select(pheno, sample_id, cluster));

qdraw(
	ggplot(d, aes(x=x, y=value, colour=cluster, group=sample_id)) +
		theme_classic() +
		geom_line(alpha=0.3) +
		facet_wrap(~ cluster, nrow=1) +
		scale_x_continuous(breaks = c(0, 0.5, 1)) +
		xlab("methylation beta value") + ylab("density") +
		theme(legend.position = "none")
	,
	width = 10,
	file = insert(pdf.fn, "density")
)

