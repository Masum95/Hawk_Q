library(foreach)


library(parallel)
library(iterators)



library(doParallel)


cl<-makeCluster(1,outfile="") #outfile="" necessary for printing inside foreach loop
registerDoParallel(cl)

print('here');
con = file("out_unique_top.kmerDiff", "r");
totalKmers <- as.numeric(system2('wc', stdout = TRUE,  args=c('-l','<','out_unique_top.kmerDiff')));
print(totalKmers);
Z <- read.table("pcs.evec",header=FALSE);

input <- read.table("gwas_eigenstratX.ind");

Y <- matrix(nrow=nrow(Z),ncol=1); # row = sample_size, 
cov <- matrix(nrow=nrow(Z),ncol=1); #row = sample size 

totals <- read.table("gwas_eigenstratX.total"); # row = sample size

counts<- vector(length=length(Y));


for(i in 1:length(Y))
{

	Y[i,1]=input[i,3]; #phenotype value 
}




# totalKmers <- read.table("total_kmers.txt",header=FALSE);
# totalKmers <- length(readLines(con));
n<-nrow(Z);

cat("", file = "pvals_top.txt")

CHUNK_SIZE=2000;

ptm <- proc.time()

cnt = 0;
while(TRUE)
{	
	if( (totalKmers - cnt) <CHUNK_SIZE )
	{
		CHUNK_SIZE = totalKmers - cnt
	}

	kmercounts <- read.table(con,nrow=CHUNK_SIZE); # significant kmers 
	cnt = (cnt + CHUNK_SIZE);
	nr=nrow(kmercounts); 

	print(nr);
	ls<-foreach(j=icount(nr), .combine=cbind) %dopar%  #loop over significant kmers  
	{



		for(i in 1:length(Y))
		{
			counts[i]=kmercounts[j,2+i]/totals[i,1]; #normalized counts for jth kmer of each sample
			#print(kmercounts[j,2+i]);
		}

		#	model1<-glm(formula = Y ~ counts, family = binomial(link = "logit"));
	

		model2<-lm(formula = Y ~ Z[,1]+Z[,2]+totals[,1]+counts); #Y -> Quant Value, Z[,1] -> components along PC1 , Z[,2] -> components 		along PC2, totals[,1] -> total kmers of individual, counts[i] -> ( normalized kmer count in j'th kmer 
		

		#summary(model1);

		#	v1<-anova(model1, test="Chisq");

		summary(model2);
		# print('lm');
		#  print(model2);
		# print('anova value');
		v2<-anova(model2, test="Chisq");
		# print(v2);
		#	rbind(v1$'Pr(>Chi)'[2],v2$'Pr(>Chi)'[13]);

		v2$'Pr(>F)'[4];
		# print(v2$'Pr(>F)'[4]);
	}
	print(cnt);
	print(totalKmers);
	print(cnt==totalKmers);
	write(ls,file='pvals_top.txt',ncolumns=1,append=TRUE,sep='\t');

	if( (totalKmers - cnt) == 0)
	{
		print(nr);
		break;
	}


	
}


close(con); 

print(proc.time() - ptm)


stopCluster(cl)

