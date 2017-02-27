##sandpi.py
##Answer for Class

##DETERMINE S and PI
##input number of samples and a list of genera (separated by ',') to include for the genus
##EXAMPLE python sandpi.py 5 Homo,Gorilla
import sys
import pysam
import numpy as np
import random

pD="/home/melinda_yang/greatapes/" ##Set this to whatever directory your data is in
samplecount,genera=int(sys.argv[1]),sys.argv[2].split(',') ##Read input options into variables
bedfile=open(pD+"Intersect_filtered_cov8_chr21_rand1000.bed",'r') ##BED FILE
vcfin=pysam.VariantFile(pD+"greatapes.fixedchr21.vcf.gz") ##VCF FILE
header = list(vcfin.header.samples)  ##Get list of all IDs in header

#####      MAKE GENOTYPE ARRAY FOR ALL INDIVIDUALS
seqary=np.zeros((1,len(header))) ##EACH ROW IS A SNP, EACH COLUMN IS AN INDIVIDUAL
num=0 ##Index to count number of SNPs that we wanted to keep for analysis 
for ind,line in enumerate(bedfile): ##Loop over every set of SNPs in BED file
    x=line.split()
    ##if ind%100==0: print ind  ##A way I check to see how long it's taking to go through bed file
    for pos in range(int(x[1]),int(x[2])):
        for entry in vcfin.fetch(x[0],pos-1,pos): ##Use .fetch (pysam method) to run through and grab the line for each position in the VCF
            if len(entry.alleles) != 2: continue ##Check if biallelic; continue if not
            ref,alt=tuple(entry.alleles) ##Not necessary for computation, but way to call alleles in vcf
            ##Grab genotypes for all individuals, marking them as 0 or 1 depending on genotype
            genos=[] 
            for geno in str(entry).split('\t')[9:]:
                vcfgeno=geno.split(':')[0]
                if vcfgeno in ['1/1']: genos.append(1)
                elif vcfgeno in ['0/0','./.']: genos.append(0) ##Treat ./. as fixed for reference allele
                elif vcfgeno in ['1/0','0/1']: genos.append(random.choice([0,1])) ##Randomly choose allele for heterozygotes
                else: print 'what is vcfgeno?',vcfgeno ##Just a check to make sure I don't miss anything
            genos=np.array(genos)
            if num==0: seqary=genos ##Create array of genotypes at first line with information
            else: seqary=np.vstack((seqary,genos)) ##Add new genotypes for next line to array of genotypes
            num+=1 ##Keeping track of every SNP we used in analysis
##Save all genotype information into an array using np.savetxt (use np.loadtxt to retrieve)
np.savetxt(pD+"greatapes_genos.txt",seqary,fmt="%.1f",comments="#".join(header)) 
#####

##If you do not want to make the genotype array every run, you can comment out everything
##included in the ##### above
seqary=np.loadtxt(pD+"greatapes_genos.txt") ##Loads genotype array if you ran the part in #####

##Find random set of indices (with count based on "samplecount") for each genera included.
##Add these indices to dictionary with keys as each genera
myapeinds={}
for ind,speciesind in enumerate(header):
    genusname=speciesind.split('_')[0] ##Grab the genus from the entire individual's ID
    if genusname not in myapeinds: myapeinds[genusname]=[ind] ##Initialize each dictionary entry
    else: myapeinds[genusname].append(ind) ##Add new indices to correct dictionary entry

##Get random sample of each genus
for mygenus in genera: myapeinds[mygenus]=sorted(random.sample(myapeinds[mygenus],samplecount))
print myapeinds  ##Dictionary of indices for each genus wanted      

##This for loop goes through the genotype array and picks out the columns 
##including the individuals I want to calculate S and pi
for myape in myapeinds: 
    myary=seqary[:,myapeinds[myape]] ##Get array of genotypes for individuals you want
    ##Find S: Add 1 for each site that is not fixed reference or fixed alternative
    S=0
    for i in np.sum(myary,axis=1):
        if i not in [0,myary.shape[1]]: S+=1
        
    #Calculate pi: For every pair of individuals, subtract their array of sites and find all 
    ##non-zero entries (where the two individuals are dissimilar). These are given a value of 1
    ##using np.count_nonzero. Thus, summing across the new array of differences, gives the total 
    ##number of differences between those two individuals, which are added to a list, for which 
    ##calculate the average.
    pairdiff=[]
    for myind1 in range(myary.shape[1]):
        for myind2 in range(myary.shape[1]):
            if myind1>=myind2: continue ##Don't accidentally include a pair more than once!
            else: pairdiff.append(np.sum(np.count_nonzero(myary[:,myind1]-myary[:,myind2])))
    pi = sum(pairdiff)/float(len(pairdiff)) #Get average for pi
    
    print myape, S, pi
