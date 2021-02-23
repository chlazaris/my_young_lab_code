#Humans
aws s3 --no-sign-request --region eu-west-1 sync s3://ngi-igenomes/igenomes/Homo_sapiens/UCSC/hg19/Sequence/STARIndex/ ./references/Homo_sapiens/UCSC/hg19/Sequence/STARIndex/
aws s3 --no-sign-request --region eu-west-1 sync s3://ngi-igenomes/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/ ./references/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/
aws s3 --no-sign-request --region eu-west-1 sync s3://ngi-igenomes/igenomes/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes/ ./references/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes/
aws s3 --no-sign-request --region eu-west-1 sync s3://ngi-igenomes/igenomes/Homo_sapiens/UCSC/hg19/Sequence/AbundantSequences/ ./references/Homo_sapiens/UCSC/hg19/Sequence/AbundantSequences/
# mm10
aws s3 --no-sign-request --region eu-west-1 sync s3://ngi-igenomes/igenomes/Mus_musculus/UCSC/mm10/Annotation/Genes/ ./references/Mus_musculus/UCSC/mm10/Annotation/Genes/ --exclude "*" --include "genes.gtf"
aws s3 --no-sign-request --region eu-west-1 sync s3://ngi-igenomes/igenomes/Mus_musculus/UCSC/mm10/Annotation/Genes/ ./references/Mus_musculus/UCSC/mm10/Annotation/Genes/ --exclude "*" --include "genes.bed"
aws s3 --no-sign-request --region eu-west-1 sync s3://ngi-igenomes/igenomes/Mus_musculus/UCSC/mm10/Sequence/BismarkIndex/ ./references/Mus_musculus/UCSC/mm10/Sequence/BismarkIndex/
aws s3 --no-sign-request --region eu-west-1 sync s3://ngi-igenomes/igenomes/Mus_musculus/UCSC/mm10/Sequence/BowtieIndex/ ./references/Mus_musculus/UCSC/mm10/Sequence/BowtieIndex/
aws s3 --no-sign-request --region eu-west-1 sync s3://ngi-igenomes/igenomes/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/ ./references/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/
aws s3 --no-sign-request --region eu-west-1 sync s3://ngi-igenomes/igenomes/Mus_musculus/UCSC/mm10/Sequence/BWAIndex/ ./references/Mus_musculus/UCSC/mm10/Sequence/BWAIndex/
aws s3 --no-sign-request --region eu-west-1 sync s3://ngi-igenomes/igenomes/Mus_musculus/UCSC/mm10/Sequence/STARIndex/ ./references/Mus_musculus/UCSC/mm10/Sequence/STARIndex/
aws s3 --no-sign-request --region eu-west-1 sync s3://ngi-igenomes/igenomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/ ./references/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/
aws s3 --no-sign-request --region eu-west-1 sync s3://ngi-igenomes/igenomes/Mus_musculus/UCSC/mm10/Sequence/Chromosomes/ ./references/Mus_musculus/UCSC/mm10/Sequence/Chromosomes/
aws s3 --no-sign-request --region eu-west-1 sync s3://ngi-igenomes/igenomes/Mus_musculus/UCSC/mm10/Sequence/AbundantSequences/ ./references/Mus_musculus/UCSC/mm10/Sequence/AbundantSequences/
# mm9
aws s3 --no-sign-request --region eu-west-1 sync s3://ngi-igenomes/igenomes/Mus_musculus/UCSC/mm9/Annotation/Genes/ ./references/Mus_musculus/UCSC/mm9/Annotation/Genes/ --exclude "*" --include "genes.gtf"
aws s3 --no-sign-request --region eu-west-1 sync s3://ngi-igenomes/igenomes/Mus_musculus/UCSC/mm9/Annotation/Genes/ ./references/Mus_musculus/UCSC/mm9/Annotation/Genes/ --exclude "*" --include "genes.bed"
aws s3 --no-sign-request --region eu-west-1 sync s3://ngi-igenomes/igenomes/Mus_musculus/UCSC/mm9/Sequence/BismarkIndex/ ./references/Mus_musculus/UCSC/mm9/Sequence/BismarkIndex/
aws s3 --no-sign-request --region eu-west-1 sync s3://ngi-igenomes/igenomes/Mus_musculus/UCSC/mm9/Sequence/BowtieIndex/ ./references/Mus_musculus/UCSC/mm9/Sequence/BowtieIndex/
aws s3 --no-sign-request --region eu-west-1 sync s3://ngi-igenomes/igenomes/Mus_musculus/UCSC/mm9/Sequence/Bowtie2Index/ ./references/Mus_musculus/UCSC/mm9/Sequence/Bowtie2Index/
aws s3 --no-sign-request --region eu-west-1 sync s3://ngi-igenomes/igenomes/Mus_musculus/UCSC/mm9/Sequence/BWAIndex/ ./references/Mus_musculus/UCSC/mm9/Sequence/BWAIndex/
aws s3 --no-sign-request --region eu-west-1 sync s3://ngi-igenomes/igenomes/Mus_musculus/UCSC/mm9/Sequence/STARIndex/ ./references/Mus_musculus/UCSC/mm9/Sequence/STARIndex/
aws s3 --no-sign-request --region eu-west-1 sync s3://ngi-igenomes/igenomes/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/ ./references/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/
aws s3 --no-sign-request --region eu-west-1 sync s3://ngi-igenomes/igenomes/Mus_musculus/UCSC/mm9/Sequence/Chromosomes/ ./references/Mus_musculus/UCSC/mm9/Sequence/Chromosomes/
aws s3 --no-sign-request --region eu-west-1 sync s3://ngi-igenomes/igenomes/Mus_musculus/UCSC/mm9/Sequence/AbundantSequences/ ./references/Mus_musculus/UCSC/mm9/Sequence/AbundantSequences/
