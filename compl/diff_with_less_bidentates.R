bidentates128 <- read.csv('bidentates128', header = FALSE)
oldBidentates <- read.csv('finalSmi/finalSmiBidentate.txt', header = FALSE, sep =' ')

hate20_smi <- (oldBidentates[!(oldBidentates$V1 %in% bidentates128$V1),])

hate20_ind <- as.numeric(rownames(hate20_smi))

seq1 <- as.character(seq(1,148)) # current labs
seq2 <- as.character(seq(1,128)) # target labs
bad_list <- as.character(hate20_ind) # bad labs
count <-0
final_mapping <- list()

for (i in seq1){
  if (i %in% bad_list)
  {
    final_mapping[i]<-'zzz'
  }
  else{
    count <- count + 1 
    final_mapping[i]<-seq2[count]
  }
}

# V1 V2 V3
# 21    [N+]-[O]-[O]-[N+]  1  4
# 22    [N+]-[S]-[S]-[N+]  1  4
# 52  [O+]-[N+]-[N+]-[O+]  1  4
# 53    [O+]-[O]-[O]-[O+]  1  4
# 54  [O+]-[P+]-[P+]-[O+]  1  4
# 55    [O+]-[S]-[S]-[O+]  1  4
# 73    [O]-[C+]=[C+]-[O]  1  4
# 74  [O]-[CH+]-[CH+]-[O]  1  6
# 77    [O]-[N+]-[N+]-[O]  1  4
# 81    [O]-[P+]-[P+]-[O]  1  4
# 85    [P+]-[O]-[O]-[P+]  1  4
# 86    [P+]-[S]-[S]-[P+]  1  4
# 116 [S+]-[N+]-[N+]-[S+]  1  4
# 117   [S+]-[O]-[O]-[S+]  1  4
# 118 [S+]-[P+]-[P+]-[S+]  1  4
# 119   [S+]-[S]-[S]-[S+]  1  4
# 137   [S]-[C+]=[C+]-[S]  1  4
# 138 [S]-[CH+]-[CH+]-[S]  1  6
# 141   [S]-[N+]-[N+]-[S]  1  4
# 145   [S]-[P+]-[P+]-[S]  1  4
# 
# 21  22  52  53  54  55  73  74  77  81  85  86 116 117 118 119 137 138 141 145
# these numbers are 1 based. to compare with the results eqn, subtract 1