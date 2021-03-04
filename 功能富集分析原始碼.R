#functional enrichment
df_vector <- df_5_out_ct %>% 
  mutate(count = NULL)

#### Configure 

rm(list=ls())

#comPATH='/media/md1200'

casePATH='/home/leotsai0127/functional'  #output path

outfolder1=file.path(casePATH,'fun.ann_35gene')  #output folder

environmentDir=file.path(comPATH,'/analysis/Script/function-annotation/environment/')

#### Choose your species 

source(file.path(environmentDir,'code','annotation_human.r'))  #Human (Cytoscape)

#source(file.path(environmentDir,'code','annotation_mouse.r')) #Mouse (Cytoscape)

#source(file.path(environmentDir,'code','annotation_rat.r')) #Rat (Cytoscape)

#### data input 

genelist = df_vector$gene
genelist = as.vector(genelist)

#### function annotation 

if (!file.exists(outfolder1)) dir.create(outfolder1)

setwd(outfolder1)

annotation(genelist,environmentDir)
