# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 13:50:39 2020

@author: etb
"""
# this script is to automatically create 
# a manifest file for QIIME2
# note that this is based on how UQAM
# (University of Quebec at Montreal)
# names their sequencing files -
# you will likely have to adjust
# the below based on your file names

def auto_manifest(mydir, out):
    import pandas as pd
    import os
    
    #creates a pandas dataframe with the columns specified
    df=pd.DataFrame(columns=['sample-id','absolute-filepath','direction'])

    #turns the provided string mydir into bytes so you can move into the
    #specified directory
    directory=os.fsencode(mydir)

    #empty lists that will contain the values for sample-id, absolute-filepath,
    #and direction
    ids=[]
    paths=[]
    directions=[]

    #iterates through each file in the directory, only considering files that
    #end with .fastq.gz
    for file in os.listdir(directory):
        filename=os.fsdecode(file)
        if filename.endswith('.fastq.gz')==True:
            #adds the full file path to the paths list
            file_loc=mydir+'/'+filename
            paths.append(file_loc)
           
            #splits the filename using underscore as the delimiter
            #then adds this name to the ids list
            myid1=filename.split('_')
            ids.append(myid1[0])
            
            #determines the direction of the read based on the filename and 
            #adds the appropriate direction to the directions list
            #if filename.endswith('_R1.fastq.gz')==True:
            if filename.endswith('_R1_001.fastq.gz')==True:
                directions.append('forward')
            #elif filename.endswith('_R2.fastq.gz')==True:
            elif filename.endswith('_R2_001.fastq.gz')==True:
                directions.append('reverse')
   
    #adds the lists as values for the columns specified
    df['sample-id']=ids
    df['absolute-filepath']=paths
    df['direction']=directions

    #turns the dataframe into a csv file
    df.to_csv(out+'.csv',index=False)
    return

#example input
# give the location of your reads for "mydir"
# give the location of where you want
# the output file in "out

#mydir='/home/16S/Reads'
#out='/home/16S'

#auto_manifest(mydir,out)