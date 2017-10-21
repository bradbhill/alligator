'''
Genomics.py
B Hill
Package to deal with the human genome.
Link to said genome: 
    ftp://ftp.ensembl.org/pub/release-85/gff3/homo_sapiens/Homo_sapiens.GRCh38.85.gff3.gz

    Features -- Initialize with Features(dataFrame,name) or optionally, Features()
        .df is a pandas DataFrame object.
        .name is a string, updated procedurally as data is transformed.
        .colNames is the list of data frame field names.
        .buildNames() adds a 'gene' field/column to the data; it extracts this info
            from the 'attributes' field.



'''
import pandas as pd
import numpy as np
import operator as op

class Features:

    def __init__(self,dataFrame=None,name='Untitled',):
        '''
        Initialize Features object.
        dataFrame :: pandas DataFrame
        '''
        self.name=name
        self.colNames=['seqid', 'source', 'type', 'start', 'end', 'score',\
                       'strand','phase', 'attributes']
        self.df = dataFrame

    def __repr__(self):
        return str(self.df.sample(10))

    def seqIDs(self):
        # List
        return self.df['seqid'].unique()

    def types(self):
        # List
        return self.df['type'].unique()

    def sources(self):
        # List
        return self.df['source'].unique()

    def buildNames(self):
        # Adds gene column to self.df
        def getGene(string):
            # Returns gene name from attributes column format
            l=string.split(';')
            for i in range(len(l)):
                try:
                    x=l[i]
                except IndexError:
                    return 'index error'
                for i in range(len(x)):
                    if x[:i+1]=='Name=':
                        return x[i+1:]
            return ''
        print(self.name+': getting gene names...')
        self.df['gene']=[getGene(s) for s in self.df['attributes']]
        self.name+='--names'

    def makeBED(self):
        '''
        Given an out-filename, it will create a .BED file from self,
        which preserves chromosome name, chromosome start, end, and gene name
        '''
        print(self.name+': writing BED file...')
        chrnames=['chr'+c for c in self.df['seqid']]
        starts=self.df['start']
        ends=self.df['end']
        names=self.df['gene']
        d=pd.DataFrame({'chrom':chrnames,'chromStart':starts,'chromEnd':ends,'name':names})
        filename=self.name+'.bed'
        if len(filename)>0:
            d.to_csv(path_or_buf='./'+filename,sep='\t',index=False,chunksize=16)
        else:
            return 0
        
    def assembled(self):
        '''
        Will return features which have been located on the chromosome.
        '''
        print(self.name+': filtering unassembled...')
        chrs = [str(_) for _ in range(1, 23)] + ['X', 'Y', 'MT']
        return Features(self.df[self.df['seqid'].isin(chrs)].copy(),self.name+'--a')

    def filtered(self,colName,val,r=op.eq):
        print(self.name+': subsetting \"'+colName+' == '+str(val)+'\"...')
        nm=self.name+'--'+colName+'-'+str(val)
        return Features(self.df[r(self.df[colName],val)].copy(),nm)

def main():

    colNames=['seqid', 'source', 'type', 'start', 'end', 'score',\
                       'strand','phase', 'attributes']
    #'Homo_sapiens.GRCh38.85.gff3.gz'
    fName=input()
    print('Loading data...')
    dataFrame = pd.read_csv(fName, compression='gzip',sep='\t',\
        comment='#', low_memory=False,header=None,names=colNames)
    
    print('Constructing Features object...')
    build=dataFrame[dataFrame.source.isin(['ensembl', 'havana', 'ensembl_havana'])].copy()
    hg=Features(build,'humanGenome').assembled()
    hg.buildNames()

    snos=hg.filtered('type','snoRNA')
    snos.makeBED()
    
    
main()
    



