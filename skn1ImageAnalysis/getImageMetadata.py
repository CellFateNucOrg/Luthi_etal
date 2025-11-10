#!/usr/bin/env python

from nd2reader import ND2Reader
import numpy as np
import re
import os
import glob
import nd2
import random
from fnmatch import fnmatch
import pandas as pd 
np.random.seed(6)
import argparse


def extractImageMetadata(work_dir,file_list):
        os.makedirs(work_dir+'/metadata', exist_ok=True)
        dicList=[]
        for file_name in file_list:
                base_name=str(re.sub('.nd2','', os.path.basename(file_name)))
                with nd2.ND2File(file_name) as ndfile:
                        #with open(work_dir+'/metadata/'+base_name+'_description.txt', "w") as text_file:
                        #        text_file.write(ndfile.text_info['description'])
                        with open(work_dir+'/metadata/'+base_name+'_capturing.txt', "w") as text_file:
                                text_file.write(ndfile.text_info['capturing'])
                        cap = pd.read_csv(work_dir+'/metadata/'+base_name+'_capturing.txt',sep=":",index_col=False,header=None,names=['attribute','value'])
                        data = {'base_name': base_name,
                                'numZ': ndfile.sizes['Z'],
                                'numXY': ndfile.sizes['Y'],
                                'channelCount': ndfile.attributes.channelCount,
                                'Zstep': ndfile.experiment[0].parameters.stepUm,
                                'bitsPercomponentInMemory': ndfile.attributes.bitsPerComponentInMemory,
                                'bitsPercomponentInSignificant': ndfile.attributes.bitsPerComponentSignificant,
                                'channelwavelength': ndfile.metadata.channels[0].channel.name,
                                'objectiveMagnification': ndfile.metadata.channels[0].microscope.objectiveMagnification,
                                'objectiveName': ndfile.metadata.channels[0].microscope.objectiveName,
                                'objectiveNumericalAperture': ndfile.metadata.channels[0].microscope.objectiveNumericalAperture,
                                'pinholediamterUm': ndfile.metadata.channels[0].microscope.pinholeDiameterUm,
                                #'MultiExcitation': cap.loc[cap['attribute']=='MultiExcitation','value'].values[0],
                                'dXPos': int(ndfile.unstructured_metadata()['ImageMetadataSeqLV|0']['SLxPictureMetadata']['dXPos']),
                                'dYPos': int(ndfile.unstructured_metadata()['ImageMetadataSeqLV|0']['SLxPictureMetadata']['dYPos']),
                                'path': ndfile.path,
                        }
                        samples = list(cap.index[cap['attribute'].str.contains('Sample')])
                        samples.append(len(cap))
                        for i in range(len(samples)-1):   
                                sampledf = cap.iloc[samples[i]:samples[i+1],:]
                                #print(sampledf)
                                capdata = {'C'+str(i+1)+'_Exposure': sampledf.loc[sampledf['attribute'].str.contains('Exposure'),'value'].values[0].strip(),
                                        'C'+str(i+1)+'_Binning': sampledf.loc[sampledf['attribute'].str.contains('Binning'),'value'].values[0].strip(),
                                        'C'+str(i+1)+'_ScanMode': sampledf.loc[sampledf['attribute'].str.contains('Scan Mode'),'value'].values[0].strip(),
                                        'C'+str(i+1)+'_Denoise.ai': sampledf.loc[sampledf['attribute'].str.contains('Denoise.ai'),'attribute'].values[0].strip().replace("Denoise.ai ",""),
                                        'C'+str(i+1)+'_Clarify.ai': sampledf.loc[sampledf['attribute'].str.contains('Clarify.ai'),'attribute'].values[0].strip().replace("Clarify.ai ","")
                                }
                                data.update(capdata)
                        dicList.append(data)
                        #with open(work_dir+'/metadata/'+base_name+'_unstructured.json', "w") as json_file: 
                        #        jsondump(ndfile.unstructured_metadata(), json_file)
                        os.remove(work_dir+'/metadata/'+base_name+'_capturing.txt')

        
        df=pd.DataFrame(dicList)
        df.to_csv(path_or_buf=work_dir+"/metadata/metadataTable.csv")
        return(df)
 

def get_args():
	parser = argparse.ArgumentParser()
	parser.add_argument("-d","--work_dir", help="Working directory in which to do analysis", required=True)
	return parser.parse_args()


def main():
    args=get_args()
    work_dir = args.work_dir
    print("setting up directories")
    pattern=work_dir+'/**/*.nd2'
    file_list=glob.glob(pattern,recursive=True)
    print("getting metadata")
    extractImageMetadata(work_dir,file_list)
    print("done")


if __name__=='__main__':
	main()

