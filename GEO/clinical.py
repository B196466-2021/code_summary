import pandas as pd
import os
import re

##batch
def clinic_info_process(dt):
    f=pd.read_excel(dt,engine='openpyxl').fillna(value='NA')
    #去掉列名中的'!Sample_'
    for x in f.columns:
        l=re.sub('!Sample_','',x)
        f=f.rename(columns={x:l})
    #提取行名中：前面的字段并修改列名为该字段
    for x in f.columns:
        if 'characteristics_ch1' in x:
            l=f[x].iloc[0].split(': ')[0]
            f=f.rename(columns={x:l})
    # 去掉列名中的'_ch1'
    for x in f.columns:
        if '_ch1' in x:
            l=re.sub('_ch1','',x)
            f=f.rename(columns={x:l})
    #清理每一行的字符
    for x in f.columns:
        if ': ' in f[x].iloc[0]:
            f[x]=f[x].str.replace(x+':','',regex=False)
    #save
    with pd.ExcelWriter ('Breast_cancer_clinical_HTA.xlsx',mode='a') as ff:
        f.to_excel(ff,sheet_name=dt.split('_')[0],index=None)

##Single
def clinic_info_process(dt):
    f=pd.read_excel(dt,engine='openpyxl').fillna(value='NA')
    #去掉列名中的'!Sample_'
    for x in f.columns:
        l=re.sub('!Sample_','',x)
        f=f.rename(columns={x:l})
    #提取行名中：前面的字段并修改列名为该字段
    for x in f.columns:
        if 'characteristics_ch1' in x:
            l=f[x].iloc[0].split(': ')[0]
            f=f.rename(columns={x:l})
    # 去掉列名中的'_ch1'
    for x in f.columns:
        if '_ch1' in x:
            l=re.sub('_ch1','',x)
            f=f.rename(columns={x:l})
    #清理每一行的字符
    for x in f.columns:
        if ': ' in f[x].iloc[0]:
            f[x]=f[x].str.replace(x+': ','',regex=False)
    #save
    f.to_excel(dt.split('.')[0]+'_info.xlsx',sheet_name=dt.split('_')[0],index=None)

#手动
#去掉列名中的'!Sample_'
for x in f.columns:
    l=re.sub('!Sample_','',x)
    f=f.rename(columns={x:l})
#提取行名中：前面的字段并修改列名为该字段
for x in f.columns:
    if 'characteristics_ch1' in x:
        l=f[x].iloc[0].split(': ')[0]
        f=f.rename(columns={x:l})

# 去掉列名中的'_ch1'
for x in f.columns:
    if '_ch1' in x:
        l=re.sub('_ch1','',x)
        f=f.rename(columns={x:l})
#清理每一行的字符
for x in f.columns:
    if ': ' in f[x].iloc[0]:
        f[x]=f[x].str.replace(x+': ','',regex=False)

f.to_excel('CRC-No-Survival-clinical-RSTA-1135.xlsx',index=None,sheet_name='GSE131418')

##处理Gene Symbol
for i in range(len(f)):
    if not f['Symbol'][i].startswith('--- // --- //'):
        f['Gene'][i]=f['Symbol'][i].split(",")[0].split('//')[1]

##2023-02-02 breast cancer merge
b=[x for x in os.listdir() if x.endswith('txt')] #read RMA normalised matrix
f=[pd.read_csv(x,sep='\t') for x in b]
ff=pd.concat(f,axis=1).T.drop_duplicates().T.rename(columns={'Unnamed: 0':'ID'})
anno=pd.read_csv('/mnt/d/23 solid tumors from GEO/platform/GPL17586-45144.txt',sep='\t',skiprows=15)[['ID','gene_assignment']]
new=pd.merge(ff,anno,on=['ID'])
gene=new.pop('Gene Symbol')
new.insert(1,'GeneSymbol',gene)
del new['ID']
new[~pd.isnull(new['GeneSymbol'])].to_csv('Merge_Breast_cancer_RSTA_platform_1090_41024.txt',sep='\t',index=None)


