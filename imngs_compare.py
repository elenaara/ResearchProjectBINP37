# -*- coding: utf-8 -*-
"""
--- IMNGS COMPARE ---

Created on Wed Apr 13 11:44:40 2022

@author: Elena Aramendía

This script takes two files including SRA numbers and checks if matches for an
SRA run in one of the files are also present in the other file.
NCBI info files can be included and in case the SRA is not found, the script
checks if codes for experiment, SRA study or Bioproject match.
I wrote the script according to my files (for example in the first SRA file
there are repeated SRA bc each experiment had many matching sequences),
with the goal of studying Breviates location and their co-occurrence with 
Arcobacter in different environments.

To study the environmental location I focused on METAGENOMES and filtered the
sampling sites found in the NCBI file (in the "ScientificName" field) 
keeping the ones specified as "[environment] metagenome".

OUTPUT:
3 plots in 3 different files
    - Pie plot representing % of SRA experiments (unique SRA numbers) found
    in both files, i.e matching bothe organisms.
    - Pie plot representing environments in which the SRAs in file 1 are found.
    - Bar plot representing environments in which the SRAs in file 1 + file 2 are found.

2 text files:
    - SRA number for file 1 co-occurring sequences
    - SRA number for file 2 co-occurring sequences
    
OPTIONS:
    - i1: First input SRA file
    - i2: Second input SRA file
    - N1: NCBI file with SRA info for the 1st file, in csv format
    - N2: NCBI file with SRA info for the 2nd file, in csv format
    - A: FIlter out experiments labeled as "AMPLICON" in the "LibraryStrategy" 
    field of the NCBI info file (both files)
    - f: Format for output plots (pdf, svg, jpeg...). Default is pdf.
    - o: Basename for output files. Default is "plot".
    
EXAMPLE USAGE:
    python imngs_compare.py -i1 SRA_FILE1 -i2 SRA_FILE2 -N1 NCBI_info1 -N2 NCBI_info2 -f svg
    
    python imngs_compare.py -i1 prob_brevs_sra.txt -i2 arco_sra_unique.txt -N1 prob_brevs_RunInfo.csv -N2 arco_RunInfo.csv
"""
#%% Modules
import argparse
import pandas as pd
import matplotlib.pyplot as plt

#%% Arg parser

parser = argparse.ArgumentParser(description=' This script takes two files including SRA numbers and checks if matches for an SRA run in one of the files are also present in the other file. Outputs different plots summing up the results.',
                                 epilog='Example usage: python imngs_compare.py -i1 SRA_FILE1 -i2 SRA_FILE2 -N1 NCBI_info1 -N2 NCBI_info2 -f svg')

parser.add_argument(                                                  # IMNGS abundance matrix
    '-i1',
    type = argparse.FileType('r'),
    dest = 'sra1',
    required = True,
    help = 'First SRA file'
    )

parser.add_argument(                                                  # IMNGS abundance matrix
    '-i2',
    type = argparse.FileType('r'),
    dest = 'sra2',
    required = True,
    help = 'Second SRA file'
    )

parser.add_argument(                                                  # NCBI file
    '-N1',
    type = argparse.FileType('r'),
    dest = 'inf1',
    required = True,
    help = 'NCBI file with SRA info for the 1st file, in csv format.'
    )

parser.add_argument(                                                  # NCBI file
    '-N2',
    type = argparse.FileType('r'),
    dest = 'inf2',
    required = True,
    help = 'NCBI file with SRA info for the 2nd file, in csv format.'
    )


parser.add_argument(                                                  # plot counts
    '--A',
    action='store_true',
    default=False,
    dest = 'no_amplicon',
    required = False,
    help = 'Filter experiments labeled as "AMPLICON".'
    )

parser.add_argument(                                                  # output 
    '--f',
    type = str,
    dest = 'format',
    required = False,
    default = 'pdf',
    help = 'Format for output plots (pdf, svg, jpeg...). Default is pdf.'
    )

parser.add_argument(                                                  # output 
    '-o',
    type = str,
    dest = 'out',
    required = False,
    default = 'plot',
    help = 'Basename for output files. Default is "plot".'
    )


args = parser.parse_args() 

#%% TEST FILES
"""
ab1=pd.read_csv("C:/Users/earam/Desktop/cole/BINP37_SpringProj/SILVA_IMNGS/prob_brevs_sra.txt",header=None)
ab2 = pd.read_csv("C:/Users/earam/Desktop/cole/BINP37_SpringProj/SILVA_IMNGS/arco_sra_unique.txt",header=None)
inf1 = pd.read_csv("C:/Users/earam/Desktop/cole/BINP37_SpringProj/SILVA_IMNGS/prob_brevs_RunInfo.csv")
inf2 = pd.read_csv("C:/Users/earam/Desktop/cole/BINP37_SpringProj/SILVA_IMNGS/arco_RunInfo.csv",dtype={'Subject_ID':'str','Body_Site':'str'})
no_amp = False
out_file = "AMP"
plot_format = "pdf"
"""
#%% Read input

# SRA number files
ab1 = pd.read_csv(args.sra1)
ab2 = pd.read_csv(args.sra2)

# NCBI RunInfo files
inf1 = pd.read_csv(args.inf1)
inf2 = pd.read_csv(args.inf2)

## OPTIONS
# Filter out amplicon experiments
no_amp = args.no_amplicon 

# Output Plot options
out_file = args.out
plot_format = args.format
    

# %% Check Brev + Arco 
# The first SRA file is read and for each SRA number checks if it is in the oher file

# Name columns in the files
ab1.columns = ["Sample_Name"]
ab2.columns = ["Sample_Name"]

# NOT FILTERING AMPLICON
if no_amp == False:
        
    # Create df for counts
    # Get unique SRAs from ab1
    samples = set(ab1['Sample_Name'])
    # Empty column for counts
    counts = [0]*len(samples)
    symb_counts = pd.DataFrame({'counts':counts},index=samples)
    # SRA/BioProject/whatever numbers
    symb_number1= set()
    symb_number2= set()
    
    # For each SRA (column 1)
    for i in ab1['Sample_Name']: # For each SRA number
        if i in list(ab2['Sample_Name']):  # If it is found in the arcobacter file, add count
            symb_counts.loc[i,'counts'] += 1
            symb_number1.add(i)
            symb_number2.add(i)
            
            # IF NOT, check BioProject, SRA Study, Experiment
        elif i not in list(ab2['Sample_Name']):
            if i in list(inf1['Run']):
                exp = inf1.loc[inf1['Run'] == i]['Experiment'].iloc[0]
                study = inf1.loc[inf1['Run'] == i]['SRAStudy'].iloc[0]
                bioproject = inf1.loc[inf1['Run'] == i]['BioProject'].iloc[0]
                
                if bioproject in list(inf2['BioProject']) and not pd.isna(bioproject):
                    symb_counts.loc[i,'counts'] += 1
                    sra = list(inf2.loc[inf2['BioProject'] == bioproject]['Run'])[0]
                    symb_number2.add(sra)
                    symb_number1.add(i)
                    
                elif exp in list(inf2['Experiment']) and not pd.isna(exp):
                    symb_counts.loc[i,'counts'] += 1
                    sra = list(inf2.loc[inf2['Experiment'] == exp]['Run'])[0]
                    symb_number2.add(sra)
                    symb_number1.add(i)
                    
                elif study in list(inf2['SRAStudy']) and not pd.isna(study):
                    symb_counts.loc[i,'counts'] += 1
                    sra = list(inf2.loc[inf2['SRAStudy'] == study]['Run'])[0]
                    symb_number2.add(sra)
                    symb_number1.add(i)

    
    # Calculate % of experiments with both organisms
    counts = len(symb_counts.loc[~(symb_counts==0).all(axis=1)])
    total = len(inf1['Run'])
    per_symb = (counts*100)/total
    print('{:.2f}% of the SRAs match both organisms'.format(per_symb))

    
#FILTER OUT AMPLICON EXPERIMENTS
if no_amp:
    inf1 = inf1.loc[inf1["LibraryStrategy"] != "AMPLICON"]
    inf2 = inf2.loc[inf2["LibraryStrategy"] != "AMPLICON"]
    
    # Create df for counts
    # Get unique SRAs from ab1
    samples = set(ab1['Sample_Name'])
    # Empty column for counts
    counts = [0]*len(samples)
    symb_counts = pd.DataFrame({'counts':counts},index=samples)
    
    # SRA/BioProject/whatever numbers
    symb_number1= set()
    symb_number2= set()
    
    # For each SRA (column 1)
    for i in ab1['Sample_Name']: # For each SRA number
        if i in list(ab2['Sample_Name']):  # If it is found in the arcobacter file, add count
            symb_counts.loc[i,'counts'] += 1
            symb_number1.add(i)
            symb_number2.add(i)
            
            # IF NOT, check BioProject, SRA Study, Experiment
        elif i not in list(ab2['Sample_Name']):
            if i in list(inf1['Run']):
                exp = inf1.loc[inf1['Run'] == i]['Experiment'].iloc[0]
                study = inf1.loc[inf1['Run'] == i]['SRAStudy'].iloc[0]
                bioproject = inf1.loc[inf1['Run'] == i]['BioProject'].iloc[0]
                
                if bioproject in list(inf2['BioProject']) and not pd.isna(bioproject):
                    symb_counts.loc[i,'counts'] += 1
                    sra = list(inf2.loc[inf2['BioProject'] == bioproject]['Run'])[0]
                    symb_number2.add(sra)
                    symb_number1.add(i)
                    
                elif exp in list(inf2['Experiment']) and not pd.isna(exp):
                    symb_counts.loc[i,'counts'] += 1
                    sra = list(inf2.loc[inf2['Experiment'] == exp]['Run'])[0]
                    symb_number2.add(sra)
                    symb_number1.add(i)
                    
                elif study in list(inf2['SRAStudy']) and not pd.isna(study):
                    symb_counts.loc[i,'counts'] += 1
                    sra = list(inf2.loc[inf2['SRAStudy'] == study]['Run'])[0]
                    symb_number2.add(sra)
                    symb_number1.add(i)
    
    counts = len(symb_counts.loc[~(symb_counts==0).all(axis=1)])
    total = len(inf1['Run'])
    per_symb = (counts*100)/total
    print('{:.2f}% of the SRAs match both organisms'.format(per_symb))

# Print SRAs that match both organisims
sra_file1 = out_file + '_sras1.txt'
sra_file2 = out_file + '_sras2.txt'

with open(sra_file1, "w") as out1, open(sra_file2, "w") as out2:
    for i in symb_number1:
        print(i, file=out1)
    for j in symb_number2:
        print(j, file=out2)
        
# %% PLOT PERCENTAGE
# Color palette
purple_colors = ["#99C4C8","#B1BCE6","#9A86A4","#97C4B8","#e3d5ca"]

# Name for output file
plot_file = out_file + '_perplot.' + plot_format # name for output file

# Get percentage of experiments with both organisms
per = [per_symb,100-per_symb]


# Title of the plot
if no_amp:
    title = "% of experiments matching both Breviates and Arcobacter sequences (excluding amplicon data)"
else:
    title = "% of experiments matching both Breviates and Arcobacter sequences"

plt.rcParams["figure.figsize"] = (10,10)
patches, texts, pcts = plt.pie(per,
                               #labels=labels,
                               autopct='%.1f%%',
                               textprops={'size': 'large'},
                               colors=[purple_colors[2],"#e8e8e4"]
                               )

pcts[1].set_alpha(0)
plt.title(title, fontsize=15, fontweight=600)
plt.savefig(plot_file,format=plot_format) # save plot

#%% Check environments of brevs SRAs

# First get all different sample sites
allsites = set()
for i in inf1['ScientificName']:
    # Filter the sites
    # Only metagenome where the origin is specified (ie, 'soil metagenome')
    site = i.split(' ')
    if len(site) == 2 and site[1] == 'metagenome':
        allsites.add(i)
for i in inf2['ScientificName']:
    # Filter the sites
    # Only metagenome where the origin is specified (ie, 'soil metagenome')
    site = i.split(' ')
    if len(site) == 2 and site[1] == 'metagenome':
        allsites.add(i)
    

# Create new df
env_counts = [0]*len(allsites)
brev_envs = pd.DataFrame({'# of SRA experiments':env_counts}, index = allsites)

# Sites were we find the SRA 
sites = set()
# 
for i in list(ab1['Sample_Name']):
    if i in list(inf1['Run']):
        site = inf1.loc[inf1['Run'] == i]['ScientificName'].iloc[0]
        sites.add(site)
        if site in allsites:
            brev_envs.loc[site,'# of SRA experiments'] += 1

brev_envs = brev_envs.loc[~(brev_envs==0).all(axis=1)]
brev_envs = brev_envs.sort_values("# of SRA experiments")

## This is sequences, not exps
# In each SRA experiment there are different matches for brev sequences
# The ab1 SRA file has all the sequences (repeated SRAs)



#%% PIE PLOT

# Name for output file
plot_file = out_file + '_envplot.' + plot_format # name for output file

# Get labels (environments capitalized)
labels = list()
for i in brev_envs.index:
    label = i.split(" ")[0].capitalize()
    labels.append(label)

# Plot
plt.rcParams["figure.figsize"] = (10,10)
patches, texts, pcts = plt.pie(brev_envs['# of SRA experiments'],
                               labels=labels,
                               autopct='%.1f%%',
                               textprops={'size': 'large'},
                               colors=purple_colors
                               )
for i, patch in enumerate(patches):
    per = str(pcts[i]).split(',')[2].strip().strip(')').strip("'")
    per = float(per.rstrip('%'))
    #if per < 3.9:
        #pcts[i].set_alpha(0)
    if per < 1.2:
        pcts[i].set_alpha(0)
        texts[i].set_alpha(0)

plt.setp(pcts, color='white', fontweight=500)
plt.setp(texts, fontweight=400)
plt.title("% of Breviate sequences in different environments", fontsize=15, fontweight=600)
plt.savefig(plot_file,format=plot_format) # save plot
#%% Check environments of SRA with both organisms

# SRA with both organisms
symb_counts0 = symb_counts.loc[~(symb_counts==0).all(axis=1)]
symb_SRA = symb_counts.loc[~(symb_counts==0).all(axis=1)].index

# New df for environments
counts = [0]*len(allsites)
symb_envs = pd.DataFrame({'# of SRA experiments':counts}, index = allsites)

# Sites were we find the SRA 
sites = set()

for i in symb_SRA:
    if i in list(inf1['Run']):
        site = inf1.loc[inf1['Run'] == i]['ScientificName'].iloc[0]
        #sites.add(site)
        if site in allsites:
            sites.add(site)
            symb_envs.loc[site,'# of SRA experiments'] += 1
    #elif i in list(inf2['Run']):
     #   site = inf2.loc[inf2['Run'] == i]['ScientificName'].iloc[0]
      #  sites.add(site)
       # if site in allsites:
        #    #sites.add(site)
         #   symb_envs.loc[site,'# of SRA experiments'] += 1

symb_envs = symb_envs.loc[~(symb_envs==0).all(axis=1)] # Remove rows with 0


#%% Plot both things
# new df with bothe Brev sites and Brev + Arcobacter
all_envs = pd.DataFrame({'Breviates':brev_envs['# of SRA experiments'], 'Breviates + Arcobacter':(0*len(brev_envs.index))}, index=brev_envs.index)
for env in symb_envs.index:
    all_envs.loc[env,'Breviates + Arcobacter'] = symb_envs.loc[env,'# of SRA experiments']

# Sort values
all_envs = all_envs.sort_values('Breviates')


#%% PLOT

# Title of the plot
if no_amp:
    title = "Number of SRA sequences matching both Breviate and Arcobacter sequences (excluding amplicon data)"
else:
    title = "Number of SRA sequences matching both Breviate and Arcobacter sequences"

# Name for output file
plot_file = out_file + '_sitesboth.' + plot_format 

# Plot
plt.rcParams["figure.figsize"] = (10,6)
plot = all_envs['Breviates'].plot.barh(color=purple_colors[0],legend=True)
plot = all_envs['Breviates + Arcobacter'].plot.barh(color=purple_colors[2],legend=True)
plt.xlabel(title)
plt.grid(color='grey', linestyle='dotted', linewidth=0.5)
plt.tight_layout()
plt.savefig(plot_file,format=plot_format) # save plot

