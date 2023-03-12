#!/usr/bin/env python3
# script: SPLASHPeaksStructures.py
# author: Daniel Desiro'

# dependencies: Python (v3.9.7), Vienna RNA Package (v2.5.0)

import sys
import os
import re
try:
    from RNA import fold_compound, cvar, CONSTRAINT_DB, CONSTRAINT_DB_DEFAULT
except:
    print("Error: The mandatory library ViennaRNA 2.5.0 is missing, please read the README.md for more information!")
    exit()




################################################################################
## settings
################################################################################

def settings():
    ############################################################################
    ## options
    opt = dict()
    opt["rst"] = "computer" # input counting option (biology or computer)
    opt["rvc"] = True       # create reverse complement for structure prediction
    opt["rvp"] = True       # reverse input positions
    opt["pkw"] = 20         # peak extension
    opt["rvo"] = True       # reverse output positions
    ############################################################################
    ## handle data
    mainDir = os.getcwd()
    inGenomesA = ["SC35MHA4x", "SC35MHA8x", "SC35Mwt", "HA5_1", "HAmuttoNAterm", "MmuttoHA5_1", "NAterm"]
    inReadsA   = ["merged_peak_tables_final"]*len(inGenomesA)
    inGenomesB = ["WyNAUdsub", "WyNAwt"]
    inReadsB   = ["PR8-Wy_merged_peak_tables"]*len(inGenomesB)
    for infasta,inSPLASH in zip(inGenomesA+inGenomesB,inReadsA+inReadsB):
        opt["pfx"] = os.path.join(mainDir,"results","SPLASHPeaksStructures",infasta) # output prefix
        opt["fsa"] = os.path.join(mainDir,"data",f"{infasta}.fa")                    # fasta genome file
        opt["stb"] = os.path.join(mainDir,"data",f"{inSPLASH}.tsv")                  # SPLASH table file
    return opt




################################################################################
## main
################################################################################

def main(opt):
    ########################################################################
    ## create output directory
    makeDir(opt)
    ########################################################################
    ## read fasta file
    fasta_dict = readFasta(opt.fsa, rvc=opt.rvc)
    ########################################################################
    ## read SPLASH reads
    tlist, header = readTable(opt.stb)
    ########################################################################
    ## extract SPLASH reads
    tlist = extractReads(tlist, fasta_dict, opt)
    ########################################################################
    ## Reverse outputs
    if opt.rvo: tlist = reverseOutput(tlist)
    ########################################################################
    ## write fasta reads
    header = ["number", "aSeq", "aType", "aStrand", "bSeq", "bType", "bStrand", "aPeak", "bPeak",
              "ai", "aj", "bi", "bj", "mfe", "RNA", "structure",
              "pai", "paj", "pbi", "pbj", "peak_mfe", "peak_RNA", "peak_structure"]
    writeTable(tlist, header, opt)




################################################################################
## functions
################################################################################

def readFasta(fasta, rvc=""):
    ## read fasta file
    fasta_dict, RNA = dict(), ""
    with open(fasta, "r") as infa:
        for line in infa:
            line = line.strip()
            if re.match(r">", line):
                if RNA: fasta_dict[name] = revComp(RNA, rvc)
                RNA, name = "", line[1:].split()[0]
            else:
                RNA += line
        fasta_dict[name] = revComp(RNA, rvc)
    return fasta_dict

def revComp(RNA, rvc):
    ## complement dictionary, or transform DNA to RNA
    RNA = RNA.upper()
    D2Rc = {"A":"U","T":"A","U":"A","C":"G","G":"C","R":"Y","Y":"R","M":"K",\
            "K":"M","S":"W","W":"S","B":"V","V":"B","D":"H","H":"D","N":"N"}
    if rvc:
        RNA = "".join(D2Rc[i] for i in RNA)
        RNA = RNA[::-1]
    else:
        RNA = RNA.replace("T","U")
    return RNA

def readTable(inTable):
    ## read tables and create links
    tlist = list()
    with open(inTable, "r") as tabin:
        header = next(tabin)
        header = header.strip().split(sep="\t")
        for i,line in enumerate(tabin):
            lnk = dict()
            line = line.strip().split(sep="\t")
            if len(header) != len(line):
                print(len(header),len(line))
                print(f"Error: Not the same number of header and list elements!")
                sys.exit()
            for name,item in zip(header,line):
                lnk[name] = item
            lk = links(**lnk)
            tlist.append(lk)
    return tlist, header

class links(object):
    def __init__(self, **data):
        self.__dict__.update((k,self.trans(v)) for k,v in data.items())
    def trans(self, s):
        if isinstance(s, str):
            try: return int(s)
            except ValueError:
                try: return float(s)
                except ValueError:
                    if s in ["True", "False"]: return s == "True"
                    else: return s
        else: return s
    def plot(self, sep):
        ldat = sep.join([f"{var}" for key,var in vars(self).items()])
        return ldat

def extractReads(tlist, fasta_dict, opt):
    ## extract SPLASH reads
    new_tlist = list()
    total = len(tlist)
    for iter,lk in enumerate(tlist):
        print(f"Status: Interactions {((iter+1)/total)*100:.2f} % ...      ", end="\r")
        lk.alen, lk.blen = len(fasta_dict[lk.aSeq]), len(fasta_dict[lk.bSeq])
        if lk.aj > lk.alen or lk.bj > lk.blen: continue
        if opt.rst == "computer":
            lk.ai, lk.bi = lk.ai+1, lk.bi+1 # change to bio for rev comp
            lk.aPeak, lk.bPeak = lk.aPeak+1, lk.bPeak+1
        if opt.rvp:
            lk.ai, lk.aj, lk.bi, lk.bj = lk.alen-lk.aj, lk.alen-lk.ai+1, lk.blen-lk.bj, lk.blen-lk.bi+1 # removed -1
            lk.aPeak, lk.bPeak = lk.alen-lk.aPeak, lk.blen-lk.bPeak
        else:
            lk.ai, lk.bi = lk.ai-1, lk.bi-1
            lk.aPeak, lk.bPeak = lk.aPeak-1, lk.bPeak-1
        if not (lk.__dict__.get("aStrand",False) and lk.__dict__.get("bStrand",False)):
            if opt.rvc:
                lk.aStrand = "-"
                lk.bStrand = "-"
            else:
                lk.aStrand = "+"
                lk.bStrand = "+"
        ########################################################################
        ## fold peak positions
        lk.pai, lk.paj = max(lk.aPeak-opt.pkw+1,0), min(lk.aPeak+opt.pkw,lk.alen)
        lk.pbi, lk.pbj = max(lk.bPeak-opt.pkw+1,0), min(lk.bPeak+opt.pkw,lk.blen)
        lk.peak_RNA = f"{fasta_dict[lk.aSeq][lk.pai:lk.paj]}&{fasta_dict[lk.bSeq][lk.pbi:lk.pbj]}"
        constraint = f"{'<'*(lk.paj-lk.pai)}{'>'*(lk.pbj-lk.pbi)}"
        lk.peak_mfe, pattern = doCofold(lk.peak_RNA, constraint, opt)
        lk.peak_structure = pattern[:len(lk.peak_RNA.split("&")[0])]+"&"+pattern[len(lk.peak_RNA.split("&")[0]):]
        ########################################################################
        ## structure
        lk.RNA = f"{fasta_dict[lk.aSeq][lk.ai:lk.aj]}&{fasta_dict[lk.bSeq][lk.bi:lk.bj]}"
        constraint = f"{'<'*(lk.aj-lk.ai)}{'>'*(lk.bj-lk.bi)}"
        lk.mfe, pattern = doCofold(lk.RNA, constraint, opt)
        lk.structure = pattern[:len(lk.RNA.split("&")[0])]+"&"+pattern[len(lk.RNA.split("&")[0]):]
    return new_tlist

def doCofold(RNA, constraint, opt):
    ## do Cofold
    cvar.noLP = int(1)
    fc = fold_compound(RNA)
    fc.constraints_add(constraint, CONSTRAINT_DB | CONSTRAINT_DB_DEFAULT)
    pattern, mfe = fc.mfe()
    return mfe, pattern

def reverseOutput(tlist):
    ## reverse output
    for lk in tlist:
        lk.aStrand, lk.bStrand = "+", "+"
        # reverse base
        lk.ai, lk.aj, lk.bi, lk.bj = lk.alen-lk.aj, lk.alen-lk.ai, lk.blen-lk.bj, lk.blen-lk.bi
        lk.RNA = lk.RNA.split("&")[0][::-1]+"&"+lk.RNA.split("&")[1][::-1]
        lk.structure = lk.structure.split("&")[0][::-1]+"&"+lk.structure.split("&")[1][::-1]
        # reverse peak
        lk.aPeak, lk.bPeak = lk.alen-lk.aPeak-1, lk.blen-lk.bPeak-1
        lk.pai, lk.paj, lk.pbi, lk.pbj = lk.alen-lk.paj, lk.alen-lk.pai, lk.blen-lk.pbj, lk.blen-lk.pbi
        lk.peak_RNA = lk.peak_RNA.split("&")[0][::-1]+"&"+lk.peak_RNA.split("&")[1][::-1]
        lk.peak_structure = lk.peak_structure.split("&")[0][::-1]+"&"+lk.peak_structure.split("&")[1][::-1]
    return tlist

def writeTable(tlist, header, opt):
    ## write interaction table
    f_name = os.path.splitext(os.path.basename(os.path.abspath(opt.fsa)))[0]
    ext = tlist[0].aStrand
    with open(os.path.join(opt.pfx, f"{f_name}_structures{ext}.tsv"), "w") as tabout:
        tabout.write("\t".join(header)+"\n")
        for lk in tlist:
            line_list = [f"{lk.__dict__[h]+1}" if h in ["ai","bi","cai","cbi","pai","pbi","aPeak","bPeak"] else f"{lk.__dict__[h]}" for h in header]
            tabout.write("\t".join(line_list)+"\n")

def makeDir(opt):
    ## create directory
    if not os.path.isdir(opt.pfx):
        os.mkdir(opt.pfx)




################################################################################
## parser
################################################################################

class options(object):
    def __init__(self, **data):
        self.__dict__.update((k,v) for k,v in data.items())

if __name__ == "__main__":

    ############################################################################
    ## main
    copt = options(**settings())
    main(copt)
