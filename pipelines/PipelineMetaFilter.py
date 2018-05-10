'''
classes and utility functions for pipeline_MetaFilter.py

'''

import os

class SortMeRNA:    
    def __init__(self,seqdat,outfile,params):
        
        self.seqdat = seqdat
        self.outfile = outfile
        self.outdir = os.getcwd()+"/"+os.path.dirname(outfile)
        self.indir = os.getcwd()+"/"
        self.params = params
        self.filelocation = ""
        self.statementlist = ["mkdir -p {}".format(self.outdir)]
        self.rminter = False
        self.checkInterleave()
        self.buildStatement()
        self.deInterleave()
        
    #add statements to generate interleaved files if necessary
    def checkInterleave(self):
        #if single end or already interleaved just use the original file location
        if self.seqdat.paired == False or self.seqdat.interleaved == True:
            self.filelocation = os.getcwd()+"/"+self.seqdat.filepath
        #else generate commands to make temp interleaved reads
        else:
            self.statementlist.append("mkdir -p {}/interleaved".format(self.outdir))
            self.filelocation = self.outdir+"/interleaved/{}".format(self.seqdat.cleanname+"."+self.seqdat.fileformat)
            self.statementlist.append("seqtk mergepe {} {} >{}".format(self.indir+self.seqdat.filepath,self.indir+self.seqdat.pairedname,self.filelocation))

    #make the main call to sortmerna and clean interleaved if necessary
    def buildStatement(self):
        sortlist = ["sortmerna"]
        sortlist.append(self.refList())
        if self.params["SortMeRNA_paired"] == "out":
            sortlist.append("--paired_out")
        else:
            sortlist.append("--paired_in")
        sortlist.append("--reads {}".format(self.filelocation))
        sortlist.append("--aligned {}".format(self.outdir+"/aligned_"+self.seqdat.cleanname))
        sortlist.append("--other {}".format(self.outdir+"/other_"+self.seqdat.cleanname))
        if self.params["SortMeRNA_fastx"] == "true":
            sortlist.append("--fastx")
        if self.params["SortMeRNA_sam"] == "true":
            sortlist.append("--sam")
        if self.params["SortMeRNA_sq"] == "true":
            sortlist.append("--SQ")
        if self.params["SortMeRNA_blast"] != "false":
            sortlist.append("--blast {}".format(self.params["SortMeRNA_blast"]))
        if self.params["SortMeRNA_log"] == "true":
            sortlist.append("--log")
        if self.params["SortMeRNA_num_alignments"] != "false":
            sortlist.append("--num_alignments {}".format(self.params["SortMeRNA_num_alignments"]))
        if self.params["SortMeRNA_best"] != "false":
            sortlist.append("--best {}".format(self.params["SortMeRNA_best"]))
            sortlist.append("--min_lis {}".format(self.params["SortMeRNA_min_lis"]))
        if self.params["SortMeRNA_print_all_reads"] == "true":
            sortlist.append("--print_all_reads")
        sortlist.append("--match {}".format(self.params["SortMeRNA_match"]))
        sortlist.append("--mismatch {}".format(self.params["SortMeRNA_mismatch"]))
        sortlist.append("--gap_open {}".format(self.params["SortMeRNA_gap_open"]))
        sortlist.append("--gap_ext {}".format(self.params["SortMeRNA_gap_ext"]))
        sortlist.append("-N {}".format(self.params["SortMeRNA_n"]))
        if self.params["SortMeRNA_f"] == "true":
            sortlist.append("-F")
        if self.params["SortMeRNA_r"] == "true":
            sortlist.append("-R")
        sortlist.append("-a {}".format(self.params["SortMeRNA_threads"]))
        sortlist.append("-e {}".format(self.params["SortMeRNA_e"]))
        sortlist.append("-m {}".format(self.params["SortMeRNA_memory"]))
        if self.params["SortMeRNA_v"] == "true":
            sortlist.append("-v")
        self.statementlist.append(" ".join(sortlist))

    #if input was interleaved remove interleaved temp file (if necessary) and deinterleave output
    def deInterleave(self):
        if self.seqdat.paired == False:
            pass
        else:
            if self.seqdat.interleaved == False:
                #remove temp interleave file if one was made
                self.statementlist.append("rm -r {}".format(self.outdir+"/interleaved/"))
            else:
                pass
            #split the interleaved file to seperate paired ends
            otherloc = self.outdir+"/other_"+os.path.basename(self.filelocation)
            binloc = os.path.dirname(os.popen('which sortmerna').read())
            self.statementlist.append("bash {}/scripts/unmerge-paired-reads.sh {} {} {}".format(binloc,
                                                                                                otherloc,
                                                                                                otherloc+".1",otherloc+".2"))
            
    #abstract out making reference command as it is long 
    def refList(self):
        reffastas = self.params["SortMeRNA_rna_refs"].split(",")
        if self.params["SortMeRNA_rna_index"] != "false":
            refindex = self.params["SortMeRNA_rna_index"].split(",")
        else:
            refindex = reffastas
            refindex = [os.path.basename(i) for i in refindex]
            refindex = [os.getcwd()+"/ref_index.dir/{}-db".format(i) for i in refindex]
        paired = [",".join([x[0],x[1]]) for x in zip(reffastas,refindex)]
        return("--ref {}".format(":".join(paired)))
            

    def build(self):
        return(" && ".join(self.statementlist))
            
