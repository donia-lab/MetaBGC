import pysam
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


class BAMStats:
    """
    A class that include common functions to process bam alignment file from Bowtie2

    """
    def __init__(self, bamfile=None, quant_df=None, errfile=None):
        if bamfile:
            self.bam = pysam.AlignmentFile(bamfile, "rb")
            self.references = self.bam.references   
        if quant_df:
            self.quant_df = pd.read_csv(quant_df, header=None, delimiter="\t")
            self.quant_df.columns = ["Reference", "len", "mapped", "unmapped"]
        if errfile:
            with open(errfile, "r") as f:
                lines = f.readlines()
            self.pairs = int(lines[1].strip().split(" ")[0])
            if len(lines) > 15:
                self.singletons = int(lines[14].strip().split(" ")[0])
            else:
                self.singletons = 0
            self.tot_reads = 2 * self.pairs + self.singletons

    def coverage_breadth(self):
        """
        Calculate coverage breadth of reference sequences

        """
        cov_breadth = []
        for reference in self.bam.references: # for each sequence
            depth_per_pos = np.sum(self.bam.count_coverage(reference), axis=0)
            totlen = len(depth_per_pos)
            mapped_len = len(np.nonzero(depth_per_pos)[0])
            cov_breadth.append((reference, mapped_len / totlen))
        return pd.DataFrame(cov_breadth, columns=["Reference", "Breadth"])

    def coverage_depth(self):
        """
        Calculate coverage depth of reference sequences in terms of RPKM

        """
        df = self.quant_df.copy()[:-1]
        df["RPKM"] = df["mapped"] / df["len"] / self.tot_reads * 1e9
        return df[["Reference", "RPKM", "mapped"]]

    def coverage(self):
        """
        Calculate coverage breadth and coverage depth (RPKM) of reference sequences

        """
        breadth_df = self.coverage_breadth()
        depth_df = self.coverage_depth()
        self.cov = pd.merge(breadth_df, depth_df, on="Reference")
        return self.cov

    def depth_profile(self, reference=None):
        """
        Calculate coverage depth of reference sequences per position

        """
        if not hasattr(self, "depth_dict"):
            self.depth_dict = {}
        if not reference:
            for ref in self.references: # for each scaffold
                depth_per_pos = np.sum(self.bam.count_coverage(ref), axis=0)
                self.depth_dict[ref] = depth_per_pos
        else:
            # given specific reference sequence
            depth_per_pos = np.sum(self.bam.count_coverage(reference), axis=0)
            self.depth_dict[reference] = depth_per_pos
        return self.depth_dict

    def plot(self, reference, out, smoothing=None):
        """
        Plot coverage depth profile of a reference sequence 
        
        """
        try:
            plot_data = pd.Series(self.depth_dict[reference])
        except:
            plot_data = pd.Series(self.depth_profile(reference)[reference])
        if smoothing:
            # calculate rolling mean
            plot_data = plot_data.rolling(smoothing).mean()
        plt.cla()
        ax = plot_data.plot(title="Coverage depth profile of {}".format(reference))
        ax.set_xlabel("position")
        ax.set_ylabel("coverage")
        return ax.get_figure().savefig(out)

    def retrieve_reads(self, reference=None):
        """
        Retrieve reads mapping to a reference sequence or all 
        references.

        Parameters
        ----------
            reference : str, optional (default=None)
                Name of the reference sequence. If None, all mapped 
                reads to all reference sequences would be fetched. 

        Returns
        -------
            reads : list
                A list of ids of reads mapped to the reference 
                sequence if reference is specified. Otherwise,
                retrieve all mapped reads.
        """
        reads = []
        for read in self.bam.fetch(reference):
            reads.append(read.qname)
        return reads
