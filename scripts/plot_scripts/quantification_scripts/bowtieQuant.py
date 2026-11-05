#!/bin/env python

'''Calculate the RPKM from a bowtie2 alignment.
'''

# Shuo Wang's RPKM quantifier
import BAMStats

def quantify(bamFilename, quantFilename, summaryFilename, readsFilename=None):

    bamStats = BAMStats.BAMStats(bamFilename, quantFilename, summaryFilename)

    if readsFilename:
        f = open(readsFilename, 'w')

        needHeader = True
        for read in bamStats.bam.fetch():
            if not read.is_unmapped:
                readDict = read.to_dict()
                if needHeader:
                    f.write('{}\n'.format('\t'.join(['reference name', 'read name', 'ref start pos', 'sequence', 'quality'])))
                    needHeader = False

                refStart = int(readDict['ref_pos'])
                f.write('{}\n'.format('\t'.join([readDict['ref_name'], readDict['name'], str(refStart), readDict['seq'], readDict['qual']])))

        f.close()
        
    return bamStats.coverage()


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='%s' % __doc__)

    parser.add_argument('--bamfile', required=True, help='BAM file from bowtie2')
    parser.add_argument('--quantfile', required=True, help='Output file from samtools idxstats')
    parser.add_argument('--summaryfile', help="The alignment summary from bowtie2 (from bowtie2's standard error")
    parser.add_argument('--sample_identifier', required=True, help="A sample identifier for labeling the output.")
    parser.add_argument('--readsfile', help="File to contain information on matching reads.")
    
    args = parser.parse_args()

    coverage = quantify(args.bamfile, args.quantfile, args.summaryfile, args.readsfile)

    for query, breadth, rpkm, mapped in coverage.iloc[:,:].values:
        print('%s\t%s\t%s\t%s\t%s' % (query, breadth, rpkm, mapped, args.sample_identifier))
