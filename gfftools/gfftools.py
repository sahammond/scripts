"""Classes and functions to work with GFF3 files
Written by S. Austin Hammond, BCGSC

GFF: main class
Gene: object to hold gene records and associated transcripts
Transcript: object to hold exons, cds, and info associated with a
    particular transcript. Includes methods to print exon and CDS
    segments, and check if terminal CDS segments have start/stop
    codons.
"""

import re

class GFF(object):
    """Custom class to hold GFF features.

    Stores each column whole in objects, and also
    breaks out common attributes into their own
    objects. Initializes slots for certain optional
    attributes too.

    feature: a raw gff line (with tab chars, newline, etc.)
    """

    def __init__(self,feature):

        self.raw = feature

        rec = feature.strip("\n").split("\t")
        attr = rec[8].split(";")

        self.seqid = rec[0]
        self.source = rec[1]
        self.type = rec[2]
        self.start = rec[3]
        self.end = rec[4]
        self.score = rec[5]
        self.strand = rec[6]
        self.offset = rec[7]
        self.attributes = rec[8]
        self.locusid = "" # slot to hold locus id
        self.product = ""
        self.name = ""
        self.note = ""
        self.incomplete = ""

        for entry in attr:
            if len(re.findall("ID=",entry)) > 0:
                self.id = re.sub("ID=","",entry)
            elif len(re.findall("Parent=",entry)) > 0:
                self.parent = re.sub("Parent=","",entry).split(",") # some features have > 1 parent, sep by ","
            elif len(re.findall("product=",entry)) > 0:
                self.product = re.sub("product=","",entry)
            elif len(re.findall("Name=",entry)) > 0:
                self.name = re.sub("Name=","",entry)
            elif len(re.findall("note=",entry)) > 0:
                self.note = re.sub("note=","",entry)


def check_start_stop(sequence):
    """
    Check if a sequence begins with ATG and ends with TAA, TGA, or TAG
    Returns tuple of start - stop as Bools: (T|F,T|F)
    """

    starts = ['ATG']
    stops = ['TAA','TGA','TAG']

    startc = False
    stopc = False

    if sequence[0:3] in starts:
        startc = True

    if sequence[-3:] in stops:
        stopc = True

    return (startc, stopc)


class Transcript(object):
    """Custom class to hold GFF transcript features.

    Stores the exons and cds in lists by position. Note that if the feature is on the
    negative strand then they will need to be printed out in reverse order.
    """

    def add_exon(self,feature):
        if len(self.exons) == 0:
            self.exons.append(feature)

        else:
            if int(feature.start) == int(self.exons[-1].start):
                print ("There is already an exon with that start in"
                        "this transcript. Please check your input.")

            elif int(feature.start) > int(self.exons[-1].start):
                self.exons.append(feature)

            else:
                for pos in range(0,len(self.exons)):
                    if int(feature.start) < int(self.exons[pos].start):
                        self.exons.insert(pos,feature)
                        break


    def add_cds(self,feature):
        if len(self.cds) == 0:
            self.cds.append(feature)

        else:
            if int(feature.start) == int(self.cds[-1].start):
                print ("There is already a cds segment with that start in"
                        "this transcript. Please check your input.")

            elif int(feature.start) > int(self.cds[-1].start):
                self.cds.append(feature)

            else:
                for pos in range(0,len(self.cds)):
                    if int(feature.start) < int(self.cds[pos].start):
                        self.cds.insert(pos,feature)
                        break


    def check_start_stop(self,sequence):
        """
        Check if a sequence begins with ATG and ends with TAA, TGA, or TAG
        Returns tuple of start - stop as Bools: (T|F,T|F)
        """
    
        starts = ['ATG']
        stops = ['TAA','TGA','TAG']
    
        startc = False
        stopc = False
    
        if sequence[0:3] in starts:
            startc = True
    
        if sequence[-3:] in stops:
            stopc = True
    
        self.start_complete = startc
        self.stop_complete = stopc


    def print_single_seg(self,feature,format='gff'):
        """Print the only exons/CDS segment in a transcript.
        Format = 'gff' or 'tbl'
        """
        if format == 'gff':
            print feature.raw

        # TODO handle incomplete first seq
        elif format == 'tbl':
            if feature.type == "exon":
                if feature.strand == "-":
                    print "\t".join([feature.end,feature.start,"mRNA"])

                else:
                    print "\t".join([feature.start,feature.end,"mRNA"])

                if self.transcript.product:
                    print "\t".join(["\t\t\tproduct",self.transcript.product])

                if self.transcript.note:
                    print "\t".join(["\t\t\tnote",self.transcript.note])

                print "\t".join(["\t\t\tprotein_id","gnl|BCGSC|"+feature.parent])
                print "\t".join(["\t\t\ttransctript_id","gnl|BCGSC|"+feature.parent])

            elif feature.type == "CDS":
                if feature.strand == "-":
                    print "\t".join([feature.end,feature.start,"CDS"])

                else:
                    print "\t".join([feature.start,feature.end,"CDS"])

                if self.transcript.product:
                    print "\t".join(["\t\t\tproduct",self.transcript.product])

                if self.transcript.note:
                    print "\t".join(["\t\t\tnote",self.transcript.note])

                print "\t\t\tcodon_start\t1"
                print "\t".join(["\t\t\tprotein_id","gnl|BCGSC|"+feature.parent])
                print "\t".join(["\t\t\ttransctript_id","gnl|BCGSC|"+feature.parent])


    def print_internal_seg(self,feature,format='gff'):
        """Print a given exon/CDS line or segment.
        Format = 'gff' or 'tbl'
        """

        if format == 'gff':
            print feature.raw

        elif format == 'tbl':
            if feature.strand == "-":
                print "\t".join([feature.end,feature.start])

            else:
                print "\t".join([feature.start,feature.end])


    def print_first_seg(self,feature,format='gff'):
        """Print a given exon line or segment.
        Would work for CDS segment too.
        Format = 'gff' or 'tbl'
        """

        if format == 'gff':
            print feature.raw

        # TODO handle incomplete first seq
        elif format == 'tbl':
            if feature.type == "exon":
                if feature.strand == "-":
                    print "\t".join([feature.end,feature.start,"mRNA"])

                else:
                    print "\t".join([feature.start,feature.end,"mRNA"])

            elif feature.type == "CDS":
                if feature.strand == "-":
                    print "\t".join([feature.end,feature.start,"CDS"])

                else:
                    print "\t".join([feature.start,feature.end,"CDS"])


    def print_last_seg(self,feature,format='gff'):
        """Print the final exon line or segment.
        Would work for CDS segment too.
        Format = 'gff' or 'tbl'
        """

        if format == 'gff':
            print feature.raw

        # TODO handle incomplete first seq
        elif format == 'tbl':
            if feature.type == "exon":
                if feature.strand == "-":
                    print "\t".join([feature.end,feature.start])

                else:
                    print "\t".join([feature.start,feature.end])

                if self.transcript.product:
                    print "\t".join(["\t\t\tproduct",self.transcript.product])

                if self.transcript.note:
                    print "\t".join(["\t\t\tnote",self.transcript.note])

                print "\t".join(["\t\t\tprotein_id","gnl|BCGSC|"+feature.parent])
                print "\t".join(["\t\t\ttransctript_id","gnl|BCGSC|"+feature.parent])

            elif feature.type == "CDS":
                if feature.strand == "-":
                    print "\t".join([feature.end,feature.start])

                else:
                    print "\t".join([feature.start,feature.end])

                if self.transcript.product:
                    print "\t".join(["\t\t\tproduct",self.transcript.product])

                if self.transcript.note:
                    print "\t".join(["\t\t\tnote",self.transcript.note])

                print "\t\t\tcodon_start\t1"
                print "\t".join(["\t\t\tprotein_id","gnl|BCGSC|"+feature.parent])
                print "\t".join(["\t\t\ttransctript_id","gnl|BCGSC|"+feature.parent])


    def print_transcript(self,format='gff'):
        """Print a transcript with exons (and cds segments).
        Format = 'gff' or 'tbl'
        """

        if format == 'gff':
                print self.transcript.raw

                for segment in self.exons:
                    print segment.raw

                if self.cds:
                    for segment in self.cds:
                        print segment.raw

        elif format == 'tbl':
                if len(self.exons) == 1:
                    print_single_seg(self.exons[0])

                elif self.strand == "-":
                    if len(self.exons) == 2:
                        print_first_seg(self.exons[1])
                        print_last_seg(self.exons[0])
    
                    else:
                        print_first_seg(self.exons[-1])
    
                        for segment in range(0,len(self.exons))[1:-1][::-1]:
                            print_internal_seg(segment)
    
                        else:
                            print_last_seg(self.exons[0])

                else:
                    if len(self.exons) == 2:
                        print_first_seg(self.exons[0])
                        print_last_seg(self.exons[1])
    
                    else:
                        print_first_seg(self.exons[0])
    
                        for segment in range(0,len(self.exons))[1:-1]:
                            print_internal_seg(segment)
    
                        else:
                            print_last_seg(self.exons[-1])

                if self.cds:
                    if len(self.exons) == 1:
                        print_single_seg(self.exons[0])
    
                    elif self.strand == "-":
                        if len(self.exons) == 2:
                            print_first_seg(self.exons[1])
                            print_last_seg(self.exons[0])
        
                        else:
                            print_first_seg(self.exons[-1])
        
                            for segment in range(0,len(self.exons))[1:-1][::-1]:
                                print_internal_seg(segment)
        
                            else:
                                print_last_seg(self.exons[0])
    
                    else:
                        if len(self.exons) == 2:
                            print_first_seg(self.exons[0])
                            print_last_seg(self.exons[1])
        
                        else:
                            print_first_seg(self.exons[0])
        
                            for segment in range(0,len(self.exons))[1:-1]:
                                print_internal_seg(segment)
        
                            else:
                                print_last_seg(self.exons[-1])

    def __init__(self,feature):
        self.strand = feature.strand
        self.transcript = feature
        self.exons = []
        self.cds = []
        self.start_complete = ''
        self.stop_complete = ''


class Gene(object):
    """Class to hold whole gene entries from gff3 file

    gene:   a GFF object with type 'gene'

    Stores strand, type (mRNA, ncRNA, etc.), and transcripts. For each transcript,
    the exons, cds, and whether the cds is complete (at either end) can be added.
    The latter info is via the check_start_stop method, which requires a sequence string.
    """

    def check_start_stop(self,sequence):
        """
        Check if a sequence begins with ATG and ends with TAA, TGA, or TAG
        Returns tuple of start - stop as Bools: (T|F,T|F)
        """
    
        starts = ['ATG']
        stops = ['TAA','TGA','TAG']
    
        startc = False
        stopc = False
    
        if sequence[0:3] in starts:
            startc = True
    
        if sequence[-3:] in stops:
            stopc = True
    
#        return (startc, stopc)
        self.start_complete = startc
        self.stop_complete = stopc


    def add_transcript(self,feature,f=0):
        """
        Add a transcript GFF object to the gene. Will be put into a dict called
        'transcript', with key = transcript.id
        """

        if feature.id in self.transcript:
            print "Transcript already associated with gene."

        else:
            self.transcript[feature.id] = Transcript(feature)
#            self.transcript[transcript.id] = transcript


    def print_gene(self,format='gff'):
        """Print a gene line.
        Format = 'gff' or 'tbl'
        """

        if format == 'gff':
                print self.gene.raw

        elif format == 'tbl':
                if self.strand == "-":
                    if self.gene.locusid:
                        print "".join([self.gene.end,"\t",self.gene.start,"\tgene\n",
                                        "\t\t\tlocus_id\t",self.gene.locusid])
                    else:
                        print "".join([self.gene.end,"\t",self.gene.start,"\tgene\n",
                                        "\t\t\tlocus_id\t",self.gene.id])
                else:
                    if self.gene.locusid:
                        print "".join([self.gene.start,"\t",self.gene.end,"\tgene\n",
                                        "\t\t\tlocus_id\t",self.gene.locusid])
                    else:
                        print "".join([self.gene.start,"\t",self.gene.end,"\tgene\n",
                                        "\t\t\tlocus_id\t",self.gene.id])


    def __init__(self,gene):
        self.gene = gene
        self.strand = gene.strand
        self.transcript = {}
#TODO mangage incompleteness

### EOF ###
