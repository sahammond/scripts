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
        self.locusid = ""
        self.product = ""
        self.name = ""
        self.note = ""
        self.incomplete = ""

        for entry in attr:
            if len(re.findall("ID=",entry)) > 0:
                self.id = re.sub("ID=","",entry)
            elif len(re.findall("Parent=",entry)) > 0:
                # some features have > 1 parent, sep by ","
                self.parent = re.sub("Parent=","",entry).split(",")
            elif len(re.findall("product=",entry)) > 0:
                self.product = re.sub("product=","",entry)
            elif len(re.findall("Name=",entry)) > 0:
                self.name = re.sub("Name=","",entry)
            elif len(re.findall("note=",entry)) > 0:
                self.note = re.sub("note=","",entry)


class Transcript(object):

    """Custom class to hold GFF transcript features.

    Stores the exons and cds in lists by position. Note that if the feature is on the
    negative strand then they will need to be printed out in reverse order.

    """

    def add_exon(self,feature):
        if len(self.exons) == 0:
            self.exons.append(feature)
        else:
            # TODO currently only compares to final exon in list for redundancy
            if int(feature.start) == int(self.exons[-1].start):
                print ("There is already an exon with that start in "
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
                print ("There is already a cds segment with that start in "
                        "this transcript. Please check your input.")
            elif int(feature.start) > int(self.cds[-1].start):
                self.cds.append(feature)
            else:
                for pos in range(0,len(self.cds)):
                    if int(feature.start) < int(self.cds[pos].start):
                        self.cds.insert(pos,feature)
                        break


    def check_start_stop(self,sequence):
        """Check for start and stop codons.

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


    def check_5prime_partial_mrna(self,feature):
        """Apply NCBI's partial marks to a segment, if applicable.

        NCBI's rules are basically that an mRNA is incomplete if it
        doesn't include a 5'/3' UTR. Therefore if the first exon
        segment's start is identical to the first CDS segment, then
        the mRNA is 5' incomplete. Likewise, if the last exon's
        end is the same as the last CDS segment's, then it's 3'
        incomplete.

        """
        # if there are no CDS, there can be no check
        if len(self.cds) < 1:
            return feature
        else:
            if int(self.exons[0].start) == int(self.cds[0].start):
                feature.start = "<"+feature.start
                return feature
            else:
                return feature


    def check_3prime_partial_mrna(self,feature):
        """Apply NCBI's partial marks to a segment, if applicable.

        NCBI's rules are basically that an mRNA is incomplete if it
        doesn't include a 5'/3' UTR. Therefore if the first exon
        segment's start is identical to the first CDS segment, then
        the mRNA is 5' incomplete. Likewise, if the last exon's
        end is the same as the last CDS segment's, then it's 3'
        incomplete.

        """
        # if there are no CDS, there can be no check
        if len(self.cds) < 1:
            return feature
        else:
            if int(self.exons[-1].start) == int(self.cds[-1].start):
                feature.start = ">"+feature.start
                return feature
            else:
                return feature


    def print_single_seg(self,feature,product_type=None,outform=None):
        """Print the only exon/CDS segment in a transcript.

        product_type = 'mRNA' or 'ncRNA'
        outform = 'gff' or 'tbl'

        """
        if outform is None:
            outform = 'tbl'
        if product_type is None:
            product_type = 'mRNA'
        if outform == 'gff':
            return feature.raw
        # TODO handle incomplete first seg
        elif outform == 'tbl':
            outbuff = []
            if feature.type == "exon":
                if feature.strand == "-":
                    if product_type == "mRNA":
                        outbuff.append("\t".join([feature.end,feature.start,"mRNA"]))
                    elif product_type == "ncRNA":
                        outbuff.append("\t".join([feature.end,feature.start,"ncRNA"]))
                if feature.strand == "+":
                    if product_type == "mRNA":
                        outbuff.append("\t".join([feature.start,feature.end,"mRNA"]))
                    elif product_type == "ncRNA":
                        outbuff.append("\t".join([feature.start,feature.end,"ncRNA"]))
                if product_type == "ncRNA":
                    outbuff.append("\t\t\tncRNA_class\tlncRNA")
                if self.transcript.product:
                    outbuff.append("\t".join(["\t\t\tproduct",self.transcript.product]))
                if self.transcript.note:
                    outbuff.append("\t".join(["\t\t\tnote",self.transcript.note]))
                if product_type == "mRNA":
                    outbuff.append("".join(["\t\t\tprotein_id\tgnl|BCGSC|",
                                            feature.parent[0]]))
                    outbuff.append("".join(["\t\t\ttranscript_id\tgnl|BCGSC|",
                                            feature.parent[0],"_mRNA\n"]))
                if product_type == "ncRNA":
                    outbuff.append("".join(["\t\t\ttranscript_id\tgnl|BCGSC|",
                                            feature.parent[0],"_ncRNA\n"]))
            elif feature.type == "CDS":
                if feature.strand == "-":
                    outbuff.append("\t".join([feature.end,feature.start,"CDS"]))
                else:
                    outbuff.append("\t".join([feature.start,feature.end,"CDS"]))

                if self.transcript.product:
                    outbuff.append("\t".join(["\t\t\tproduct",self.transcript.product]))
                if self.transcript.note:
                    outbuff.append("\t".join(["\t\t\tnote",self.transcript.note]))

                outbuff.append("\t\t\tcodon_start\t1")
                outbuff.append("".join(["\t\t\tprotein_id\tgnl|BCGSC|",
                                        feature.parent[0]]))
                outbuff.append("".join(["\t\t\ttranscript_id\tgnl|BCGSC|",
                                        feature.parent[0],"_mRNA\n"]))

            return "\n".join(outbuff)


    def print_internal_seg(self,feature,outform=None):
        """Print a given exon/CDS line or segment.

        Format = 'gff' or 'tbl'

        """
        if outform is None:
            outform = 'tbl'
        if outform == 'gff':
            return feature.raw
        elif outform == 'tbl':
            outbuff = []
            if feature.strand == "-":
                outbuff.append("\t".join([feature.end,feature.start]))
            else:
                outbuff.append("\t".join([feature.start,feature.end]))

            return "\n".join(outbuff)

    def print_first_seg(self,feature,product_type=None,outform=None):
        """Print the first exon/CDS line or segment.

        product_type = 'mRNA' or 'ncRNA'
        outform = 'gff' or 'tbl'

        """
        if outform is None:
            outform = 'tbl'
        if product_type is None:
            product_type = 'mRNA'
        if outform == 'gff':
            return feature.raw
        # TODO handle incomplete first seg
        elif outform == 'tbl':
            outbuff = []
            if feature.type == "exon":
                if feature.strand == "-":
                    if product_type == "mRNA":
                        outbuff.append("\t".join([feature.end,feature.start,"mRNA"]))
                    if product_type == "ncRNA":
                        outbuff.append("\t".join([feature.end,feature.start,"ncRNA"]))
                if feature.strand == "+":
                    if product_type == "mRNA":
                        outbuff.append("\t".join([feature.start,feature.end,"mRNA"]))
                    if product_type == "ncRNA":
                        outbuff.append("\t".join([feature.start,feature.end,"ncRNA"]))
            elif feature.type == "CDS":
                if feature.strand == "-":
                    outbuff.append("\t".join([feature.end,feature.start,"CDS"]))
                if feature.strand == "+": ###
                    outbuff.append("\t".join([feature.start,feature.end,"CDS"]))

            return "\n".join(outbuff)


    def print_last_seg(self,feature,product_type=None,outform=None):
        """Print the final exon/CDS line or segment.

        product_type = 'mRNA' or 'ncRNA'
        outform = 'gff' or 'tbl'

        """
        if outform is None:
            outform = 'tbl'
        if product_type is None:
            product_type = 'mRNA'
        if outform == 'gff':
            return feature.raw
        # TODO handle incomplete first seq
        elif outform == 'tbl':
            outbuff = []
            if feature.type == "exon":
                if feature.strand == "-":
                    outbuff.append("\t".join([feature.end,feature.start]))
                if feature.strand == "+":
                    outbuff.append("\t".join([feature.start,feature.end]))
                if self.transcript.product:
                    outbuff.append("\t".join(["\t\t\tproduct",self.transcript.product]))
                if self.transcript.note:
                    outbuff.append("\t".join(["\t\t\tnote",self.transcript.note]))
                if product_type == "mRNA":
                    outbuff.append("".join(["\t\t\tprotein_id\tgnl|BCGSC|",
                                            feature.parent[0]]))
                    outbuff.append("".join(["\t\t\ttranscript_id\tgnl|BCGSC|",
                                            feature.parent[0],"_mRNA\n"]))
                if product_type == "ncRNA":
                    outbuff.append("\t\t\tncRNA_class\tlncRNA")
                    outbuff.append("".join(["\t\t\ttranscript_id\tgnl|BCGSC|",
                                            feature.parent[0],"_ncRNA\n"]))
            elif feature.type == "CDS":
                if feature.strand == "-":
                    outbuff.append("\t".join([feature.end,feature.start]))
                if feature.strand == "+":
                    outbuff.append("\t".join([feature.start,feature.end]))
                if self.transcript.product:
                    outbuff.append("\t".join(["\t\t\tproduct",self.transcript.product]))
                if self.transcript.note:
                    outbuff.append("\t".join(["\t\t\tnote",self.transcript.note]))

                outbuff.append("\t\t\tcodon_start\t1")
                outbuff.append("".join(["\t\t\tprotein_id\tgnl|BCGSC|",
                                        feature.parent[0]]))
                outbuff.append("".join(["\t\t\ttranscript_id\tgnl|BCGSC|",
                                        feature.parent[0],"_mRNA\n"]))

            return "\n".join(outbuff)

    def print_transcript(self,product_type=None,outform=None):
        """Print a transcript with exons (and cds segments).

        product_type = 'mRNA' or 'ncRNA'
        outform = 'gff' or 'tbl'

        """
        if outform is None:
            outform = 'tbl'
        if product_type is None:
            product_type = 'mRNA'
        if outform == 'gff':
            outbuff = []
            outbuff.append(self.transcript.raw)

            for segment in self.exons:
                outbuff.append(segment.raw)

            if self.cds:
                for segment in self.cds:
                    outbuff.append(segment.raw)

            return "\n".join(outbuff)
        elif outform == 'tbl':
            outbuff = []
            if len(self.exons) == 1:
                return self.print_single_seg(self.exons[0],product_type,outform)
            if len(self.exons) > 1:
                if self.strand == "-":
                    if len(self.exons) == 2:
                        outbuff.append(self.print_first_seg(self.exons[1],product_type,
                                        outform))
                        outbuff.append(self.print_last_seg(self.exons[0],product_type,
                                        outform))
                    if len(self.exons) > 2:
                        outbuff.append(self.print_first_seg(self.exons[-1],product_type,
                                        outform))
                        for segment in range(0,len(self.exons))[1:-1][::-1]:
                           outbuff.append(self.print_internal_seg(self.exons[segment],
                                            outform))
                        else:
                            outbuff.append(self.print_last_seg(self.exons[0],
                                            product_type,outform))
                if self.strand == "+":
                    if len(self.exons) == 2:
                        outbuff.append(self.print_first_seg(self.exons[0],product_type,
                                        outform))
                        outbuff.append(self.print_last_seg(self.exons[1],product_type,
                                        outform))
                    if len(self.exons) > 2:
                        outbuff.append(self.print_first_seg(self.exons[0],product_type,
                                        outform))
                        for segment in range(0,len(self.exons))[1:-1]:
                            outbuff.append(self.print_internal_seg(self.exons[segment],
                                            outform))
                        else:
                            outbuff.append(self.print_last_seg(self.exons[-1],
                                            product_type,outform))
            if self.cds:
                if len(self.cds) == 1:
                    outbuff.append(self.print_single_seg(self.cds[0],product_type,
                                    outform))
                if len(self.cds) > 1:
                    if self.strand == "-":
                        if len(self.cds) == 2:
                            outbuff.append(self.print_first_seg(self.cds[1],
                                            product_type,outform))
                            outbuff.append(self.print_last_seg(self.cds[0],product_type,
                                            outform))
                        if len(self.cds) > 2:
                            outbuff.append(self.print_first_seg(self.cds[-1],
                                            product_type,outform))
                            for segment in range(0,len(self.cds))[1:-1][::-1]:
                                outbuff.append(self.print_internal_seg(self.cds[segment],
                                                outform))
                            else:
                                outbuff.append(self.print_last_seg(self.cds[0],
                                            product_type,outform))
                    if self.strand == "+":
                        if len(self.cds) == 2:
                            outbuff.append(self.print_first_seg(self.cds[0],product_type,
                                            outform))
                            outbuff.append(self.print_last_seg(self.cds[1],product_type,
                                            outform))
                        if len(self.cds) > 2:
                            outbuff.append(self.print_first_seg(self.cds[0],product_type,
                                            outform))
                            for segment in range(0,len(self.cds))[1:-1]:
                                outbuff.append(self.print_internal_seg(self.cds[segment],
                                                outform))
                            else:
                                outbuff.append(self.print_last_seg(self.cds[-1],
                                                product_type,outform))

            #print "\n".join(outbuff)
            return "\n".join(outbuff)


    def __init__(self,feature):
        self.strand = feature.strand
        self.transcript = feature
        self.exons = []
        self.cds = []
        self.start_complete = ''
        self.stop_complete = ''


class Gene(object):

    """Class to hold whole gene entries from gff3 file.

    gene:   a GFF object with type 'gene'

    Stores strand, type (mRNA, ncRNA, etc.), and transcripts. For each transcript,
    the exons, cds, and whether the cds is complete (at either end) can be added.
    The latter info is via the check_start_stop method, which requires a sequence string.

    """

    def check_start_stop(self,sequence):
        """Check for start and stop codons.

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


    def add_transcript(self,feature,f=0):
        """Add a transcript GFF object to the gene.

        Will be put into a dict called 'transcript', with key = transcript.id

        """
        if feature.id in self.transcript:
            print "Transcript already associated with gene."
        else:
            self.transcript[feature.id] = Transcript(feature)


    def print_gene(self,outform=None):
        """Print a gene line.

        Format = 'gff' or 'tbl'

        """
        if outform is None:
            outform = 'tbl'
        if outform == 'gff':
            return self.gene.raw
        elif outform == 'tbl':
            outbuff = []
            if self.gene.strand == "-":
                if self.gene.locusid:
                    outbuff.append("".join([self.gene.end,"\t",self.gene.start,"\tgene\n",
                                    "\t\t\tlocus_id\t",self.gene.locusid,"\n"]))
                else:
                    outbuff.append("".join([self.gene.end,"\t",self.gene.start,"\tgene\n",
                                    "\t\t\tlocus_id\t",self.gene.id,"\n"]))
            if self.gene.strand == "+":
                if self.gene.locusid:
                    outbuff.append("".join([self.gene.start,"\t",self.gene.end,"\tgene\n",
                                    "\t\t\tlocus_id\t",self.gene.locusid,"\n"]))
                else:
                    outbuff.append("".join([self.gene.start,"\t",self.gene.end,"\tgene\n",
                                    "\t\t\tlocus_id\t",self.gene.id,"\n"]))

            #print "\n".join(outbuff)
            return "\n".join(outbuff)

    def __init__(self,gene):
        self.gene = gene
        self.id = gene.id
        self.transcript = {}

# TODO mangage incompleteness

### EOF ###
