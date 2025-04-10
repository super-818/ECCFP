import re
import Bio.Seq
from .config import args, genome, fastq
class Base(object):
    def __init__(self):
        self.nt = {}
        self.last = ''

    def match(self, nt):
        if nt in self.nt:
            self.nt[nt] += 1
        else:
            self.nt[nt] = 1
        self.last = nt

    def insert(self, nt):
        self.nt[self.last] -= 1
        if self.last != '-':
            nt = self.last + nt
        if nt in self.nt:
            self.nt[nt] += 1
        else:
            self.nt[nt] = 1

    def call(self):
        alt = None
        cnt = 0
        depth = 0
        for nt in self.nt:
            depth += self.nt[nt]
            if self.nt[nt] > cnt:
                alt = nt
                cnt = self.nt[nt]
        return alt, cnt, depth

class Variant(object):
    def __init__(self, chr, pos, ref, alt, count, depth, site, strand):
        self.chr = chr
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.count = count
        self.depth = depth
        self.site = site
        self.strand = strand
        self.type = ''
        self.maxPos = pos

    def varType(self):
        if self.alt == '-':
            if len(self.ref) <= 50:
                if len(self.ref) == 1:
                    self.type = 'single_Base_InDel'
                else:
                    self.type = 'InDel'
            else:
                self.type = 'Deletion'
        else:
            if len(self.alt) == 1 and len(self.ref) == 1:
                if (self.ref == 'A' and self.alt == 'G') or (self.ref == 'G' and self.alt == 'A') or (
                        self.ref == 'C' and self.alt == 'T') or (self.ref == 'T' and self.alt == 'C'):
                    self.type = 'transition'
                elif (self.ref == 'A' and self.alt in ('C', 'T')) or (self.ref == 'G' and self.alt in ('C', 'T')) or (
                        self.ref == 'C' and self.alt in ('A', 'G')) or (self.ref == 'T' and self.alt in ('A', 'G')):
                    self.type = 'transversion'
                else:
                    self.type = 'single_Base_mutation'
            else:
                if len(self.alt) <= 50:
                    if len(self.alt) == 2 and self.alt[0] == self.ref and self.alt[1] == self.ref:
                        self.type = 'single_Base_InDel'
                    else:
                        self.type = 'InDel'
                else:
                    self.type = 'Insertion'
        return self

class consensusSequence(object):
    def __init__(self):
        self.genome = genome
        self.fastq = fastq
        self.sequence = []
        self.varient = []
        self.minDP = args.minDP
        self.minAF = args.minAF
        self.reference = ''

    def alignment(self, ref, frag, cigar_, seq):
        refChr, refStart, refEnd, refStrand = ref
        fragChr, fragStart, fragEnd, fragStrand = frag
        if fragChr != refChr or fragStart > refEnd or refStart > fragEnd:
            return False
        cigar_re = re.compile(r"(\d+)([MDI])")
        ces = [(int(m.group(1)), m.group(2)) for m in cigar_re.finditer(cigar_)]
        read_cur = 0
        if refStrand == '+':
            ref_cur = fragStart - refStart
        else:
            ref_cur = refEnd - fragEnd
            ces = ces[::-1]
        for step, op in ces:
            if op == 'M':
                for _ in range(step):
                    if ref_cur >= 0 and ref_cur < len(self.sequence):
                        self.sequence[ref_cur].match(seq[read_cur])
                    read_cur += 1
                    ref_cur += 1
            elif op == 'D':
                for _ in range(step):
                    if ref_cur >= 0 and ref_cur < len(self.sequence):
                        self.sequence[ref_cur].match('-')
                    ref_cur += 1
            elif op == 'I':
                if ref_cur > 0 and ref_cur <= len(self.sequence):
                    self.sequence[ref_cur - 1].insert(seq[read_cur:read_cur + step])
                read_cur += step
        return True

    def getConsensus(self, ref, readsID, unit, cigar, readsSite):
        chr_, start, end, strand = ref
        self.reference = self.genome[chr_][start - 1:end].seq
        if strand == '-':
            self.reference = str(Bio.Seq.Seq(self.reference).reverse_complement())
        self.sequence = [Base() for _ in range(len(self.reference))]
        for readname, frag, cigarUnit, site in zip(readsID, unit, cigar, readsSite):
            unitStart, unitEnd = site
            seq = self.fastq[readname].seq[unitStart - 1:unitEnd]
            _ = self.alignment(ref, frag, cigarUnit, seq)
        return self

    def qualified(self, var, pos):
        if var.depth < self.minDP:
            return False
        if float(var.count) / var.depth < self.minAF:
            return False
        # Filter del and ins in homopolymer region
        homoins = False
        if len(var.alt) > 1:
            homoins = True
            for i in range(1, len(var.alt)):
                if var.alt[i] != var.ref:
                    homoins = False
                    break
        if var.alt == '-' or homoins:
            if pos + 2 < len(self.reference) and self.reference[pos + 1] == var.ref and self.reference[
                pos + 2] == var.ref:
                return False
            if pos - 1 >= 0 and pos + 1 < len(self.reference) and self.reference[pos - 1] == var.ref and self.reference[
                pos + 1] == var.ref:
                return False
            if pos - 2 >= 0 and self.reference[pos - 1] == var.ref and self.reference[pos - 2] == var.ref:
                return False
        return True

    def getSeqVar(self, refSite, circle):
        seq = [''] * len(self.reference)
        var = []
        for i in range(len(self.reference)):
            alt, cnt, depth = self.sequence[i].call()
            ref = self.reference[i]

            if alt is None or alt == ref:
                seq[i] = ref
            else:
                chr_, start, end, strand = refSite
                if strand == '+':
                    pos = start + i
                else:
                    pos = end - i
                v = Variant(chr_, pos, ref, alt, cnt, depth, circle, strand)
                if self.qualified(v, i):
                    if alt != '-':
                        seq[i] = alt
                        var.append(v)
                    else:
                        if len(var) > 0 and var[-1].alt == '-' and abs(pos - var[-1].maxPos) == len(var[-1].ref):
                            if pos < var[-1].pos:
                                var[-1].pos = pos
                            var[-1].ref += ref
                            var[-1].count = min(var[-1].count, cnt)
                            var[-1].depth = min(var[-1].depth, depth)
                            var[-1] = var[-1].varType()
                        else:
                            var.append(v)
                else:
                    seq[i] = ref
        seq = ''.join(seq)
        return seq, var