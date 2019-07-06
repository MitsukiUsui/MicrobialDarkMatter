from logging import getLogger

logger = getLogger(__name__)


class GffRecord:
    def __init__(self, line):
        line = line.strip()
        line_split = line.split('\t')
        assert len(line_split) == 9

        self.seqid = line_split[0]
        self.source = line_split[1]
        self.type = line_split[2]
        self.start = int(line_split[3])
        self.end = int(line_split[4])
        self.score = float(line_split[5])
        self.strand = line_split[6]
        self.phase = line_split[7]
        self.attributes = self.decode_attributes(line_split[8])

    def __str__(self):
        return "\t".join([self.seqid, self.source, self.type, str(self.start), str(self.end), str(self.score),
                          self.strand, self.phase, self.encode_attributes(self.attributes)])

    def __repr__(self):
        return "<GffRecord({} {}:{})>".format(self.seqid, self.start, self.end)

    @staticmethod
    def decode_attributes(text):
        atts = {}
        for subtext in text.split(';'):
            if len(subtext) > 0:
                try:
                    key, val = subtext.split('=')
                    atts[key] = val
                except ValueError as e:
                    logger.debug("found unparsable attribute subtext {}".format(subtext), e)
        return atts

    @staticmethod
    def encode_attributes(atts):
        text = ""
        for key, val in atts.items():
            text += "{}={};".format(key, val)
        return text


def parse_gff(gff_fp):
    records = []
    with open(gff_fp, 'r') as f:
        for line in f:
            if len(line) > 0 and line[0] != '#':
                records.append(GffRecord(line))
    return records
