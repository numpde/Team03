# RA, 2020-11-05


import re
import typing

SEP = '\t'


def dot_is_none(s: str) -> str:
    return (None if (s == ".") else str(s))


class RawDataline:
    def __init__(self, line: str):
        line = line.split(SEP)
        self.chrom = str(line[0])
        self.pos = int(line[1])
        self.id = dot_is_none(line[2])
        self.ref = str(line[3])
        self.alt = dot_is_none(line[4])
        self.qual = float(line[5])
        self.filter = str(line[6])
        self.info = str(line[7])
        self.format = (str(line[8]) if (len(line) >= 8) else None)
        self.samples = line[9:] or None

    def __str__(self):
        fields = [
            self.chrom, self.pos, self.id, self.ref, self.alt,
            (F"{self.qual:f}").rstrip('0').rstrip('.'),
            self.filter, self.info,
            self.format
        ]
        return "\t".join(map(str, fields + self.samples))


def proxy(fd: typing.TextIO):
    """
    Separates the meta-info, the header and the data lines.

    Returns a structure with fields
        meta, a list of meta-info lines
        header, the header as string
        datalines, an iterable over the remaining lines

    The file must be open while using `datalines`.
    """

    from idiva.io import Oneliner

    class _:
        meta: typing.List[str] = []
        header: str = None
        datalines: typing.Iterable[RawDataline] = None

    oneliner = Oneliner(fd)

    # First read in the comments
    for line in oneliner:
        if line.startswith("##"):
            _.meta.append(line[2:].strip())
        else:
            break

    # Assume now comes the header
    assert oneliner.last.startswith("#")
    _.header = oneliner.last[1:]

    # The remainder are datalines that can be iterated
    _.datalines = map(RawDataline, oneliner)

    return _


class ReadVCF:
    def __init__(self, fd: typing.TextIO):
        self.meta = {}
        self.datalines = None

        self._parse_proxy(proxy(fd))

    def _parse_proxy(self, proxy):
        # Parse the meta lines

        for line in proxy.meta:
            try:
                (k, v) = line.split("=", maxsplit=1)
            except ValueError:
                pass
            else:
                try:
                    (id, more) = re.fullmatch(r"^<ID=([^,]+),(.*)>$", v).groups()
                except AttributeError:
                    self.meta[k] = v
                else:
                    assignments = re.findall(r'([^,]+)[=]("[^"]+"|[^"][^,]*)', more)
                    assert (more == ",".join([F"{a}={b}" for (a, b) in assignments]))

                    more = dict(assignments)

                    K = 'Number'
                    if (K in more):
                        V = more[K]
                        # https://gist.github.com/inutano/f0a2f5c219ab4920c5b5
                        more[K] = (None if (V == ".") else V if (V in "ARG") else int(V))

                    # print(self.meta[k], more)
                    self.meta.setdefault(k, {})[id] = more

        # Parse the header line

        header = str(proxy.header)

        default_header = "CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO".split(SEP)
        (*mandatory, extra) = header.split(SEP, maxsplit=len(default_header))
        assert (mandatory == default_header)

        extra = extra.split(SEP)
        assert (extra[0] == "FORMAT")

        # https://www.internationalgenome.org/wiki/Analysis/vcf4.0/
        # [..] followed by a FORMAT column header,
        # then an arbitrary number of sample IDs.
        self.sample_ids = extra[1:]

        # Datalines forwarding

        self.datalines = proxy.datalines

    def __iter__(self):
        return iter(self.datalines)

