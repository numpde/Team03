# RA, 2020-11-05

from idiva import log

import contextlib
import re
import io
import typing

SEP = '\t'
TCGA = {"T", "C", "G", "A"}


# KEY_COLS = ["CHROM", "POS", "REF", "ALT"]


def parse_gt(gt: str) -> typing.Tuple[int, int]:
    (a, b) = gt.split("|")
    return (int(a), int(b))


def is_genomic_string(s: str):
    return set(s).issubset(TCGA)


class RawDataline:
    def __init__(self, line: str):
        def dot_is_none(s: str) -> str:
            return (None if (s == ".") else str(s))

        try:
            line = line.split(SEP)
            self.chrom = str(line[0])
            self.pos = int(line[1])
            self.id = dot_is_none(line[2])
            self.ref = str(line[3])
            self.alt = dot_is_none(line[4])
            self.qual = (None if (line[5] == ".") else float(line[5]))
            self.filter = str(line[6])
            self.info = str(line[7])
            self.format = (str(line[8]) if (len(line) >= 9) else None)
            self.samples = line[9:] or None
        except:
            log.error(F"Can't process {line}.")
            raise

    def __str__(self):
        fields = [
            self.chrom, self.pos, self.id, self.ref, self.alt,
            (
                (F"{self.qual:f}").rstrip('0').rstrip('.')  # This is not the VCF dot
                if (self.qual is not None)
                else "."  # This is the VCF dot
            ),
            self.filter, self.info,
            self.format
        ]
        if self.samples:
            fields += self.samples
        while fields[-1] is None:
            fields.pop()
        return "\t".join(map(str, fields))


def proxy(fd: io.TextIOBase):
    """
    Separates the meta-info, the header and the data lines.

    Returns a structure with fields
        meta, a list of meta-info lines
        header, the header as string
        datalines, an iterable over the remaining lines
        dataline_start_pos, start location in the stream

    The file must be open while using `datalines`.
    """

    from idiva.io import Oneliner

    class VCFProxy:
        meta: typing.List[str] = []
        header: str = None
        datalines: typing.Iterable[RawDataline] = None
        dataline_start_pos: int = None

    oneliner = Oneliner(fd, buffer_lines=False)

    # First read in the comments
    for line in oneliner:
        if line.startswith("##"):
            VCFProxy.meta.append(line[2:].strip())
        else:
            break

    # Assume now comes the header
    assert oneliner.last.startswith("#")
    VCFProxy.header = oneliner.last[1:]

    VCFProxy.dataline_start_pos = oneliner.fd.tell()
    oneliner.buffer_lines = True  # Ok after `tell`

    # The remainder are datalines that can be iterated
    VCFProxy.datalines = map(RawDataline, oneliner)

    return VCFProxy


class ReadVCF:
    default_header = "CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO".split(',')
    fd_tracer = set()

    def __init__(self, fd: io.TextIOBase, rewind=True):
        assert (fd.tell() == 0) or (not rewind), "Specify rewind=False for a used file descriptor."

        self._md5 = None

        self._fd = fd
        self.meta = {}
        self.header: list = []
        self.datalines: typing.Iterable[RawDataline] = []
        self.dataline_start_pos: int = None

        self._parse_proxy(proxy(fd))

    @classmethod
    @contextlib.contextmanager
    def open(ReadVCF, file, rewind=True):
        assert rewind
        if isinstance(file, ReadVCF):
            log.warning(F"Attempt to reopen VCF {file.fd}.")
            with file.rewind_when_done:
                yield file
        else:
            from idiva.io import open_maybe_gz
            with open_maybe_gz(file) as fd:
                vcf = ReadVCF(fd)
                with vcf.rewind_when_done:
                    yield vcf

    @property
    def md5(self) -> str:
        from idiva.utils.files import checksum_md5
        if not self._md5:
            self._md5 = checksum_md5(self._fd)
        return self._md5

    @property
    def fd(self) -> io.TextIOBase:
        return self._fd

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

        self.header = str(proxy.header).split(SEP)
        ndefault = len(self.default_header)
        (mandatory, extra) = (self.header[0:ndefault], self.header[ndefault:])
        assert (mandatory == self.default_header), F"Unexpected header: {mandatory} vs {self.default_header}"

        if extra:
            assert (extra[0] == "FORMAT")

            # https://www.internationalgenome.org/wiki/Analysis/vcf4.0/
            # [..] followed by a FORMAT column header,
            # then an arbitrary number of sample IDs.
            self.sample_ids = extra[1:]

        # Datalines forwarding

        self.datalines = proxy.datalines
        self.dataline_start_pos = proxy.dataline_start_pos

    def preload_all(self):
        from idiva.utils import seek_then_rewind
        with seek_then_rewind(self.fd, seek=None):
            self.datalines = list(iter(self.datalines))

        return self

    def __iter__(self) -> typing.Iterator[RawDataline]:
        return iter(self.datalines)

    @property
    @contextlib.contextmanager
    def rewind_when_done(self):
        from idiva.utils import seek_then_rewind
        with seek_then_rewind(self.fd, seek=None):
            yield


def parse_info(info: str):
    """
    HK, 2020-12-02
    RA, 2020-12-04
    """

    # See
    # https://github.com/numpde/Team03/blob/5de2c7d/project2/solution/idiva/db/dbSNP.py#L40
    raise NotImplementedError


def align(*, case: ReadVCF, ctrl: ReadVCF):
    raise NotImplementedError("See https://github.com/numpde/Team03/blob/fae828d/project2/solution/idiva/io/vcf.py#L202")

