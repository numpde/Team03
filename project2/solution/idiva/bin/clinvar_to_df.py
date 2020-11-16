import os
from pathlib import Path

from idiva.db.clinvar import clinvar_to_df

"""
Script that reads a clinvar vcf file and saves it as structured data in a csv file. Can then be used for future queries.
It is not recommended to run this script as it takes too much time. 
The clinvar file can be downloaded from https://www.dropbox.com/s/u8wttgtkvv4mpde/clinvar.gzip
"""

BASE = Path(__file__).parent.parent.parent / "data"
output_path = BASE / "clinvar.gzip"

if not os.path.exists(BASE):
    os.mkdir(BASE)

clinvar_to_df(output_path, make_checkpoints=True, resume_checkpoint=BASE / 'clinvarchckpt0.gzip')
