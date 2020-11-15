from idiva.db.clinvar import clinvar_to_df
from pathlib import Path
import os

"""
Script that reads a clinvar vcf file and saves it as structured data in a csv file. Can then be used for future queries.
"""

BASE = Path(__file__).parent.parent.parent / "data"
output_path = BASE / "clinvar.csv"

if not os.path.exists(BASE):
    os.mkdir(BASE)

df = clinvar_to_df()
df.to_csv(output_path)
