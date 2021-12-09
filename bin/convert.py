#!/usr/bin/env python
import click
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, CompoundLocation
from pathlib import Path


@click.command()
@click.argument("gbk-file")
def fix_feature_type(gbk_file):
    assert Path(gbk_file).exists(), f"File {gbk_file} is not existed!"
    recs = [rec for rec in SeqIO.parse(gbk_file, "genbank")]
    if len(recs) > 0:
        for rec in recs:
            for feat in rec.features:
                feat.type = "CDS"
        with open(f"new_{Path(gbk_file).name.replace('.gbk', '')}.gbk", "w") as oh:
            SeqIO.write(rec, oh, "genbank")


@click.command()
@click.argument("gbk-file")
@click.option("--ori-gene", default=["ori", "p15A ori", "pSC101 ori"], help="Ori genes", show_default=True)
def rotate(gbk_file, ori_gene):
    assert Path(gbk_file).exists(), f"File {gbk_file} is not existed!"
    file_name = Path(gbk_file).name.split(".")[0]
    recs = [rec for rec in SeqIO.parse(gbk_file, "genbank")]
    suffix = "final"
    print(ori_gene)
    if len(recs) > 0:
        for rec in recs:
            # print(rec)
            rot_pos = 0
            for feat in rec.features:
                if feat.type == "rep_origin" and feat.qualifiers["label"]:
                    if feat.qualifiers["label"][0] in ori_gene:
                        print(feat.qualifiers["label"][0])
                        if isinstance(feat.location, CompoundLocation) and len(feat.location.parts) == 2:
                            if feat.location.strand == 1:
                                rot_pos = feat.location.parts[0].start
                            if feat.location.strand == -1:
                                rot_pos = feat.location.parts[1].start
                        elif isinstance(feat.location, FeatureLocation):
                            rot_pos = feat.location.start
                    
            header = f">{file_name}"
            new_seq = rec.seq[rot_pos:] + rec.seq[:rot_pos]
            print(rec.seq[:10])
            print(new_seq[:10])
            with open(f"{file_name}.{suffix}.fasta", "w") as fa:
                fa.write(header)
                fa.write("\n")
                fa.write(str(new_seq))
                fa.write("\n")


@click.group()
def cli():
    pass


cli.add_command(fix_feature_type)
cli.add_command(rotate)

if __name__ == "__main__":
    cli()
