from pathlib import Path
import random
import os
import logging

from Bio import SeqIO
from tqdm.auto import tqdm
import click
import numpy as np

from agora import AgoraDatabase, AgoraIntermediateSearch, AgoraSearch, AgoraModeling, AgoraResultSet


logging.basicConfig(level=logging.INFO)


def get_hh_scop40_hie():
    hie = {}
    for record in SeqIO.parse('data/hh_scop40.fasta', 'fasta'):
        sf = '.'.join(record.description.split()[1].split('.')[:3])
        hie.setdefault(sf, [])
        hie[sf].append(record)
    return hie


@click.group()
def main():
    pass


@main.command()
@click.option('--in-fasta', type=click.Path(exists=True), required=True, help='Query sequence')
@click.option('--out-dir', type=click.Path(), required=True, help='Directory for results and intermediate files')
@click.option('--last-max-hits', type=click.INT, default=500, show_default=True, help='Maximum number of hits')
@click.option('--last-label', type=click.STRING, default='PDB100', show_default=True, help='Node label for last layer')
@click.option('--last-db', type=click.STRING, default='pdbaa', show_default=True, help='Database name for last layer')
@click.option('--num-align', type=click.INT, default=100, show_default=True, help='Number of alignments generated in results')
@click.option('--overwrite', is_flag=True, show_default=True, help='Overwrite results or not')
@click.option('--num-threads', type=click.INT, default=int(os.cpu_count()), show_default=True, help='Number of used CPU cores')
def iss(in_fasta, out_dir, overwrite, num_threads, last_label, last_db, num_align, last_max_hits):
    cmd_iss(in_fasta, out_dir, overwrite, num_threads, last_label, last_db, last_max_hits)
    cmd_search(in_fasta, out_dir, overwrite, num_align, last_label, last_db)
    cmd_show_results(in_fasta, out_dir, last_db)


def cmd_iss(in_fasta, out_dir, overwrite, num_threads, last_label, last_db, last_max_hits):
    domains = list(SeqIO.parse(in_fasta, 'fasta'))
    # random.shuffle(domains)
    logging.info('Running intermediate sequence search')
    for q in tqdm(domains):
        ofile = Path(out_dir)/f'{q.id}.graphml'
        if ofile.exists() and not overwrite:
            continue
        AgoraIntermediateSearch(q, 'QUERY')\
            .deltablast('uniref50', 'cdd_delta', 'UniRef50', max_hits=50, num_threads=num_threads)\
            .deltablast('uniref50', 'cdd_delta', 'UniRef50', max_hits=50, num_threads=num_threads)\
            .deltablast(last_db, 'cdd_delta', last_label, max_hits=last_max_hits, num_threads=num_threads, last=True)\
            .out(ofile)


def cmd_search(in_fasta, graphml_dir, overwrite, num_align, last_label, last_db):
    domains = list(SeqIO.parse(in_fasta, 'fasta'))
    # random.shuffle(domains)
    logging.info('Generating alignments')
    for q in tqdm(domains):
        f_graphml = Path(graphml_dir)/f'{q.id}.graphml'
        f_npy = Path(graphml_dir)/f'{q.id}.npy'
        if f_npy.exists() and not overwrite:
            continue
        results = AgoraSearch().search_graph(q, f_graphml, num_align, last_label, last_db)
        np.save(f_npy, results)


def cmd_show_results(in_fasta, out_dir, last_db):
    for q in list(SeqIO.parse(in_fasta, 'fasta')):
        f_npy = Path(out_dir)/f'{q.id}.npy'
        result_set = AgoraResultSet()
        result_set.from_numpy_file(f_npy)
        for i, r in enumerate(result_set.items[1:]):
            record = AgoraSearch()._get_seq(last_db, r.template_id)
            print(f'No.{i+1} {record.id} | {record.name} {record.description}')
            print(format(r.alignment, 'clustal'))


@main.command()
@click.argument('in_fasta', type=click.Path(exists=True))
@click.argument('graphml_dir', type=click.Path(exists=True))
@click.argument('atom_files_dir', type=click.Path(exists=True))
@click.option('--num-align', type=click.INT, default=100)
def modeller_automodel(in_fasta, graphml_dir, atom_files_dir, num_align):
    domains = list(SeqIO.parse(in_fasta, 'fasta'))
    random.shuffle(domains)
    for q in tqdm(domains):
        AgoraModeling().modeller_automodel(q, Path(graphml_dir)/f'{q.id}.npy', num_align, Path(atom_files_dir))


# @main.command()
# @click.argument('in_file', type=click.Path(exists=True))
# @click.argument('blastdb', type=click.Path())
# @click.argument('cdd_delta', type=click.Path())
# @click.argument('structuredb', type=click.Path(exists=True))
# @click.argument('label', type=click.STRING)
# @click.argument('out_dir', type=click.Path())
# @click.option('--sub-dir-key-from', type=click.INT)
# @click.option('--sub-dir-key-to', type=click.INT)
def make_database(in_file, blastdb, cdd_delta, structuredb, label, out_dir, sub_dir_key_from, sub_dir_key_to):
    seq_index = SeqIO.index(in_file, 'fasta')
    seq_id_list = list(seq_index)
    random.shuffle(seq_id_list)
    for seq_id in tqdm(seq_id_list):
        p = Path(out_dir)/seq_id[sub_dir_key_from:sub_dir_key_to+1]/f'{seq_id}.graphml'
        if p.exists():
            continue
        AgoraDatabase(
            Path(blastdb), Path(cdd_delta), Path(structuredb), (sub_dir_key_from, sub_dir_key_to), label, None
        ).make_database(seq_index[seq_id], p)


# @main.command()
# @click.argument('in_file', type=click.Path(exists=True))
# @click.argument('graphml_dir', type=click.Path(exists=True))
# @click.option('--sub-dir-key-from', type=click.INT)
# @click.option('--sub-dir-key-to', type=click.INT)
def import_database(in_file, graphml_dir, sub_dir_key_from, sub_dir_key_to):
    agoradb = AgoraDatabase(Path(), Path(), Path(), (sub_dir_key_from, sub_dir_key_to), '')
    agoradb.clear_database()
    for seq_id in tqdm(list(SeqIO.index(in_file, 'fasta'))):
        p = Path(graphml_dir)/seq_id[sub_dir_key_from:sub_dir_key_to+1]/f'{seq_id}.graphml'
        if not p.exists():
            continue
        agoradb.import_graphml(p)


if __name__ == '__main__':
    main()
