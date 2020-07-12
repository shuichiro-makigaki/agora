from io import StringIO
import os
from pathlib import Path
import tempfile
import shutil

from Bio import SeqIO, SearchIO, SeqRecord, AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import generic_protein
from Bio.Blast.Applications import NcbipsiblastCommandline, NcbideltablastCommandline, NcbimakeblastdbCommandline
from tqdm.auto import tqdm
import numpy as np
import pandas


# delta_results = {}
# for record in tqdm(list(SeqIO.parse('evaluation.fasta', 'fasta'))):
#     stdout, _ = NcbideltablastCommandline(
#         rpsdb=Path(os.getenv('TMPDIR'))/'blastdb'/'cdd_delta',
#         db=Path(os.getenv('TMPDIR'))/'blastdb'/'scop95', evalue=9999, outfmt=5, num_threads=8
#     )(stdin=record.format('fasta'))
#     delta_results[record.id] = stdout
# np.save('data/delta_scop95.npy', np.array(delta_results))


class AgoraModeling:
    def _remove_gaps(self, record_a: SeqRecord, record_b: SeqRecord):
        seq_a = record_a.seq.tomutable()
        seq_b = record_b.seq.tomutable()
        i = 0
        while i < len(seq_a):
            if seq_a[i] in ['-', '.'] and seq_b[i] in ['-', '.']:
                seq_a[i] = '@'
                seq_b[i] = '@'
            i += 1
        while True:
            try:
                seq_a.remove('@')
            except ValueError:
                break
        while True:
            try:
                seq_b.remove('@')
            except ValueError:
                break
        i = 0
        while i < len(seq_a):
            if seq_a[i] == '.':
                seq_a[i] = '-'
            if seq_b[i] == '.':
                seq_b[i] = '-'
            i += 1
        seq_a = seq_a.toseq()
        seq_a.alphabet = generic_protein
        seq_b = seq_b.toseq()
        seq_b.alphabet = generic_protein
        assert len(seq_a) == len(seq_b)
        record_a.seq = seq_a
        record_b.seq = seq_b
        record_a = record_a.upper()
        record_b = record_b.upper()
        return record_a, record_b

    def modeller_automodel(self, query: SeqRecord, results: Path, num_align: int, atom_files_dir: Path):
        from modeller import environ, log
        from modeller.automodel import automodel
        raw_df = pandas.read_csv('data/delta_new_hits_thresh10_xdug10_xdg10_xdgf10.csv')
        for row in tqdm(raw_df.itertuples(), total=raw_df.shape[0]):
            try:
                aln = SearchIO.read(StringIO(row.XML), 'blast-xml').hsps[0].aln
            except IndexError:
                continue
            assert aln[0].id == row.Query and aln[1].id == row.Hit
            q_rec, t_rec = self._remove_gaps(aln[0], aln[-1])
            q_rec.name = ''
            t_rec.name = ''
            q_rec.description = f'sequence:{row.Query}::::::::'
            t_rec.description = f'structureX:{row.Hit}::{row.Hit[-2].upper()}::{row.Hit[-2].upper()}::::'
            aln = MultipleSeqAlignment([q_rec, t_rec])
            out_d = results.resolve()
            cwd = os.getcwd()
            with tempfile.TemporaryDirectory() as tmpdir:
                try:
                    os.chdir(tmpdir)
                    AlignIO.write(aln, 'aln.pir', 'pir')
                    log.none()
                    env = environ()
                    env.io.atom_files_directory = [(atom_files_dir/aln[1].id[2:4]).resolve().as_posix()]
                    mod = automodel(env, 'aln.pir', knowns=[aln[1].id], sequence=aln[0].id)
                    mod.make()
                    shutil.copy(list(Path().glob('*.B*.pdb'))[0], out_d/f'{aln[0].id}_{aln[1].id}.pdb')
                except:
                    pass
                finally:
                    os.chdir(cwd)


AgoraModeling().modeller_automodel(None, Path('data/evaluation/delta_new_hits_thresh10_xdug10_xdg10_xdgf10'), None, Path('/data/DB/ASTRAL/pdbstyle-2.07'))
