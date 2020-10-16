import os
from pathlib import Path
from io import StringIO
import subprocess
from subprocess import PIPE
import re
import math
import tempfile
import shutil
import logging
import warnings

from Bio import BiopythonExperimentalWarning
warnings.simplefilter('ignore', BiopythonExperimentalWarning)

from Bio.Blast.Applications import NcbideltablastCommandline
from Bio import SearchIO, SeqIO, AlignIO, Align
from Bio.Align import substitution_matrices, PairwiseAligner
from Bio.SearchIO._model import HSP
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_protein
from Bio.AlignIO import MultipleSeqAlignment
from Bio.PDB import PDBParser, CaPPBuilder
import networkx
from neo4j import GraphDatabase
import numpy as np


class TMalignCommandLine:
    def __init__(self, protein_a=None, protein_b=None, binary='TMalign'):
        self.protein_A = protein_a
        self.protein_B = protein_b
        self.cmd = binary
        self.tmscore = None
        self.stdout = None
        self.stderr = None
        self.alignment = None

    def _parse(self):
        for i, l in enumerate(self.stdout):
            if re.match('^TM-score=', l):
                self.tmscore = (float(self.stdout[i].split(' ')[1]), float(self.stdout[i + 1].split(' ')[1]))
                break
        for i, l in enumerate(self.stdout):
            if re.match('^\(":" denotes', l):
                a = SeqRecord(Seq(self.stdout[i + 1], alphabet=generic_protein),
                              id=Path(self.protein_A).stem + '&' + Path(self.protein_B).stem,
                              description=f'TM-score={self.tmscore[0]}')
                b = SeqRecord(Seq(self.stdout[i + 3], alphabet=generic_protein),
                              id=Path(self.protein_B).stem + '&' + Path(self.protein_A).stem,
                              description=f'TM-score={self.tmscore[1]}')
                self.alignment = MultipleSeqAlignment([a, b])
                break

    def run(self):
        arg = [self.cmd, self.protein_A, self.protein_B]
        p = subprocess.run(arg, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        self.stdout = p.stdout.splitlines()
        self.stderr = p.stderr.splitlines()
        self._parse()


class AgoraIntermediateSearch:
    def __init__(self, query: SeqRecord, label: str):
        self.query = query
        self.intermediates = [self.query]
        self.graph = networkx.DiGraph()
        self.graph.add_node(self.query.id, labels=label)

    def deltablast(self, blastdb: str, rpsdb: str, label: str, max_hits: int, evalue=9999, num_threads=os.cpu_count(),
                   last=False):
        networkx.set_node_attributes(self.graph, {
            self.query.id: {'labels': self.graph.nodes[self.query.id]['labels'] + f':{label}'}})
        if last:
            stdout, _ = NcbideltablastCommandline(
                db=blastdb, rpsdb=rpsdb, num_threads=num_threads, evalue=evalue, outfmt=5
            )(stdin=self.query.format('fasta'))
            for hit in SearchIO.read(StringIO(stdout), 'blast-xml').hits[:max_hits]:
                if hit.id in self.graph:
                    networkx.set_node_attributes(self.graph,
                                                 {hit.id: {'labels': self.graph.nodes[hit.id]['labels'] + f':{label}'}})
                else:
                    self.graph.add_node(hit.id, labels=label)
                self.graph.add_edge(self.query.id, hit.id, evalue=float(hit.hsps[0].evalue))
        else:
            self.graph.add_edge(self.query.id, self.query.id, evalue=0.0)
        stdout, _ = NcbideltablastCommandline(
            db=blastdb, rpsdb=rpsdb, num_threads=num_threads, evalue=evalue, outfmt=5
        )(stdin='\n'.join([_.format('fasta') for _ in self.intermediates]))
        self.intermediates = []
        for results in SearchIO.parse(StringIO(stdout), 'blast-xml'):
            for hit in results.hits[:max_hits]:
                if hit.id in self.graph:
                    networkx.set_node_attributes(self.graph,
                                                 {hit.id: {'labels': self.graph.nodes[hit.id]['labels'] + f':{label}'}})
                else:
                    self.graph.add_node(hit.id, labels=label)
                self.graph.add_edge(results.id, hit.id, evalue=float(hit.hsps[0].evalue))
                self.intermediates.append(hit.hsps[0].hit)
        return self

    def out(self, path: Path):
        path.parent.mkdir(parents=True, exist_ok=True)
        networkx.write_graphml(self.graph, path.as_posix())


class AgoraDatabase:
    def __init__(self, blastdb: Path, rpsdb: Path, structuredb: Path, sub_dir_key_range: (int, int),
                 label: str, auth=('neo4j', 'password')):
        self.rpsdb = rpsdb
        self.blastdb = blastdb
        self.label = label
        self.structuredb = structuredb
        self.key_range = sub_dir_key_range
        if auth is None:
            self.driver = None
        else:
            self.driver = GraphDatabase.driver('bolt://localhost:7687', auth=auth,
                                               encrypted=False, max_retry_time=15)

    def deltablast(self, sequence: Seq):
        stdout, _ = NcbideltablastCommandline(
            rpsdb=self.rpsdb.as_posix(), db=self.blastdb.as_posix(),
            num_threads=os.cpu_count(), outfmt=5)(stdin=str(sequence))
        return SearchIO.read(StringIO(stdout), 'blast-xml').hits

    def make_database(self, seq_record: SeqRecord, graphml_path: Path):
        hits = [_ for _ in self.deltablast(seq_record.seq) if _.id != seq_record.id]
        graph = networkx.DiGraph()
        for hit in hits:
            tmalign = TMalignCommandLine(
                self.structuredb / seq_record.id[self.key_range[0]:self.key_range[1] + 1] / f'{seq_record.id}.ent',
                self.structuredb / hit.id[self.key_range[0]:self.key_range[1] + 1] / f'{hit.id}.ent')
            tmalign.run()
            if tmalign.tmscore is None:
                tmscore = -1
            else:
                tmscore = tmalign.tmscore[0]
            # ToDo: スコア正規化の方向は合ってる？
            graph.add_node(seq_record.id)
            graph.add_node(hit.id)
            graph.add_edge(seq_record.id, hit.id,
                           label='SIMILAR_TO', evalue=min([float(_.evalue) for _ in hit.hsps]), tmscore=tmscore)
        graphml_path.parent.mkdir(exist_ok=True, parents=True)
        networkx.write_graphml(graph, graphml_path)

    def clear_database(self):
        with self.driver.session() as session:
            session.run('MATCH (n) DETACH DELETE n')

    def import_graphml(self, graphml_path: Path):
        with self.driver.session() as session:
            session.run(
                f'CALL apoc.import.graphml("{graphml_path.resolve()}", {{storeNodeIds: true}})')


class AgoraSearch:
    def __init__(self, blastdb=None, rpsdb=None, neo4j_auth=None):
        self.rpsdb = rpsdb
        self.blastdb = blastdb
        if neo4j_auth is not None:
            self._driver = GraphDatabase.driver('bolt://localhost:7687', auth=neo4j_auth,
                                                encrypted=False, max_retry_time=3)

    def get_initial_nodes(self, seq_records: [SeqRecord]):
        stdout, _ = NcbideltablastCommandline(
            rpsdb=self.rpsdb, db=self.blastdb, outfmt=5, num_threads=os.cpu_count()
        )(stdin='\n'.join([_.format('fasta') for _ in seq_records]))
        result = {}
        for results in SearchIO.parse(StringIO(stdout), 'blast-xml'):
            result.setdefault(results.query.id, [])
            result[results.query.id] = [
                {'name': _.id, 'evalue': min([float(__.evalue) for __ in _.hsps]), 'tmscore': 0.0}
                for _ in results.hits if _.id != results.query.id]
        return result

    def exec_cypher(self, seq_records: [SeqRecord], initial_nodes):
        results = []
        with self._driver.session() as session:
            try:
                for seq_record in seq_records:
                    session.run('MATCH (n:Query{name:$qid}) DETACH DELETE n', qid=seq_record.id)
                    for n in initial_nodes:
                        session.run('''
                            MATCH (n{name:$name}) MERGE (q:Query{name:$qid}) MERGE (q)-[:RELATED{evalue:$e,tmscore:$t}]->(n)
                        ''', name=n['name'], e=n['weight'], t=n['tmscore'], qid=seq_record.id)
                    rows = session.run('''
                        MATCH path = (n:Query{name:$qid})-[*..4]->(hit)
                        WHERE ALL(x IN nodes(path)[1..] WHERE x.name <> $qid)
                        RETURN hit, reduce(w=1.0, r IN relationships(path) | log(w)+log(r.weight)) AS score, path
                        ORDER BY score
                    ''', qid=seq_record.id)
                    for row in rows:
                        if row['hit']['description'] in [_['hit'] for _ in results]:
                            continue
                        path = {
                            'nodes': [_['name'] for _ in row['path'].nodes],
                            'relationships': [float(_['weight']) for _ in row['path'].relationships]
                        }
                        results.append({'hit': row['hit']['description'], 'score': float(row['score']), 'path': path})
                        if len(results) > 200:
                            break
            except Exception as e:
                print(e)
        return results

    def _get_seq(self, db_name: str, seq_id: str) -> SeqRecord:
        res = subprocess.run(['blastdbcmd', '-db', db_name, '-entry', seq_id],
                             stdout=PIPE, stderr=PIPE, universal_newlines=True, check=True)
        return SeqIO.read(StringIO(res.stdout), 'fasta')

    def _local_align(self, record_a: SeqRecord, record_b: SeqRecord, open_gap_score: int):
        aligner = Align.PairwiseAligner()
        aligner.mode = 'local'
        aligner.substitution_matrix = substitution_matrices.load('BLOSUM62')
        aligner.open_gap_score = open_gap_score
        aligner.extend_gap_score = -1
        aln = aligner.align(record_a.seq.ungap('-').upper(), record_b.seq.ungap('-').upper())[0]
        seq_a = Seq(str(aln).splitlines()[0].replace(' ', '-'), generic_protein)
        seq_b = Seq(str(aln).splitlines()[2].replace(' ', '-'), generic_protein)
        return MultipleSeqAlignment(
            [SeqRecord(seq_a, id=record_a.id), SeqRecord(seq_b, id=record_b.id)],
            annotations={'score': aln.score, 'path': aln.path, 'aligned': aln.aligned}
        )

    def _deltablast_alignment(self, record_a: SeqRecord, record_b: SeqRecord):
        with tempfile.TemporaryDirectory() as t:
            tmpdir = Path(t)
            SeqIO.write(record_a, tmpdir/'query.fasta', 'fasta')
            SeqIO.write(record_b, tmpdir/'subject.fasta', 'fasta')
            out, _ = NcbideltablastCommandline(query=(tmpdir/'query.fasta').as_posix(),
                                               subject=(tmpdir/'subject.fasta').as_posix(),
                                               use_sw_tback=True, outfmt=5, rpsdb='cdd_delta')()
        return SearchIO.read(StringIO(out), 'blast-xml')

    def _extend_region(self, query: SeqRecord, hit: SeqRecord, hsp: HSP, padding: int):
        q_s = hsp.query_start - padding if hsp.query_start - padding >= 0 else 0
        q_e = hsp.query_end + padding if hsp.query_end + padding <= len(query) else len(query)
        h_s = hsp.hit_start - padding if hsp.hit_start - padding >= 0 else 0
        h_e = hsp.hit_end + padding if hsp.hit_end + padding <= len(hit) else len(hit)
        return SeqRecord(query.seq[q_s:q_e], id=query.id), SeqRecord(hit.seq[h_s:h_e], id=hit.id)

    def _mafft_merge2(self, sequences: [SeqRecord]):
        with tempfile.TemporaryDirectory() as t:
            tmpdir = Path(t)
            SeqIO.write(sequences, tmpdir/'input.fasta', 'fasta')
            (tmpdir/'table').write_text('1 2\n3 4\n5 6')
            res = subprocess.run(['mafft', '--clustalout', '--merge', 'table', 'input.fasta'],
                                 cwd=tmpdir, stdout=PIPE, stderr=PIPE, universal_newlines=True, check=True)
        return res.stdout, res.stderr

    def _mafft_merge1(self, sequences: [SeqRecord]):
        with tempfile.TemporaryDirectory() as t:
            tmpdir = Path(t)
            SeqIO.write(sequences, tmpdir/'input.fasta', 'fasta')
            (tmpdir/'table').write_text('1 2\n3 4')
            res = subprocess.run(['mafft', '--clustalout', '--merge', 'table', 'input.fasta'],
                                 cwd=tmpdir, stdout=PIPE, stderr=PIPE, universal_newlines=True, check=True)
        return res.stdout, res.stderr

    def search_layer2_evalue_sum_i_local(self, query: SeqRecord, graphml: Path, num_align: int):
        graph = networkx.read_graphml(graphml)
        path = []
        for nei1 in graph.neighbors(query.id):
            if 'UniRef50' not in graph.nodes[nei1]['labels'].split(':'):
                continue
            for nei2 in graph.neighbors(nei1):
                if 'UniRef50' not in graph.nodes[nei2]['labels'].split(':'):
                    continue
                for nei3 in graph.neighbors(nei2):
                    if 'SCOP95' not in graph.nodes[nei3]['labels'].split(':'):
                        continue
                    score = graph.get_edge_data(query.id, nei1)['evalue'] \
                        + graph.get_edge_data(nei1, nei2)['evalue'] + graph.get_edge_data(nei2, nei3)['evalue']
                    path.append([query.id, None, nei1, None, nei2, None, nei3, None, score])
        # dedup
        results, seen = [], []
        for p in sorted(path, key=lambda _: _[-1]):
            if p[6] in seen:
                continue
            seen.append(p[6])
            results.append(p)
        for r in results[:num_align]:
            n1_seq = query if r[2] == query.id else self._get_seq('uniref50', r[2])
            n2_seq = query if r[4] == query.id else self._get_seq('uniref50', r[4])
            n3_seq = query if r[6] == query.id else self._get_seq('scop95', r[6])
            # 1st
            try:
                hsp1 = self._deltablast_alignment(query, n1_seq)[0][0]
            except IndexError:
                logging.warning(f'No HSP: {query.id} -X-> {n1_seq.id} ---> {n2_seq.id} ---> {n3_seq.id}')
                continue
            i_query, i_n1_seq = self._extend_region(query, n1_seq, hsp1, 5)
            aln1 = self._local_align(i_query, i_n1_seq, -11)
            r[1] = (aln1.format('clustal'), aln1.annotations)
            # 2nd
            try:
                hsp2 = self._deltablast_alignment(i_n1_seq, n2_seq)[0][0]
            except IndexError:
                logging.warning(f'No HSP: {query.id} ---> {n1_seq.id} -X-> {n2_seq.id} ---> {n3_seq.id}')
                continue
            i_n1_seq, i_n2_seq = self._extend_region(i_n1_seq, n2_seq, hsp2, 5)
            aln2 = self._local_align(i_n1_seq, i_n2_seq, -11)
            r[3] = (aln2.format('clustal'), aln2.annotations)
            # 3rd
            try:
                hsp3 = self._deltablast_alignment(i_n2_seq, n3_seq)[0][0]
            except IndexError:
                logging.warning(f'No HSP: {query.id} ---> {n1_seq.id} ---> {n2_seq.id} -X-> {n3_seq.id}')
                continue
            i_n2_seq, i_n3_seq = self._extend_region(i_n2_seq, n3_seq, hsp3, 20)
            aln3 = self._local_align(i_n2_seq, i_n3_seq, -11)
            r[5] = (aln3.format('clustal'), aln3.annotations)
            aln1 = AlignIO.read(StringIO(r[1][0]), 'clustal')
            aln2 = AlignIO.read(StringIO(r[3][0]), 'clustal')
            aln3 = AlignIO.read(StringIO(r[5][0]), 'clustal')
            r[7] = self._mafft_merge2([aln1[0], aln1[1], aln2[0], aln2[1], aln3[0], aln3[1]])

        return results

    def search_layer2_evalue_sum_local_i_local(self, query: SeqRecord, graphml: Path, num_align: int, last_label, last_db):
        graph = networkx.read_graphml(graphml)
        path = []
        for nei1 in graph.neighbors(query.id):
            if 'UniRef50' not in graph.nodes[nei1]['labels'].split(':'):
                continue
            for nei2 in graph.neighbors(nei1):
                if 'UniRef50' not in graph.nodes[nei2]['labels'].split(':'):
                    continue
                for nei3 in graph.neighbors(nei2):
                    if last_label not in graph.nodes[nei3]['labels'].split(':'):
                        continue
                    score = graph.get_edge_data(query.id, nei1)['evalue'] \
                        + graph.get_edge_data(nei1, nei2)['evalue'] + graph.get_edge_data(nei2, nei3)['evalue']
                    path.append([query.id, None, nei1, None, nei2, None, nei3, None, score])
        # dedup
        results, seen = [], []
        for p in sorted(path, key=lambda _: _[-1]):
            if p[6] in seen:
                continue
            seen.append(p[6])
            results.append(p)
        for r in results[:num_align+1]:
            # Query
            # intermediate 1: n1
            # intermediate 2: n2
            # Hit           : n3
            try:
                n1_seq = query if r[2] == query.id else self._get_seq('uniref50', r[2])
                n2_seq = query if r[4] == query.id else self._get_seq('uniref50', r[4])
                n3_seq = query if r[6] == query.id else self._get_seq(last_db, r[6])
            except Exception as e:
                print(e)
                return results
            # Smith-Waterman
            aln_local = self._local_align(query, n3_seq, -11)
            r[7] = (aln_local.format('clustal'), 'local', aln_local.annotations, None, None)
            # 1st
            try:
                hsp1 = self._deltablast_alignment(query, n1_seq)[0][0]
            except IndexError:
                continue
            i_query, i_n1_seq = self._extend_region(query, n1_seq, hsp1, 5)
            aln1 = self._local_align(i_query, i_n1_seq, -11)
            r[1] = (aln1.format('clustal'), aln1.annotations)
            # 2nd
            try:
                hsp2 = self._deltablast_alignment(i_n1_seq, n2_seq)[0][0]
            except IndexError:
                continue
            i_n1_seq, i_n2_seq = self._extend_region(i_n1_seq, n2_seq, hsp2, 5)
            aln2 = self._local_align(i_n1_seq, i_n2_seq, -11)
            r[3] = (aln2.format('clustal'), aln2.annotations)
            # 3rd
            try:
                hsp3 = self._deltablast_alignment(i_n2_seq, n3_seq)[0][0]
            except IndexError:
                continue
            i_n2_seq, i_n3_seq = self._extend_region(i_n2_seq, n3_seq, hsp3, 20)
            aln3 = self._local_align(i_n2_seq, i_n3_seq, -11)
            r[5] = (aln3.format('clustal'), aln3.annotations)
            ### Merge
            # Q->n1->n2->n3 tag=ilocal
            aln1 = AlignIO.read(StringIO(r[1][0]), 'clustal')
            aln2 = AlignIO.read(StringIO(r[3][0]), 'clustal')
            aln3 = AlignIO.read(StringIO(r[5][0]), 'clustal')
            out, _ = self._mafft_merge2([aln1[0], aln1[1], aln2[0], aln2[1], aln3[0], aln3[1]])
            aln_ilocal = AlignIO.read(StringIO(out), 'clustal')
            # Q->n3: aln_local tag=local
            # Q->n1->n3 aln1+aln4 tag=ilocal1
            aln4 = self._local_align(i_n1_seq, i_n3_seq, -11)
            out, _ = self._mafft_merge1([aln1[0], aln1[1], aln4[0], aln4[1]])
            aln_ilocal1 = AlignIO.read(StringIO(out), 'clustal')
            # Q->n2->n3 aln5+aln3 tag=ilocal2
            aln5 = self._local_align(i_query, i_n2_seq, -11)
            out, _ = self._mafft_merge1([aln5[0], aln5[1], aln3[0], aln3[1]])
            aln_ilocal2 = AlignIO.read(StringIO(out), 'clustal')
            # ToDo: Select from aln_ilocal, aln_local, aln_ilocal1, aln_ilocal2
            best = sorted([
                {'key': 'local', 'value': self._aligned_number(aln_local)},
                {'key': 'ilocal', 'value': self._aligned_number(aln_ilocal)},
                {'key': 'ilocal1', 'value': self._aligned_number(aln_ilocal1)},
                {'key': 'ilocal2', 'value': self._aligned_number(aln_ilocal2)}
            ], key=lambda _: _['value'], reverse=True)[0]
            if best['key'] == 'local':
                r[7] = (aln_local.format('clustal'), 'local', aln_local.annotations, aln_ilocal.format('clustal'), aln_ilocal1.format('clustal'), aln_ilocal2.format('clustal'))
            elif best['key'] == 'ilocal1':
                r[7] = (aln_ilocal1.format('clustal'), 'ilocal1', None, aln_local.format('clustal'), aln_local.annotations, aln_ilocal.format('clustal'), aln_ilocal2.format('clustal'))
            elif best['key'] == 'ilocal2':
                r[7] = (aln_ilocal2.format('clustal'), 'ilocal2', None, aln_local.format('clustal'), aln_local.annotations, aln_ilocal.format('clustal'), aln_ilocal1.format('clustal'))
            else:
                r[7] = (aln_ilocal.format('clustal'), 'ilocal', None, aln_local.format('clustal'), aln_local.annotations, aln_ilocal1.format('clustal'), aln_ilocal2.format('clustal'))

        return results

    def _aligned_number(self, aln: MultipleSeqAlignment):
        return len([i for i in range(aln.get_alignment_length()) if '-' not in aln[:, i]])

    def search_layer1_evalue_sum_i_local(self, query: SeqRecord, graphml: Path, num_align: int, padding: int):
        graph = networkx.read_graphml(graphml)
        path = []
        for nei1 in graph.neighbors(query.id):
            if 'UniRef50' not in graph.nodes[nei1]['labels'].split(':'):
                continue
            for nei2 in graph.neighbors(nei1):
                if 'SCOP95' not in graph.nodes[nei2]['labels'].split(':'):
                    continue
                score = graph.get_edge_data(query.id, nei1)['evalue'] + graph.get_edge_data(nei1, nei2)['evalue']
                path.append([query.id, None, nei1, None, None, None, nei2, None, score])
        path = sorted(path, key=lambda _: _[-1])
        # dedup
        results, seen = [], []
        for p in path:
            if p[6] in seen:
                continue
            seen.append(p[6])
            results.append(p)
        for r in results[:num_align]:
            n1_seq = query if r[2] == query.id else self._get_seq('uniref50', r[2])
            n2_seq = query if r[6] == query.id else self._get_seq('scop95', r[6])
            # 1st
            try:
                hsp1 = self._deltablast_alignment(query, n1_seq)[0][0]
            except IndexError:
                logging.warning(f'No HSP: {query.id} -X-> {n1_seq.id} ---> {n2_seq.id}')
                continue
            i_query, i_n1_seq = self._extend_region(query, n1_seq, hsp1, padding)
            aln1 = self._local_align(i_query, i_n1_seq, -11)
            r[1] = (aln1.format('clustal'), aln1.annotations)
            # 2nd
            try:
                hsp2 = self._deltablast_alignment(i_n1_seq, n2_seq)[0][0]
            except IndexError:
                logging.warning(f'No HSP: {query.id} ---> {n1_seq.id} -X-> {n2_seq.id}')
                continue
            i_n1_seq, i_n2_seq = self._extend_region(i_n1_seq, n2_seq, hsp2, padding)
            aln2 = self._local_align(i_n1_seq, i_n2_seq, -11)
            r[3] = (aln2.format('clustal'), aln2.annotations)
            aln1 = AlignIO.read(StringIO(r[1][0]), 'clustal')
            aln2 = AlignIO.read(StringIO(r[3][0]), 'clustal')
            r[7] = self._mafft_merge1([aln1[0], aln1[1], aln2[0], aln2[1]])

        return results

    def search_layer2_evalue_sum_i_blast(self, query: SeqRecord, graphml: Path, num_align: int):
        graph = networkx.read_graphml(graphml)
        path = []
        for nei1 in graph.neighbors(query.id):
            if 'UniRef50' not in graph.nodes[nei1]['labels'].split(':'):
                continue
            for nei2 in graph.neighbors(nei1):
                if 'UniRef50' not in graph.nodes[nei2]['labels'].split(':'):
                    continue
                for nei3 in graph.neighbors(nei2):
                    if 'SCOP95' not in graph.nodes[nei3]['labels'].split(':'):
                        continue
                    score = graph.get_edge_data(query.id, nei1)['evalue'] \
                            + graph.get_edge_data(nei1, nei2)['evalue'] + graph.get_edge_data(nei2, nei3)['evalue']
                    path.append([query.id, None, nei1, None, nei2, None, nei3, None, score])
        path = sorted(path, key=lambda _: _[-1])
        # dedup
        results, seen = [], []
        for p in path:
            if p[-3] in seen:
                continue
            seen.append(p[-3])
            results.append(p)
        for r in results[:num_align]:
            n1_seq = query if r[2] == query.id else self._get_seq('uniref50', r[2])
            n2_seq = query if r[4] == query.id else self._get_seq('uniref50', r[4])
            n3_seq = query if r[6] == query.id else self._get_seq('scop95', r[6])
            with tempfile.TemporaryDirectory() as t:
                tmpdir = Path(t)
                SeqIO.write(query, tmpdir / 'query.fasta', 'fasta')
                SeqIO.write(n1_seq, tmpdir / 'n1.fasta', 'fasta')
                r[1], _ = NcbideltablastCommandline(query=(tmpdir / 'query.fasta').as_posix(),
                                                    subject=(tmpdir / 'n1.fasta').as_posix(),
                                                    use_sw_tback=True, outfmt=5, rpsdb='cdd_delta')()
                aln1 = SearchIO.read(StringIO(r[1]), 'blast-xml')[0][0].aln
                SeqIO.write(aln1, tmpdir / '1.aln', 'fasta')
                SeqIO.write(aln1[1], tmpdir / 'n1.fasta', 'fasta')
                SeqIO.write(n2_seq, tmpdir / 'n2.fasta', 'fasta')
                r[3], _ = NcbideltablastCommandline(query=(tmpdir / 'n1.fasta').as_posix(),
                                                    subject=(tmpdir / 'n2.fasta').as_posix(),
                                                    use_sw_tback=True, outfmt=5, rpsdb='cdd_delta')()
                aln2 = SearchIO.read(StringIO(r[3]), 'blast-xml')[0][0].aln
                SeqIO.write(aln2, tmpdir / '2.aln', 'fasta')
                SeqIO.write(aln2[1], tmpdir / 'n2.fasta', 'fasta')
                SeqIO.write(n3_seq, tmpdir / 'n3.fasta', 'fasta')
                r[5], _ = NcbideltablastCommandline(query=(tmpdir / 'n2.fasta').as_posix(),
                                                    subject=(tmpdir / 'n3.fasta').as_posix(),
                                                    use_sw_tback=True, outfmt=5, rpsdb='cdd_delta')()
                try:
                    aln3 = SearchIO.read(StringIO(r[5]), 'blast-xml')[0][0].aln
                    SeqIO.write(aln3, tmpdir / '3.aln', 'fasta')
                    SeqIO.write([aln1[0], aln1[1], aln2[0], aln2[1], aln3[0], aln3[1]], tmpdir / 'input.fasta', 'fasta')
                    (tmpdir / 'table').write_text('1 2\n3 4\n5 6')
                    res = subprocess.run(['mafft', '--clustalout', '--merge', 'table', 'input.fasta'],
                                         cwd=tmpdir, stdout=PIPE, stderr=PIPE, universal_newlines=True, check=True)
                    r[7] = (res.stdout, res.stderr)
                except:
                    r[7] = (None, None)
        return results

    def search_layer2_evalue_sum_i_blast_merged(self, query: SeqRecord, graphml: Path, num_align: int):
        graph = networkx.read_graphml(graphml)
        path = []
        for nei1 in graph.neighbors(query.id):
            if 'UniRef50' not in graph.nodes[nei1]['labels'].split(':'):
                continue
            for nei2 in graph.neighbors(nei1):
                if 'UniRef50' not in graph.nodes[nei2]['labels'].split(':'):
                    continue
                for nei3 in graph.neighbors(nei2):
                    if 'SCOP95' not in graph.nodes[nei3]['labels'].split(':'):
                        continue
                    score = graph.get_edge_data(query.id, nei1)['evalue'] \
                            + graph.get_edge_data(nei1, nei2)['evalue'] + graph.get_edge_data(nei2, nei3)['evalue']
                    path.append([query.id, None, nei1, None, nei2, None, nei3, None, score])
        path = sorted(path, key=lambda _: _[-1])
        # dedup
        results, seen = [], []
        for p in path:
            if p[-3] in seen:
                continue
            seen.append(p[-3])
            results.append(p)
        for r in results[:num_align]:
            n1_seq = query if r[2] == query.id else self._get_seq('uniref50', r[2])
            n2_seq = query if r[4] == query.id else self._get_seq('uniref50', r[4])
            n3_seq = query if r[6] == query.id else self._get_seq('scop95', r[6])
            with tempfile.TemporaryDirectory() as t:
                tmpdir = Path(t)
                SeqIO.write(query, tmpdir / 'query.fasta', 'fasta')
                SeqIO.write(n1_seq, tmpdir / 'n1.fasta', 'fasta')
                r[1], _ = NcbideltablastCommandline(query=(tmpdir / 'query.fasta').as_posix(),
                                                    subject=(tmpdir / 'n1.fasta').as_posix(),
                                                    use_sw_tback=True, outfmt=5, rpsdb='cdd_delta')()
                hsp1 = SearchIO.read(StringIO(r[1]), 'blast-xml')[0][0]
                rec1, rec2 = hsp1.aln[0], hsp1.aln[1]
                seq1, seq2 = rec1.seq.tomutable(), rec2.seq.tomutable()
                aln1 = hsp1.aln
                SeqIO.write(aln1[1], tmpdir / 'n1.fasta', 'fasta')
                SeqIO.write(n2_seq, tmpdir / 'n2.fasta', 'fasta')
                r[3], _ = NcbideltablastCommandline(query=(tmpdir / 'n1.fasta').as_posix(),
                                                    subject=(tmpdir / 'n2.fasta').as_posix(),
                                                    use_sw_tback=True, outfmt=5, rpsdb='cdd_delta')()
                hsp2 = SearchIO.read(StringIO(r[3]), 'blast-xml')[0][0]
                rec3, rec4 = hsp2.aln[0], hsp2.aln[1]
                seq3, seq4 = rec3.seq.tomutable(), rec4.seq.tomutable()
                aln2 = hsp2.aln
                SeqIO.write(aln2[1], tmpdir / 'n2.fasta', 'fasta')
                SeqIO.write(n3_seq, tmpdir / 'n3.fasta', 'fasta')
                r[5], _ = NcbideltablastCommandline(query=(tmpdir / 'n2.fasta').as_posix(),
                                                    subject=(tmpdir / 'n3.fasta').as_posix(),
                                                    use_sw_tback=True, outfmt=5, rpsdb='cdd_delta')()
                try:
                    hsp3 = SearchIO.read(StringIO(r[5]), 'blast-xml')[0][0]
                except IndexError:
                    print(hsp1)
                    r[7] = (None, None)
                    continue
                rec5, rec6 = hsp3.aln[0], hsp3.aln[1]
                seq5, seq6 = rec5.seq.tomutable(), rec6.seq.tomutable()

                seq3 = '-' * hsp2.query_start + seq3
                seq4 = '-' * hsp2.query_start + seq4
                seq5 = '-' * (hsp2.query_start + hsp3.query_start) + seq5
                seq6 = '-' * (hsp2.query_start + hsp3.query_start) + seq6

                for i in range(hsp2.query_start, len(seq2)):
                    if i >= len(seq3):
                        seq3.append('-')
                        seq4.append('-')
                    elif seq2[i] == '-' and seq3[i] != '-':
                        seq3.insert(i, '-')
                        seq4.insert(i, '-')
                    elif seq3[i] == '-' and seq2[i] != '-':
                        seq2.insert(i, '-')
                        seq1.insert(i, '-')

                append = len(seq2) - len(seq3)
                seq3 = seq3 + '-' * append
                seq4 = seq4 + '-' * append

                for i in range(hsp2.query_start + hsp3.query_start, len(seq4)):
                    if i >= len(seq5):
                        seq5.append('-')
                        seq6.append('-')
                    elif seq4[i] == '-' and seq5[i] != '-':
                        seq5.insert(i, '-')
                        seq6.insert(i, '-')
                    elif seq5[i] == '-' and seq4[i] != '-':
                        seq4.insert(i, '-')
                        seq3.insert(i, '-')
                        seq2.insert(i, '-')
                        seq1.insert(i, '-')

                append = len(seq4) - len(seq5)
                seq5 = seq5 + '-' * append
                seq6 = seq6 + '-' * append

                rec1.seq = seq1
                rec2.seq = seq2
                rec3.seq = seq3
                rec4.seq = seq4
                rec5.seq = seq5
                rec6.seq = seq6
                try:
                    r[7] = (MultipleSeqAlignment([rec1, rec2, rec3, rec4, rec5, rec6]).format('clustal'), None)
                except ValueError:
                    print(hsp1)
                    r[7] = (None, None)
                    continue

        return results

    def search_layer1_evalue_sum_i_blast(self, query: SeqRecord, graphml: Path, num_align: int):
        graph = networkx.read_graphml(graphml)
        path = []
        for nei1 in graph.neighbors(query.id):
            if 'UniRef50' not in graph.nodes[nei1]['labels'].split(':'):
                continue
            for nei2 in graph.neighbors(nei1):
                if 'SCOP95' not in graph.nodes[nei2]['labels'].split(':'):
                    continue
                score = graph.get_edge_data(query.id, nei1)['evalue'] + graph.get_edge_data(nei1, nei2)['evalue']
                path.append([query.id, None, nei1, None, None, None, nei2, None, score])
        path = sorted(path, key=lambda _: _[-1])
        # dedup
        results, seen = [], []
        for p in path:
            if p[6] in seen:
                continue
            seen.append(p[6])
            results.append(p)
        for r in results[:num_align]:
            n1_seq = query if r[2] == query.id else self._get_seq('uniref50', r[2])
            n2_seq = query if r[6] == query.id else self._get_seq('scop95', r[6])
            with tempfile.TemporaryDirectory() as t:
                tmpdir = Path(t)
                SeqIO.write(query, tmpdir / 'query.fasta', 'fasta')
                SeqIO.write(n1_seq, tmpdir / 'n1.fasta', 'fasta')
                r[1], _ = NcbideltablastCommandline(query=(tmpdir / 'query.fasta').as_posix(),
                                                    subject=(tmpdir / 'n1.fasta').as_posix(),
                                                    use_sw_tback=True, outfmt=5, rpsdb='cdd_delta')()
                aln1 = SearchIO.read(StringIO(r[1]), 'blast-xml')[0][0].aln
                SeqIO.write(aln1[1], tmpdir / 'n1.fasta', 'fasta')
                SeqIO.write(n2_seq, tmpdir / 'n2.fasta', 'fasta')
                r[3], _ = NcbideltablastCommandline(query=(tmpdir / 'n1.fasta').as_posix(),
                                                    subject=(tmpdir / 'n2.fasta').as_posix(),
                                                    use_sw_tback=True, outfmt=5, rpsdb='cdd_delta')()
                try:
                    aln2 = SearchIO.read(StringIO(r[3]), 'blast-xml')[0][0].aln
                    SeqIO.write([aln1[0], aln1[1], aln2[0], aln2[1]], tmpdir / 'input.fasta', 'fasta')
                    (tmpdir / 'table').write_text('1 2\n3 4')
                    res = subprocess.run(['mafft', '--clustalout', '--merge', 'table', 'input.fasta'],
                                         cwd=tmpdir, stdout=PIPE, stderr=PIPE, universal_newlines=True, check=True)
                    r[7] = (res.stdout, res.stderr)
                except Exception as error:
                    print(f'{n1_seq}, {n2_seq}, {error}')
                    r[7] = (None, None)

        return results

    def search_layer2_evalue_sum_local(self, query: SeqRecord, graphml: Path, num_align: int):
        graph = networkx.read_graphml(graphml)
        path = []
        for nei1 in graph.neighbors(query.id):
            if 'UniRef50' not in graph.nodes[nei1]['labels'].split(':'):
                continue
            for nei2 in graph.neighbors(nei1):
                if 'UniRef50' not in graph.nodes[nei2]['labels'].split(':'):
                    continue
                for nei3 in graph.neighbors(nei2):
                    if 'SCOP95' not in graph.nodes[nei3]['labels'].split(':'):
                        continue
                    score = graph.get_edge_data(query.id, nei1)['evalue'] \
                            + graph.get_edge_data(nei1, nei2)['evalue'] + graph.get_edge_data(nei2, nei3)['evalue']
                    path.append([query.id, None, nei1, None, nei2, None, nei3, None, score])
        path = sorted(path, key=lambda _: _[-1])
        # dedup
        results, seen = [], []
        for p in path:
            if p[-3] in seen:
                continue
            seen.append(p[-3])
            results.append(p)
        for r in results[:num_align]:
            n3_seq = query if r[6] == query.id else self._get_seq('scop95', r[6])
            aln1 = self._local_align(query, n3_seq, -11)
            r[7] = (aln1.format('clustal'), aln1.annotations)
        return results

    def search_graphml_pagerank(self, query: str, graphml: Path):
        graph = networkx.read_graphml(graphml)
        for e, d in graph.edges.items():
            graph.edges[e]['r_evalue'] = math.exp(-d['evalue'] / 100)
        pagerank = networkx.pagerank_numpy(graph, weight='r_evalue', alpha=1.0)
        path = [(_, pagerank[_]) for _ in pagerank if 'SCOP95' in graph.nodes[_]['labels'].split(':')]
        path = sorted(path, key=lambda _: _[1])
        return path

    def search_graphml_clocent(self, query: str, graphml: Path):
        graph = networkx.read_graphml(graphml)
        path = []
        for n in [_ for _ in graph.nodes() if 'SCOP95' in graph.nodes[_]['labels'].split(':')]:
            cent = networkx.algorithms.centrality.closeness_centrality(graph, distance='evalue', u=n)
            path.append((n, cent))
        path = sorted(path, key=lambda _: _[1], reverse=True)
        return path

    def search_graph(self, query: SeqRecord, graphml: Path, num_align: int, last_label, last_db):
        return self.search_layer2_evalue_sum_local_i_local(query, graphml, num_align, last_label, last_db)


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

    def _remove_missing_res(self, record: SeqRecord, pdb: Path):
        structure = PDBParser().get_structure(record.id, pdb)
        sequence = ''.join([str(_.get_sequence()) for _ in CaPPBuilder().build_peptides(structure, aa_only=False)])
        path = PairwiseAligner().align(record.seq.ungap('-'), sequence)[0].path
        gaps = []
        for i, _ in enumerate(path[:-1]):
            if path[i][1] == path[i+1][1]:
                gaps.append((path[i][0], path[i+1][0]))
        gaps = list(reversed(gaps))
        mut = record.seq.tomutable()
        for gap in gaps:
            i = 0
            for k, res in enumerate(mut):
                if res == '-':
                    continue
                if gap[0] <= i < gap[1]:
                    mut[k] = '-'
                i += 1
        record.seq = mut.toseq()
        return record

    def modeller_automodel(self, query: SeqRecord, results: Path, num_align: int, atom_files_dir: Path):
        from modeller import environ
        from modeller.automodel import automodel
        for model_index, r in enumerate(np.load(results, allow_pickle=True)[:num_align]):
            try:
                aln = AlignIO.read(StringIO(r[-2][0]), 'clustal')
            except:
                logging.error(f'Failed to parse alignment: {r[0]} -> {r[2]} -> {r[4]} -> {r[6]}')
                continue
            assert query.id == aln[0].id and aln[-1].id == r[-3]
            q_rec, t_rec = self._remove_gaps(aln[0], aln[-1])
            try:
                t_rec = self._remove_missing_res(
                    t_rec, (atom_files_dir/aln[-1].id[2:4]/f'{aln[-1].id}.ent').resolve().as_posix())
            except FileNotFoundError as e:
                logging.exception(e)
                continue
            q_rec.name, t_rec.name = '', ''
            q_rec.description = f'sequence:{q_rec.id}::::::::'
            t_rec.description = f'structureX:{t_rec.id}::{t_rec.id[-2].upper()}::{t_rec.id[-2].upper()}::::'
            aln = MultipleSeqAlignment([q_rec, t_rec])
            out_d = results.resolve().parent
            if (out_d/f'{aln[0].id}_{model_index+1}.pdb').exists():
                continue
            cwd = os.getcwd()
            with tempfile.TemporaryDirectory() as tmpdir:
                try:
                    os.chdir(tmpdir)
                    AlignIO.write(aln, 'aln.pir', 'pir')
                    env = environ()
                    env.io.atom_files_directory = [(atom_files_dir/aln[1].id[2:4]).resolve().as_posix()]
                    mod = automodel(env, 'aln.pir', knowns=[aln[1].id], sequence=aln[0].id)
                    mod.make()
                    shutil.copy(list(Path().glob('*.B*.pdb'))[0], out_d/f'{aln[0].id}_{model_index+1}.pdb')
                except Exception as e:
                    logging.error(f'knowns=[{aln[1].id}], sequence={aln[0].id}')
                    logging.exception(e)
                finally:
                    os.chdir(cwd)


class AgoraResult:
    def __init__(self, template_id, aln_clustal):
        self.template_id = template_id
        msa = AlignIO.read(StringIO(aln_clustal), 'clustal', alphabet=generic_protein)
        start = 0
        for i, v in enumerate(msa[0]):
            if v != '-':
                start = i
                break
        end = None
        for i, v in enumerate(msa[0][::-1]):
            if v != '-':
                if i == 0:
                    end = None
                else:
                    end = -i
                break
        msa = msa[:, start:end]
        self.alignment = msa


class AgoraResultSet:
    def __init__(self):
        self.items = []

    def from_numpy_file(self, filename: Path):
        resultset = np.load(filename, allow_pickle=True)
        for r in resultset[:100]:
            self.items.append(AgoraResult(r[6], r[7][0]))
