from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

import pandas
from Bio import SeqIO, SeqRecord
from tqdm.auto import tqdm

from agora import TMalignCommandLine


def get_hie(seq_index):
    hie = {}
    for name in seq_index:
        record = seq_index[name]
        sf = '.'.join(record.description.split()[1].split('.')[:3])
        hie.setdefault(sf, [])
        hie[sf].append(record)
    return hie


def _tmalign(data_dir: str, query: SeqRecord, rank: int):
    tmalign = TMalignCommandLine(f'data/evaluation/targets/{query.id}.ent', f'{data_dir}/{query.id}_{rank+1}.pdb')
    tmalign.run()
    return (tmalign.tmscore, rank, query)


def _process(data_dir: str, key: str):
    test_data = get_hie(SeqIO.index('evaluation.fasta', 'fasta'))
    with ProcessPoolExecutor() as executor:
        futures = []
        for rank in range(100):
            for sf_sccs in test_data:
                futures.append(executor.submit(_tmalign, data_dir, test_data[sf_sccs][0], rank))
        tmscore_d = {'Aligner': [], 'Rank': [], 'TM-score': [], 'PDB': []}
        for future in tqdm(as_completed(futures), total=len(futures)):
            data = future.result()
            if data[0] is not None:
                tmscore_d['Aligner'].append(key)
                tmscore_d['Rank'].append(data[1]+1)
                tmscore_d['PDB'].append(f'{data[2].id}_{data[1]+1}.pdb')
                tmscore_d['TM-score'].append(data[0][0])
        pandas.DataFrame.from_dict(tmscore_d).to_csv(f'{data_dir}/tmscore.csv')


def main():
    # data_dir = 'data/evaluation/delta_u50_50_u50_50_s95_500_evalue_sum_local'
    # if not Path(f'{data_dir}/tmscore.csv').exists():
    #     _process(data_dir, 'SW')
    # data_dir = 'data/evaluation/delta_u50_50_u50_50_s95_500_evalue_sum_i_local_p20'
    # if not Path(f'{data_dir}/tmscore.csv').exists():
    #     _process(data_dir, 'iSW_2_p20')
    # data_dir = 'data/evaluation/delta_u50_50_u50_50_s95_500_evalue_sum_i_local_p5'
    # if not Path(f'{data_dir}/tmscore.csv').exists():
    #     _process(data_dir, 'iSW_2_p5')
    # data_dir = 'data/evaluation/delta_u50_50_u50_50_s95_500_evalue_sum_i_local_p5p5p20'
    # if not Path(f'{data_dir}/tmscore.csv').exists():
    #     _process(data_dir, 'iSW_2_p5p5p20')
    # data_dir = 'data/evaluation/delta_u50_50_u50_50_s95_500_evalue_sum_local_i_local_p5p5p20'
    # if not Path(f'{data_dir}/tmscore.csv').exists():
    #     _process(data_dir, 'SWiSW_2_p5p5p20')
    # data_dir = 'data/evaluation/delta_u50_50_u50_50_s95_500_evalue_sum_i_blast'
    # if not Path(f'{data_dir}/tmscore.csv').exists():
    #     _process(data_dir, 'iBLAST_2')
    data_dir = 'data/evaluation/delta_u50_50_u50_50_s95_500_evalue_sum_i_blast'
    if not Path(f'{data_dir}/tmscore.csv').exists():
        _process(data_dir, 'iBLAST_2')


if __name__ == '__main__':
    main()
