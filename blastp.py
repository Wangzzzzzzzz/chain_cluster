from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
from io import StringIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import SeqIO
from Bio.PDB import PDBList, PDBParser
from Bio.PDB.Polypeptide import PPBuilder, three_to_one
import os
import numpy as np 
import pandas as pd
import json

# extract the sequences based on the dataset we have 
def obtain_seq(score, sequences):
    score  = './dataFile/' + score 
    sequences = './dataFile/' + sequences
    score_file = pd.read_csv(score, sep='\t')
    wt_1 = score_file.iloc[:,0]
    wt_2 = score_file.iloc[:,1]
    # concatinate the wt_1 and wt_2 file 
    wild_type = list(wt_1) + list(wt_2)
    # get only the unique items and then sort them
    wild_type = set(wild_type)
    
    # read in the sequence file and compile into a dictionary
    seq_dict = {}
    with open(sequences, 'r') as s:
        for l in s:
            # consider only the wildtype
            content = l.split('\t')
            sequence = content[1].replace(' ','').replace('\n','')
            header = content[0]
            header_info = header.split('_')
            if (len(header_info) == 2) and (header in wild_type):
                seq_dict[header] = sequence

    return seq_dict

def write_to_fasta(seq_dict, fasta_name):
    # generate the SeqRecord 
    seqRec_list = [SeqRecord(Seq(seq_str, IUPAC.protein), id=seq_key)
                   for seq_key, seq_str in seq_dict.items()]
    with open(fasta_name, 'w') as output_hl:
        SeqIO.write(seqRec_list, output_hl, "fasta")

def run_blast(query, subject):
    output = NcbiblastpCommandline(
        query=query, subject=subject, outfmt=5, evalue=0.05)()[0]
    # read all parser result
    blast_result_record = list(NCBIXML.parse(StringIO(output)))
    return blast_result_record

def output_json(blast_res, output_name):
    #output the result in JSON file, compile the results into a conplex dictionary
    dict_blast_res = {}
    # each object in the list is the alignments find for one
    # one sequence in the query set
    for query_seq in blast_res:
        # query seq's name
        query_seq_name = query_seq.query[0:6]
        if len(query_seq.alignments):
            # create an empty dictionary for that query seq & add this query seq to dictionary
            dict_blast_res[query_seq_name] = dict()
        else:
            continue
        # go through all alignments of all one sequence in query
        for alignment in query_seq.alignments:
            # count of hsp, since it does not have unique identifier
            hsp_id = 1
            # subject seq (aligned sequence in database)'s name
            subjt_seq_name = alignment.title[0:6]
            if len(alignment.hsps):
                # add the subjt name to the sub dict, and create an empty dict for hsp
                dict_blast_res[query_seq_name][subjt_seq_name] = dict()
            else:
                continue
            # go through all of the high score pairs (represent the regions in the hit seq
            # that contains significant alignment, usually there is only one)
            for hsp in alignment.hsps:
                hsp_discriptor = {
                    "score": hsp.score,
                    "E_Value": hsp.expect,
                    "Query Cover": hsp.align_length/query_seq.query_length*100,
                    "Perc_Identity": hsp.identities/hsp.align_length*100,
                    "Num_of_Identity": hsp.identities,
                    "Aligned_length": hsp.align_length,
                    "Align_Seq_Length": alignment.length,
                    "Query_Seq_Length": query_seq.query_length,
                    "query_start_pos": hsp.query_start,
                    "aligned_query_seq": hsp.query,
                    "aligned_match_pts": hsp.match,
                    "aligned_sbjct_seq": hsp.sbjct,
                    "sbjct_start_pos":hsp.sbjct_start,
                }
                # filter of the Perc_Identity, change if needed
                if hsp_discriptor['Perc_Identity'] < 25:
                    continue
                # add the hsp_id to dict & attach the discriptor of hsp to the dict
                dict_blast_res[query_seq_name][subjt_seq_name][hsp_id] = hsp_discriptor
                hsp_id += 1

            # remove empty dicts
            if not dict_blast_res[query_seq_name][subjt_seq_name]:
                del dict_blast_res[query_seq_name][subjt_seq_name]

        if not dict_blast_res[query_seq_name]:
            del dict_blast_res[query_seq_name]

    #if there are something to output
    if dict_blast_res:
        with open('./blast_rs/' + output_name, 'w') as output_hl:
            json.dump(dict_blast_res, output_hl,indent=4)
        return 1
    
    # if there is nothing to output:
    else:
        return 0


def generate_seq_file(score_file, save_file):
    score_file = './dataFile/' + score_file
    sf = pd.read_csv(score_file, sep='\t')
    mut_chains = sf.iloc[:,0]

    mut_dict = dict()
    mut_track = set()
    pdb_track = set()
    for chain in mut_chains:
        info = chain.split('_')
        pdb_id = info[0]
        chain_id = info[1]
        wt_aa = info[2][0:3]
        mu_aa = info[2][-3:]
        mu_pos = int(''.join(filter(lambda x: x.isdigit(), info[2])))
        if not chain in mut_track:
            mut_track.add(chain)
            if pdb_id in pdb_track:
                mut_dict[pdb_id].append({'chain_id':chain_id,
                                         'wt_aa': wt_aa,
                                         'mu_aa': mu_aa,
                                         'mu_pos': mu_pos,
                                         'name': chain})
            else:
                mut_dict[pdb_id] = [{'chain_id': chain_id,
                                     'wt_aa': wt_aa,
                                     'mu_aa': mu_aa,
                                     'mu_pos': mu_pos,
                                     'name': chain}]
                pdb_track.add(pdb_id)
    del mut_track
    del pdb_track
                
    parser = PDBParser()
    seq_builder = PPBuilder()
    pdb_dl_handle = PDBList()
    PDB_DIR = './dataFile/PDB_dl'
    # check if pdb file exists
    mut_collect = dict()
    for pdb_id in mut_dict.keys():
        if not os.path.exists(PDB_DIR+'/pdb'+pdb_id.lower()+'.ent'):
            pdb_dl_handle.retrieve_pdb_file(pdb_code=pdb_id, file_format='pdb', overwrite=False, pdir=PDB_DIR)
        pdb_file = PDB_DIR+'/pdb'+pdb_id.lower()+'.ent'
        model = parser.get_structure(pdb_id, pdb_file)[0]

        for mutation in mut_dict[pdb_id]:
            protein_chain = model[mutation['chain_id']]
            sequence = "".join([str(pp.get_sequence())
                                for pp in seq_builder.build_peptides(protein_chain)])
            sequence = sequence.replace('\n', '').replace(' ', '')
            assert sequence[mutation['mu_pos']-1] == three_to_one(mutation['wt_aa']), 'Wt amino acid failed to match'
            mut_Seq_list = list(sequence)
            mut_Seq_list[mutation['mu_pos']-1] = three_to_one(mutation['mu_aa'])
            mut_Seq = ''.join(mut_Seq_list)
            mut_collect[mutation['name']] = mut_Seq
    
    with open(save_file, 'w') as output_hl:
        for k, v in mut_collect.items():
            output_hl.write(k+'\t'+v+'\n')

        
# obtian seq_without seq file
def obtian_seq_wo_seq_file(score_file):
    score_file = './dataFile/' + score_file
    sf = pd.read_csv(score_file,sep='\t')
    chains_involved = sf.iloc[:,0]
    pdb = dict()
    pdb_track = set()
    for chain in chains_involved:
        chain_name = chain[0:6]
        pdb_name = chain[0:4]
        # if we encounter a old pdb
        if pdb_name in pdb_track:
            pdb[pdb_name].add(chain_name)
        # else, we have a new pdb
        else:
            # update the track file
            pdb_track.add(pdb_name)
            pdb[pdb_name] = {chain_name}

    # create the link to the PDB database and retrive all the file 
    # related to the files, store them locally under ./dataFile/PDB_dl/
    PDB_DIR = './dataFile/PDB_dl'
    if not os.path.exists(PDB_DIR):
        os.mkdir(PDB_DIR)
    # create the download handle
    pdb_dl_handle = PDBList()
    # download all of the pdb files
    for item in pdb.keys():
        if not os.path.exists(PDB_DIR+'/pdb'+item.lower()+'.ent'):
            pdb_dl_handle.retrieve_pdb_file(pdb_code=item, file_format='pdb', overwrite=False,pdir=PDB_DIR)
    
    # for each pdb, we will construct the sequence
    seq_dict = dict()
    parser = PDBParser()
    seq_builder = PPBuilder()
    # key is the pdb_id, value is the chain in a 
    for pdb_id, chain_names in pdb.items():
        pdb_file = PDB_DIR+'/pdb'+pdb_id.lower()+'.ent'
        model = parser.get_structure(pdb_id, pdb_file)[0]
        for chain in chain_names:
            # extract the last letter, which is the chain name
            chain_id = chain[-1]
            protein_chain = model[chain_id]
            sequence = "".join([str(pp.get_sequence()) for pp in seq_builder.build_peptides(protein_chain)])
            sequence = sequence.replace('\n','').replace(' ','') # clean the bad chars
            seq_dict[chain] = sequence

    return seq_dict


def main():
    if not os.path.exists('./fasta_db'):
        os.mkdir('./fasta_db')

    if not os.path.exists('./blast_rs'):
        os.mkdir('./blast_rs')

    # run blast on skempi_v1 and skempi_v2
    skempi_v1 = obtain_seq('./SKP1402m.ddg.txt','./SKP1402m.seq.txt')
    write_to_fasta(skempi_v1, './fasta_db/skempi_v1.fasta')
    skempi_v2 = obtain_seq('./3G_S487_test_dataset_top1.txt','./skempi_v2.singlemut.mut4.seq.txt')
    write_to_fasta(skempi_v2, './fasta_db/skempi_v2.fasta')
    blast_skempi_v1_skempi_v2 = run_blast(query='./fasta_db/skempi_v1.fasta',
                                          subject='./fasta_db/skempi_v2.fasta')
    blast_skempi_v1_skempi_v1 = run_blast(query='./fasta_db/skempi_v1.fasta',
                                          subject='./fasta_db/skempi_v1.fasta')
    blast_skempi_v2_skempi_v2 = run_blast(query='./fasta_db/skempi_v2.fasta',
                                          subject='./fasta_db/skempi_v2.fasta')

    # run blast on skempi_v1 and NM
    NM = obtain_seq('./NM_test.scores.txt','NM_test.seq.txt')
    write_to_fasta(NM, './fasta_db/NM.fasta')
    blast_skempi_v1_NM = run_blast(query='./fasta_db/skempi_v1.fasta',subject='./fasta_db/NM.fasta')

    # run blast on skepmi_v1 vs MDM2
    MDM2 = obtian_seq_wo_seq_file('features_MDM2-p53_test_dataset_top1.tsv')
    write_to_fasta(MDM2, './fasta_db/MDM2.fasta')
    blast_skempi_v1_MDM2 = run_blast(query='./fasta_db/skempi_v1.fasta', subject='./fasta_db/MDM2.fasta')
    generate_seq_file('features_MDM2-p53_test_dataset_top1.tsv', './dataFile/MDM2_test.seq.txt')


    # output the json file for all of the above results
    if not output_json(blast_skempi_v1_skempi_v1, 'skempi_v1_skempi_v1.json'):
        print('No alignments for skempi_v1 & skempi_v1')
    if not output_json(blast_skempi_v1_skempi_v2, 'skempi_v1_skempi_v2.json'):
        print('No alignments for skempi_v1 & skempi_v2')
    if not output_json(blast_skempi_v2_skempi_v2, 'skempi_v2_skempi_v2.json'):
        print('No alignments for skempi_v2 & skempi_v2')
    if not output_json(blast_skempi_v1_NM, 'skempi_v1_NM.json'):
        print('No alignments for skempi_v1 & NM')
    if not output_json(blast_skempi_v1_MDM2, 'skempi_v1_MDM2.json'):
        print('No alignments for skempi_v1 & MDM2')

if __name__ == "__main__":
    main()

