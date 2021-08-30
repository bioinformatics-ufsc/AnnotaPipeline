# Criei os arquivos teste_celegans com poucas prot pra testar
# Substitui as anotaoes do .aa pra ser igual o do trangeli (anotacao completa)
# pq é isso que vai ta la qndo pega
# -----------------
# Cria um dicionario do .aa com a id (g1.t1) como key e a linha inteira como value (organism | ....)
# parseia o codingseq em fasta
# ve se a key (g10.g1) ta no id do .codingseq (>NC_003279.8.g10.t1)
# se tiver, busca o valor da key (a linha inteira) e troca o record.id
# salva e deu pra bola


# file into array of lines
from Bio import SeqIO

def annotate_codingseq(aa_fasta, codingseq_fasta, basename):
    # Linha de anotacao completa
    anno_all = [line.strip() for line in open(aa_fasta) if ">" in line]
    anno_all = [sub.replace(">", "") for sub in anno_all]
    # IDs simplificadas
    anno_ids = [name.split("|")[0].strip() for name in anno_all]

    corrrected_ids = dict(zip(anno_ids, anno_all))

    id_dict  = SeqIO.to_dict(SeqIO.parse(codingseq_fasta, "fasta"))

    corrected_fasta = str(f"AnnotaPipeline_{basename}_Transcripts.fasta")

    with open(corrected_fasta, "w") as corrected:
        for key, record in id_dict.items():
            for id_key in corrrected_ids.keys():
                if id_key in key:
                    record.id = corrrected_ids.get(id_key)
                    record.description = ""
                    record.seq = record.seq.upper()
                    SeqIO.write(record, corrected, "fasta-2line")

annotate_codingseq("teste_celegans.aa", "teste_celegans.codingseq", "teste_normal")

