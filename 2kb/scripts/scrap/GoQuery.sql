use anopheles_gambiae_core_1312_73_3;

select distinct
    GENE_ID.stable_id as "ensembl.gene",
    RNA_ID.stable_id as "ensembl.transcript",
    PROT_ID.stable_id as "ensembl.translation",
    GO.acc as "go.acc",
    GO.name as "go.name",
    GOXREF.linkage_type as "evidence"
from
    ensembl_go_54.term as GO,
    external_db as EXTDB0,
    external_db as EXTDB1,
    object_xref as OX0,
    object_xref as OX1,
    xref as XREF0,
    xref as XREF1,
    transcript as RNA,
    transcript_stable_id as RNA_ID,
    gene as GENE,
    gene_stable_id as GENE_ID,
    translation as PROT,
    translation_stable_id as PROT_ID,
    go_xref as GOXREF
where
    XREF0.dbprimary_acc="NM_030621" and
    XREF0.external_db_id=EXTDB0.external_db_id and
    EXTDB0.db_name="RefSeq_dna" and
    OX0.xref_id=XREF0.xref_id and
    RNA.gene_id=GENE.gene_id and
    GENE.gene_id= GENE_ID.gene_id and
    RNA.transcript_id=OX0.ensembl_id and
    RNA_ID.transcript_id=RNA.transcript_id and
    PROT.transcript_id = RNA.transcript_id and
    OX1.ensembl_id=PROT.translation_id and
    PROT.translation_id=PROT_ID.translation_id and
    OX1.ensembl_object_type='Translation' and
    OX1.xref_id=XREF1.xref_id and
    GOXREF.object_xref_id=OX1.object_xref_id and 
    XREF1.external_db_id=EXTDB1.external_db_id and
    EXTDB1.db_name="GO" and
    GO.acc=XREF1.dbprimary_acc

order by GO.acc;