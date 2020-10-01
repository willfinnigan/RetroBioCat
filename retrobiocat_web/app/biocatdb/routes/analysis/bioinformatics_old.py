

def create_alignment_jobs(enzyme_type):
    current_app.app_context().push()

    enz_type_obj = EnzymeType.objects(enzyme_type=enzyme_type)[0]

    make_db_job = current_app.alignment_queue.enqueue(all_by_all_blast.make_blast_db_for_enzyme_type, enzyme_type)

    seqs = Sequence.objects(db.Q(enzyme_type=enzyme_type) & db.Q(bioinformatics_ignore__ne=True))
    seqs_for_alignment = []
    for seq in seqs:
        if seq.sequence != '' and seq.sequence is not None and seq.blast is None:
            seqs_for_alignment.append(seq)
        else:
            seq.alignments_made = datetime.datetime.now()
            seq.save()

    last_alignment_job = None
    for seq in seqs_for_alignment:
        last_alignment_job = current_app.alignment_queue.enqueue(all_by_all_blast.single_enz_all_by_all_blast, seq.enzyme_name,
                                                                 'biocatdb', enzyme_type, depends_on=make_db_job)
    for seq in UniRef90.objects(enzyme_type=enz_type_obj):
        last_alignment_job = current_app.alignment_queue.enqueue(all_by_all_blast.single_enz_all_by_all_blast, seq.enzyme_name,
                                                                 'uniref', enzyme_type, depends_on=make_db_job)

    current_app.task_queue.enqueue(set_bioinformatics_status, enzyme_type, 'Idle', depends_on=last_alignment_job)